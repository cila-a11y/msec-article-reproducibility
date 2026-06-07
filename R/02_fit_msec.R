# -----------------------------------------------------------------------------
# Maximum-likelihood fitting and likelihood-ratio tools for the MSEC model
# The optimizer uses closed constraints for delta_i >= 0; no softplus link.
# -----------------------------------------------------------------------------


msec_pack_par <- function(mu, sigma, rho, delta, alpha,
                          fix_delta0 = FALSE, fix_alpha0 = FALSE) {
  mu <- unname(as.numeric(mu))
  sigma <- unname(as.numeric(sigma))
  rho <- unname(as.numeric(rho))[1]
  delta <- unname(as.numeric(delta))
  alpha <- unname(as.numeric(alpha))

  x <- c(mu1 = mu[1], mu2 = mu[2],
         eta1 = log(sigma[1]), eta2 = log(sigma[2]),
         zeta = atanh(msec_clip(rho, -0.999, 0.999)))
  if (!fix_delta0) x <- c(x, delta1 = delta[1], delta2 = delta[2])
  if (!fix_alpha0) x <- c(x, alpha1 = alpha[1], alpha2 = alpha[2])
  x
}


msec_unpack_par <- function(par, fix_delta0 = FALSE, fix_alpha0 = FALSE) {
  nm <- names(par)
  par <- as.numeric(par)

  expected <- c('mu1', 'mu2', 'eta1', 'eta2', 'zeta')
  if (!fix_delta0) expected <- c(expected, 'delta1', 'delta2')
  if (!fix_alpha0) expected <- c(expected, 'alpha1', 'alpha2')

  if (is.null(nm)) {
    if (length(par) != length(expected)) {
      stop('Parameter vector has no names and has unexpected length')
    }
    names(par) <- expected
  } else {
    names(par) <- nm
  }

  getn <- function(name) {
    val <- par[name]
    if (length(val) != 1L || is.na(val)) stop('Missing parameter: ', name)
    unname(val)
  }

  mu <- c(getn('mu1'), getn('mu2'))
  sigma <- exp(c(getn('eta1'), getn('eta2')))
  rho <- tanh(getn('zeta'))
  delta <- if (fix_delta0) c(0, 0) else c(getn('delta1'), getn('delta2'))
  alpha <- if (fix_alpha0) c(0, 0) else c(getn('alpha1'), getn('alpha2'))
  delta <- pmax(delta, 0)
  Sigma <- msec_make_sigma(sigma, rho)
  lambda <- msec_lambda_from_alpha(alpha, Sigma)
  list(mu = mu, sigma = sigma, rho = rho, Sigma = Sigma,
       delta = delta, alpha = alpha, lambda = lambda)
}

msec_negloglik_par <- function(par, W, generator = 'gaussian', nu = 6,
                               fix_delta0 = FALSE, fix_alpha0 = FALSE,
                               eps = 1e-12) {
  obj <- try({
    pp <- msec_unpack_par(par, fix_delta0, fix_alpha0)
    ll <- msec_loglik(W, pp$mu, pp$Sigma, pp$delta, pp$alpha,
                      generator = generator, nu = nu, eps = eps)
    if (!is.finite(ll)) return(1e100)
    -ll
  }, silent = TRUE)
  if (inherits(obj, 'try-error') || !is.finite(obj)) 1e100 else obj
}

msec_default_starts <- function(W, fix_delta0 = FALSE, fix_alpha0 = FALSE,
                                max_starts = 4L) {
  W <- as.matrix(W)
  mu0 <- apply(W, 2L, median)
  deltas <- if (fix_delta0) {
    list(c(0, 0))
  } else {
    list(c(0.15, 0.15), c(0.5, 0.3), c(0.8, 0.5), c(0, 0), c(1.1, 0.7))
  }
  starts <- list()
  for (d in deltas) {
    T0 <- msec_q(sweep(W, 2L, mu0, '-'), d)
    S0 <- try(cov(T0), silent = TRUE)
    if (inherits(S0, 'try-error') || any(!is.finite(S0))) S0 <- diag(apply(T0, 2L, var, na.rm = TRUE), 2L)
    if (any(!is.finite(S0))) S0 <- diag(c(1, 1), 2L)
    S0 <- msec_near_pd_2x2(S0 + diag(1e-6, 2L))
    comps <- msec_sigma_to_components(S0)
    sigma0 <- pmax(c(comps['sigma1'], comps['sigma2']), 1e-4)
    rho0 <- msec_clip(comps['rho'], -0.95, 0.95)
    alpha_starts <- list(c(0, 0))
    if (!fix_alpha0) {
      Z <- scale(T0)
      sk <- suppressWarnings(colMeans(Z^3, na.rm = TRUE))
      if (all(is.finite(sk)) && sqrt(sum(sk^2)) > 1e-8) {
        dir <- sk / sqrt(sum(sk^2))
        alpha_starts <- c(alpha_starts, list(0.5 * dir, -0.5 * dir, 1.0 * dir))
      } else {
        alpha_starts <- c(alpha_starts, list(c(0.5, 0), c(0, 0.5), c(-0.5, 0.5)))
      }
    }
    for (a in alpha_starts) {
      starts[[length(starts) + 1L]] <- msec_pack_par(mu0, sigma0, rho0, d, a, fix_delta0, fix_alpha0)
    }
  }
  starts[seq_len(min(length(starts), max_starts))]
}

msec_fit <- function(W, generator = c('gaussian', 't'), nu = 6,
                     fix_delta0 = FALSE, fix_alpha0 = FALSE,
                     starts = NULL, max_starts = 4L, maxit = 600L,
                     max_delta = 5, max_abs_alpha = 20,
                     compute_hessian = FALSE, compute_grad = TRUE,
                     eps = 1e-12, verbose = FALSE) {
  generator <- match.arg(generator)
  W <- as.matrix(W)
  if (is.null(starts)) starts <- msec_default_starts(W, fix_delta0, fix_alpha0, max_starts)

  best <- NULL
  for (s in seq_along(starts)) {
    p0 <- starts[[s]]
    lower <- rep(-Inf, length(p0)); upper <- rep(Inf, length(p0))
    names(lower) <- names(upper) <- names(p0)
    lower[c('eta1', 'eta2')] <- log(1e-5)
    upper[c('eta1', 'eta2')] <- log(1e5)
    lower['zeta'] <- atanh(-0.999)
    upper['zeta'] <- atanh(0.999)
    if (!fix_delta0) {
      lower[c('delta1', 'delta2')] <- 0
      upper[c('delta1', 'delta2')] <- max_delta
    }
    if (!fix_alpha0) {
      lower[c('alpha1', 'alpha2')] <- -max_abs_alpha
      upper[c('alpha1', 'alpha2')] <- max_abs_alpha
    }
    p0 <- pmin(pmax(p0, lower), upper)
    opt <- try(optim(p0, msec_negloglik_par, W = W, generator = generator, nu = nu,
                     fix_delta0 = fix_delta0, fix_alpha0 = fix_alpha0, eps = eps,
                     method = 'L-BFGS-B', lower = lower, upper = upper,
                     control = list(maxit = maxit, factr = 1e7)), silent = TRUE)
    if (inherits(opt, 'try-error')) next
    if (is.null(best) || opt$value < best$value) best <- opt
    if (verbose) message('start ', s, ': nll=', signif(opt$value, 7), ' conv=', opt$convergence)
  }

  if (is.null(best)) {
    return(list(ok = FALSE, convergence = 999L, message = 'all starts failed', logLik = -Inf))
  }

  pp <- msec_unpack_par(best$par, fix_delta0, fix_alpha0)
  ll <- -best$value
  grad <- rep(NA_real_, length(best$par))
  grad_norm <- alpha_grad_norm <- NA_real_
  if (compute_grad) {
    grad <- try(msec_num_grad(msec_negloglik_par, best$par, W = W, generator = generator,
                              nu = nu, fix_delta0 = fix_delta0, fix_alpha0 = fix_alpha0, eps = eps),
                silent = TRUE)
    if (inherits(grad, 'try-error')) grad <- rep(NA_real_, length(best$par))
    names(grad) <- names(best$par)
    grad_norm <- sqrt(sum(grad^2, na.rm = TRUE))
    if (!fix_alpha0 && all(c('alpha1', 'alpha2') %in% names(grad))) {
      alpha_grad_norm <- sqrt(sum(grad[c('alpha1', 'alpha2')]^2, na.rm = TRUE))
    }
  }

  Hess <- NULL
  if (compute_hessian) {
    Hess <- try(optimHess(best$par, msec_negloglik_par, W = W, generator = generator,
                          nu = nu, fix_delta0 = fix_delta0, fix_alpha0 = fix_alpha0, eps = eps),
                silent = TRUE)
    if (inherits(Hess, 'try-error')) Hess <- NULL
    if (!is.null(Hess)) dimnames(Hess) <- list(names(best$par), names(best$par))
  }

  interest <- msec_interest(pp$mu, pp$Sigma, pp$delta, pp$alpha)
  list(ok = TRUE, convergence = best$convergence, message = best$message,
       value = best$value, logLik = ll, par = best$par,
       mu = pp$mu, sigma = pp$sigma, rho = pp$rho, Sigma = pp$Sigma,
       delta = pp$delta, alpha = pp$alpha, lambda = pp$lambda,
       interest = interest, grad = grad, grad_norm = grad_norm,
       alpha_grad_norm = alpha_grad_norm, Hessian = Hess,
       fix_delta0 = fix_delta0, fix_alpha0 = fix_alpha0,
       generator = generator, nu = nu)
}

msec_fit_success <- function(fit, grad_tol = 5e-2, alpha_grad_tol = 5e-3) {
  isTRUE(fit$ok) && identical(as.integer(fit$convergence), 0L) &&
    is.finite(fit$logLik) &&
    (is.na(fit$grad_norm) || fit$grad_norm < grad_tol) &&
    (is.na(fit$alpha_grad_norm) || fit$alpha_grad_norm < alpha_grad_tol)
}


msec_lr_lambda <- function(W, generator = 'gaussian', nu = 6,
                           max_starts = 4L, maxit = 600L, ...) {
  null <- msec_fit(W, generator, nu, fix_delta0 = FALSE, fix_alpha0 = TRUE,
                   max_starts = max_starts, maxit = maxit, ...)

  starts_full <- NULL
  if (isTRUE(null$ok)) {
    null_as_full <- msec_pack_par(
      mu = null$mu,
      sigma = null$sigma,
      rho = null$rho,
      delta = null$delta,
      alpha = c(0, 0),
      fix_delta0 = FALSE,
      fix_alpha0 = FALSE
    )
    default_starts <- msec_default_starts(
      W,
      fix_delta0 = FALSE,
      fix_alpha0 = FALSE,
      max_starts = max_starts
    )
    starts_full <- c(list(null_as_full), default_starts)
  }

  full <- msec_fit(W, generator, nu, fix_delta0 = FALSE, fix_alpha0 = FALSE,
                   starts = starts_full,
                   max_starts = max_starts, maxit = maxit, ...)

  raw_lr <- if (isTRUE(full$ok) && isTRUE(null$ok)) {
    2 * (full$logLik - null$logLik)
  } else {
    NA_real_
  }

  lr <- if (is.finite(raw_lr)) max(0, raw_lr) else NA_real_
  pval <- if (is.finite(lr)) pchisq(lr, df = 2, lower.tail = FALSE) else NA_real_

  list(full = full, null = null, LR = lr, LR_raw = raw_lr,
       nested_gap = if (isTRUE(full$ok) && isTRUE(null$ok)) full$logLik - null$logLik else NA_real_,
       p_chisq = pval)
}


msec_lr_delta <- function(W, generator = 'gaussian', nu = 6,
                          max_starts = 4L, maxit = 600L, ...) {
  null <- msec_fit(W, generator, nu, fix_delta0 = TRUE, fix_alpha0 = FALSE,
                   max_starts = max_starts, maxit = maxit, ...)

  starts_full <- NULL
  if (isTRUE(null$ok)) {
    null_as_full <- msec_pack_par(
      mu = null$mu,
      sigma = null$sigma,
      rho = null$rho,
      delta = c(0, 0),
      alpha = null$alpha,
      fix_delta0 = FALSE,
      fix_alpha0 = FALSE
    )
    default_starts <- msec_default_starts(
      W,
      fix_delta0 = FALSE,
      fix_alpha0 = FALSE,
      max_starts = max_starts
    )
    starts_full <- c(list(null_as_full), default_starts)
  }

  full <- msec_fit(W, generator, nu, fix_delta0 = FALSE, fix_alpha0 = FALSE,
                   starts = starts_full,
                   max_starts = max_starts, maxit = maxit, ...)

  raw_lr <- if (isTRUE(full$ok) && isTRUE(null$ok)) {
    2 * (full$logLik - null$logLik)
  } else {
    NA_real_
  }

  lr <- if (is.finite(raw_lr)) max(0, raw_lr) else NA_real_
  pval <- if (is.finite(lr)) pchisq(lr, df = 2, lower.tail = FALSE) else NA_real_

  list(full = full, null = null, LR = lr, LR_raw = raw_lr,
       nested_gap = if (isTRUE(full$ok) && isTRUE(null$ok)) full$logLik - null$logLik else NA_real_,
       p_chisq = pval)
}

msec_boot_lr_delta <- function(W, generator = 'gaussian', nu = 6, B_boot = 199L,
                               seed = NULL, max_starts = 3L, maxit = 500L,
                               verbose = FALSE, ...) {
  if (!is.null(seed)) set.seed(seed)
  obs <- msec_lr_delta(W, generator, nu, max_starts = max_starts, maxit = maxit, ...)
  if (!is.finite(obs$LR) || !isTRUE(obs$null$ok)) {
    return(list(obs = obs, boot_stats = rep(NA_real_, B_boot), p_boot = NA_real_))
  }
  th0 <- obs$null
  boot_stats <- rep(NA_real_, B_boot)
  for (b in seq_len(B_boot)) {
    Wb <- msec_r_msec(nrow(W), th0$mu, th0$Sigma, th0$delta, th0$alpha,
                      generator = generator, nu = nu)
    fb <- msec_lr_delta(Wb, generator, nu, max_starts = max_starts, maxit = maxit, ...)
    boot_stats[b] <- fb$LR
    if (verbose && b %% 25L == 0L) message('bootstrap ', b, '/', B_boot)
  }
  list(obs = obs, boot_stats = boot_stats,
       p_boot = msec_plus_one_pvalue(obs$LR, boot_stats))
}

msec_interest_from_internal <- function(par, fix_delta0 = FALSE, fix_alpha0 = FALSE) {
  pp <- msec_unpack_par(par, fix_delta0, fix_alpha0)
  msec_interest(pp$mu, pp$Sigma, pp$delta, pp$alpha)
}

msec_per_obs_scores_numeric <- function(W, fit, step = NULL, eps = 1e-12) {
  W <- as.matrix(W)
  par <- fit$par
  p <- length(par)
  n <- nrow(W)
  if (is.null(step)) step <- pmax(1e-5, abs(par) * 1e-5)
  S <- matrix(NA_real_, nrow = n, ncol = p)
  colnames(S) <- names(par)
  fvec <- function(pa) {
    pp <- msec_unpack_par(pa, fit$fix_delta0, fit$fix_alpha0)
    msec_loglik_vec(W, pp$mu, pp$Sigma, pp$delta, pp$alpha,
                    generator = fit$generator, nu = fit$nu, eps = eps)
  }
  for (j in seq_len(p)) {
    h <- step[j]
    p1 <- par; p2 <- par
    p1[j] <- p1[j] + h
    p2[j] <- p2[j] - h
    S[, j] <- (fvec(p1) - fvec(p2)) / (2 * h)
  }
  S
}

msec_wald_intervals_interest <- function(W, fit, level = 0.95,
                                         sandwich = FALSE, eps = 1e-12) {
  if (is.null(fit$Hessian)) {
    H <- try(optimHess(fit$par, msec_negloglik_par, W = W,
                       generator = fit$generator, nu = fit$nu,
                       fix_delta0 = fit$fix_delta0, fix_alpha0 = fit$fix_alpha0, eps = eps),
             silent = TRUE)
  } else {
    H <- fit$Hessian
  }
  if (inherits(H, 'try-error') || is.null(H) || any(!is.finite(H))) return(NULL)
  H <- (H + t(H)) / 2
  Binv <- try(solve(H), silent = TRUE)
  if (inherits(Binv, 'try-error')) Binv <- try(solve(H + diag(1e-6, nrow(H))), silent = TRUE)
  if (inherits(Binv, 'try-error')) return(NULL)
  Vpar <- Binv
  if (sandwich) {
    S <- try(msec_per_obs_scores_numeric(W, fit, eps = eps), silent = TRUE)
    if (!inherits(S, 'try-error')) {
      Meat <- t(S) %*% S
      Vpar <- Binv %*% Meat %*% Binv
    }
  }
  J <- try(msec_num_jacobian(msec_interest_from_internal, fit$par,
                             fix_delta0 = fit$fix_delta0, fix_alpha0 = fit$fix_alpha0),
           silent = TRUE)
  if (inherits(J, 'try-error')) return(NULL)
  Vint <- J %*% Vpar %*% t(J)
  est <- msec_interest_from_internal(fit$par, fit$fix_delta0, fit$fix_alpha0)
  se <- sqrt(pmax(diag(Vint), 0))
  z <- qnorm(1 - (1 - level) / 2)
  data.frame(parameter = names(est), estimate = as.numeric(est), se = se,
             lower = as.numeric(est) - z * se,
             upper = as.numeric(est) + z * se,
             row.names = NULL)
}

msec_t_weights <- function(W, fit) {
  pp <- fit
  Sinv <- solve(pp$Sigma)
  Y <- sweep(as.matrix(W), 2L, pp$mu, '-')
  Tm <- msec_q(Y, pp$delta)
  U <- rowSums((Tm %*% Sinv) * Tm)
  (pp$nu + 2) / (pp$nu + U)
}


msec_t_scale_update_once <- function(W, fit, denominator = c('n', 'sumw')) {
  denominator <- match.arg(denominator)
  if (fit$generator != 't') stop('Scale update ablation is only for t generator')
  Y <- sweep(as.matrix(W), 2L, fit$mu, '-')
  Tm <- msec_q(Y, fit$delta)
  ww <- msec_t_weights(W, fit)
  den <- if (denominator == 'n') nrow(W) else sum(ww)
  S <- crossprod(Tm * sqrt(ww)) / den
  S <- msec_near_pd_2x2(S)
  comps <- msec_sigma_to_components(S)
  c(sigma1 = unname(comps['sigma1']),
    sigma2 = unname(comps['sigma2']),
    rho = unname(comps['rho']))
}

