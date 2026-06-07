# -----------------------------------------------------------------------------
# Core probability functions for the bivariate MSEC model
# Corrected specification: both Gaussian and Student-t use the same
# skew-symmetric density 2 f0(t; Sigma, g) F_c(lambda' t).
# -----------------------------------------------------------------------------

msec_make_sigma <- function(sigma = c(1, 1), rho = 0) {
  sigma <- as.numeric(sigma)
  if (length(sigma) != 2L) stop('sigma must have length 2')
  if (any(!is.finite(sigma)) || any(sigma <= 0)) stop('sigma must be positive')
  rho <- as.numeric(rho)
  if (!is.finite(rho) || abs(rho) >= 1) stop('rho must be in (-1,1)')
  matrix(c(sigma[1]^2, rho * sigma[1] * sigma[2],
           rho * sigma[1] * sigma[2], sigma[2]^2), nrow = 2L, byrow = TRUE)
}

msec_sigma_to_components <- function(Sigma) {
  Sigma <- (Sigma + t(Sigma)) / 2
  s1 <- sqrt(Sigma[1, 1])
  s2 <- sqrt(Sigma[2, 2])
  rho <- Sigma[1, 2] / (s1 * s2)
  c(sigma1 = s1, sigma2 = s2, rho = msec_clip(rho, -0.999999, 0.999999))
}

msec_chol_lower <- function(Sigma) {
  t(chol(msec_near_pd_2x2(Sigma)))
}

msec_lambda_from_alpha <- function(alpha, Sigma) {
  A <- msec_chol_lower(Sigma)
  as.numeric(solve(t(A), as.numeric(alpha)))
}

msec_alpha_from_lambda <- function(lambda, Sigma) {
  A <- msec_chol_lower(Sigma)
  as.numeric(t(A) %*% as.numeric(lambda))
}

msec_q <- function(y, delta) {
  y <- as.matrix(y)
  delta <- as.numeric(delta)
  if (length(delta) != 2L) stop('delta must have length 2')
  out <- y
  out[, 1] <- sign(y[, 1]) * abs(y[, 1])^(delta[1] + 1)
  out[, 2] <- sign(y[, 2]) * abs(y[, 2])^(delta[2] + 1)
  out
}

msec_q_inv <- function(t, delta) {
  t <- as.matrix(t)
  delta <- as.numeric(delta)
  out <- t
  out[, 1] <- sign(t[, 1]) * abs(t[, 1])^(1 / (delta[1] + 1))
  out[, 2] <- sign(t[, 2]) * abs(t[, 2])^(1 / (delta[2] + 1))
  out
}

msec_log_jacobian <- function(y, delta, eps = 1e-12) {
  y <- as.matrix(y)
  delta <- as.numeric(delta)
  ay <- pmax(abs(y), eps)
  log(delta[1] + 1) + delta[1] * log(ay[, 1]) +
    log(delta[2] + 1) + delta[2] * log(ay[, 2])
}

msec_generator_log_g <- function(u, generator = c('gaussian', 't'), nu = 6) {
  generator <- match.arg(generator)
  if (generator == 'gaussian') return(-0.5 * u)
  -0.5 * (nu + 2) * log1p(u / nu)
}

msec_generator_log_Fc <- function(x, generator = c('gaussian', 't'), nu = 6) {
  generator <- match.arg(generator)
  if (generator == 'gaussian') return(pnorm(x, log.p = TRUE))
  pt(x, df = nu, log.p = TRUE)
}

msec_generator_Fc <- function(x, generator = c('gaussian', 't'), nu = 6) {
  generator <- match.arg(generator)
  if (generator == 'gaussian') return(pnorm(x))
  pt(x, df = nu)
}

msec_loglik_vec <- function(W, mu, Sigma, delta, alpha = c(0, 0),
                            generator = c('gaussian', 't'), nu = 6,
                            eps = 1e-12) {
  generator <- match.arg(generator)
  W <- as.matrix(W)
  if (ncol(W) != 2L) stop('W must be an n by 2 matrix')
  mu <- as.numeric(mu)
  delta <- as.numeric(delta)
  alpha <- as.numeric(alpha)
  Sigma <- msec_near_pd_2x2(Sigma)
  Sinv <- solve(Sigma)
  logdet <- as.numeric(determinant(Sigma, logarithm = TRUE)$modulus)
  lambda <- msec_lambda_from_alpha(alpha, Sigma)
  Y <- sweep(W, 2L, mu, '-')
  Tm <- msec_q(Y, delta)
  U <- rowSums((Tm %*% Sinv) * Tm)
  eta <- as.numeric(Tm %*% lambda)
  # c_g = 2*pi for both exp(-u/2) and (1+u/nu)^(-(nu+2)/2) in dimension 2.
  const <- -log(pi) - 0.5 * logdet
  out <- const + msec_generator_log_g(U, generator, nu) +
    msec_generator_log_Fc(eta, generator, nu) +
    msec_log_jacobian(Y, delta, eps = eps)
  out
}

msec_loglik <- function(W, mu, Sigma, delta, alpha = c(0, 0),
                        generator = c('gaussian', 't'), nu = 6,
                        eps = 1e-12) {
  v <- msec_loglik_vec(W, mu, Sigma, delta, alpha, generator, nu, eps)
  if (any(!is.finite(v))) return(-Inf)
  sum(v)
}

msec_r_ec <- function(n, Sigma, generator = c('gaussian', 't'), nu = 6) {
  generator <- match.arg(generator)
  A <- msec_chol_lower(Sigma)
  Z <- matrix(rnorm(2L * n), ncol = 2L)
  if (generator == 't') {
    scale <- sqrt(rchisq(n, df = nu) / nu)
    Z <- Z / scale
  }
  Z %*% t(A)
}

msec_r_se <- function(n, Sigma, alpha = c(0, 0), generator = c('gaussian', 't'),
                      nu = 6, batch_factor = 2.4, max_batches = 10000L) {
  generator <- match.arg(generator)
  lambda <- msec_lambda_from_alpha(alpha, Sigma)
  out <- matrix(NA_real_, nrow = n, ncol = 2L)
  filled <- 0L
  batches <- 0L
  while (filled < n && batches < max_batches) {
    batches <- batches + 1L
    need <- n - filled
    m <- max(64L, ceiling(batch_factor * need))
    Tprop <- msec_r_ec(m, Sigma, generator, nu)
    pacc <- msec_generator_Fc(as.numeric(Tprop %*% lambda), generator, nu)
    keep <- runif(m) <= pacc
    if (any(keep)) {
      accepted <- Tprop[keep, , drop = FALSE]
      take <- min(nrow(accepted), need)
      out[(filled + 1L):(filled + take), ] <- accepted[seq_len(take), , drop = FALSE]
      filled <- filled + take
    }
  }
  if (filled < n) stop('Accept-reject sampler failed to generate enough observations')
  colnames(out) <- c('T1', 'T2')
  out
}

msec_r_msec <- function(n, mu, Sigma, delta, alpha = c(0, 0),
                        generator = c('gaussian', 't'), nu = 6) {
  Tm <- msec_r_se(n, Sigma, alpha, generator, nu)
  W0 <- msec_q_inv(Tm, delta)
  W <- sweep(W0, 2L, as.numeric(mu), '+')
  colnames(W) <- c('W1', 'W2')
  W
}


msec_standard_params <- function(n = 300, generator = 'gaussian', nu = 6,
                                 rho = 0.35, delta = c(0.7, 0.4),
                                 alpha = c(1.2, -0.8)) {
  n <- as.integer(n)[1]
  generator <- as.character(generator)[1]
  nu <- as.numeric(nu)[1]
  rho <- as.numeric(rho)[1]
  delta <- unname(as.numeric(delta))
  alpha <- unname(as.numeric(alpha))

  mu <- c(0.3, -0.2)
  sigma <- c(1.1, 0.9)
  Sigma <- msec_make_sigma(sigma, rho)
  lambda <- msec_lambda_from_alpha(alpha, Sigma)
  list(n = n, generator = generator, nu = nu, mu = mu,
       sigma = sigma, rho = rho, Sigma = Sigma, delta = delta,
       alpha = alpha, lambda = lambda)
}


msec_interest <- function(mu, Sigma, delta, alpha) {
  comps <- msec_sigma_to_components(Sigma)
  lambda <- msec_lambda_from_alpha(alpha, Sigma)
  c(mu1 = unname(mu[1]),
    mu2 = unname(mu[2]),
    sigma1 = unname(comps['sigma1']),
    sigma2 = unname(comps['sigma2']),
    rho = unname(comps['rho']),
    delta1 = unname(delta[1]),
    delta2 = unname(delta[2]),
    alpha1 = unname(alpha[1]),
    alpha2 = unname(alpha[2]),
    lambda1 = unname(lambda[1]),
    lambda2 = unname(lambda[2]),
    alpha_norm = sqrt(sum(alpha^2)),
    alpha_angle = msec_angle_deg(alpha))
}

