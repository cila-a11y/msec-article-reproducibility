# -----------------------------------------------------------------------------
# Monte Carlo designs for the MSEC paper
# -----------------------------------------------------------------------------

msec_true_row <- function(params) {
  interest <- msec_interest(params$mu, params$Sigma, params$delta, params$alpha)
  as.list(interest)
}

msec_fit_row <- function(fit, prefix = '') {
  if (!isTRUE(fit$ok)) {
    vals <- rep(NA_real_, 13L)
    names(vals) <- c('mu1','mu2','sigma1','sigma2','rho','delta1','delta2',
                     'alpha1','alpha2','lambda1','lambda2','alpha_norm','alpha_angle')
    return(as.list(setNames(vals, paste0(prefix, names(vals)))))
  }
  vals <- fit$interest
  as.list(setNames(as.numeric(vals), paste0(prefix, names(vals))))
}

msec_sim_one_estimation <- function(row, rep_id, seed, max_starts = 4L,
                                    maxit = 600L, compute_coverage = FALSE) {
  set.seed(seed)
  params <- msec_standard_params(n = row$n, generator = row$generator, nu = row$nu,
                                 rho = row$rho,
                                 delta = c(row$delta1, row$delta2),
                                 alpha = c(row$alpha1, row$alpha2))
  W <- msec_r_msec(params$n, params$mu, params$Sigma, params$delta, params$alpha,
                   generator = params$generator, nu = params$nu)
  fit <- msec_fit(W, generator = params$generator, nu = params$nu,
                  max_starts = max_starts, maxit = maxit,
                  compute_hessian = compute_coverage, compute_grad = TRUE)
  ok_strict <- msec_fit_success(fit)
  ans <- c(list(scenario_id = row$scenario_id, study = row$study, rep = rep_id,
                generator = params$generator, n = params$n, nu = params$nu,
                ok = isTRUE(fit$ok), convergence = ifelse(is.null(fit$convergence), NA, fit$convergence),
                strict_success = ok_strict, logLik = fit$logLik,
                grad_norm = fit$grad_norm, alpha_grad_norm = fit$alpha_grad_norm),
           msec_fit_row(fit, 'hat_'),
           as.list(setNames(as.numeric(msec_interest(params$mu, params$Sigma, params$delta, params$alpha)),
                            paste0('true_', names(msec_interest(params$mu, params$Sigma, params$delta, params$alpha))))))

  if (isTRUE(compute_coverage) && isTRUE(fit$ok)) {
    obs <- msec_wald_intervals_interest(W, fit, sandwich = FALSE)
    sw <- msec_wald_intervals_interest(W, fit, sandwich = TRUE)
    for (par in c('rho','delta1','delta2','alpha1','alpha2')) {
      tv <- ans[[paste0('true_', par)]]
      if (!is.null(obs) && par %in% obs$parameter) {
        oo <- obs[obs$parameter == par, ]
        ans[[paste0('cover_obs_', par)]] <- is.finite(oo$lower) && oo$lower <= tv && tv <= oo$upper
      } else ans[[paste0('cover_obs_', par)]] <- NA
      if (!is.null(sw) && par %in% sw$parameter) {
        ss <- sw[sw$parameter == par, ]
        ans[[paste0('cover_sand_', par)]] <- is.finite(ss$lower) && ss$lower <= tv && tv <= ss$upper
      } else ans[[paste0('cover_sand_', par)]] <- NA
    }
  }
  as.data.frame(ans, stringsAsFactors = FALSE)
}

msec_sim_one_lr_lambda <- function(row, rep_id, seed, max_starts = 4L, maxit = 600L) {
  set.seed(seed)
  params <- msec_standard_params(n = row$n, generator = row$generator, nu = row$nu,
                                 rho = row$rho,
                                 delta = c(row$delta1, row$delta2),
                                 alpha = c(row$alpha1, row$alpha2))
  W <- msec_r_msec(params$n, params$mu, params$Sigma, params$delta, params$alpha,
                   generator = params$generator, nu = params$nu)
  lr <- msec_lr_lambda(W, generator = params$generator, nu = params$nu,
                       max_starts = max_starts, maxit = maxit)
  cbind(data.frame(scenario_id = row$scenario_id, study = row$study, rep = rep_id,
                   generator = params$generator, n = params$n, nu = params$nu,
                   LR_lambda = lr$LR, p_lambda = lr$p_chisq,
                   full_ok = isTRUE(lr$full$ok), null_ok = isTRUE(lr$null$ok),
                   full_logLik = lr$full$logLik, null_logLik = lr$null$logLik,
                   full_strict_success = msec_fit_success(lr$full),
                   stringsAsFactors = FALSE),
        as.data.frame(msec_fit_row(lr$full, 'hat_'), stringsAsFactors = FALSE),
        as.data.frame(as.list(setNames(as.numeric(msec_interest(params$mu, params$Sigma, params$delta, params$alpha)),
                                       paste0('true_', names(msec_interest(params$mu, params$Sigma, params$delta, params$alpha))))),
                      stringsAsFactors = FALSE))
}

msec_sim_one_lr_delta <- function(row, rep_id, seed, max_starts = 3L, maxit = 500L) {
  set.seed(seed)
  params <- msec_standard_params(n = row$n, generator = row$generator, nu = row$nu,
                                 rho = row$rho,
                                 delta = c(row$delta1, row$delta2),
                                 alpha = c(row$alpha1, row$alpha2))
  W <- msec_r_msec(params$n, params$mu, params$Sigma, params$delta, params$alpha,
                   generator = params$generator, nu = params$nu)
  boot <- msec_boot_lr_delta(W, generator = params$generator, nu = params$nu,
                             B_boot = row$B_boot, seed = seed + 500000L,
                             max_starts = max_starts, maxit = maxit)
  lr <- boot$obs
  cbind(data.frame(scenario_id = row$scenario_id, study = row$study, rep = rep_id,
                   generator = params$generator, n = params$n, nu = params$nu,
                   LR_delta = lr$LR, p_delta_chisq = lr$p_chisq,
                   p_delta_boot = boot$p_boot,
                   boot_ok = sum(is.finite(boot$boot_stats)),
                   full_ok = isTRUE(lr$full$ok), null_ok = isTRUE(lr$null$ok),
                   full_logLik = lr$full$logLik, null_logLik = lr$null$logLik,
                   stringsAsFactors = FALSE),
        as.data.frame(msec_fit_row(lr$full, 'hat_'), stringsAsFactors = FALSE),
        as.data.frame(as.list(setNames(as.numeric(msec_interest(params$mu, params$Sigma, params$delta, params$alpha)),
                                       paste0('true_', names(msec_interest(params$mu, params$Sigma, params$delta, params$alpha))))),
                      stringsAsFactors = FALSE))
}

msec_sim_one_ablation_t <- function(row, rep_id, seed, max_starts = 4L, maxit = 600L) {
  set.seed(seed)
  params <- msec_standard_params(n = row$n, generator = 't', nu = row$nu,
                                 rho = row$rho,
                                 delta = c(row$delta1, row$delta2),
                                 alpha = c(row$alpha1, row$alpha2))
  W <- msec_r_msec(params$n, params$mu, params$Sigma, params$delta, params$alpha,
                   generator = 't', nu = params$nu)
  fit <- msec_fit(W, generator = 't', nu = params$nu,
                  max_starts = max_starts, maxit = maxit)
  corr <- wrong <- c(sigma1 = NA_real_, sigma2 = NA_real_, rho = NA_real_)
  if (isTRUE(fit$ok)) {
    corr <- msec_t_scale_update_once(W, fit, denominator = 'n')
    wrong <- msec_t_scale_update_once(W, fit, denominator = 'sumw')
  }
  data.frame(scenario_id = row$scenario_id, study = row$study, rep = rep_id,
             generator = 't', n = params$n, nu = params$nu,
             ok = isTRUE(fit$ok), logLik = fit$logLik,
             sigma1_correct = corr['sigma1'], sigma2_correct = corr['sigma2'], rho_correct = corr['rho'],
             sigma1_wrong = wrong['sigma1'], sigma2_wrong = wrong['sigma2'], rho_wrong = wrong['rho'],
             delta_logLik = NA_real_,
             stringsAsFactors = FALSE)
}

msec_run_scenario_rep <- function(row, rep_id, seed, max_starts = 4L, maxit = 600L) {
  study <- as.character(row$study)
  if (study %in% c('baseline_accuracy', 'robust_rmse', 'sample_size', 'grid')) {
    return(msec_sim_one_estimation(row, rep_id, seed, max_starts, maxit, compute_coverage = FALSE))
  }
  if (study == 'coverage') {
    return(msec_sim_one_estimation(row, rep_id, seed, max_starts, maxit, compute_coverage = TRUE))
  }
  if (study == 'lr_lambda') {
    return(msec_sim_one_lr_lambda(row, rep_id, seed, max_starts, maxit))
  }
  if (study == 'lr_delta_boot') {
    return(msec_sim_one_lr_delta(row, rep_id, seed, max_starts = min(max_starts, 3L), maxit = maxit))
  }
  if (study == 'ablation_t') {
    return(msec_sim_one_ablation_t(row, rep_id, seed, max_starts, maxit))
  }
  stop('Unknown study: ', study)
}
