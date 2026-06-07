# -----------------------------------------------------------------------------
# Summary functions for Monte Carlo output
# -----------------------------------------------------------------------------

msec_param_names <- function() {
  c('mu1','mu2','sigma1','sigma2','rho','delta1','delta2','alpha1','alpha2','lambda1','lambda2')
}

msec_summarise_estimation <- function(d) {
  pars <- msec_param_names()
  out <- list()
  for (sc in unique(d$scenario_id)) {
    ds <- d[d$scenario_id == sc, , drop = FALSE]
    for (p in pars) {
      hp <- paste0('hat_', p)
      tp <- paste0('true_', p)
      if (!all(c(hp, tp) %in% names(ds))) next
      ok <- is.finite(ds[[hp]]) & is.finite(ds[[tp]])
      if ('strict_success' %in% names(ds)) ok <- ok & (isTRUE(ds$strict_success) | ds$strict_success %in% TRUE)
      if (!any(ok)) ok <- is.finite(ds[[hp]]) & is.finite(ds[[tp]])
      err <- ds[[hp]][ok] - ds[[tp]][ok]
      out[[length(out) + 1L]] <- data.frame(
        scenario_id = sc,
        study = unique(ds$study)[1],
        generator = unique(ds$generator)[1],
        n = unique(ds$n)[1],
        parameter = p,
        true = mean(ds[[tp]][ok], na.rm = TRUE),
        mean = mean(ds[[hp]][ok], na.rm = TRUE),
        sd = stats::sd(ds[[hp]][ok], na.rm = TRUE),
        bias = mean(err, na.rm = TRUE),
        rmse = sqrt(mean(err^2, na.rm = TRUE)),
        n_used = sum(ok),
        n_total = nrow(ds),
        success_rate = if ('strict_success' %in% names(ds)) mean(ds$strict_success %in% TRUE, na.rm = TRUE) else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(out)) return(data.frame())
  do.call(rbind, out)
}

msec_summarise_lr_lambda <- function(d) {
  d <- d[is.finite(d$LR_lambda), , drop = FALSE]
  if (!nrow(d)) return(data.frame())
  levs <- c('10%' = 0.10, '5%' = 0.05, '1%' = 0.01)
  out <- list()
  for (sc in unique(d$scenario_id)) {
    ds <- d[d$scenario_id == sc, , drop = FALSE]
    x <- data.frame(scenario_id = sc,
                    generator = unique(ds$generator)[1],
                    n = unique(ds$n)[1],
                    mean_LR = mean(ds$LR_lambda, na.rm = TRUE),
                    median_LR = stats::median(ds$LR_lambda, na.rm = TRUE),
                    n_used = nrow(ds),
                    full_ok_rate = mean(ds$full_ok %in% TRUE, na.rm = TRUE),
                    null_ok_rate = mean(ds$null_ok %in% TRUE, na.rm = TRUE),
                    stringsAsFactors = FALSE)
    for (nm in names(levs)) x[[paste0('reject_', nm)]] <- mean(ds$p_lambda < levs[[nm]], na.rm = TRUE)
    out[[length(out) + 1L]] <- x
  }
  do.call(rbind, out)
}

msec_summarise_lr_delta <- function(d) {
  d <- d[is.finite(d$LR_delta), , drop = FALSE]
  if (!nrow(d)) return(data.frame())
  levs <- c('10%' = 0.10, '5%' = 0.05, '1%' = 0.01)
  out <- list()
  for (sc in unique(d$scenario_id)) {
    ds <- d[d$scenario_id == sc, , drop = FALSE]
    x <- data.frame(scenario_id = sc,
                    generator = unique(ds$generator)[1],
                    n = unique(ds$n)[1],
                    mean_LR = mean(ds$LR_delta, na.rm = TRUE),
                    median_LR = stats::median(ds$LR_delta, na.rm = TRUE),
                    n_used = nrow(ds),
                    mean_boot_ok = mean(ds$boot_ok, na.rm = TRUE),
                    stringsAsFactors = FALSE)
    for (nm in names(levs)) {
      x[[paste0('reject_chisq_', nm)]] <- mean(ds$p_delta_chisq < levs[[nm]], na.rm = TRUE)
      x[[paste0('reject_boot_', nm)]] <- mean(ds$p_delta_boot < levs[[nm]], na.rm = TRUE)
    }
    out[[length(out) + 1L]] <- x
  }
  do.call(rbind, out)
}

msec_summarise_coverage <- function(d) {
  pars <- c('rho','delta1','delta2','alpha1','alpha2')
  out <- list()
  for (sc in unique(d$scenario_id)) {
    ds <- d[d$scenario_id == sc, , drop = FALSE]
    for (p in pars) {
      co <- paste0('cover_obs_', p)
      cs <- paste0('cover_sand_', p)
      if (!all(c(co, cs) %in% names(ds))) next
      out[[length(out) + 1L]] <- data.frame(
        scenario_id = sc, generator = unique(ds$generator)[1], n = unique(ds$n)[1],
        parameter = p,
        observed = mean(ds[[co]] %in% TRUE, na.rm = TRUE),
        sandwich = mean(ds[[cs]] %in% TRUE, na.rm = TRUE),
        n_used = sum(!is.na(ds[[co]]) | !is.na(ds[[cs]])),
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(out)) return(data.frame())
  do.call(rbind, out)
}

msec_summarise_ablation <- function(d) {
  if (!nrow(d)) return(data.frame())
  out <- list()
  for (sc in unique(d$scenario_id)) {
    ds <- d[d$scenario_id == sc, , drop = FALSE]
    out[[length(out) + 1L]] <- data.frame(
      scenario_id = sc, generator = 't', n = unique(ds$n)[1],
      sigma1_correct_mean = mean(ds$sigma1_correct, na.rm = TRUE),
      sigma1_correct_sd = stats::sd(ds$sigma1_correct, na.rm = TRUE),
      sigma2_correct_mean = mean(ds$sigma2_correct, na.rm = TRUE),
      sigma2_correct_sd = stats::sd(ds$sigma2_correct, na.rm = TRUE),
      rho_correct_mean = mean(ds$rho_correct, na.rm = TRUE),
      rho_correct_sd = stats::sd(ds$rho_correct, na.rm = TRUE),
      sigma1_wrong_mean = mean(ds$sigma1_wrong, na.rm = TRUE),
      sigma1_wrong_sd = stats::sd(ds$sigma1_wrong, na.rm = TRUE),
      sigma2_wrong_mean = mean(ds$sigma2_wrong, na.rm = TRUE),
      sigma2_wrong_sd = stats::sd(ds$sigma2_wrong, na.rm = TRUE),
      rho_wrong_mean = mean(ds$rho_wrong, na.rm = TRUE),
      rho_wrong_sd = stats::sd(ds$rho_wrong, na.rm = TRUE),
      n_used = sum(ds$ok %in% TRUE, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

msec_write_latex_simple <- function(d, file, digits = 3) {
  msec_safe_mkdir(dirname(file))
  con <- file(file, open = 'wt')
  on.exit(close(con), add = TRUE)
  cat('% Automatically generated by combine_results.R\n', file = con)
  if (!nrow(d)) {
    cat('% Empty table\n', file = con, append = TRUE)
    return(invisible(file))
  }
  dd <- d
  for (j in seq_along(dd)) if (is.numeric(dd[[j]])) dd[[j]] <- round(dd[[j]], digits)
  cat('\\begin{tabular}{', paste(rep('l', ncol(dd)), collapse = ''), '}\n', sep = '', file = con, append = TRUE)
  cat('\\toprule\n', file = con, append = TRUE)
  cat(paste(names(dd), collapse = ' & '), ' \\\\\n', sep = '', file = con, append = TRUE)
  cat('\\midrule\n', file = con, append = TRUE)
  for (i in seq_len(nrow(dd))) {
    cat(paste(as.character(dd[i, ]), collapse = ' & '), ' \\\\\n', sep = '', file = con, append = TRUE)
  }
  cat('\\bottomrule\n\\end{tabular}\n', file = con, append = TRUE)
  invisible(file)
}
