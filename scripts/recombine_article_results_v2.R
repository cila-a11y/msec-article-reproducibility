args <- commandArgs(trailingOnly = TRUE)

infile <- if (length(args) >= 1) args[1] else "results/article_summary_raw_v2/all_replications.csv"
outdir <- if (length(args) >= 2) args[2] else "results/article_summary_v2"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

d <- read.csv(infile, stringsAsFactors = FALSE)

write.csv(d, file.path(outdir, "all_replications_article.csv"), row.names = FALSE)

est_studies <- c("baseline_accuracy", "robust_rmse", "sample_size", "coverage")

param_names <- c(
  "mu1", "mu2",
  "sigma1", "sigma2", "rho",
  "delta1", "delta2",
  "alpha1", "alpha2",
  "lambda1", "lambda2",
  "alpha_norm", "alpha_angle"
)

summarise_estimation <- function(d) {
  d <- d[d$study %in% est_studies, , drop = FALSE]
  out <- list()

  for (sc in unique(d$scenario_id)) {
    ds <- d[d$scenario_id == sc, , drop = FALSE]

    base_ok <- rep(TRUE, nrow(ds))
    if ("ok" %in% names(ds)) base_ok <- base_ok & (ds$ok %in% TRUE)
    if ("logLik" %in% names(ds)) base_ok <- base_ok & is.finite(ds$logLik)

    for (p in param_names) {
      hp <- paste0("hat_", p)
      tp <- paste0("true_", p)
      if (!all(c(hp, tp) %in% names(ds))) next

      ok <- base_ok & is.finite(ds[[hp]]) & is.finite(ds[[tp]])
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
        ok_rate = mean(base_ok, na.rm = TRUE),
        conv0_rate = if ("convergence" %in% names(ds)) mean(ds$convergence == 0, na.rm = TRUE) else NA_real_,
        strict_rate = if ("strict_success" %in% names(ds)) mean(ds$strict_success %in% TRUE, na.rm = TRUE) else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }

  if (!length(out)) return(data.frame())
  do.call(rbind, out)
}

summarise_fit_quality <- function(d) {
  d <- d[d$study %in% est_studies, , drop = FALSE]
  out <- list()

  for (sc in unique(d$scenario_id)) {
    ds <- d[d$scenario_id == sc, , drop = FALSE]

    alpha_boundary <- rep(NA, nrow(ds))
    if (all(c("hat_alpha1", "hat_alpha2") %in% names(ds))) {
      alpha_boundary <- abs(ds$hat_alpha1) > 19.5 | abs(ds$hat_alpha2) > 19.5
    }

    out[[length(out) + 1L]] <- data.frame(
      scenario_id = sc,
      study = unique(ds$study)[1],
      generator = unique(ds$generator)[1],
      n = unique(ds$n)[1],
      n_total = nrow(ds),
      ok_rate = if ("ok" %in% names(ds)) mean(ds$ok %in% TRUE, na.rm = TRUE) else NA_real_,
      conv0_rate = if ("convergence" %in% names(ds)) mean(ds$convergence == 0, na.rm = TRUE) else NA_real_,
      strict_rate = if ("strict_success" %in% names(ds)) mean(ds$strict_success %in% TRUE, na.rm = TRUE) else NA_real_,
      alpha_boundary_rate = mean(alpha_boundary, na.rm = TRUE),
      median_grad_norm = if ("grad_norm" %in% names(ds)) median(ds$grad_norm, na.rm = TRUE) else NA_real_,
      median_alpha_grad_norm = if ("alpha_grad_norm" %in% names(ds)) median(ds$alpha_grad_norm, na.rm = TRUE) else NA_real_,
      stringsAsFactors = FALSE
    )
  }

  if (!length(out)) return(data.frame())
  do.call(rbind, out)
}

summarise_wald_coverage <- function(d) {
  pars <- c("rho", "delta1", "delta2", "alpha1", "alpha2")
  d <- d[d$study == "coverage", , drop = FALSE]
  if (!nrow(d)) return(data.frame())

  out <- list()

  for (sc in unique(d$scenario_id)) {
    ds <- d[d$scenario_id == sc, , drop = FALSE]

    for (p in pars) {
      co <- paste0("cover_obs_", p)
      cs <- paste0("cover_sand_", p)
      if (!all(c(co, cs) %in% names(ds))) next

      out[[length(out) + 1L]] <- data.frame(
        scenario_id = sc,
        generator = unique(ds$generator)[1],
        n = unique(ds$n)[1],
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

write.csv(summarise_estimation(d),
          file.path(outdir, "summary_estimation_article.csv"),
          row.names = FALSE)

write.csv(summarise_fit_quality(d),
          file.path(outdir, "summary_fit_quality_article.csv"),
          row.names = FALSE)

write.csv(summarise_wald_coverage(d),
          file.path(outdir, "summary_wald_coverage_diagnostic.csv"),
          row.names = FALSE)

cat("wrote article summaries to", outdir, "\n")
