#!/usr/bin/env Rscript

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

est <- read.csv("results/article_summary_v2/summary_estimation_article.csv", stringsAsFactors = FALSE)
lr_lambda <- read.csv("results/article_summary_v2/summary_lr_lambda_article.csv", stringsAsFactors = FALSE)
lr_delta <- read.csv("results/article_summary_v2/summary_lr_delta_article.csv", stringsAsFactors = FALSE)
coverage <- read.csv("results/article_summary_v2/summary_wald_coverage_diagnostic.csv", stringsAsFactors = FALSE)

fig_pdf <- "figures/pdf"
dir.create(fig_pdf, recursive = TRUE, showWarnings = FALSE)

save_pdf <- function(name, plot_fun, width = 6.5, height = 4.5) {
  grDevices::pdf(file.path(fig_pdf, paste0(name, ".pdf")),
                 width = width, height = height, useDingbats = FALSE)
  par(mar = c(4.2, 4.2, 1.0, 1.0))
  plot_fun()
  grDevices::dev.off()
}

get_col <- function(d, candidates) {
  for (cc in candidates) {
    if (cc %in% names(d)) return(d[[cc]])
  }
  stop("none of these columns exist: ", paste(candidates, collapse = ", "))
}

plot_rmse_sample_size <- function() {
  pars <- c("rho", "delta1", "delta2", "alpha1", "alpha2")
  ss <- est[est$study == "sample_size" & est$parameter %in% pars, ]

  plot(NA, NA,
       xlim = range(ss$n, na.rm = TRUE),
       ylim = range(ss$rmse, na.rm = TRUE),
       xlab = "sample size",
       ylab = "rmse",
       xaxt = "n")

  axis(1, at = sort(unique(ss$n)))

  pch_vals <- c(16, 17, 15, 1, 2)
  lty_vals <- c(1, 2, 3, 4, 5)

  for (j in seq_along(pars)) {
    z <- ss[ss$parameter == pars[j], ]
    z <- z[order(z$n), ]
    lines(z$n, z$rmse, lty = lty_vals[j])
    points(z$n, z$rmse, pch = pch_vals[j])
  }

  legend("topright", legend = pars, lty = lty_vals, pch = pch_vals, bty = "n", cex = 0.85)
}

plot_bias_sample_size <- function() {
  pars <- c("rho", "delta1", "delta2", "alpha1", "alpha2")
  ss <- est[est$study == "sample_size" & est$parameter %in% pars, ]

  plot(NA, NA,
       xlim = range(ss$n, na.rm = TRUE),
       ylim = range(ss$bias, na.rm = TRUE),
       xlab = "sample size",
       ylab = "bias",
       xaxt = "n")

  axis(1, at = sort(unique(ss$n)))
  abline(h = 0, lty = 2)

  pch_vals <- c(16, 17, 15, 1, 2)
  lty_vals <- c(1, 2, 3, 4, 5)

  for (j in seq_along(pars)) {
    z <- ss[ss$parameter == pars[j], ]
    z <- z[order(z$n), ]
    lines(z$n, z$bias, lty = lty_vals[j])
    points(z$n, z$bias, pch = pch_vals[j])
  }

  legend("bottomright", legend = pars, lty = lty_vals, pch = pch_vals, bty = "n", cex = 0.85)
}

plot_baseline_rmse <- function() {
  pars <- c("mu1", "mu2", "sigma1", "sigma2", "rho", "delta1", "delta2", "alpha1", "alpha2")
  z <- est[est$scenario_id == "baseline_gaussian_n300" & est$parameter %in% pars, ]
  z$parameter <- factor(z$parameter, levels = pars)
  z <- z[order(z$parameter), ]

  barplot(z$rmse,
          names.arg = z$parameter,
          las = 2,
          ylab = "rmse",
          xlab = "parameter")
}

plot_lr_lambda <- function() {
  vals <- get_col(lr_lambda, c("reject_5.", "reject_5%", "reject_5pct", "reject_5"))
  labs <- c("size gaussian", "power gaussian", "size student-t", "power student-t")

  barplot(vals,
          names.arg = labs,
          las = 2,
          ylim = c(0, 1),
          ylab = "rejection rate at 5%",
          xlab = "")
  abline(h = 0.05, lty = 2)
}

plot_lr_delta <- function() {
  chisq <- get_col(lr_delta, c("reject_chisq_5.", "reject_chisq_5%", "reject_chisq_5pct", "reject_chisq_5"))
  boot <- get_col(lr_delta, c("reject_boot_5.", "reject_boot_5%", "reject_boot_5pct", "reject_boot_5"))

  vals <- rbind(chisq, boot)
  colnames(vals) <- c("size gaussian", "power gaussian")

  barplot(vals,
          beside = TRUE,
          ylim = c(0, 1),
          ylab = "rejection rate at 5%",
          xlab = "",
          legend.text = c("chi-square", "bootstrap"),
          args.legend = list(bty = "n", x = "topleft"))
  abline(h = 0.05, lty = 2)
}

plot_wald_coverage <- function() {
  z <- coverage[coverage$scenario_id %in% c("coverage_gaussian_n300", "coverage_t6_n300"), ]

  z$label <- paste(
    z$parameter,
    ifelse(z$generator == "t", "student-t", "gaussian"),
    sep = ", "
  )

  vals <- rbind(z$observed, z$sandwich)
  colnames(vals) <- z$label

  barplot(vals,
          beside = TRUE,
          las = 2,
          ylim = c(0, 1),
          ylab = "coverage",
          xlab = "",
          legend.text = c("observed", "sandwich"),
          args.legend = list(bty = "n", x = "bottomleft", cex = 0.8))
  abline(h = 0.95, lty = 2)
}

save_pdf("mc_rmse_by_sample_size", plot_rmse_sample_size)
save_pdf("mc_bias_by_sample_size", plot_bias_sample_size)
save_pdf("mc_baseline_rmse", plot_baseline_rmse, width = 7.0, height = 4.8)
save_pdf("lr_skewness_rejection_rate", plot_lr_lambda, width = 7.0, height = 4.8)
save_pdf("lr_warping_rejection_rate", plot_lr_delta, width = 6.5, height = 4.5)
save_pdf("wald_coverage_diagnostic", plot_wald_coverage, width = 8.0, height = 5.0)

cat("monte carlo pdf figures written to:\n")
cat("  ", fig_pdf, "\n")
