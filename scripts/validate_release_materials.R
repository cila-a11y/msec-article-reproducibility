required_files <- c(
  "results/article_summary_v2/summary_estimation_article.csv",
  "results/article_summary_v2/summary_fit_quality_article.csv",
  "results/article_summary_v2/summary_lr_lambda_article.csv",
  "results/article_summary_v2/summary_lr_delta_article.csv",
  "results/article_summary_v2/summary_wald_coverage_diagnostic.csv",

  "results/faithful_v2/faithful_data_summary.csv",
  "results/faithful_v2/faithful_data_extra.csv",
  "results/faithful_v2/faithful_model_comparison.csv",
  "results/faithful_v2/faithful_parameter_estimates.csv",
  "results/faithful_v2/faithful_lr_delta_bootstrap.csv",

  "figures/pdf/mc_rmse_by_sample_size.pdf",
  "figures/pdf/mc_bias_by_sample_size.pdf",
  "figures/pdf/mc_baseline_rmse.pdf",
  "figures/pdf/lr_skewness_rejection_rate.pdf",
  "figures/pdf/lr_warping_rejection_rate.pdf",
  "figures/pdf/wald_coverage_diagnostic.pdf",

  "figures/faithful/pdf/faithful_empirical_density.pdf",
  "figures/faithful/pdf/faithful_fitted_density.pdf",
  "figures/faithful/pdf/faithful_whitened_coordinates.pdf",
  "figures/faithful/pdf/faithful_radial_qq.pdf"
)

missing <- required_files[!file.exists(required_files)]

if (length(missing)) {
  cat("missing required files:\n")
  print(missing)
  stop("release materials are incomplete")
}

est <- read.csv("results/article_summary_v2/summary_estimation_article.csv", stringsAsFactors = FALSE)
bad_n <- est[est$n_used / est$n_total < 0.95, ]

if (nrow(bad_n)) {
  print(unique(bad_n[, c("scenario_id", "study", "n_used", "n_total")]))
  stop("some estimation summaries use too few replications")
}

fitq <- read.csv("results/article_summary_v2/summary_fit_quality_article.csv", stringsAsFactors = FALSE)
bad_alpha <- fitq[is.finite(fitq$alpha_boundary_rate) & fitq$alpha_boundary_rate > 0.02, ]

if (nrow(bad_alpha)) {
  print(bad_alpha)
  stop("some scenarios have too many alpha-boundary fits")
}

lr_l <- read.csv("results/article_summary_v2/summary_lr_lambda_article.csv", stringsAsFactors = FALSE)
lr_d <- read.csv("results/article_summary_v2/summary_lr_delta_article.csv", stringsAsFactors = FALSE)

cat("required files present:", length(required_files), "\n")
cat("estimation summaries validated\n")
cat("fit quality validated\n")
cat("lr lambda rows:", nrow(lr_l), "\n")
cat("lr delta rows:", nrow(lr_d), "\n")
cat("OK: release materials are complete\n")
