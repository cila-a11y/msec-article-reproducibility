est <- read.csv("results/article_summary_v2/summary_estimation_article.csv", stringsAsFactors = FALSE)
fitq <- read.csv("results/article_summary_v2/summary_fit_quality_article.csv", stringsAsFactors = FALSE)

cat("fit quality:\n")
print(fitq)

cat("\nminimum n_used by scenario:\n")
print(aggregate(n_used ~ scenario_id + study + generator + n + n_total,
                data = est,
                FUN = min))

bad_n <- est[est$n_used / est$n_total < 0.95, ]
if (nrow(bad_n)) {
  cat("\nScenarios with n_used/n_total < 0.95:\n")
  print(unique(bad_n[, c("scenario_id", "study", "generator", "n", "n_used", "n_total")]))
  stop("Some estimation summaries still use too few replications.")
}

bad_alpha <- fitq[is.finite(fitq$alpha_boundary_rate) & fitq$alpha_boundary_rate > 0.02, ]
if (nrow(bad_alpha)) {
  cat("\nScenarios with alpha boundary rate > 0.02:\n")
  print(bad_alpha)
  stop("Some scenarios have too many alpha-boundary fits.")
}

key <- est[est$scenario_id %in% c("baseline_gaussian_n300", "robust_rmse_gaussian", "robust_rmse_t6") &
             est$parameter %in% c("mu1","mu2","sigma1","sigma2","rho","delta1","delta2","alpha1","alpha2"),
           c("scenario_id","parameter","true","mean","bias","rmse","n_used","n_total")]

cat("\nkey estimation summaries:\n")
print(key, row.names = FALSE)

base <- est[est$scenario_id == "baseline_gaussian_n300", ]

check <- base[base$parameter %in% c("mu1","mu2","sigma1","sigma2","rho","delta1","delta2","alpha1","alpha2"), ]

cat("\nbaseline absolute biases:\n")
print(check[, c("parameter","true","mean","bias","rmse")], row.names = FALSE)

if (abs(check$bias[check$parameter == "alpha1"]) > 0.25) {
  stop("alpha1 bias in baseline is still too large.")
}

if (abs(check$bias[check$parameter == "alpha2"]) > 0.20) {
  stop("alpha2 bias in baseline is still too large.")
}

if (abs(check$bias[check$parameter == "rho"]) > 0.05) {
  stop("rho bias in baseline is still too large.")
}

cat("\nOK: article v2 estimation summaries pass validation checks.\n")
