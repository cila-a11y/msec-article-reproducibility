#!/usr/bin/env Rscript

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

dir.create("figures/monte_carlo", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/pdf", recursive = TRUE, showWarnings = FALSE)
dir.create("results/explanatory_figures", recursive = TRUE, showWarnings = FALSE)

make_sigma <- function(rho, sigma = c(1, 1)) {
  matrix(
    c(
      sigma[1]^2,
      rho * sigma[1] * sigma[2],
      rho * sigma[1] * sigma[2],
      sigma[2]^2
    ),
    nrow = 2,
    byrow = TRUE
  )
}

qwarp <- function(Y, delta) {
  Y <- as.matrix(Y)
  out <- Y
  out[, 1] <- sign(Y[, 1]) * abs(Y[, 1])^(delta[1] + 1)
  out[, 2] <- sign(Y[, 2]) * abs(Y[, 2])^(delta[2] + 1)
  out
}

lambda_from_alpha <- function(alpha, Sigma) {
  A <- t(chol(Sigma))
  as.vector(solve(t(A), alpha))
}

msec_logdens_gaussian <- function(W, rho, delta, alpha, eps = 1e-12) {
  W <- as.matrix(W)
  Sigma <- make_sigma(rho)
  Tm <- qwarp(W, delta)
  Sinv <- solve(Sigma)
  u <- rowSums((Tm %*% Sinv) * Tm)
  lambda <- lambda_from_alpha(alpha, Sigma)
  eta <- as.vector(Tm %*% lambda)
  logdet <- as.numeric(determinant(Sigma, logarithm = TRUE)$modulus)

  logJ <- rowSums(
    matrix(log(delta + 1), nrow = nrow(W), ncol = 2, byrow = TRUE) +
      sweep(log(pmax(abs(W), eps)), 2, delta, "*")
  )

  -log(pi) - 0.5 * logdet - 0.5 * u + pnorm(eta, log.p = TRUE) + logJ
}

density_grid <- function(rho, delta, alpha, lim = 2.6, n = 180) {
  x <- seq(-lim, lim, length.out = n)
  y <- seq(-lim, lim, length.out = n)
  grid <- expand.grid(x, y)
  z <- exp(msec_logdens_gaussian(grid, rho = rho, delta = delta, alpha = alpha))
  list(x = x, y = y, z = matrix(z, nrow = n, ncol = n))
}

draw_panel <- function(g, label) {
  plot(
    NA, NA,
    xlim = range(g$x),
    ylim = range(g$y),
    xlab = "coordinate 1",
    ylab = "coordinate 2",
    asp = 1
  )
  contour(g$x, g$y, g$z, add = TRUE, drawlabels = FALSE, nlevels = 9)
  usr <- par("usr")
  text(
    usr[1] + 0.05 * diff(usr[1:2]),
    usr[4] - 0.06 * diff(usr[3:4]),
    labels = label,
    adj = c(0, 1),
    cex = 0.85
  )
}

save_pdf <- function(file, fun, width, height) {
  pdf(file, width = width, height = height, useDingbats = FALSE)
  on.exit(dev.off(), add = TRUE)
  fun()
}

## ------------------------------------------------------------
## first addition: delta sensitivity
## ------------------------------------------------------------

delta_cases <- data.frame(
  panel = c("a", "b", "c", "d"),
  delta1 = c(0.00, 0.30, 0.70, 1.10),
  delta2 = c(0.00, 0.30, 0.70, 1.10),
  rho = 0.70,
  alpha1 = 0,
  alpha2 = 0,
  description = c(
    "unwarped gaussian core",
    "weak symmetric warping",
    "moderate symmetric warping",
    "strong symmetric warping"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  delta_cases,
  "results/explanatory_figures/delta_sensitivity_parameters.csv",
  row.names = FALSE
)

save_pdf(
  "figures/monte_carlo/delta_sensitivity_density.pdf",
  function() {
    old <- par(no.readonly = TRUE)
    on.exit(par(old), add = TRUE)
    par(mfrow = c(2, 2), mar = c(4.0, 4.0, 1.0, 0.8))

    for (i in seq_len(nrow(delta_cases))) {
      delta <- c(delta_cases$delta1[i], delta_cases$delta2[i])
      g <- density_grid(rho = 0.70, delta = delta, alpha = c(0, 0), lim = 2.4)
      draw_panel(
        g,
        paste0("delta = (", sprintf("%.1f", delta[1]), ", ", sprintf("%.1f", delta[2]), ")")
      )
    }
  },
  width = 7.2,
  height = 7.0
)

file.copy(
  "figures/monte_carlo/delta_sensitivity_density.pdf",
  "figures/pdf/delta_sensitivity_density.pdf",
  overwrite = TRUE
)

## ------------------------------------------------------------
## second addition: weak, moderate, strong scenarios
## ------------------------------------------------------------

shape_cases <- data.frame(
  scenario = c("weak", "moderate baseline", "strong"),
  rho = c(0.35, 0.35, 0.70),
  delta1 = c(0.30, 0.70, 1.00),
  delta2 = c(0.30, 0.40, 0.70),
  alpha1 = c(0.60, 1.20, 1.20),
  alpha2 = c(-0.40, -0.80, -0.80),
  description = c(
    "shallow troughs and weak modal separation",
    "baseline geometry used in the main monte carlo design",
    "pronounced troughs and stronger diagonal separation"
  ),
  stringsAsFactors = FALSE
)

shape_cases$alpha_norm <- sqrt(shape_cases$alpha1^2 + shape_cases$alpha2^2)
shape_cases$alpha_angle <- atan2(shape_cases$alpha2, shape_cases$alpha1) * 180 / pi

write.csv(
  shape_cases,
  "results/explanatory_figures/mc_shape_scenarios.csv",
  row.names = FALSE
)

save_pdf(
  "figures/monte_carlo/mc_multimodality_scenarios.pdf",
  function() {
    old <- par(no.readonly = TRUE)
    on.exit(par(old), add = TRUE)
    par(mfrow = c(1, 3), mar = c(4.0, 4.0, 1.0, 0.8))

    for (i in seq_len(nrow(shape_cases))) {
      delta <- c(shape_cases$delta1[i], shape_cases$delta2[i])
      alpha <- c(shape_cases$alpha1[i], shape_cases$alpha2[i])
      g <- density_grid(rho = shape_cases$rho[i], delta = delta, alpha = alpha, lim = 2.5)
      draw_panel(g, shape_cases$scenario[i])
    }
  },
  width = 10.0,
  height = 3.8
)

file.copy(
  "figures/monte_carlo/mc_multimodality_scenarios.pdf",
  "figures/pdf/mc_multimodality_scenarios.pdf",
  overwrite = TRUE
)

cat("written figures:\n")
cat("  figures/monte_carlo/delta_sensitivity_density.pdf\n")
cat("  figures/monte_carlo/mc_multimodality_scenarios.pdf\n")
cat("written numerical files:\n")
cat("  results/explanatory_figures/delta_sensitivity_parameters.csv\n")
cat("  results/explanatory_figures/mc_shape_scenarios.csv\n")
