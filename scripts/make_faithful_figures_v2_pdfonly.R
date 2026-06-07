#!/usr/bin/env Rscript

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

source("R/00_utils.R")
msec_source_all(getwd())

out_dir <- "results/faithful_v2"
fig_pdf <- "figures/faithful/pdf"

dir.create(fig_pdf, recursive = TRUE, showWarnings = FALSE)

fit_file <- file.path(out_dir, "faithful_application_fits.rds")
if (!file.exists(fit_file)) {
  stop("missing fitted application object: ", fit_file)
}

obj <- readRDS(fit_file)

W <- obj$W
fit_g <- obj$fit_gaussian

save_pdf <- function(name, plot_fun, width = 5.8, height = 5.2) {
  pdf_file <- file.path(fig_pdf, paste0(name, ".pdf"))
  grDevices::pdf(pdf_file, width = width, height = height, useDingbats = FALSE)
  par(mar = c(4.2, 4.2, 1.0, 1.0))
  plot_fun()
  grDevices::dev.off()
  invisible(pdf_file)
}

kde2d_base <- function(x, y, n = 150, lims = NULL, h = NULL) {
  x <- as.numeric(x)
  y <- as.numeric(y)

  if (is.null(lims)) {
    lims <- c(range(x), range(y))
  }

  gx <- seq(lims[1], lims[2], length.out = n)
  gy <- seq(lims[3], lims[4], length.out = n)

  if (is.null(h)) {
    hx <- 1.06 * stats::sd(x) * length(x)^(-1/5)
    hy <- 1.06 * stats::sd(y) * length(y)^(-1/5)
    h <- c(hx, hy)
  }

  zx <- outer(gx, x, function(a, b) stats::dnorm((a - b) / h[1]) / h[1])
  zy <- outer(gy, y, function(a, b) stats::dnorm((a - b) / h[2]) / h[2])

  z <- zx %*% t(zy) / length(x)

  list(x = gx, y = gy, z = z)
}

msec_log_density_points <- function(W, mu, Sigma, delta, alpha,
                                    generator = "gaussian", nu = 6,
                                    eps = 1e-12) {
  W <- as.matrix(W)
  Y <- sweep(W, 2, mu, "-")
  Tm <- msec_q(Y, delta)

  Sinv <- solve(Sigma)
  u <- rowSums((Tm %*% Sinv) * Tm)

  lambda <- msec_lambda_from_alpha(alpha, Sigma)
  eta <- as.vector(Tm %*% lambda)

  logdet <- as.numeric(determinant(Sigma, logarithm = TRUE)$modulus)

  if (generator == "gaussian") {
    logg <- -0.5 * u
    logFc <- stats::pnorm(eta, log.p = TRUE)
  } else if (generator == "t") {
    logg <- -0.5 * (nu + 2) * log1p(u / nu)
    logFc <- stats::pt(eta, df = nu, log.p = TRUE)
  } else {
    stop("unknown generator")
  }

  logJ <- rowSums(
    matrix(log(delta + 1), nrow = nrow(Y), ncol = 2, byrow = TRUE) +
      sweep(log(pmax(abs(Y), eps)), 2, delta, "*")
  )

  -log(pi) - 0.5 * logdet + logg + logFc + logJ
}

plot_faithful_empirical_density <- function() {
  kd <- kde2d_base(
    W[, 1],
    W[, 2],
    n = 150,
    lims = c(
      min(W[, 1]) - 0.4, max(W[, 1]) + 0.4,
      min(W[, 2]) - 0.4, max(W[, 2]) + 0.4
    )
  )

  plot(
    W[, 1], W[, 2],
    xlab = "eruption duration (robustly standardized)",
    ylab = "waiting time (robustly standardized)",
    pch = 16,
    cex = 0.45
  )
  contour(kd$x, kd$y, kd$z, add = TRUE, drawlabels = FALSE, nlevels = 8)
}

plot_faithful_fitted_density <- function() {
  xg <- seq(min(W[, 1]) - 0.4, max(W[, 1]) + 0.4, length.out = 150)
  yg <- seq(min(W[, 2]) - 0.4, max(W[, 2]) + 0.4, length.out = 150)
  grd <- expand.grid(xg, yg)

  ld <- msec_log_density_points(
    grd,
    mu = fit_g$mu,
    Sigma = fit_g$Sigma,
    delta = fit_g$delta,
    alpha = fit_g$alpha,
    generator = "gaussian",
    nu = 6
  )

  zz <- matrix(exp(ld), nrow = length(xg), ncol = length(yg))

  plot(
    W[, 1], W[, 2],
    xlab = "eruption duration (robustly standardized)",
    ylab = "waiting time (robustly standardized)",
    pch = 16,
    cex = 0.45
  )
  contour(xg, yg, zz, add = TRUE, drawlabels = FALSE, nlevels = 8)
}

plot_faithful_whitened <- function() {
  Y <- sweep(W, 2, fit_g$mu, "-")
  Tm <- msec_q(Y, fit_g$delta)

  A <- t(chol(fit_g$Sigma))
  Z <- t(solve(A, t(Tm)))

  plot(
    Z[, 1], Z[, 2],
    xlab = "whitened coordinate 1",
    ylab = "whitened coordinate 2",
    pch = 16,
    cex = 0.45
  )

  abline(h = 0, v = 0, lty = 2)

  a <- fit_g$alpha
  nr <- sqrt(sum(a^2))

  if (is.finite(nr) && nr > 0) {
    lim <- max(abs(Z), na.rm = TRUE)
    arrows(
      0, 0,
      0.65 * lim * a[1] / nr,
      0.65 * lim * a[2] / nr,
      length = 0.10
    )
    text(
      0.70 * lim * a[1] / nr,
      0.70 * lim * a[2] / nr,
      labels = "skewness direction",
      pos = 4,
      cex = 0.8
    )
  }
}

plot_faithful_radial_qq <- function() {
  Y <- sweep(W, 2, fit_g$mu, "-")
  Tm <- msec_q(Y, fit_g$delta)

  A <- t(chol(fit_g$Sigma))
  Z <- t(solve(A, t(Tm)))

  r2 <- rowSums(Z^2)

  theo <- stats::qchisq(stats::ppoints(length(r2)), df = 2)
  samp <- sort(r2)

  plot(
    theo, samp,
    xlab = "theoretical quantiles",
    ylab = "sample radial quantiles",
    pch = 16,
    cex = 0.55
  )
  abline(0, 1, lty = 2)
}

save_pdf("faithful_empirical_density", plot_faithful_empirical_density)
save_pdf("faithful_fitted_density", plot_faithful_fitted_density)
save_pdf("faithful_whitened_coordinates", plot_faithful_whitened)
save_pdf("faithful_radial_qq", plot_faithful_radial_qq)

cat("faithful pdf figures written to:\n")
cat("  ", fig_pdf, "\n")
