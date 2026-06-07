#!/usr/bin/env Rscript

source("R/00_utils.R")
msec_source_all(getwd())

options(bitmapType = "cairo")

out_dir <- "results/faithful_v2"
fig_pdf <- "figures/faithful/pdf"
fig_png <- "figures/faithful/png"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_pdf, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_png, recursive = TRUE, showWarnings = FALSE)

## ------------------------------------------------------------
## graphics helpers
## ------------------------------------------------------------

open_png_device <- function(filename, width, height, res = 300) {
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(filename, width = width, height = height,
                  units = "in", res = res)
    return(invisible(TRUE))
  }

  if (capabilities("cairo")) {
    grDevices::png(filename, width = width, height = height,
                   units = "in", res = res, type = "cairo")
    return(invisible(TRUE))
  }

  stop("no headless png device available; pdf files will still be produced")
}

save_figure <- function(name, plot_fun, width = 5.8, height = 5.2) {
  pdf_file <- file.path(fig_pdf, paste0(name, ".pdf"))
  png_file <- file.path(fig_png, paste0(name, ".png"))

  grDevices::pdf(pdf_file, width = width, height = height, useDingbats = FALSE)
  par(mar = c(4.2, 4.2, 1.0, 1.0))
  plot_fun()
  grDevices::dev.off()

  ok_png <- TRUE
  tryCatch({
    open_png_device(png_file, width = width, height = height, res = 300)
    par(mar = c(4.2, 4.2, 1.0, 1.0))
    plot_fun()
    grDevices::dev.off()
  }, error = function(e) {
    ok_png <<- FALSE
    message("png skipped for ", name, ": ", conditionMessage(e))
  })

  invisible(list(pdf = pdf_file, png = if (ok_png) png_file else NA_character_))
}

## ------------------------------------------------------------
## data and robust standardization
## ------------------------------------------------------------

faith <- datasets::faithful
W_raw <- as.matrix(faith[, c("eruptions", "waiting")])
colnames(W_raw) <- c("eruption duration", "waiting time")

center <- apply(W_raw, 2L, stats::median)
scalev <- apply(W_raw, 2L, stats::mad, constant = 1.4826)

W <- sweep(sweep(W_raw, 2L, center, "-"), 2L, scalev, "/")
W <- as.matrix(W)
colnames(W) <- c("eruption duration", "waiting time")

data_summary <- data.frame(
  variable = colnames(W_raw),
  mean = colMeans(W_raw),
  sd = apply(W_raw, 2L, stats::sd),
  median = apply(W_raw, 2L, stats::median),
  iqr = apply(W_raw, 2L, stats::IQR),
  mad_scaled = scalev,
  minimum = apply(W_raw, 2L, min),
  maximum = apply(W_raw, 2L, max),
  standardized_sd = apply(W, 2L, stats::sd),
  stringsAsFactors = FALSE
)

data_extra <- data.frame(
  n = nrow(W),
  robust_correlation = stats::cor(W[, 1], W[, 2]),
  stringsAsFactors = FALSE
)

write.csv(data_summary, file.path(out_dir, "faithful_data_summary.csv"), row.names = FALSE)
write.csv(data_extra, file.path(out_dir, "faithful_data_extra.csv"), row.names = FALSE)

## ------------------------------------------------------------
## log-density helper for fitted MSEC model
## ------------------------------------------------------------

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

## ------------------------------------------------------------
## model fitting
## ------------------------------------------------------------

set.seed(20260612)

fit_g <- msec_fit(
  W,
  generator = "gaussian",
  nu = 6,
  max_starts = 40,
  maxit = 2500,
  compute_grad = TRUE
)

fit_t <- msec_fit(
  W,
  generator = "t",
  nu = 6,
  max_starts = 40,
  maxit = 2500,
  compute_grad = TRUE
)

lr_g <- msec_lr_lambda(
  W,
  generator = "gaussian",
  nu = 6,
  max_starts = 30,
  maxit = 2000,
  compute_grad = TRUE
)

lr_t <- msec_lr_lambda(
  W,
  generator = "t",
  nu = 6,
  max_starts = 30,
  maxit = 2000,
  compute_grad = TRUE
)

## ------------------------------------------------------------
## bootstrap LR for warping under selected Gaussian model
## ------------------------------------------------------------

msec_boot_lr_delta_parallel <- function(W, generator = "gaussian", nu = 6,
                                        B_boot = 199L, seed = 20260613L,
                                        max_starts_obs = 40L,
                                        maxit_obs = 2500L,
                                        max_starts_boot = 20L,
                                        maxit_boot = 1500L) {
  set.seed(seed)

  obs <- msec_lr_delta(
    W,
    generator = generator,
    nu = nu,
    max_starts = max_starts_obs,
    maxit = maxit_obs,
    compute_grad = TRUE
  )

  if (!is.finite(obs$LR) || !isTRUE(obs$null$ok)) {
    return(list(obs = obs, boot_stats = rep(NA_real_, B_boot),
                p_boot = NA_real_, boot_ok = 0L))
  }

  th0 <- obs$null
  seeds <- seed + seq_len(B_boot)

  ncores <- as.integer(Sys.getenv("SLURM_NTASKS", "1"))
  ncores <- max(1L, min(ncores, B_boot))

  boot_stats <- parallel::mclapply(
    seq_len(B_boot),
    function(b) {
      set.seed(seeds[b])
      Wb <- msec_r_msec(
        nrow(W),
        th0$mu,
        th0$Sigma,
        th0$delta,
        th0$alpha,
        generator = generator,
        nu = nu
      )

      fb <- try(
        msec_lr_delta(
          Wb,
          generator = generator,
          nu = nu,
          max_starts = max_starts_boot,
          maxit = maxit_boot,
          compute_grad = FALSE
        ),
        silent = TRUE
      )

      if (inherits(fb, "try-error")) return(NA_real_)
      fb$LR
    },
    mc.cores = ncores
  )

  boot_stats <- unlist(boot_stats)
  p_boot <- msec_plus_one_pvalue(obs$LR, boot_stats)

  list(
    obs = obs,
    boot_stats = boot_stats,
    p_boot = p_boot,
    boot_ok = sum(is.finite(boot_stats))
  )
}

lr_delta_g <- msec_boot_lr_delta_parallel(
  W,
  generator = "gaussian",
  nu = 6,
  B_boot = 199,
  seed = 20260613,
  max_starts_obs = 40,
  maxit_obs = 2500,
  max_starts_boot = 20,
  maxit_boot = 1500
)

## ------------------------------------------------------------
## tables
## ------------------------------------------------------------

fit_row <- function(fit, lr, generator_label) {
  k <- length(fit$par)
  a <- fit$alpha
  alpha_norm <- sqrt(sum(a^2))
  alpha_angle <- atan2(a[2], a[1]) * 180 / pi

  data.frame(
    generator = generator_label,
    logLik = fit$logLik,
    k = k,
    AIC = -2 * fit$logLik + 2 * k,
    BIC = -2 * fit$logLik + log(nrow(W)) * k,
    LR_lambda = lr$LR,
    p_lambda = lr$p_chisq,
    alpha_norm = alpha_norm,
    alpha_angle = alpha_angle,
    ok = isTRUE(fit$ok),
    convergence = fit$convergence,
    grad_norm = fit$grad_norm,
    alpha_grad_norm = fit$alpha_grad_norm,
    stringsAsFactors = FALSE
  )
}

comparison <- rbind(
  fit_row(fit_g, lr_g, "gaussian"),
  fit_row(fit_t, lr_t, "student-t(6)")
)

param_table <- data.frame(
  parameter = c(
    "mu1", "mu2", "sigma1", "sigma2", "rho",
    "delta1", "delta2", "alpha1", "alpha2",
    "lambda1", "lambda2"
  ),
  gaussian = c(
    fit_g$mu[1], fit_g$mu[2],
    fit_g$sigma[1], fit_g$sigma[2],
    fit_g$rho,
    fit_g$delta[1], fit_g$delta[2],
    fit_g$alpha[1], fit_g$alpha[2],
    fit_g$lambda[1], fit_g$lambda[2]
  ),
  student_t6 = c(
    fit_t$mu[1], fit_t$mu[2],
    fit_t$sigma[1], fit_t$sigma[2],
    fit_t$rho,
    fit_t$delta[1], fit_t$delta[2],
    fit_t$alpha[1], fit_t$alpha[2],
    fit_t$lambda[1], fit_t$lambda[2]
  ),
  stringsAsFactors = FALSE
)

lr_delta_table <- data.frame(
  generator = "gaussian",
  LR_delta = lr_delta_g$obs$LR,
  p_chisq = lr_delta_g$obs$p_chisq,
  p_boot = lr_delta_g$p_boot,
  B_boot = length(lr_delta_g$boot_stats),
  boot_ok = lr_delta_g$boot_ok,
  full_logLik = lr_delta_g$obs$full$logLik,
  null_logLik = lr_delta_g$obs$null$logLik,
  stringsAsFactors = FALSE
)

write.csv(comparison, file.path(out_dir, "faithful_model_comparison.csv"), row.names = FALSE)
write.csv(param_table, file.path(out_dir, "faithful_parameter_estimates.csv"), row.names = FALSE)
write.csv(lr_delta_table, file.path(out_dir, "faithful_lr_delta_bootstrap.csv"), row.names = FALSE)

saveRDS(
  list(
    fit_gaussian = fit_g,
    fit_t6 = fit_t,
    lr_gaussian = lr_g,
    lr_t6 = lr_t,
    lr_delta_gaussian = lr_delta_g,
    center = center,
    scale = scalev,
    W = W,
    W_raw = W_raw,
    data_summary = data_summary,
    data_extra = data_extra
  ),
  file.path(out_dir, "faithful_application_fits.rds")
)

## LaTeX fragments
fmt <- function(x, digits = 3) formatC(x, format = "f", digits = digits)

tex_comp <- file.path(out_dir, "faithful_model_comparison.tex")
cat("\\begin{tabular}{lrrrrrrr}\n", file = tex_comp)
cat("\\toprule\n", file = tex_comp, append = TRUE)
cat("Generator & logLik & AIC & BIC & LR$_\\lambda$ & $p_\\lambda$ & $\\|\\widehat{\\bm\\alpha}\\|$ & angle \\\\\n", file = tex_comp, append = TRUE)
cat("\\midrule\n", file = tex_comp, append = TRUE)
for (i in seq_len(nrow(comparison))) {
  cat(
    comparison$generator[i], " & ",
    fmt(comparison$logLik[i], 3), " & ",
    fmt(comparison$AIC[i], 3), " & ",
    fmt(comparison$BIC[i], 3), " & ",
    fmt(comparison$LR_lambda[i], 3), " & ",
    ifelse(comparison$p_lambda[i] < 0.001, "$<0.001$", fmt(comparison$p_lambda[i], 3)), " & ",
    fmt(comparison$alpha_norm[i], 3), " & ",
    fmt(comparison$alpha_angle[i], 2), "$^\\circ$ \\\\\n",
    sep = "",
    file = tex_comp,
    append = TRUE
  )
}
cat("\\bottomrule\n\\end{tabular}\n", file = tex_comp, append = TRUE)

tex_par <- file.path(out_dir, "faithful_parameter_estimates.tex")
cat("\\begin{tabular}{lrr}\n", file = tex_par)
cat("\\toprule\n", file = tex_par, append = TRUE)
cat("Parameter & Gaussian & Student-$t(6)$ \\\\\n", file = tex_par, append = TRUE)
cat("\\midrule\n", file = tex_par, append = TRUE)
for (i in seq_len(nrow(param_table))) {
  cat(
    param_table$parameter[i], " & ",
    fmt(param_table$gaussian[i], 4), " & ",
    fmt(param_table$student_t6[i], 4), " \\\\\n",
    sep = "",
    file = tex_par,
    append = TRUE
  )
}
cat("\\bottomrule\n\\end{tabular}\n", file = tex_par, append = TRUE)

tex_delta <- file.path(out_dir, "faithful_lr_delta_bootstrap.tex")
cat("\\begin{tabular}{lrrrr}\n", file = tex_delta)
cat("\\toprule\n", file = tex_delta, append = TRUE)
cat("Generator & LR$_\\delta$ & bootstrap $p$ & $B_{\\rm boot}$ & valid bootstraps \\\\\n", file = tex_delta, append = TRUE)
cat("\\midrule\n", file = tex_delta, append = TRUE)
cat(
  "Gaussian & ",
  fmt(lr_delta_table$LR_delta, 3), " & ",
  ifelse(lr_delta_table$p_boot < 0.001, "$<0.001$", fmt(lr_delta_table$p_boot, 3)), " & ",
  lr_delta_table$B_boot, " & ",
  lr_delta_table$boot_ok, " \\\\\n",
  sep = "",
  file = tex_delta,
  append = TRUE
)
cat("\\bottomrule\n\\end{tabular}\n", file = tex_delta, append = TRUE)

## ------------------------------------------------------------
## figures
## ------------------------------------------------------------

plot_faithful_empirical_density <- function() {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("MASS package is required for kde2d")
  }

  kd <- MASS::kde2d(W[, 1], W[, 2], n = 150)

  plot(W[, 1], W[, 2],
       xlab = "eruption duration (robustly standardized)",
       ylab = "waiting time (robustly standardized)",
       pch = 16,
       cex = 0.45)
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

  plot(W[, 1], W[, 2],
       xlab = "eruption duration (robustly standardized)",
       ylab = "waiting time (robustly standardized)",
       pch = 16,
       cex = 0.45)
  contour(xg, yg, zz, add = TRUE, drawlabels = FALSE, nlevels = 8)
}

plot_faithful_whitened <- function() {
  Y <- sweep(W, 2, fit_g$mu, "-")
  Tm <- msec_q(Y, fit_g$delta)

  A <- t(chol(fit_g$Sigma))
  Z <- t(solve(A, t(Tm)))

  plot(Z[, 1], Z[, 2],
       xlab = "whitened coordinate 1",
       ylab = "whitened coordinate 2",
       pch = 16,
       cex = 0.45)

  abline(h = 0, v = 0, lty = 2)

  a <- fit_g$alpha
  nr <- sqrt(sum(a^2))
  if (is.finite(nr) && nr > 0) {
    lim <- max(abs(Z), na.rm = TRUE)
    arrows(0, 0,
           0.65 * lim * a[1] / nr,
           0.65 * lim * a[2] / nr,
           length = 0.10)
    text(0.70 * lim * a[1] / nr,
         0.70 * lim * a[2] / nr,
         labels = "skewness direction",
         pos = 4,
         cex = 0.8)
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

  plot(theo, samp,
       xlab = "theoretical quantiles",
       ylab = "sample radial quantiles",
       pch = 16,
       cex = 0.55)
  abline(0, 1, lty = 2)
}

save_figure("faithful_empirical_density", plot_faithful_empirical_density)
save_figure("faithful_fitted_density", plot_faithful_fitted_density)
save_figure("faithful_whitened_coordinates", plot_faithful_whitened)
save_figure("faithful_radial_qq", plot_faithful_radial_qq)

## manuscript-compatible copies
if (file.exists(file.path(fig_png, "faithful_fitted_density.png"))) {
  file.copy(file.path(fig_png, "faithful_fitted_density.png"),
            "faithful_density2d.png",
            overwrite = TRUE)
}

if (file.exists(file.path(fig_png, "faithful_whitened_coordinates.png"))) {
  file.copy(file.path(fig_png, "faithful_whitened_coordinates.png"),
            "whitened_gaussian_faithful.png",
            overwrite = TRUE)
}

cat("application outputs written to:\n")
cat("  ", out_dir, "\n")
cat("figures written to:\n")
cat("  ", fig_pdf, "\n")
cat("  ", fig_png, "\n")
