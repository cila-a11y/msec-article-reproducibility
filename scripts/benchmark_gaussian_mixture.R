#!/usr/bin/env Rscript
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
           VECLIB_MAXIMUM_THREADS="1", NUMEXPR_NUM_THREADS="1")
source("R/00_utils.R")
msec_source_all(getwd())

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default=NULL) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) return(default)
  args[idx + 1]
}
mode <- get_arg("--mode", "sim")
chunk <- as.integer(get_arg("--chunk", "0"))
n_chunks <- as.integer(get_arg("--n-chunks", "1"))
seed <- as.integer(get_arg("--seed", "20260608"))
B <- as.integer(get_arg("--B", "200"))
n_train <- as.integer(get_arg("--n", "300"))
n_test <- as.integer(get_arg("--ntest", "5000"))
grid_n <- as.integer(get_arg("--grid-n", "120"))
gmm_starts <- as.integer(get_arg("--gmm-starts", "20"))
msec_starts <- as.integer(get_arg("--msec-starts", "20"))
msec_maxit <- as.integer(get_arg("--msec-maxit", "1500"))

out_dir <- "results/gaussian_mixture_benchmark"
fig_dir <- "figures/monte_carlo"
fig_legacy <- "figures/pdf"
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(fig_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(fig_legacy, recursive=TRUE, showWarnings=FALSE)

logsumexp <- function(x) { m <- max(x); m + log(sum(exp(x - m))) }
near_pd2 <- function(S, eps=1e-7) {
  S <- (S + t(S))/2
  ee <- eigen(S, symmetric=TRUE)
  ee$vectors %*% diag(pmax(ee$values, eps), length(ee$values)) %*% t(ee$vectors)
}
dmvnorm2_log <- function(X, mean, Sigma) {
  X <- as.matrix(X)
  R <- try(chol(Sigma), silent=TRUE)
  if (inherits(R, "try-error")) R <- chol(near_pd2(Sigma))
  Xc <- sweep(X, 2, mean, "-")
  Z <- t(backsolve(R, t(Xc), transpose=TRUE))
  q <- rowSums(Z^2)
  logdet <- 2 * sum(log(diag(R)))
  -0.5 * (2 * log(2*pi) + logdet + q)
}
msec_log_density_points <- function(W, mu, Sigma, delta, alpha, eps=1e-12) {
  W <- as.matrix(W)
  Y <- sweep(W, 2, mu, "-")
  Tm <- msec_q(Y, delta)
  Sinv <- solve(Sigma)
  u <- rowSums((Tm %*% Sinv) * Tm)
  lambda <- msec_lambda_from_alpha(alpha, Sigma)
  eta <- as.vector(Tm %*% lambda)
  logdet <- as.numeric(determinant(Sigma, logarithm=TRUE)$modulus)
  logJ <- rowSums(matrix(log(delta + 1), nrow=nrow(Y), ncol=2, byrow=TRUE) +
                    sweep(log(pmax(abs(Y), eps)), 2, delta, "*"))
  -log(pi) - 0.5*logdet - 0.5*u + pnorm(eta, log.p=TRUE) + logJ
}
msec_density_points <- function(W, mu, Sigma, delta, alpha) {
  exp(msec_log_density_points(W, mu, Sigma, delta, alpha))
}


gmm_k <- function(G) (G - 1) + 2*G + 3*G
gmm_logdens <- function(X, fit) {
  X <- as.matrix(X)
  G <- length(fit$pi)
  logcomp <- matrix(NA_real_, nrow=nrow(X), ncol=G)
  for (g in seq_len(G)) {
    logcomp[, g] <- log(fit$pi[g]) + dmvnorm2_log(X, fit$mu[g, ], fit$Sigma[, , g])
  }
  apply(logcomp, 1, logsumexp)
}
gmm_density <- function(X, fit) exp(gmm_logdens(X, fit))

gmm_fit_one_start <- function(X, G, init_class=NULL, maxit=500, tol=1e-7) {
  X <- as.matrix(X); n <- nrow(X); d <- ncol(X)
  if (G == 1) {
    mu_g <- matrix(colMeans(X), nrow=1)
    Sigma_g <- array(near_pd2(stats::cov(X) + diag(1e-6, d)), dim=c(d,d,1))
    ll <- sum(dmvnorm2_log(X, mu_g[1,], Sigma_g[,,1]))
    return(list(G=1, pi=1, mu=mu_g, Sigma=Sigma_g, logLik=ll, iter=1))
  }
  if (is.null(init_class)) init_class <- sample(rep(seq_len(G), length.out=n))
  pi_g <- pmax(as.numeric(tabulate(init_class, nbins=G))/n, 1e-4)
  pi_g <- pi_g / sum(pi_g)
  mu_g <- matrix(0, G, d)
  Sigma_g <- array(0, dim=c(d,d,G))
  S_all <- near_pd2(stats::cov(X) + diag(1e-6, d))
  for (g in seq_len(G)) {
    idx <- which(init_class == g)
    if (length(idx) < d + 1) {
      mu_g[g,] <- X[sample(seq_len(n), 1),]
      Sigma_g[,,g] <- S_all
    } else {
      mu_g[g,] <- colMeans(X[idx,,drop=FALSE])
      Sigma_g[,,g] <- near_pd2(stats::cov(X[idx,,drop=FALSE]) + diag(1e-6, d))
    }
  }
  ll_old <- -Inf; ll <- -Inf
  for (iter in seq_len(maxit)) {
    logcomp <- matrix(NA_real_, nrow=n, ncol=G)
    for (g in seq_len(G)) logcomp[,g] <- log(pi_g[g]) + dmvnorm2_log(X, mu_g[g,], Sigma_g[,,g])
    lse <- apply(logcomp, 1, logsumexp)
    ll <- sum(lse)
    tau <- exp(logcomp - lse)
    nk <- pmax(colSums(tau), 1e-8)
    pi_g <- nk / n
    for (g in seq_len(G)) {
      wg <- tau[,g]
      mu_g[g,] <- colSums(X * wg) / nk[g]
      Xc <- sweep(X, 2, mu_g[g,], "-")
      Sigma_g[,,g] <- near_pd2(crossprod(Xc * sqrt(wg))/nk[g] + diag(1e-6, d))
    }
    if (is.finite(ll_old) && abs(ll - ll_old) < tol * (1 + abs(ll_old))) break
    ll_old <- ll
  }
  list(G=G, pi=pi_g, mu=mu_g, Sigma=Sigma_g, logLik=ll, iter=iter)
}
gmm_fit <- function(X, G, n_starts=20, maxit=500) {
  X <- as.matrix(X); n <- nrow(X)
  starts <- list()
  if (G == 1) starts[[1]] <- rep(1, n) else {
    km <- try(kmeans(X, centers=G, nstart=10), silent=TRUE)
    if (!inherits(km, "try-error")) starts[[1]] <- km$cluster
  }
  while (length(starts) < n_starts) starts[[length(starts)+1]] <- sample(rep(seq_len(G), length.out=n))
  fits <- lapply(starts, function(cl) try(gmm_fit_one_start(X, G, cl, maxit=maxit), silent=TRUE))
  fits <- Filter(function(z) !inherits(z, "try-error") && is.finite(z$logLik), fits)
  if (!length(fits)) stop("all Gaussian mixture starts failed for G=", G)
  fits[[which.max(vapply(fits, function(z) z$logLik, numeric(1)))]]
}
gmm_fit_select_bic <- function(X, Gmax=4, n_starts=20) {
  fits <- lapply(seq_len(Gmax), function(G) gmm_fit(X, G, n_starts=n_starts))
  n <- nrow(X)
  crit <- data.frame(G=seq_len(Gmax),
                     logLik=vapply(fits, function(z) z$logLik, numeric(1)),
                     k=vapply(seq_len(Gmax), gmm_k, numeric(1)))
  crit$AIC <- -2*crit$logLik + 2*crit$k
  crit$BIC <- -2*crit$logLik + log(n)*crit$k
  best <- which.min(crit$BIC)
  list(best=fits[[best]], all=fits, criteria=crit, selected_G=crit$G[best])
}


grid_for_ise <- function(lim=3.2, n=120) {
  x <- seq(-lim, lim, length.out=n); y <- seq(-lim, lim, length.out=n)
  list(points=as.matrix(expand.grid(x, y)), dx=diff(x)[1], dy=diff(y)[1])
}
ise_grid <- function(fhat, ftrue, dx, dy) sum((fhat - ftrue)^2) * dx * dy
fit_msec_gaussian <- function(W, max_starts=msec_starts, maxit=msec_maxit) {
  msec_fit(W, generator="gaussian", nu=6, max_starts=max_starts, maxit=maxit, compute_grad=FALSE)
}
scenario_params <- function(scenario, n=300) {
  delta <- switch(scenario,
                  weak=c(0.15, 0.10),
                  moderate=c(0.70, 0.40),
                  strong=c(1.30, 0.90),
                  stop("unknown scenario"))
  msec_standard_params(n=n, generator="gaussian", nu=6, rho=0.35,
                       delta=delta, alpha=c(1.2, -0.8))
}
msec_fit_row <- function(fit, train, test, true_params, grid) {
  k <- length(fit$par); ll <- fit$logLik; n <- nrow(train)
  test_log_score <- mean(msec_log_density_points(test, fit$mu, fit$Sigma, fit$delta, fit$alpha))
  fhat <- msec_density_points(grid$points, fit$mu, fit$Sigma, fit$delta, fit$alpha)
  ftrue <- msec_density_points(grid$points, true_params$mu, true_params$Sigma, true_params$delta, true_params$alpha)
  data.frame(model="gaussian msec", selected_G=1, logLik=ll,
             AIC=-2*ll+2*k, BIC=-2*ll+log(n)*k,
             test_log_score=test_log_score,
             grid_ISE=ise_grid(fhat, ftrue, grid$dx, grid$dy),
             stringsAsFactors=FALSE)
}
gmm_fit_row <- function(fit_select, train, test, true_params, grid) {
  fit <- fit_select$best; G <- fit_select$selected_G
  k <- gmm_k(G); ll <- fit$logLik; n <- nrow(train)
  test_log_score <- mean(gmm_logdens(test, fit))
  fhat <- gmm_density(grid$points, fit)
  ftrue <- msec_density_points(grid$points, true_params$mu, true_params$Sigma, true_params$delta, true_params$alpha)
  data.frame(model="gaussian mixture", selected_G=G, logLik=ll,
             AIC=-2*ll+2*k, BIC=-2*ll+log(n)*k,
             test_log_score=test_log_score,
             grid_ISE=ise_grid(fhat, ftrue, grid$dx, grid$dy),
             stringsAsFactors=FALSE)
}
run_one_sim <- function(scenario, rep_id, n=n_train, ntest=n_test) {
  set.seed(seed + 100000 * match(scenario, c("weak","moderate","strong")) + rep_id)
  pars <- scenario_params(scenario, n=n)
  train <- msec_r_msec(n, pars$mu, pars$Sigma, pars$delta, pars$alpha, generator="gaussian", nu=6)
  test <- msec_r_msec(ntest, pars$mu, pars$Sigma, pars$delta, pars$alpha, generator="gaussian", nu=6)
  grid <- grid_for_ise(lim=3.2, n=grid_n)
  fm <- try(fit_msec_gaussian(train), silent=TRUE)
  fg <- try(gmm_fit_select_bic(train, Gmax=4, n_starts=gmm_starts), silent=TRUE)
  rows <- list()
  if (!inherits(fm, "try-error") && isTRUE(fm$ok)) {
    rows[[1]] <- msec_fit_row(fm, train, test, pars, grid)
  } else {
    rows[[1]] <- data.frame(model="gaussian msec", selected_G=1, logLik=NA, AIC=NA, BIC=NA, test_log_score=NA, grid_ISE=NA)
  }
  if (!inherits(fg, "try-error")) {
    rows[[2]] <- gmm_fit_row(fg, train, test, pars, grid)
  } else {
    rows[[2]] <- data.frame(model="gaussian mixture", selected_G=NA, logLik=NA, AIC=NA, BIC=NA, test_log_score=NA, grid_ISE=NA)
  }
  out <- do.call(rbind, rows)
  out$scenario <- scenario; out$rep <- rep_id; out$n <- n; out$ntest <- ntest
  out
}
se <- function(x) { x <- x[is.finite(x)]; if (length(x)<=1) NA_real_ else sd(x)/sqrt(length(x)) }
mode_int <- function(x) { x <- x[is.finite(x)]; if (!length(x)) return(NA_real_); ux <- sort(unique(x)); ux[which.max(tabulate(match(x, ux)))] }
summarise_sim <- function(d) {
  key <- unique(d[, c("scenario","model")]); out <- list()
  for (i in seq_len(nrow(key))) {
    ds <- d[d$scenario==key$scenario[i] & d$model==key$model[i],]
    out[[length(out)+1]] <- data.frame(
      scenario=key$scenario[i], model=key$model[i],
      selected_G_mean=mean(ds$selected_G, na.rm=TRUE),
      selected_G_mode=mode_int(ds$selected_G),
      BIC_mean=mean(ds$BIC, na.rm=TRUE), BIC_se=se(ds$BIC),
      test_log_score_mean=mean(ds$test_log_score, na.rm=TRUE), test_log_score_se=se(ds$test_log_score),
      grid_ISE_mean=mean(ds$grid_ISE, na.rm=TRUE), grid_ISE_se=se(ds$grid_ISE),
      n_used=sum(is.finite(ds$BIC)), stringsAsFactors=FALSE)
  }
  do.call(rbind, out)
}


write_sim_tex <- function(s, file) {
  fmt <- function(x, d=3) ifelse(is.finite(x), formatC(x, format="f", digits=d), "--")
  fmtse <- function(m, z, d=3) if (!is.finite(m)) "--" else paste0(fmt(m,d), " (", fmt(z,d), ")")
  cat("\\begin{tabular}{llrrrr}\n", file=file)
  cat("\\toprule\n", file=file, append=TRUE)
  cat("Scenario & Model & Selected $G$ & BIC & Test log score & Grid ISE \\\\\n", file=file, append=TRUE)
  cat("\\midrule\n", file=file, append=TRUE)
  for (sc in c("weak","moderate","strong")) {
    for (mo in c("gaussian msec","gaussian mixture")) {
      z <- s[s$scenario==sc & s$model==mo,]
      if (!nrow(z)) next
      cat(sc, " & ", mo, " & ",
          ifelse(mo=="gaussian msec", "1", fmt(z$selected_G_mode,0)), " & ",
          fmtse(z$BIC_mean,z$BIC_se,1), " & ",
          fmtse(z$test_log_score_mean,z$test_log_score_se,3), " & ",
          fmtse(z$grid_ISE_mean,z$grid_ISE_se,4), " \\\\\n",
          sep="", file=file, append=TRUE)
    }
    if (sc != "strong") cat("\\addlinespace\n", file=file, append=TRUE)
  }
  cat("\\bottomrule\n\\end{tabular}\n", file=file, append=TRUE)
}
plot_metric <- function(s, metric, outfile, ylab) {
  pdf(outfile, width=7.2, height=4.6, useDingbats=FALSE)
  old <- par(no.readonly=TRUE); on.exit({par(old); dev.off()}, add=TRUE)
  par(mar=c(4.2,4.2,1.0,1.0))
  scen <- c("weak","moderate","strong")
  mods <- c("gaussian msec","gaussian mixture")
  mat <- matrix(NA_real_, nrow=2, ncol=3)
  rownames(mat) <- c("msec","mixture"); colnames(mat) <- scen
  for (i in seq_along(mods)) {
    for (j in seq_along(scen)) {
      z <- s[s$model==mods[i] & s$scenario==scen[j],]
      if (nrow(z)) mat[i,j] <- z[[metric]]
    }
  }
  barplot(mat, beside=TRUE, ylab=ylab, xlab="scenario",
          legend.text=rownames(mat), args.legend=list(bty="n", x="topleft"))
}
make_benchmark_figures <- function(s) {
  plot_metric(s, "grid_ISE_mean", "figures/monte_carlo/gaussian_mixture_benchmark_ise.pdf", "grid ise")
  file.copy("figures/monte_carlo/gaussian_mixture_benchmark_ise.pdf",
            "figures/pdf/gaussian_mixture_benchmark_ise.pdf", overwrite=TRUE)
  plot_metric(s, "test_log_score_mean", "figures/monte_carlo/gaussian_mixture_benchmark_logscore.pdf", "test log score")
  file.copy("figures/monte_carlo/gaussian_mixture_benchmark_logscore.pdf",
            "figures/pdf/gaussian_mixture_benchmark_logscore.pdf", overwrite=TRUE)
}
faithful_data <- function() {
  W_raw <- as.matrix(datasets::faithful[, c("eruptions","waiting")])
  center <- apply(W_raw, 2L, median)
  scalev <- apply(W_raw, 2L, mad, constant=1.4826)
  as.matrix(sweep(sweep(W_raw, 2L, center, "-"), 2L, scalev, "/"))
}
cv_logscore_msec <- function(W, fold, K) {
  vals <- numeric(K)
  for (k in seq_len(K)) {
    train <- W[fold != k,,drop=FALSE]; test <- W[fold == k,,drop=FALSE]
    fit <- fit_msec_gaussian(train, max_starts=min(msec_starts, 20), maxit=msec_maxit)
    vals[k] <- mean(msec_log_density_points(test, fit$mu, fit$Sigma, fit$delta, fit$alpha))
  }
  mean(vals)
}
cv_logscore_gmm <- function(W, fold, K) {
  vals <- numeric(K); selected <- integer(K)
  for (k in seq_len(K)) {
    train <- W[fold != k,,drop=FALSE]; test <- W[fold == k,,drop=FALSE]
    fit <- gmm_fit_select_bic(train, Gmax=4, n_starts=gmm_starts)
    vals[k] <- mean(gmm_logdens(test, fit$best)); selected[k] <- fit$selected_G
  }
  list(score=mean(vals), selected=selected)
}
run_faithful_benchmark <- function() {
  W <- faithful_data(); n <- nrow(W)
  fit_file <- "results/faithful_v2/faithful_application_fits.rds"
  if (file.exists(fit_file)) {
    fit_m <- readRDS(fit_file)$fit_gaussian
  } else {
    fit_m <- fit_msec_gaussian(W, max_starts=40, maxit=2500)
  }
  gm <- gmm_fit_select_bic(W, Gmax=4, n_starts=max(gmm_starts, 30))
  set.seed(seed); K <- 10; fold <- sample(rep(seq_len(K), length.out=n))
  cv_m <- cv_logscore_msec(W, fold, K); cv_g <- cv_logscore_gmm(W, fold, K)
  km <- length(fit_m$par); kg <- gmm_k(gm$selected_G)
  rows <- rbind(
    data.frame(model="gaussian msec", selected_G=1, logLik=fit_m$logLik,
               AIC=-2*fit_m$logLik+2*km, BIC=-2*fit_m$logLik+log(n)*km,
               cv_log_score=cv_m, stringsAsFactors=FALSE),
    data.frame(model="gaussian mixture", selected_G=gm$selected_G, logLik=gm$best$logLik,
               AIC=-2*gm$best$logLik+2*kg, BIC=-2*gm$best$logLik+log(n)*kg,
               cv_log_score=cv_g$score, stringsAsFactors=FALSE))
  write.csv(rows, file.path(out_dir, "faithful_benchmark.csv"), row.names=FALSE)
  fmt <- function(x, d=3) formatC(x, format="f", digits=d)
  tex <- file.path(out_dir, "faithful_benchmark.tex")
  cat("\\begin{tabular}{lrrrrr}\n\\toprule\n", file=tex)
  cat("Model & Selected $G$ & LL & AIC & BIC & CV log score \\\\\n\\midrule\n", file=tex, append=TRUE)
  for (i in seq_len(nrow(rows))) {
    cat(rows$model[i], " & ", rows$selected_G[i], " & ", fmt(rows$logLik[i],3), " & ",
        fmt(rows$AIC[i],3), " & ", fmt(rows$BIC[i],3), " & ", fmt(rows$cv_log_score[i],3), " \\\\\n",
        sep="", file=tex, append=TRUE)
  }
  cat("\\bottomrule\n\\end{tabular}\n", file=tex, append=TRUE)
  rows
}


if (mode == "sim") {
  scenarios <- c("weak","moderate","strong")
  jobs <- expand.grid(scenario=scenarios, rep=seq_len(B), stringsAsFactors=FALSE)
  jobs <- jobs[seq_len(nrow(jobs)) %% n_chunks == chunk, , drop=FALSE]
  if (!nrow(jobs)) {
    out <- data.frame()
  } else {
    rows <- vector("list", nrow(jobs))
    for (i in seq_len(nrow(jobs))) {
      cat("chunk", chunk, "scenario", jobs$scenario[i], "rep", jobs$rep[i], "\n", file=stderr())
      rows[[i]] <- tryCatch(run_one_sim(jobs$scenario[i], jobs$rep[i]),
        error=function(e) data.frame(
          model=c("gaussian msec","gaussian mixture"), selected_G=c(1,NA),
          logLik=NA, AIC=NA, BIC=NA, test_log_score=NA, grid_ISE=NA,
          scenario=jobs$scenario[i], rep=jobs$rep[i], n=n_train, ntest=n_test,
          error=conditionMessage(e), stringsAsFactors=FALSE))
    }
    out <- do.call(rbind, rows)
  }
  file <- file.path(out_dir, sprintf("sim_chunk_%04d.csv", chunk))
  write.csv(out, file, row.names=FALSE)
  cat("wrote", file, "rows", nrow(out), "\n", file=stderr())
} else if (mode == "combine") {
  files <- sort(list.files(out_dir, pattern="^sim_chunk_.*\\.csv$", full.names=TRUE))
  if (!length(files)) stop("no simulation chunk files found")
  d <- do.call(rbind, lapply(files, read.csv, stringsAsFactors=FALSE))
  write.csv(d, file.path(out_dir, "simulated_benchmark_all.csv"), row.names=FALSE)
  s <- summarise_sim(d)
  write.csv(s, file.path(out_dir, "simulated_benchmark_summary.csv"), row.names=FALSE)
  write_sim_tex(s, file.path(out_dir, "simulated_benchmark_summary.tex"))
  make_benchmark_figures(s)
  cat("combined", length(files), "files\n")
  cat("wrote simulated benchmark summaries\n")
} else if (mode == "faithful") {
  rows <- run_faithful_benchmark()
  print(rows)
} else {
  stop("unknown mode: ", mode)
}
