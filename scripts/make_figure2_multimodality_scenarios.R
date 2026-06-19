#!/usr/bin/env Rscript

## ------------------------------------------------------------------
## Figure 2: no-warping, weak, moderate, and strong warping scenarios
##
## This script must be run from the root of the repository:
##
##   Rscript scripts/make_figure2_multimodality_scenarios.R
##
## It uses the baseline parameters stated in the manuscript:
##
##   mu    = (0.3, -0.2)
##   sigma = (1.1, 0.9)
##   rho   = 0.35
##   alpha = (1.2, -0.8)
##
## Only delta varies across the four panels.
## ------------------------------------------------------------------

local({

  ## Avoid unintended multithreading in numerical libraries.
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    NUMEXPR_NUM_THREADS = "1"
  )

  ## --------------------------------------------------------------
  ## Check working directory
  ## --------------------------------------------------------------

  required_paths <- c(
    "scripts",
    "figures",
    "results"
  )

  missing_paths <- required_paths[!file.exists(required_paths)]

  if (length(missing_paths) > 0L) {
    stop(
      paste0(
        "Run this script from the root of the repository.\n",
        "Missing paths: ",
        paste(missing_paths, collapse = ", ")
      )
    )
  }

  ## --------------------------------------------------------------
  ## Output directories
  ## --------------------------------------------------------------

  dir.create(
    "figures/monte_carlo",
    recursive = TRUE,
    showWarnings = FALSE
  )

  dir.create(
    "figures/pdf",
    recursive = TRUE,
    showWarnings = FALSE
  )

  dir.create(
    "results/explanatory_figures",
    recursive = TRUE,
    showWarnings = FALSE
  )

  ## --------------------------------------------------------------
  ## Model functions
  ## --------------------------------------------------------------

  fig2_make_sigma <- function(rho, sigma) {

    if (
      length(sigma) != 2L ||
      any(!is.finite(sigma)) ||
      any(sigma <= 0) ||
      !is.finite(rho) ||
      abs(rho) >= 1
    ) {
      stop("Invalid scale or correlation parameters.")
    }

    matrix(
      c(
        sigma[1]^2,
        rho * sigma[1] * sigma[2],
        rho * sigma[1] * sigma[2],
        sigma[2]^2
      ),
      nrow = 2L,
      byrow = TRUE
    )
  }


  fig2_qwarp <- function(Y, delta) {

    Y <- as.matrix(Y)

    if (
      ncol(Y) != 2L ||
      length(delta) != 2L ||
      any(!is.finite(delta)) ||
      any(delta < 0)
    ) {
      stop("Invalid coordinates or warp exponents.")
    }

    out <- Y

    out[, 1] <-
      sign(Y[, 1]) *
      abs(Y[, 1])^(1 + delta[1])

    out[, 2] <-
      sign(Y[, 2]) *
      abs(Y[, 2])^(1 + delta[2])

    out
  }


  fig2_lambda_from_alpha <- function(alpha, Sigma) {

    if (
      length(alpha) != 2L ||
      any(!is.finite(alpha)) ||
      !all(dim(Sigma) == c(2L, 2L))
    ) {
      stop("Invalid whitened skewness vector or scale matrix.")
    }

    ## Base R returns an upper-triangular factor R such that
    ## t(R) %*% R = Sigma. Hence A = t(R) satisfies
    ## A %*% t(A) = Sigma.
    A <- t(chol(Sigma))

    ## The manuscript uses alpha = t(A) %*% lambda.
    as.vector(
      solve(
        t(A),
        alpha
      )
    )
  }


  fig2_log_jacobian <- function(Y, delta) {

    Y <- as.matrix(Y)

    if (
      ncol(Y) != 2L ||
      length(delta) != 2L ||
      any(delta < 0)
    ) {
      stop("Invalid coordinates or warp exponents.")
    }

    n <- nrow(Y)
    log_jacobian <- rep(0, n)

    for (j in 1:2) {

      ## log(a_j), with a_j = 1 + delta_j.
      log_jacobian <-
        log_jacobian +
        log1p(delta[j])

      ## When delta_j = 0, the contribution
      ## delta_j * log|y_j| is identically zero.
      if (delta[j] > 0) {

        abs_y <- abs(Y[, j])
        positive <- abs_y > 0

        contribution <- rep(-Inf, n)

        contribution[positive] <-
          delta[j] *
          log(abs_y[positive])

        log_jacobian <-
          log_jacobian +
          contribution
      }
    }

    log_jacobian
  }


  fig2_msec_logdensity_gaussian <- function(
      W,
      mu,
      sigma,
      rho,
      delta,
      alpha
  ) {

    W <- as.matrix(W)

    if (
      ncol(W) != 2L ||
      length(mu) != 2L ||
      length(sigma) != 2L ||
      length(delta) != 2L ||
      length(alpha) != 2L ||
      any(!is.finite(mu)) ||
      any(!is.finite(sigma)) ||
      any(!is.finite(delta)) ||
      any(!is.finite(alpha)) ||
      any(sigma <= 0) ||
      any(delta < 0) ||
      !is.finite(rho) ||
      abs(rho) >= 1
    ) {
      stop("Invalid Gaussian MSEC parameters.")
    }

    ## Observed coordinates centred at the warp origin.
    Y <- sweep(
      W,
      MARGIN = 2L,
      STATS = mu,
      FUN = "-"
    )

    ## Odd power transformation.
    Tmat <- fig2_qwarp(
      Y = Y,
      delta = delta
    )

    ## Transformed-core scale matrix.
    Sigma <- fig2_make_sigma(
      rho = rho,
      sigma = sigma
    )

    Sigma_inv <- solve(Sigma)

    ## Squared transformed Mahalanobis radius.
    u <- rowSums(
      (Tmat %*% Sigma_inv) *
        Tmat
    )

    ## Native skewness corresponding to the fixed whitened vector.
    lambda <- fig2_lambda_from_alpha(
      alpha = alpha,
      Sigma = Sigma
    )

    skew_index <- as.vector(
      Tmat %*% lambda
    )

    log_det_sigma <- as.numeric(
      determinant(
        Sigma,
        logarithm = TRUE
      )$modulus
    )

    log_jacobian <- fig2_log_jacobian(
      Y = Y,
      delta = delta
    )

    ## For the Gaussian generator:
    ##
    ##   g(u) = exp(-u/2),
    ##   F_c  = Phi,
    ##   c_g  = 2*pi.
    ##
    ## Therefore:
    ##
    ##   2 / {sqrt(|Sigma|) * c_g}
    ##   =
    ##   1 / {pi * sqrt(|Sigma|)}.
    -log(pi) -
      0.5 * log_det_sigma -
      0.5 * u +
      pnorm(
        skew_index,
        log.p = TRUE
      ) +
      log_jacobian
  }


  fig2_density_grid <- function(
      mu,
      sigma,
      rho,
      delta,
      alpha,
      xlim,
      ylim,
      n = 321L
  ) {

    if (
      length(xlim) != 2L ||
      length(ylim) != 2L ||
      n < 50L
    ) {
      stop("Invalid grid specification.")
    }

    x <- seq(
      from = xlim[1],
      to = xlim[2],
      length.out = n
    )

    y <- seq(
      from = ylim[1],
      to = ylim[2],
      length.out = n
    )

    grid <- expand.grid(
      w1 = x,
      w2 = y
    )

    log_density <- fig2_msec_logdensity_gaussian(
      W = grid,
      mu = mu,
      sigma = sigma,
      rho = rho,
      delta = delta,
      alpha = alpha
    )

    density <- exp(log_density)

    list(
      x = x,
      y = y,
      z = matrix(
        density,
        nrow = length(x),
        ncol = length(y),
        byrow = FALSE
      )
    )
  }


  fig2_draw_panel <- function(
      density_grid,
      panel,
      scenario,
      delta
  ) {

    finite_density <-
      density_grid$z[
        is.finite(density_grid$z)
      ]

    if (length(finite_density) == 0L) {
      stop("The density grid contains no finite values.")
    }

    maximum_density <- max(finite_density)

    ## Common relative-density contour levels are used in every panel
    ## to emphasize the change in modal geometry.
    relative_levels <- c(
      0.03,
      0.06,
      0.10,
      0.16,
      0.25,
      0.38,
      0.52,
      0.68,
      0.82,
      0.93
    )

    contour_levels <-
      maximum_density *
      relative_levels

    panel_title <- sprintf(
      "%s %s\ndelta = (%.2f, %.2f)",
      panel,
      scenario,
      delta[1],
      delta[2]
    )

    plot(
      NA,
      NA,
      xlim = range(density_grid$x),
      ylim = range(density_grid$y),
      xlab = expression(w[1]),
      ylab = expression(w[2]),
      main = panel_title,
      asp = 1,
      xaxs = "i",
      yaxs = "i",
      cex.main = 0.88
    )

    contour(
      x = density_grid$x,
      y = density_grid$y,
      z = density_grid$z,
      levels = contour_levels,
      add = TRUE,
      drawlabels = FALSE,
      lwd = 0.95
    )

    box()
  }


  fig2_save_pdf <- function(
      filename,
      width,
      height,
      plot_function
  ) {

    grDevices::pdf(
      file = filename,
      width = width,
      height = height,
      useDingbats = FALSE,
      onefile = TRUE
    )

    on.exit(
      grDevices::dev.off(),
      add = TRUE
    )

    plot_function()

    invisible(filename)
  }

  ## --------------------------------------------------------------
  ## Parameters used in the manuscript
  ## --------------------------------------------------------------

  baseline_mu <- c(
    0.3,
    -0.2
  )

  baseline_sigma <- c(
    1.1,
    0.9
  )

  baseline_rho <- 0.35

  baseline_alpha <- c(
    1.2,
    -0.8
  )

  baseline_Sigma <- fig2_make_sigma(
    rho = baseline_rho,
    sigma = baseline_sigma
  )

  baseline_lambda <- fig2_lambda_from_alpha(
    alpha = baseline_alpha,
    Sigma = baseline_Sigma
  )

  ## --------------------------------------------------------------
  ## Four warping scenarios
  ## --------------------------------------------------------------

  shape_cases <- data.frame(
    panel = c(
      "(a)",
      "(b)",
      "(c)",
      "(d)"
    ),
    scenario = c(
      "No warping",
      "Weak",
      "Moderate",
      "Strong"
    ),
    mu1 = rep(
      baseline_mu[1],
      4L
    ),
    mu2 = rep(
      baseline_mu[2],
      4L
    ),
    sigma1 = rep(
      baseline_sigma[1],
      4L
    ),
    sigma2 = rep(
      baseline_sigma[2],
      4L
    ),
    rho = rep(
      baseline_rho,
      4L
    ),
    delta1 = c(
      0.00,
      0.15,
      0.70,
      1.30
    ),
    delta2 = c(
      0.00,
      0.10,
      0.40,
      0.90
    ),
    alpha1 = rep(
      baseline_alpha[1],
      4L
    ),
    alpha2 = rep(
      baseline_alpha[2],
      4L
    ),
    lambda1 = rep(
      baseline_lambda[1],
      4L
    ),
    lambda2 = rep(
      baseline_lambda[2],
      4L
    ),
    description = c(
      "identity warp; the model reduces to the skew-elliptical core",
      "shallow internal troughs and small modal displacement",
      "baseline setting; visible but not extreme modal separation",
      "pronounced troughs and larger separation of modal regions"
    ),
    stringsAsFactors = FALSE
  )

  shape_cases$alpha_norm <-
    sqrt(
      shape_cases$alpha1^2 +
        shape_cases$alpha2^2
    )

  shape_cases$alpha_angle_degrees <-
    atan2(
      shape_cases$alpha2,
      shape_cases$alpha1
    ) *
    180 / pi

  ## --------------------------------------------------------------
  ## Internal consistency checks
  ## --------------------------------------------------------------

  fixed_columns <- c(
    "mu1",
    "mu2",
    "sigma1",
    "sigma2",
    "rho",
    "alpha1",
    "alpha2",
    "lambda1",
    "lambda2"
  )

  for (column in fixed_columns) {

    if (
      length(
        unique(shape_cases[[column]])
      ) != 1L
    ) {
      stop(
        paste0(
          "Parameter '",
          column,
          "' must remain fixed across all four panels."
        )
      )
    }
  }

  expected_scenarios <- c(
    "No warping",
    "Weak",
    "Moderate",
    "Strong"
  )

  if (
    nrow(shape_cases) != 4L ||
    !identical(
      shape_cases$scenario,
      expected_scenarios
    )
  ) {
    stop("The four Figure 2 scenarios are not correctly specified.")
  }

  expected_delta1 <- c(
    0.00,
    0.15,
    0.70,
    1.30
  )

  expected_delta2 <- c(
    0.00,
    0.10,
    0.40,
    0.90
  )

  if (
    !isTRUE(
      all.equal(
        shape_cases$delta1,
        expected_delta1
      )
    ) ||
    !isTRUE(
      all.equal(
        shape_cases$delta2,
        expected_delta2
      )
    )
  ) {
    stop("The warp exponents do not match the manuscript table.")
  }

  ## --------------------------------------------------------------
  ## Save the exact parameter configuration
  ## --------------------------------------------------------------

  parameter_file <-
    "results/explanatory_figures/mc_shape_scenarios.csv"

  write.csv(
    shape_cases,
    file = parameter_file,
    row.names = FALSE
  )

  ## --------------------------------------------------------------
  ## Common observed-coordinate grid
  ## --------------------------------------------------------------

  xlim <-
    baseline_mu[1] +
    c(-3.0, 3.0) *
    baseline_sigma[1]

  ylim <-
    baseline_mu[2] +
    c(-3.0, 3.0) *
    baseline_sigma[2]

  density_panels <- vector(
    mode = "list",
    length = nrow(shape_cases)
  )

  for (
    i in seq_len(
      nrow(shape_cases)
    )
  ) {

    delta_i <- c(
      shape_cases$delta1[i],
      shape_cases$delta2[i]
    )

    density_panels[[i]] <-
      fig2_density_grid(
        mu = baseline_mu,
        sigma = baseline_sigma,
        rho = baseline_rho,
        delta = delta_i,
        alpha = baseline_alpha,
        xlim = xlim,
        ylim = ylim,
        n = 321L
      )
  }

  ## --------------------------------------------------------------
  ## Generate Figure 2
  ## --------------------------------------------------------------

  main_figure <-
    "figures/monte_carlo/mc_multimodality_scenarios.pdf"

  mirrored_figure <-
    "figures/pdf/mc_multimodality_scenarios.pdf"

  fig2_save_pdf(
    filename = main_figure,
    width = 7.6,
    height = 7.2,
    plot_function = function() {

      old_par <- par(
        no.readonly = TRUE
      )

      on.exit(
        par(old_par),
        add = TRUE
      )

      par(
        mfrow = c(2, 2),
        mar = c(4.1, 4.2, 3.0, 0.8),
        mgp = c(2.3, 0.7, 0),
        oma = c(0.2, 0.2, 0.2, 0.2)
      )

      for (
        i in seq_len(
          nrow(shape_cases)
        )
      ) {

        delta_i <- c(
          shape_cases$delta1[i],
          shape_cases$delta2[i]
        )

        fig2_draw_panel(
          density_grid = density_panels[[i]],
          panel = shape_cases$panel[i],
          scenario = shape_cases$scenario[i],
          delta = delta_i
        )
      }
    }
  )

  copied <- file.copy(
    from = main_figure,
    to = mirrored_figure,
    overwrite = TRUE
  )

  if (!isTRUE(copied)) {
    stop(
      paste0(
        "Figure 2 was generated, but it could not be copied to:\n",
        mirrored_figure,
        "\nClose any open PDF viewer and run the script again."
      )
    )
  }

  ## Confirm that the two stored PDF copies are identical.
  hashes <- tools::md5sum(
    c(
      main_figure,
      mirrored_figure
    )
  )

  if (
    length(unique(unname(hashes))) != 1L
  ) {
    stop("The two Figure 2 PDF copies are not identical.")
  }

  ## --------------------------------------------------------------
  ## Completion message
  ## --------------------------------------------------------------

  cat(
    paste0(
      "\nFigure 2 generated successfully.\n\n",
      "Main figure:\n  ",
      main_figure,
      "\n\n",
      "Mirrored figure:\n  ",
      mirrored_figure,
      "\n\n",
      "Parameter file:\n  ",
      parameter_file,
      "\n\n",
      "Scenarios:\n",
      "  (a) No warping: delta = (0.00, 0.00)\n",
      "  (b) Weak:       delta = (0.15, 0.10)\n",
      "  (c) Moderate:   delta = (0.70, 0.40)\n",
      "  (d) Strong:     delta = (1.30, 0.90)\n\n",
      "The location, scales, correlation, and whitened skewness\n",
      "parameters are fixed across all four panels.\n"
    )
  )
})
