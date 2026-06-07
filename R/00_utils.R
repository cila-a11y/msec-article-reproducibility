# -----------------------------------------------------------------------------
# Utility functions for the MSEC replication package
# -----------------------------------------------------------------------------

msec_parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (grepl('^--', key)) {
      key <- sub('^--', '', key)
      if (i == length(args) || grepl('^--', args[[i + 1L]])) {
        out[[key]] <- TRUE
        i <- i + 1L
      } else {
        out[[key]] <- args[[i + 1L]]
        i <- i + 2L
      }
    } else {
      i <- i + 1L
    }
  }
  out
}

msec_as_logical <- function(x, default = FALSE) {
  if (is.null(x) || length(x) == 0L || is.na(x)) return(default)
  if (is.logical(x)) return(isTRUE(x))
  tolower(as.character(x)) %in% c('1', 'true', 't', 'yes', 'y')
}

msec_as_numeric <- function(x, default = NA_real_) {
  if (is.null(x) || length(x) == 0L || is.na(x)) return(default)
  as.numeric(x)
}

msec_as_integer <- function(x, default = NA_integer_) {
  if (is.null(x) || length(x) == 0L || is.na(x)) return(default)
  as.integer(x)
}

msec_project_root <- function() {
  env <- Sys.getenv('MSEC_ROOT', unset = '')
  if (nzchar(env)) return(normalizePath(env, mustWork = FALSE))
  wd <- getwd()
  if (dir.exists(file.path(wd, 'R')) && dir.exists(file.path(wd, 'scripts'))) {
    return(normalizePath(wd, mustWork = FALSE))
  }
  parent <- normalizePath(file.path(wd, '..'), mustWork = FALSE)
  if (dir.exists(file.path(parent, 'R')) && dir.exists(file.path(parent, 'scripts'))) {
    return(parent)
  }
  normalizePath(wd, mustWork = FALSE)
}

msec_source_all <- function(root = msec_project_root()) {
  files <- c('R/00_utils.R', 'R/01_msec_core.R', 'R/02_fit_msec.R', 'R/03_simulation_designs.R', 'R/04_summaries.R')
  for (f in files) {
    path <- file.path(root, f)
    if (!file.exists(path)) stop('Missing source file: ', path)
    source(path)
  }
  invisible(TRUE)
}

msec_safe_mkdir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

msec_write_csv <- function(x, file) {
  msec_safe_mkdir(dirname(file))
  write.csv(x, file = file, row.names = FALSE, na = '')
  invisible(file)
}

msec_read_scenarios <- function(file) {
  if (!file.exists(file)) stop('Scenario file not found: ', file)
  df <- read.csv(file, stringsAsFactors = FALSE)
  required <- c('scenario_id', 'study', 'generator', 'n', 'nu', 'rho',
                'delta1', 'delta2', 'alpha1', 'alpha2', 'B', 'B_boot')
  missing <- setdiff(required, names(df))
  if (length(missing)) stop('Scenario file is missing columns: ', paste(missing, collapse = ', '))
  df
}

msec_rep_indices <- function(B, chunk_id, n_chunks) {
  if (B <= 0L) return(integer(0L))
  which((seq_len(B) - 1L) %% n_chunks == chunk_id)
}


msec_bind_rows_fill <- function(xs) {
  xs <- Filter(function(z) {
    !is.null(z) && is.data.frame(z) && nrow(z) > 0L
  }, xs)

  if (!length(xs)) return(data.frame())

  all_names <- Reduce(union, lapply(xs, names))

  xs <- lapply(xs, function(d) {
    missing <- setdiff(all_names, names(d))
    for (m in missing) d[[m]] <- NA
    d[all_names]
  })

  do.call(rbind, xs)
}

msec_combine_csvs <- function(files) {
  files <- files[file.exists(files)]
  if (!length(files)) return(data.frame())
  xs <- lapply(files, function(f) {
    z <- try(read.csv(f, stringsAsFactors = FALSE), silent = TRUE)
    if (inherits(z, 'try-error')) {
      warning('Could not read ', f)
      return(NULL)
    }
    z$source_file <- basename(f)
    z
  })
  xs <- Filter(Negate(is.null), xs)
  if (!length(xs)) return(data.frame())
  common <- Reduce(union, lapply(xs, names))
  xs <- lapply(xs, function(d) {
    missing <- setdiff(common, names(d))
    for (m in missing) d[[m]] <- NA
    d[common]
  })
  do.call(rbind, xs)
}


msec_num_grad <- function(fn, par, step = NULL, ...) {
  nm <- names(par)
  par <- as.numeric(par)
  if (!is.null(nm)) names(par) <- nm
  n <- length(par)
  if (is.null(step)) step <- pmax(1e-5, abs(par) * 1e-5)
  g <- numeric(n)
  for (j in seq_len(n)) {
    h <- step[j]
    p1 <- par; p2 <- par
    p1[j] <- p1[j] + h
    p2[j] <- p2[j] - h
    f1 <- fn(p1, ...)
    f2 <- fn(p2, ...)
    if (!is.finite(f1) || !is.finite(f2)) {
      h <- pmax(1e-6, abs(par[j]) * 1e-4)
      p1 <- par; p2 <- par
      p1[j] <- p1[j] + h
      p2[j] <- p2[j] - h
      f1 <- fn(p1, ...)
      f2 <- fn(p2, ...)
    }
    g[j] <- (f1 - f2) / (2 * h)
  }
  if (!is.null(nm)) names(g) <- nm
  g
}


msec_num_jacobian <- function(fn, par, step = NULL, ...) {
  nm <- names(par)
  par <- as.numeric(par)
  if (!is.null(nm)) names(par) <- nm
  f0 <- as.numeric(fn(par, ...))
  p <- length(par)
  q <- length(f0)
  if (is.null(step)) step <- pmax(1e-5, abs(par) * 1e-5)
  J <- matrix(NA_real_, nrow = q, ncol = p)
  if (!is.null(nm)) colnames(J) <- nm
  for (j in seq_len(p)) {
    h <- step[j]
    p1 <- par; p2 <- par
    p1[j] <- p1[j] + h
    p2[j] <- p2[j] - h
    f1 <- as.numeric(fn(p1, ...))
    f2 <- as.numeric(fn(p2, ...))
    J[, j] <- (f1 - f2) / (2 * h)
  }
  J
}

msec_near_pd_2x2 <- function(S, eps = 1e-8) {
  S <- (S + t(S)) / 2
  ee <- eigen(S, symmetric = TRUE)
  vals <- pmax(ee$values, eps)
  ee$vectors %*% diag(vals, 2L) %*% t(ee$vectors)
}

msec_clip <- function(x, lower, upper) pmin(pmax(x, lower), upper)

msec_angle_deg <- function(x) {
  if (length(x) != 2L || any(!is.finite(x))) return(NA_real_)
  atan2(x[2], x[1]) * 180 / pi
}

msec_plus_one_pvalue <- function(stat, boot_stats) {
  boot_stats <- boot_stats[is.finite(boot_stats)]
  (1 + sum(boot_stats >= stat)) / (length(boot_stats) + 1)
}
