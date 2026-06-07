#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# Run a chunk of the MSEC Monte Carlo simulations.
# Designed for Slurm multi-node execution: one R process per Slurm task.
# -----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
# Source utilities first without relying on MSEC_ROOT.
root_env <- Sys.getenv('MSEC_ROOT', unset = '')
root <- if (nzchar(root_env)) normalizePath(root_env, mustWork = FALSE) else normalizePath(file.path(getwd()), mustWork = FALSE)
if (!file.exists(file.path(root, 'R/00_utils.R'))) {
  root <- normalizePath(file.path(getwd(), '..'), mustWork = FALSE)
}
source(file.path(root, 'R/00_utils.R'))
msec_source_all(root)

opts <- msec_parse_args(args)
scenario_file <- if (!is.null(opts$scenarios)) opts$scenarios else file.path(root, 'configs/scenarios_main.csv')
out_dir <- if (!is.null(opts$out)) opts$out else file.path(root, 'results/chunks')
seed0 <- msec_as_integer(opts$seed, 7302026L)
max_starts <- msec_as_integer(opts[['max-starts']], msec_as_integer(opts[['max_starts']], 4L))
maxit <- msec_as_integer(opts$maxit, 600L)

slurm_id <- Sys.getenv('SLURM_PROCID', unset = '')
slurm_ntasks <- Sys.getenv('SLURM_NTASKS', unset = '')
chunk_id <- if (nzchar(slurm_id)) as.integer(slurm_id) else msec_as_integer(opts$chunk, 0L)
n_chunks <- if (nzchar(slurm_ntasks)) as.integer(slurm_ntasks) else msec_as_integer(opts[['n-chunks']], msec_as_integer(opts[['n_chunks']], 1L))
if (is.na(chunk_id) || chunk_id < 0L) chunk_id <- 0L
if (is.na(n_chunks) || n_chunks < 1L) n_chunks <- 1L

msec_safe_mkdir(out_dir)
sc <- msec_read_scenarios(scenario_file)

message('MSEC chunk ', chunk_id, '/', n_chunks - 1L,
        ' scenarios=', nrow(sc), ' max_starts=', max_starts, ' maxit=', maxit)

rows <- list()
count <- 0L
for (s in seq_len(nrow(sc))) {
  row <- sc[s, ]
  B <- as.integer(row$B)
  reps <- msec_rep_indices(B, chunk_id, n_chunks)
  if (!length(reps)) next
  message('scenario ', row$scenario_id, ': ', length(reps), ' replications')
  for (rep_id in reps) {
    rep_seed <- seed0 + 1000000L * s + rep_id
    z <- try(msec_run_scenario_rep(row, rep_id, rep_seed,
                                   max_starts = max_starts, maxit = maxit), silent = TRUE)
    if (inherits(z, 'try-error')) {
      z <- data.frame(scenario_id = row$scenario_id, study = row$study, rep = rep_id,
                      generator = row$generator, n = row$n, nu = row$nu,
                      error = as.character(z), stringsAsFactors = FALSE)
      warning('replication failed: scenario=', row$scenario_id, ' rep=', rep_id)
    }
    rows[[length(rows) + 1L]] <- z
    count <- count + 1L
    if (count %% 10L == 0L) {
      tmp <- msec_bind_rows_fill(rows)
      msec_write_csv(tmp, file.path(out_dir, sprintf('chunk_%04d_partial.csv', chunk_id)))
    }
  }
}

res <- if (length(rows)) msec_bind_rows_fill(rows) else data.frame()
outfile <- file.path(out_dir, sprintf('chunk_%04d.csv', chunk_id))
msec_write_csv(res, outfile)
message('wrote ', outfile, ' rows=', nrow(res))
