#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# Combine chunk outputs and generate CSV/LaTeX summaries.
# -----------------------------------------------------------------------------

root_env <- Sys.getenv('MSEC_ROOT', unset = '')
root <- if (nzchar(root_env)) normalizePath(root_env, mustWork = FALSE) else normalizePath(getwd(), mustWork = FALSE)
if (!file.exists(file.path(root, 'R/00_utils.R'))) root <- normalizePath(file.path(getwd(), '..'), mustWork = FALSE)
source(file.path(root, 'R/00_utils.R'))
msec_source_all(root)
opts <- msec_parse_args()
chunk_dir <- if (!is.null(opts$chunks)) opts$chunks else file.path(root, 'results/chunks')
out_dir <- if (!is.null(opts$out)) opts$out else file.path(root, 'results/summary')
msec_safe_mkdir(out_dir)
files <- list.files(chunk_dir, pattern = '^chunk_[0-9]+\\.csv$', full.names = TRUE)
message('Combining ', length(files), ' chunk files from ', chunk_dir)
raw <- msec_combine_csvs(files)
msec_write_csv(raw, file.path(out_dir, 'all_replications.csv'))

if (!nrow(raw)) stop('No rows to summarize')

est_studies <- raw[raw$study %in% c('baseline_accuracy','robust_rmse','sample_size','grid','coverage'), , drop = FALSE]
if (nrow(est_studies)) {
  est <- msec_summarise_estimation(est_studies)
  msec_write_csv(est, file.path(out_dir, 'summary_estimation.csv'))
  msec_write_latex_simple(est, file.path(out_dir, 'summary_estimation.tex'))
}

if ('lr_lambda' %in% raw$study) {
  lr_l <- msec_summarise_lr_lambda(raw[raw$study == 'lr_lambda', , drop = FALSE])
  msec_write_csv(lr_l, file.path(out_dir, 'summary_lr_lambda.csv'))
  msec_write_latex_simple(lr_l, file.path(out_dir, 'summary_lr_lambda.tex'))
}

if ('lr_delta_boot' %in% raw$study) {
  lr_d <- msec_summarise_lr_delta(raw[raw$study == 'lr_delta_boot', , drop = FALSE])
  msec_write_csv(lr_d, file.path(out_dir, 'summary_lr_delta.csv'))
  msec_write_latex_simple(lr_d, file.path(out_dir, 'summary_lr_delta.tex'))
}

if ('coverage' %in% raw$study) {
  cov <- msec_summarise_coverage(raw[raw$study == 'coverage', , drop = FALSE])
  msec_write_csv(cov, file.path(out_dir, 'summary_coverage.csv'))
  msec_write_latex_simple(cov, file.path(out_dir, 'summary_coverage.tex'))
}

if ('ablation_t' %in% raw$study) {
  abl <- msec_summarise_ablation(raw[raw$study == 'ablation_t', , drop = FALSE])
  msec_write_csv(abl, file.path(out_dir, 'summary_ablation_t.csv'))
  msec_write_latex_simple(abl, file.path(out_dir, 'summary_ablation_t.tex'))
}

message('Summaries written to ', out_dir)
