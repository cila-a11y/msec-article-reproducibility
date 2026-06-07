#!/bin/bash
#SBATCH --job-name=msec_art_c2
#SBATCH --account=f202500010hpcvlabuminhox
#SBATCH --qos=normal
#SBATCH --partition=normal-x86
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/msec_article_combine_v2_%j.out
#SBATCH --error=logs/msec_article_combine_v2_%j.err

set -euo pipefail

module purge
module load R

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export R_DEFAULT_PACKAGES="stats,graphics,grDevices,utils,datasets,methods,base"

cd "$SLURM_SUBMIT_DIR"

mkdir -p results/article_summary_raw_v2 results/article_summary_v2

Rscript scripts/combine_results.R \
  --chunks results/article_chunks_v2 \
  --out results/article_summary_raw_v2

Rscript scripts/recombine_article_results_v2.R \
  results/article_summary_raw_v2/all_replications.csv \
  results/article_summary_v2

cp results/summary/summary_lr_lambda.csv results/article_summary_v2/summary_lr_lambda_article.csv
cp results/summary/summary_lr_delta.csv results/article_summary_v2/summary_lr_delta_article.csv

echo "Article v2 summaries written to results/article_summary_v2"
