#!/bin/bash
#SBATCH --job-name=msec_art_v2
#SBATCH --account=f202500010hpcvlabuminhox
#SBATCH --qos=normal
#SBATCH --partition=normal-x86
#SBATCH --time=2-00:00:00
#SBATCH --nodes=4
#SBATCH --ntasks=128
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/msec_article_v2_%j.out
#SBATCH --error=logs/msec_article_v2_%j.err

set -euo pipefail

module purge
module load R

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export R_DEFAULT_PACKAGES="stats,graphics,grDevices,utils,datasets,methods,base"

cd "$SLURM_SUBMIT_DIR"
export MSEC_ROOT="$SLURM_SUBMIT_DIR"

rm -rf results/article_chunks_v2 results/article_summary_raw_v2 results/article_summary_v2
mkdir -p logs results/article_chunks_v2 results/article_summary_raw_v2 results/article_summary_v2

echo "Running article estimation/RMSE/coverage simulations v2"
echo "SLURM_JOB_ID=${SLURM_JOB_ID}"
echo "Date: $(date)"
echo "Working directory: $(pwd)"

srun Rscript scripts/run_chunk.R \
  --scenarios configs/scenarios_article_estim_coverage_v2.csv \
  --out results/article_chunks_v2 \
  --chunk "${SLURM_PROCID}" \
  --n-chunks 128 \
  --seed 20260610 \
  --max-starts 30 \
  --maxit 2000

echo "Article estimation/RMSE/coverage v2 finished at $(date)"
