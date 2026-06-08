#!/bin/bash
#SBATCH --job-name=msec_fgmix
#SBATCH --account=f202500010hpcvlabuminhox
#SBATCH --qos=normal
#SBATCH --partition=dev-x86
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/msec_faithful_gmix_%j.out
#SBATCH --error=logs/msec_faithful_gmix_%j.err

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

mkdir -p logs results/gaussian_mixture_benchmark

Rscript scripts/benchmark_gaussian_mixture.R \
  --mode faithful \
  --seed 20260608 \
  --gmm-starts 30 \
  --msec-starts 20 \
  --msec-maxit 1500
