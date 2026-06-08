#!/bin/bash
#SBATCH --job-name=msec_gmix
#SBATCH --account=f202500010hpcvlabuminhox
#SBATCH --qos=normal
#SBATCH --partition=normal-x86
#SBATCH --time=2-00:00:00
#SBATCH --nodes=4
#SBATCH --ntasks=128
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/msec_gmix_%j.out
#SBATCH --error=logs/msec_gmix_%j.err

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

mkdir -p logs results/gaussian_mixture_benchmark figures/monte_carlo figures/pdf

srun bash -c '
  Rscript scripts/benchmark_gaussian_mixture.R \
    --mode sim \
    --chunk "${SLURM_PROCID}" \
    --n-chunks 128 \
    --seed 20260608
'
