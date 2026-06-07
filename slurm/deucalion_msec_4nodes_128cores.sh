#!/bin/bash
#SBATCH --job-name=msec_mc
#SBATCH --account=f202500010hpcvlabuminhox
#SBATCH --qos=normal
#SBATCH --partition=normal-x86
#SBATCH --time=2-00:00:00
#SBATCH --nodes=4
#SBATCH --ntasks=128
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/msec_mc_%j.out
#SBATCH --error=logs/msec_mc_%j.err

set -euo pipefail

# Deucalion uses environment modules; check the exact R module with: module spider R
module purge
module load R

# One BLAS/OpenMP thread per R process. The parallelism here is across 128 Slurm tasks.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export R_DEFAULT_PACKAGES="stats,graphics,grDevices,utils,datasets,methods,base"

cd "$SLURM_SUBMIT_DIR"
export MSEC_ROOT="$SLURM_SUBMIT_DIR"
mkdir -p logs results/chunks results/summary

# Run the main simulation design. Each Slurm task receives a disjoint subset
# of Monte Carlo replications through SLURM_PROCID.
srun --cpu-bind=cores Rscript scripts/run_chunk.R \
  --scenarios configs/scenarios_main.csv \
  --out results/chunks \
  --seed 7302026 \
  --max-starts 4 \
  --maxit 600
