#!/bin/bash
#SBATCH --job-name=msec_mcfig
#SBATCH --account=f202500010hpcvlabuminhox
#SBATCH --qos=normal
#SBATCH --partition=dev-x86
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/msec_mc_figures_%j.out
#SBATCH --error=logs/msec_mc_figures_%j.err

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

mkdir -p logs figures/pdf

Rscript scripts/make_mc_figures_v2_pdfonly.R
