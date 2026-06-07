#!/bin/bash
#SBATCH --job-name=msec_combine
#SBATCH --account=f202500010hpcvlabuminhox
#SBATCH --qos=normal
#SBATCH --partition=normal-x86
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/msec_combine_%j.out
#SBATCH --error=logs/msec_combine_%j.err

set -euo pipefail
module purge
module load R
cd "$SLURM_SUBMIT_DIR"
export MSEC_ROOT="$SLURM_SUBMIT_DIR"
mkdir -p results/summary
Rscript scripts/combine_results.R --chunks results/chunks --out results/summary
