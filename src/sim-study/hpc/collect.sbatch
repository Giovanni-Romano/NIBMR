#!/bin/bash -l
#SBATCH --job-name=mode-collect
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem=8G
#SBATCH --time=01:00:00
set -euo pipefail
source src/sim-study/hpc/environment.sh
srun Rscript --no-save --no-restore src/sim-study/hpc/collect.R "$1"
