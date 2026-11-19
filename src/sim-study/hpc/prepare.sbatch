#!/bin/bash -l
#SBATCH --job-name=mode-prepare
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem=4G
#SBATCH --time=00:30:00
set -euo pipefail
# Submit from R/. Arguments: config.csv settings.R study-directory
source src/sim-study/hpc/environment.sh
srun Rscript --no-save --no-restore src/sim-study/hpc/prepare.R "$1" "$2" "$3"
