#!/bin/bash -l
#SBATCH --job-name=mode-study
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem=4G
#SBATCH --time=01:00:00
set -euo pipefail
# Arguments: study-directory config.csv row-offset
source src/sim-study/hpc/environment.sh
row=$((SLURM_ARRAY_TASK_ID + ${3:-0}))
srun Rscript --no-save --no-restore src/sim-study/hpc/run-task.R "$1" "$2" "$row"
