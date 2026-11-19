#!/bin/bash
# From R/: bash src/sim-study/hpc/submit.sh STUDY_DIR [CONFIG.csv] [CONCURRENCY]
# Optional: DEPENDENCY=<prepare-job-id> ARRAY_SIZE_LIMIT=1000 bash ...
set -euo pipefail
study_dir=$(cd "$1" && pwd)
config=${2:-"$study_dir/config.csv"}
config_dir=$(cd "$(dirname "$config")" && pwd)
config="$config_dir/$(basename "$config")"
concurrency=${3:-50}
rows=$(awk 'END {print NR-1}' "$config")
if (( rows < 1 )); then echo "No runs in $config" >&2; exit 1; fi
mkdir -p "$study_dir/logs"
# MaxArraySize limits the highest index, not just the number of tasks.
# Some sites retain a 1001 limit; split and chain arrays there. Chaining
# preserves the 50-run cap across ALL submitted chunks, even after failures.
array_size=${ARRAY_SIZE_LIMIT:-$(scontrol show config | awk '$1 == "MaxArraySize" {print $3-1}')}
array_size=${array_size:-1000}
if (( array_size < 1 || concurrency < 1 )); then echo "Invalid array/concurrency limit" >&2; exit 1; fi
dependency=${DEPENDENCY:-}
dependency_type=afterok # A preparation failure must not release the array.
for ((offset=0; offset<rows; offset+=array_size)); do
  count=$((rows-offset)); if ((count>array_size)); then count=$array_size; fi
  args=(--parsable --array="1-$count%$concurrency" --chdir="$PWD"
        --output="$study_dir/logs/%A_%a.out" --error="$study_dir/logs/%A_%a.err")
  if [[ -n "$dependency" ]]; then args+=(--dependency="$dependency_type:$dependency"); fi
  job=$(sbatch "${args[@]}" src/sim-study/hpc/array.sbatch "$study_dir" "$config" "$offset")
  dependency=${job%%;*}
  dependency_type=afterany
  echo "Submitted $dependency: config rows $((offset+1))-$((offset+count)), at most $concurrency simultaneous runs"
done
# Collect even when some tasks failed; inventory.csv and retry.csv show gaps.
collect_args=(--parsable --chdir="$PWD" --dependency="afterany:$dependency"
              --output="$study_dir/logs/collect-%j.out" --error="$study_dir/logs/collect-%j.err")
job=$(sbatch "${collect_args[@]}" src/sim-study/hpc/collect.sbatch "$study_dir")
echo "Collection job: ${job%%;*}"
