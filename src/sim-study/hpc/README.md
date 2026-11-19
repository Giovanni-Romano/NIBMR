# Simulation study on bwUniCluster 3.0

Run all commands from the **R project directory**. No jobs are submitted by
creating these files. The launcher uses Slurm's `cpu` partition and one physical
CPU per task. `environment.sh` limits BLAS/OpenMP/Stan threads to one.

## Design and runtime

There are 2 × 4 × 3 × 100 = **2,400 datasets**, each fitted with M0–M5. The CSV
`runs.csv` contains n, DGP, predictor and replicate, one run per row. Regenerate
it with `Rscript src/sim-study/hpc/make-config.R` if needed. The first 24 rows
cover every setting; seeds depend on the scenario, not its row/order.

`settings.R` is separate from the interactive worker defaults. It specifies:

- 60,000 total iterations, 10,000 burn-in, thinning 10: **5,000 saved draws per
  chain**; one M5 chain on one CPU by default.
- D_SPLINE=8; 100 M3 bootstrap refits; 95% pointwise intervals.
- The agreed 50-point curves and 50×50 surfaces, including region-specific
  surfaces; one fixed independent test set of 500 responses per DGP.
- Original priors, loss policy and calibration. Optimizer warnings are muted;
  optimizer diagnostics and native sampler warnings remain available.

At 50 simultaneous tasks, 2,400/50 = 48 is the equal-duration batch equivalent.
Slurm fills a free slot immediately; it does not wait for 50 tasks to finish.
The reported Mac benchmark (3.48 minutes) was D3/X2/n=250, total=25,000,
burn-in=10,000, thin=25, 600 saved draws and the extractor's default 100-point
rather than 50-point grid. It cannot establish the runtime of the new schedule,
n=1000, X3, or cluster CPUs. Five minutes per run is therefore not a conservative
upper bound. Use the pilot below to measure this before committing resources.
The 50-run cap is easily increased after measuring memory and cluster throughput.

## 1. Install once on Linux

Copy the R project to a cluster workspace, including `.Rprofile`, `renv/activate.R`,
`src/`, and **`sim-study/tanzania.gra`**. Do not copy the Mac `renv/library/`,
compiled Stan binaries, or previous model caches. Use a workspace with enough
capacity for saved fits and regional evaluations.

Inspect `module avail math/R`, then edit/pin `R_MODULE` in `environment.sh` or
export it before submitting. The supplied default is `math/R/4.5.1`, listed by
the cluster's software catalogue. The local lock records R 4.5.2; the available
cluster module must be checked on site. Keep R, compiler and package versions
fixed throughout one study; exact draws across platforms are not guaranteed.

The original project lockfile omitted several study packages. The separate
`hpc/renv.lock` includes the currently installed study packages and dependencies;
it does not modify the original lockfile. Restore it **once before submission**:

```bash
module load math/R/4.5.1                 # use the version you selected
Rscript -e 'renv::restore(lockfile="src/sim-study/hpc/renv.lock", prompt=FALSE)'
```

RStan requires a working C++ toolchain and the compiled GUD dependency. Resolve
installation/compiler issues here, not in 2,400 array tasks. Rscript in the batch
files deliberately reads `.Rprofile` to activate the project library.

## 2. Prepare shared inputs in a compute job

Choose a **new directory for each study configuration**; paths below are examples:

```bash
STUDY="$PWD/sim-study/runs/study-5000"
mkdir -p "$STUDY/logs"
prep=$(sbatch --parsable --output="$STUDY/logs/prepare-%j.out" \
  --error="$STUDY/logs/prepare-%j.err" \
  src/sim-study/hpc/prepare.sbatch \
  "$PWD/src/sim-study/hpc/runs.csv" "$PWD/src/sim-study/hpc/settings.R" "$STUDY")
echo "$prep"
```

Preparation copies the exact CSV/settings into `manifest.rds` and `config.csv`,
saves the four shared test sets, records `sessionInfo()`, and compiles M5 **once**
to `stan-model.rds`. This is a cluster/Linux cache, shared read-only by tasks.
Compile and sample on the same partition/toolchain; avoid architecture-specific
`-march=native` flags if the binary will run on different CPU types.

Code fingerprints prevent an edited statistical routine from silently mixing
with earlier results. To change code or settings after preparation, create a
new study directory. Repeating preparation with identical inputs reuses the
existing test sets/model cache.

## 3. Pilot all 24 settings

After preparation finishes successfully:

```bash
head -n 25 "$STUDY/config.csv" > "$STUDY/pilot.csv"
bash src/sim-study/hpc/submit.sh "$STUDY" "$STUDY/pilot.csv" 8
```

The pilot has the **full requested chain/bootstrap schedule**, not development
settings. Its results are reused in the eventual full run. Start with eight
simultaneous tasks to inspect resource usage. Initial per-task requests are
4 GiB and one hour; these are configurable ceilings, not measured requirements.
Increase them in `array.sbatch` if the pilot needs more, or lower them once
measured. The collector requests 8 GiB and reads one 100-replicate group at a
time. Preparation requests 4 GiB/30 minutes. No GPU or exclusive node is needed.

Inspect logs and Slurm accounting (replace JOBID with the printed array ID):

```bash
squeue -u "$USER"
sacct -j JOBID --format=JobID,State,Elapsed,AllocCPUS,MaxRSS,ExitCode --units=G
```

Look at the `.0` srun steps for peak memory. `summary/timings.csv` records
fitting, extraction/bootstrap, structural evaluation, test scoring and save
stages. `total` excludes startup/library loading and the final result write;
Slurm Elapsed is the end-to-end task time. Scale observed total CPU work by
50 for a rough study time, then allow for queueing and slow tasks. The 24-run
pilot gives one observation per setting, so more replicates improve this estimate.
Check actual MCMC/NUTS diagnostics before interpreting the simulation results.

## 4. Submit the study

Wait for the pilot and its collector before submitting into the same directory:

```bash
bash src/sim-study/hpc/submit.sh "$STUDY" "$STUDY/config.csv" 50
```

Alternatively, when skipping the pilot, queue the array while preparation runs:

```bash
DEPENDENCY="${prep%%;*}" bash src/sim-study/hpc/submit.sh \
  "$STUDY" "$PWD/src/sim-study/hpc/runs.csv" 50
```

Supply the original CSV here because preparation may not yet have written its
saved copy. That first dependency uses `afterok`, so a failed preparation does not release
all tasks. Full submission uses `%50`: at most 50 simultaneous runs, each with
one allocated CPU. More M5 chains, if selected, run sequentially on that CPU.

The helper reads Slurm's `MaxArraySize`. If 2,400 indices are supported it uses
one array; otherwise it splits the CSV range into arrays and chains them. Each
chunk's index plus its saved offset selects the correct row. Chunks are chained
with `afterany`, so later chunks still run after a failed replicate while the
50-task cap holds across chunks. Chunk boundaries can add some waiting for
stragglers. `ARRAY_SIZE_LIMIT=1000` overrides the detected limit. User/account
limits on submitted jobs and resource availability still apply.

A collection job is submitted with an `afterany` dependency on the last array.
It runs even if some tasks failed. R exceptions are saved in `failures/` and
return a nonzero exit code. Slurm kills/timeouts may leave no failure file;
collection still reports the missing result. No task appends to a shared RDS
or CSV: each task writes unique files and completes them by atomic rename.

## Files, resumption and collection

For each key such as `n1000_D3_X3_r015`:

- `fits/KEY.rds`: full `run_cell()` output, including beta/tau, priors,
  calibration, settings and the native objects needed for re-evaluation.
- `extracted/KEY.rds`: point estimates and intervals, if `save_extracted=TRUE`.
- `results/KEY.rds`: structural errors/coverage/width, test scores, timings
  and execution information. Its presence marks a completed run.
- `failures/KEY.rds`: latest R error and stage, if an attempt failed. A later
  successful result takes precedence; the earlier failure remains for inspection.

A completed run is skipped. A fit checkpoint avoids refitting after an
extraction failure; a saved extraction avoids repeating the bootstrap after a
later evaluation failure. Identical settings/code are required. Fit/evaluation
objects can be large, especially for X3 and unthinned native Stan output.
Measure disk use during the pilot (`du -sh "$STUDY"/*`); `save_extracted=FALSE`
omits that extra copy but preserves fits and evaluation metrics. Decide before
preparation. Never delete checkpoint files needed for a planned resumption.

`summary/` contains:

- One `n_DGP_predictor.rds` per setting: existing `aggregate_evaluations()`
  output, including bias/MAD/RMSE, coverage/width, counts, and per-replicate
  training errors. Regional metrics use only their own training observations.
- `test-scores.csv`: one row per expected run/model; RMSE, MAE and log score.
- `test-summary.csv`: their means across replicates and valid metric counts.
- `inventory.csv`: expected runs with OK/MISSING/FAILED/UNREADABLE status.
- `retry.csv`: only incomplete runs, in the original four-column format.
- `timings.csv`: stage timings of completed runs.

Missing runs remain in failure counts and training rows; incomplete intervals
are excluded from coverage denominators, not counted as noncoverage. A completed
pipeline can still contain an unavailable M4 winner; its model failure is retained
in the structural counts and score status. Summaries are conditional on the
available estimates: inspect both run and model counts.

After the preceding collector has finished, retry only incomplete runs:

```bash
bash src/sim-study/hpc/submit.sh "$STUDY" "$STUDY/summary/retry.csv" 50
# Or collect manually inside an appropriate interactive allocation:
Rscript src/sim-study/hpc/collect.R "$STUDY"
```

Do not run duplicate arrays/collectors against the same directory concurrently.
For a run that failed in Slurm, check the Slurm log/accounting before retrying
and adjust resource requests if necessary.

## Local verification

From R/, using the project R library:

```bash
Rscript --vanilla src/sim-study/tests/test-hpc-pipeline.R
python3 src/sim-study/tests/test-hpc-submit.py
Rscript --vanilla src/sim-study/tests/test-hpc-engines.R
```

The real-engine check uses explicit short test settings in a temporary directory.
It checks serialization, cached model reuse, all six engines, bootstrap,
X3/n=250 and X2/n=1000, extraction, scoring and collection. It does not establish
convergence, HPC runtime, scheduler access or the suitability of memory limits.

References: [bwUniCluster jobs/partitions](https://wiki.bwhpc.de/e/BwUniCluster3.0/Running_Jobs),
[R environment](https://wiki.bwhpc.de/e/BwUniCluster3.0/Software/R),
[software catalogue](https://www.bwhpc.de/software.php),
[Slurm arrays and dependencies](https://slurm.schedmd.com/job_array.html).
