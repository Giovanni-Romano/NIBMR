# Four-chain compilation fix

The earlier simulation adapter compiled with `save_dso = FALSE`. Such a model
can run in its original R process, but cannot be reloaded from its serialized
object in a separate R process. RStan uses such worker processes in some
parallel execution paths, including the one that reports `checkForRemoteErrors`.

Both current compilers now use **`save_dso = TRUE`**. Older in-memory models
need recompilation. This retains the compiled
binary inside the model object so worker processes can load it. Keeping
`auto_write = FALSE` disables the adjacent `.rds` model cache; retaining the
binary does not require writing a cache beside the Stan source.

From the `R` project directory, fix an existing `study` created with
`chains = 4, cores = 4`:

```r
source("src/sim-study/02-worker-stan-parallel.R")
enable_sim_study_stan_parallel(study)

# Rerun your fit. The first call recompiles; subsequent calls reuse that model.
fit <- study$fit_nils_gbi(scenario$dat, "X3", "M0", seed = scenario$seed)
```

This preserves the existing study's chain/core settings, warmup, iteration
schedule, priors, and calibration code. It replaces only its compiler and
compilation cache. Existing files are not edited. The Stan model is unchanged.
Simply changing an RStan option does not repair the previously compiled object.

For a fresh study, use the new entry point:

```r
source("src/sim-study/02-worker-stan-parallel.R")
study <- load_sim_study_stan_parallel(chains = 4, cores = 4,
                                      init_method = "dispersed")
study$MCMC_TOTAL_ITER <- 11000L
study$MCMC_BURNIN <- 1000L
study$MCMC_THIN <- 10L
s <- study$run_nils_cell("D3", "X3", 250, 15)
```

The updated loader accepts `init_method = "dispersed"` (the default) or
`init_method = "stan"`; see README.md for distributions and seed handling.
Initialization can also be switched with `study$STAN_INIT_METHOD` before a fit.

Four chains retain four times the per-chain draw count, as in the existing
adapter. Set the iteration schedule before starting your fit. For a temporary
workaround with the original loader, `chains = 4, cores = 1` runs the chains
sequentially in one R process, but does not repair the compiled object for
parallel workers or later serialization.

Run `Rscript --vanilla src/sim-study/stan/test-parallel.R` to check compilation,
transfer to four independent socket workers, all four chain IDs, and the
adapter's combined draw layout. The short test samples validate execution and
data layout, not convergence.
