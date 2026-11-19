# M5: M0 sampled with Stan

`run_cell()` fits M0–M5 by default. M5 uses Stan/NUTS and receives M0's exact
stored design, loss parameters, calibrated weight and priors. It does not run
another calibration. This comparison changes the sampler, while keeping the
posterior target fixed. `include_m5 = FALSE` keeps the original five-model run.

From the R project directory, in an interactive session:

```r
source("src/sim-study/02-worker.R")
source("src/sim-study/03-extract.R")
source("src/sim-study/07-test-evaluate.R")

# Choose an adequate warmup before an inferential run. The unchanged default
# development profile has only 10 burn-in iterations and is not adequate here.
MCMC_TOTAL_ITER <- 11000L
MCMC_BURNIN <- 1000L
MCMC_THIN <- 10L
M5_CHAINS <- 4L
M5_CORES <- 4L

s <- run_cell("D3", "X3", 250, 15)
training <- results_on_train(s)
tests <- make_test_sets()
unseen <- evaluate_test_sample(s, tests$D3)
rstan::check_hmc_diagnostics(s$fits$m5$fit)
print(s$fits$m5$fit, pars = c("beta", "tau"))
```

The existing worker CLI also includes M5. Evaluation and study-wide collection
remain separate. For just M5 after fitting M0:

```r
m5 <- fit_nils_stan(s$scenario$dat, s$scenario$predictor,
                    seed = s$scenario$seed + 7000,
                    design = s$designs$m0, reference_fit = s$fits$m0)
```

## Target and numerical checks

`nils-gbi.stan` implements the simulation model: identity link, two linear
continuous effects, nonlinear DR blocks, and the optional X3 MMR spatial block.
There are no application random intercepts. The target uses the g=1 nested
interval loss. Its near-zero evaluation now uses the same rationalization and
`expm1` as R; the mathematical loss is unchanged.

The proper Gaussian intercept prior and effective IG priors on linear and
nonlinear variances are copied from M0. Calibration uses `(P/2)/mean_cv` and
the shared BFGS routine, including its gradients, controls, folds and retries.
The old standalone adapter's `1/2` numerator and separate L-BFGS-B calibration
were inconsistent with M0 and have been replaced.

For X3, `V_tilde = Gamma_pen %*% diag(1/sqrt(lambda))`, with
`t(V_tilde) %*% K %*% V_tilde = I`. MMR coefficients therefore have independent
Gaussian priors in these coordinates. Their regional effects remain spatially
structured after multiplication by the basis. The spatial variance has
`weibull(0.5, theta_GMRF)`, matching R's SDP. Stan's automatic positive-parameter
transformation includes the log-variance Jacobian; no additional Jacobian term
belongs in the Stan model. The current debug value `theta_GMRF = 1` is preserved.

## Saved draws and settings

M5 returns `draws`, `tau`, `fhat`, `calibration`, `loss`, `intercept_prior` and
`tau_prior` in M0's formats. `tau` means standard deviation, not variance.
Training/grid effects, slope and regional intervals, structural evaluation and
unseen-data predictions accept `m5`. Its B2 density score uses the same residual
KDE procedure as M0. `acc` is NULL: IWLS block acceptance has no NUTS equivalent.

Native `fit`, `sampler_diagnostics`, `chain`, `retained_iteration`, `schedule`,
`stan_settings` and `initialization` are retained. Diagnostics contain all
post-warmup transitions; `draws`/`tau` are thinned. Retained positions match R:
`burnin+1, burnin+1+thin, ...`. Rows are concatenated by chain. Each chain uses
the full schedule, so four chains produce four times the single-chain draw
count. Defaults are one chain, one core, `adapt_delta=.95`, `max_treedepth=12`.
Always assess convergence and NUTS diagnostics before interpreting results.

The M5 sampling seed is `scenario$seed + 7000`, within the existing 10,000-seed
block. M0's calibration seed remains unchanged. Compilation is lazy and reused
within the session, with no compiled cache written beside source files.

Default starts are independently dispersed by chain: beta is Uniform(-.5,.5),
linear/spatial tau is Exp(1), and nonlinear tau is Exp(2). These are starting
distributions, not priors. Separate L'Ecuyer streams preserve reproducibility
and restore the ambient R RNG state. `M5_INIT_METHOD <- "stan"` selects Stan's
own random unconstrained initialization instead.

## Existing private Stan experiments

`load_sim_study_stan()` remains available for the earlier M0/M2-only
experiments. It loads a private worker and replaces its M0/M2 engines with
Stan; use `study$run_nils_cell()` there. It is separate from the production
M0-versus-M5 comparison above. Its `chains`, `cores`, `control` and
`STAN_INIT_METHOD` settings still apply to those private fits.
`load_sim_study_stan_parallel()` retains its PSOCK-compatible compiler.

## Verification

Run from R/ using the project's R library:

```sh
Rscript --vanilla src/sim-study/tests/test-m5-stan.R
Rscript --vanilla src/sim-study/tests/test-m5-preservation.R
Rscript --vanilla src/sim-study/tests/test-m5-engines.R
Rscript --vanilla src/sim-study/stan/test-initialization.R
```

Checks cover all X1–X3/M0–M2 layouts, coefficient and IG/SDP conditional-density
changes, transformed-target gradients, tiny residuals, real NUTS sampling,
beta/tau ordering, thinning and pooling, target reuse, extraction, evaluation
and seeds. The preservation test accepts a pre-change baseline directory for
exact R fit/RNG comparisons. The short NUTS tests can emit divergence,
tree-depth and mixing warnings; they test the implementation, not convergence.
