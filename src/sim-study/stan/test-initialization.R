# Run from R/: Rscript --vanilla src/sim-study/stan/test-initialization.R
source("src/sim-study/02-worker-stan.R")
source("src/sim-study/02-worker-stan-parallel.R")
study <- load_sim_study_stan_parallel(chains = 4L, cores = 1L,
                                     control = list(adapt_delta = 0.8, max_treedepth = 6L))
close_to <- function(x, y) stopifnot(isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = 1e-10)))
prepared_cases <- list()
for (predictor in c("X1", "X2", "X3")) {
  scenario <- study$generate_data("D3", predictor, 250, 15)
  for (model in c("M0", "M2")) {
    prepared <- study$prepare_nils_stan(scenario$dat, predictor, model, seed = scenario$seed, w = 1)
    set.seed(819)
    old_seed <- .Random.seed
    old_kind <- RNGkind()
    starts <- study$make_nils_stan_inits(prepared, seed = scenario$seed)
    stopifnot(identical(old_seed, .Random.seed), identical(old_kind, RNGkind()),
              identical(starts, study$make_nils_stan_inits(prepared, seed = scenario$seed)),
              !identical(starts, study$make_nils_stan_inits(prepared, seed = scenario$seed + 1L)),
              identical(starts[1L], make_sim_nils_stan_inits(prepared, chains = 1L, seed = scenario$seed)))
    stopifnot(length(starts) == 4L, length(unique(lapply(starts, serialize, connection = NULL))) == 4L)
    for (start in starts) {
      stopifnot(identical(dim(start$beta), prepared$data$P),
                identical(dim(start$tau2), prepared$data$n_tau),
                all(start$beta > -0.5 & start$beta < 0.5),
                all(is.finite(start$tau2)), all(start$tau2 > 0))
    }
    stopifnot(identical(study$make_nils_stan_inits(prepared, seed = scenario$seed, method = "stan"), "random"))
    prepared_cases[[paste(predictor, model)]] <- list(prepared = prepared, starts = starts)
  }
}
# Check the distributions on the requested SD scale, not on the variance scale.
p <- prepared_cases[["X3 M0"]]$prepared
draws <- make_sim_nils_stan_inits(p, chains = 2000L, seed = 137)
tau <- vapply(draws, function(x) sqrt(x$tau2), numeric(p$data$n_tau))
stopifnot(abs(mean(tau[1:2, ]) - 1) < 0.07,
          abs(mean(tau[3:4, ]) - 0.5) < 0.04,
          abs(mean(tau[5, ]) - 1) < 0.09)

# The no-existing-RNG-state case must also be restored exactly.
saved <- .Random.seed
rm(.Random.seed, envir = .GlobalEnv)
invisible(make_sim_nils_stan_inits(p, chains = 4L, seed = 137))
stopifnot(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
assign(".Random.seed", saved, envir = .GlobalEnv)
cat("All six layouts, dispersed starts, rates, reproducibility and RNG isolation passed.\n")

sm <- study$compile_nils_stan()
# Supply all four distinct starts to Stan and verify the actual constrained starts.
for (case in prepared_cases) {
  fit <- rstan::sampling(sm, data = case$prepared$data, init = case$starts,
                         chains = 4L, cores = 1L, iter = 2L, warmup = 0L,
                         algorithm = "Fixed_param", seed = 71L, refresh = 0L)
  actual <- rstan::get_inits(fit)
  for (ch in seq_len(4L)) {
    close_to(actual[[ch]]$beta, case$starts[[ch]]$beta)
    close_to(actual[[ch]]$tau2, case$starts[[ch]]$tau2)
  }
}

# Verify Stan's own unconstrained initialization, including positive tau2.
p <- prepared_cases[["X3 M0"]]$prepared
fit <- rstan::sampling(sm, data = p$data, init = "random", chains = 4L, cores = 1L,
                       iter = 2L, warmup = 0L, algorithm = "Fixed_param", seed = 71L, refresh = 0L)
actual <- rstan::get_inits(fit)
stopifnot(length(unique(lapply(actual, serialize, connection = NULL))) == 4L)
for (start in actual) {
  stopifnot(all(abs(start$beta) < 2), all(abs(log(start$tau2)) < 2))
}
again <- rstan::sampling(sm, data = p$data, init = "random", chains = 4L, cores = 1L,
                         iter = 2L, warmup = 0L, algorithm = "Fixed_param", seed = 71L, refresh = 0L)
stopifnot(identical(actual, rstan::get_inits(again)))

# Actual adapter fits exercise the selectable methods and recorded metadata.
study$MCMC_TOTAL_ITER <- 24L
study$MCMC_BURNIN <- 12L
study$MCMC_THIN <- 3L
scenario <- study$generate_data("D1", "X1", 250, 1)
for (method in c("dispersed", "stan")) {
  study$STAN_INIT_METHOD <- method
  out <- study$fit_nils_gbi(scenario$dat, "X1", "M2", seed = scenario$seed)
  stopifnot(identical(out$initialization$method, method),
            out$initialization$seed == scenario$seed,
            length(out$initialization$actual) == 4L,
            length(unique(lapply(out$initialization$actual, serialize, connection = NULL))) == 4L,
            nrow(out$draws) == 16L, all(is.finite(out$draws)))
  if (method == "stan") stopifnot(identical(out$initialization$supplied, "random"))
  else for (ch in seq_len(4L)) {
    close_to(out$initialization$actual[[ch]]$beta, out$initialization$supplied[[ch]]$beta)
    close_to(out$initialization$actual[[ch]]$tau2, out$initialization$supplied[[ch]]$tau2)
  }
}
cat("Both methods passed actual Stan initialization and four-chain adapter checks.\n")
