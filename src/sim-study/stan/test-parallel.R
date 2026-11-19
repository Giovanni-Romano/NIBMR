# Run from R/: Rscript --vanilla src/sim-study/stan/test-parallel.R
source("src/sim-study/02-worker-stan-parallel.R")
study <- load_sim_study_stan_parallel(chains = 4L, cores = 4L,
                                     control = list(adapt_delta = 0.8,
                                                    max_treedepth = 6L))
study$MCMC_TOTAL_ITER <- 24L
study$MCMC_BURNIN <- 12L
study$MCMC_THIN <- 3L
scenario <- study$generate_data("D1", "X1", 250, replicate = 1)
prepared <- study$prepare_nils_stan(scenario$dat, "X1", "M2", seed = scenario$seed)
prepared$inits <- study$make_nils_stan_inits(prepared, seed = scenario$seed)
sm <- study$compile_nils_stan()
stopifnot(sm@dso@dso_saved, identical(sm, study$compile_nils_stan()))

check_socket_workers <- function(sm, data, inits) {
  # PSOCK explicitly serializes the model to new R sessions. Unix fork-only
  # tests can miss the save_dso=FALSE defect seen in interactive RStan sessions.
  cluster <- parallel::makePSOCKcluster(4L)
  on.exit(parallel::stopCluster(cluster), add = TRUE)
  paths <- .libPaths()
  parallel::clusterCall(cluster, function(paths) {
    .libPaths(paths)
    suppressPackageStartupMessages(library(rstan))
    NULL
  }, paths)
  parallel::clusterExport(cluster, c("sm", "data", "inits"), envir = environment())
  result <- parallel::parLapply(cluster, seq_len(4L), function(chain) {
    fit <- rstan::sampling(sm, data = data, chains = 1L, cores = 1L,
                           iter = 2L, warmup = 0L, algorithm = "Fixed_param",
                           chain_id = chain, init = inits[chain], seed = 101L,
                           refresh = 0L)
    draws <- rstan::extract(fit, pars = "beta", permuted = FALSE)
    list(chain = fit@stan_args[[1L]]$chain_id, dimensions = dim(draws),
         finite = all(is.finite(draws)), binary_saved = fit@stanmodel@dso@dso_saved,
         actual = rstan::get_inits(fit)[[1L]])
  })
  stopifnot(identical(as.integer(vapply(result, `[[`, numeric(1), "chain")), seq_len(4L)),
            all(vapply(result, function(x) x$finite && x$binary_saved &&
                         identical(x$dimensions, c(2L, 1L, data$P)), logical(1))))
  for (chain in seq_len(4L)) {
    for (parameter in c("beta", "tau2"))
      stopifnot(isTRUE(all.equal(as.numeric(result[[chain]]$actual[[parameter]]),
                                 as.numeric(inits[[chain]][[parameter]]))))
  }
}
check_socket_workers(sm, prepared$data, prepared$inits)
cat("Four independent socket workers used their distinct supplied starts.\n")

# Exercise the actual public adapter's four-chain sampling and assembly.
fit <- study$fit_nils_gbi(scenario$dat, "X1", "M2", seed = scenario$seed)
raw <- rstan::extract(fit$fit, pars = "beta", permuted = FALSE, inc_warmup = FALSE)
stopifnot(fit$fit@sim$chains == 4L, identical(dim(raw), c(12L, 4L, 3L)),
          nrow(fit$draws) == 16L, length(fit$sampler_diagnostics) == 4L,
          identical(fit$chain, rep(seq_len(4L), each = 4L)),
          all(is.finite(fit$draws)))
for (chain in seq_len(4L)) {
  stopifnot(isTRUE(all.equal(unname(fit$draws[fit$chain == chain, , drop = FALSE]),
                             unname(raw[seq.int(1L, 12L, 3L), chain, ]))))
  for (parameter in c("beta", "tau2"))
    stopifnot(isTRUE(all.equal(as.numeric(fit$initialization$actual[[chain]][[parameter]]),
                               as.numeric(prepared$inits[[chain]][[parameter]]))))
}
stopifnot(isTRUE(all.equal(fit$fhat, as.vector(prepared$design$X %*% colMeans(fit$draws)))))
cat("Four-chain adapter sampling, draw counts, chain order and thinning passed.\n")
