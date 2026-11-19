# Cluster-only settings. The interactive worker defaults are not changed.
study_settings <- list(
  total_iter = 60000L, burnin = 10000L, thin = 10L, # 5,000 retained draws/chain
  # M1 defaults: k=10 -> 9 centered coefficients/effect.
  # DR uses 1 linear + D_SPLINE nonlinear coefficients, so D_SPLINE=8 matches.
  D_SPLINE = 8L, SF = 1,
  M5_CHAINS = 1L, # Additional chains run sequentially within the one-CPU task.
  M5_CONTROL = list(adapt_delta = .95, max_treedepth = 12L),
  M5_INIT_METHOD = "stan", # Stan's default random unconstrained initialization.
  OPTIM_WARN = FALSE,
  B_M3 = 100L, grid_N = 50L, probs = c(.025, .975),
  grid_chunk_size = 512L, kde_chunk_size = 256L,
  seed_base = 0, test_replicate = 101L, n_test = 500L,
  save_extracted = TRUE # Fits and evaluation metrics are always saved.
)
