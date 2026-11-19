# M5 production backend, plus the private M0/M2 Stan experiment loader.
# Source definitions only; the Stan model is compiled on the first fit.

# Generate all custom starts in the parent process, before launching workers.
# Each chain gets a distinct L'Ecuyer stream. Reusing the same model seed is
# reproducible; adding chains preserves the starts of the earlier chains.
make_sim_nils_stan_inits <- function(prepared, chains, seed = NULL,
                                     method = c("dispersed", "stan")) {
  method <- match.arg(method)
  if (!is.numeric(chains) || length(chains) != 1L || !is.finite(chains) ||
      chains < 1 || chains != floor(chains)) stop("chains must be a positive integer.")
  if (method == "stan") return("random")
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1L)
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed < 0 || seed > .Machine$integer.max || seed != floor(seed))
    stop("The initialization seed must be an integer between 0 and .Machine$integer.max.")

  # Initialization must not change the R stream used for CV or other models.
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(list = ".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed, kind = "L'Ecuyer-CMRG", normal.kind = "Inversion", sample.kind = "Rejection")
  stream <- get(".Random.seed", envir = .GlobalEnv)
  dat <- prepared$data
  n_linear <- prepared$design$q + prepared$design$p
  linear <- seq_len(n_linear)
  nonlinear <- if (dat$n_ig > n_linear) seq.int(n_linear + 1L, dat$n_ig) else integer(0)
  spatial <- if (dat$s > 0L) dat$n_tau else integer(0)
  starts <- vector("list", chains)
  for (ch in seq_len(chains)) {
    assign(".Random.seed", stream, envir = .GlobalEnv)
    beta <- runif(dat$P, min = -0.5, max = 0.5) # Includes intercept and MMR beta.
    tau <- numeric(dat$n_tau)
    tau[linear] <- rexp(length(linear), rate = 1)
    tau[nonlinear] <- rexp(length(nonlinear), rate = 2)
    tau[spatial] <- rexp(length(spatial), rate = 1)
    # Exclude a possible floating-point boundary draw; Stan requires tau2 > 0.
    while (any(tau <= 0)) {
      bad <- which(tau <= 0)
      tau[bad] <- rexp(length(bad), rate = ifelse(bad %in% nonlinear, 2, 1))
    }
    starts[[ch]] <- list(beta = array(beta, dim = dat$P),
                         tau2 = array(tau^2, dim = dat$n_tau))
    stream <- parallel::nextRNGStream(stream)
  }
  starts
}

# The production M5 fit receives M0's stored inputs. Standalone Stan experiments
# use the same shared setup, prior recipe and calibration routine as M0/M2.
prepare_sim_nils_stan <- function(worker, dat, predictor, model_type, seed = NULL,
                                  w = NULL, design = NULL, reference_fit = NULL) {
  predictor <- match.arg(predictor, c("X1", "X2", "X3"))
  model_type <- match.arg(model_type, c("M0", "M2"))
  dd <- if (is.null(design)) worker$make_design_M0(dat, predictor,
    splines = model_type == "M0") else design
  if (dd$s > 0L && !isTRUE(all.equal(unname(dd$K_GMRF), diag(1, dd$s),
                                    tolerance = 1e-10, check.attributes = FALSE)))
    stop("The simulation Stan model requires the identity MMR coefficient penalty.")
  ntau <- dd$q + dd$p + (model_type == "M0") * dd$p + as.integer(dd$s > 0)
  n_ig <- ntau - as.integer(dd$s > 0)
  tau_names <- c(paste0(colnames(dd$X_cont), ".linear"),
                 if (model_type == "M0") paste0(colnames(dd$X_cont), ".nonlinear"),
                 if (dd$s > 0L) "GMRF")
  if (is.null(reference_fit)) {
    setup <- worker$nils_setup(dd$X, dat$y)
    prior <- worker$nils_prior(dd, splines = model_type == "M0")
    loss <- list(k = setup$k, c = setup$c, s = setup$s)
    calibration <- if (is.null(w)) worker$calibrate_nils(
      dd$X, dat$y, unname(setup$start), setup$k, setup$c, K_fold = 10L, seed = seed,
      optim_warn = worker$OPTIM_WARN, optim_retry = worker$OPTIM_RETRY,
      optim_retry_factor = worker$OPTIM_RETRY_FACTOR,
      cv_optim_control = worker$CV_OPTIM_CONTROL, optim_ntry_cv = worker$OPTIM_NTRY_CV,
      cv_analytic_gradient = worker$CV_ANALYTIC_GRADIENT) else
        list(w = w, numerator = ncol(dd$X) / 2, mean_cv = NA_real_)
    intercept_prior <- prior$intercept
    tau_prior <- list(a_tau = setNames(prior$a_tau[seq_len(n_ig)], tau_names[seq_len(n_ig)]),
                      b_tau = setNames(prior$b_tau[seq_len(n_ig)], tau_names[seq_len(n_ig)]),
                      theta_GMRF = dd$theta_GMRF)
  } else {
    # Reuse rather than recalibrate: a different optimizer, fold allocation or
    # numerator would change the posterior and confound the sampler comparison.
    loss <- reference_fit$loss
    calibration <- reference_fit$calibration
    intercept_prior <- reference_fit$intercept_prior
    tau_prior <- reference_fit$tau_prior
    if (!is.null(w) && !isTRUE(all.equal(w, calibration$w)))
      stop("w must equal the reference M0 calibration.")
  }
  tau_idx <- c(seq_len(dd$q + dd$p),
    if (model_type == "M0") rep(dd$q + dd$p + seq_len(dd$p), times = dd$d),
    if (dd$s > 0L) rep(ntau, dd$s))
  coefficient_names <- c("Intercept", "x1", "x2",
    if (model_type == "M0") c(paste0("x1.s", seq_len(dd$d[1L])),
                              paste0("x2.s", seq_len(dd$d[2L]))),
    if (dd$s > 0L) paste0("gmrf", seq_len(dd$s)))
  stan_data <- list(n = nrow(dat), P = ncol(dd$X), X = dd$X, y = dat$y,
    k = as.vector(loss$k), c = loss$c, w = calibration$w, n_tau = ntau, n_ig = n_ig,
    tau_idx = as.integer(tau_idx), a_tau = array(tau_prior$a_tau, dim = n_ig),
    b_tau = array(tau_prior$b_tau, dim = n_ig), s = dd$s,
    theta_GMRF = if (dd$s > 0L) tau_prior$theta_GMRF else 1,
    intercept_mean = intercept_prior$mean, intercept_sd = intercept_prior$sd)
  list(data = stan_data, design = dd, calibration = calibration,
       fold_id = calibration$fold_id, coefficient_names = coefficient_names,
       tau_names = tau_names, loss = loss, intercept_prior = intercept_prior,
       tau_prior = tau_prior, calibration_reused = !is.null(reference_fit))
}

# Compile once per R session; retain the DSO for existing PSOCK experiments.
compile_nils_stan <- local({
  model <- NULL
  function() {
    if (is.null(model)) model <<- rstan::stan_model(
      file = "src/sim-study/stan/nils-gbi.stan", auto_write = FALSE, save_dso = TRUE)
    model
  }
})

sample_sim_nils_stan <- function(worker, prepared, seed, chains, cores, control,
                                 init_method) {
  schedule <- worker$mcmc_schedule(worker$MCMC_TOTAL_ITER, worker$MCMC_BURNIN, worker$MCMC_THIN)
  sampling_seed <- if (is.null(seed)) sample.int(.Machine$integer.max, 1L) else seed
  initial_values <- make_sim_nils_stan_inits(prepared, chains, sampling_seed, init_method)
  fit <- rstan::sampling(worker$compile_nils_stan(), data = prepared$data,
    chains = chains, cores = cores, iter = schedule$total_iter, warmup = schedule$burnin,
    thin = 1L, control = control, init = initial_values, seed = sampling_seed, refresh = 0)
  # Match R's retained positions: burnin+1, burnin+1+thin, ... . Stan's own
  # thinning starts at another offset. Pool chains only after selecting rows.
  keep <- schedule$keep - schedule$burnin
  extract_draws <- function(parameter, names, axis) {
    raw <- rstan::extract(fit, pars = parameter, permuted = FALSE, inc_warmup = FALSE)
    keys <- paste0(parameter, "[", seq_along(names), "]")
    out <- do.call(rbind, lapply(seq_len(chains), function(ch)
      matrix(raw[keep, ch, keys, drop = FALSE], nrow = length(keep), ncol = length(names))))
    worker$assert_mcmc_draws(out, chains * schedule$ndraws, paste("M5", parameter))
    dimnames(out) <- setNames(list(seq_len(nrow(out)), names), c("iteration", axis))
    out
  }
  draws <- extract_draws("beta", prepared$coefficient_names, "variable")
  tau <- extract_draws("tau", prepared$tau_names, "parameter")
  list(draws = draws, tau = tau,
       fhat = as.vector(prepared$design$X %*% colMeans(draws)),
       calibration = prepared$calibration, acc = NULL, loss = prepared$loss,
       intercept_prior = prepared$intercept_prior, tau_prior = prepared$tau_prior,
       fit = fit, backend = "stan", chain = rep(seq_len(chains), each = schedule$ndraws),
       retained_iteration = rep(schedule$keep, times = chains), schedule = schedule,
       calibration_reused = prepared$calibration_reused,
       stan_settings = list(chains = chains, cores = cores, control = control, init_method = init_method),
       initialization = list(method = init_method, seed = sampling_seed,
                             supplied = initial_values, actual = rstan::get_inits(fit)),
       sampler_diagnostics = rstan::get_sampler_params(fit, inc_warmup = FALSE))
}

# M5 is M0 sampled with Stan/NUTS. run_cell supplies M0's design and target.
fit_nils_stan <- function(dat, predictor, seed = NULL, design = NULL, reference_fit = NULL) {
  worker <- environment(sys.function())
  prepared <- prepare_sim_nils_stan(worker, dat, predictor, "M0", seed,
                                   design = design, reference_fit = reference_fit)
  sample_sim_nils_stan(worker, prepared, seed, worker$M5_CHAINS, worker$M5_CORES,
                      worker$M5_CONTROL, worker$M5_INIT_METHOD)
}

load_sim_study_stan <- function(project_dir = ".", benchmarks = FALSE,
                                chains = 1L, cores = chains,
                                control = list(adapt_delta = 0.95, max_treedepth = 12L),
                                init_method = c("dispersed", "stan")) {
  init_method <- match.arg(init_method)
  project_dir <- normalizePath(project_dir, mustWork = TRUE)
  worker_file <- file.path(project_dir, "src/sim-study/02-worker.R")
  if (!file.exists(worker_file)) stop("project_dir must point to the R project directory.")
  for (nm in c("chains", "cores")) {
    x <- get(nm)
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 1 || x != floor(x))
      stop(nm, " must be a positive integer.")
  }
  for (pkg in c("rstan", "mgcv", "geigen", "quantreg", "R2BayesX")) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Required package is missing: ", pkg)
  }
  worker <- new.env(parent = parent.frame())
  old_dir <- getwd()
  on.exit(setwd(old_dir), add = TRUE)
  setwd(project_dir)
  # Load the existing definitions/configuration without invoking the original
  # if (!interactive()) main(). With benchmarks=FALSE, defer benchmark libraries.
  load_definitions <- function(path) {
    for (expr in parse(path)) {
      if (identical(expr, quote(if (!interactive()) main())) ||
          identical(expr, quote(if (sys.nframe() == 0L) main()))) next
      if (is.call(expr) && as.character(expr[[1L]]) %in% c("library", "require")) {
        if (benchmarks) eval(expr, worker)
      } else if (is.call(expr) && identical(expr[[1L]], as.name("source"))) {
        load_definitions(as.character(expr[[2L]]))
      } else {
        eval(expr, worker)
      }
    }
  }
  load_definitions(worker_file)
  # The original spline helpers use these mgcv functions without a namespace.
  worker$s <- mgcv::s
  worker$smoothCon <- mgcv::smoothCon
  worker$PredictMat <- mgcv::PredictMat
  worker$Predict.matrix <- mgcv::Predict.matrix
  sys.source(file.path(project_dir, "src/stan/prepare-data.R"), envir = worker)
  worker$fit_nils_gbi_iwls <- worker$fit_nils_gbi
  prepare <- prepare_sim_nils_stan
  worker$prepare_nils_stan <- function(dat, predictor, model_type, seed = NULL, w = NULL,
                                           design = NULL, reference_fit = NULL)
    prepare(worker, dat, predictor, model_type, seed, w, design, reference_fit)
  make_inits <- make_sim_nils_stan_inits
  worker$STAN_INIT_METHOD <- init_method
  worker$make_nils_stan_inits <- function(prepared, seed = NULL,
                                         method = worker$STAN_INIT_METHOD)
    make_inits(prepared, chains = chains, seed = seed, method = method)

  model_file <- file.path(project_dir, "src/sim-study/stan/nils-gbi.stan")
  cache <- new.env(parent = emptyenv())
  cache$model <- NULL
  worker$compile_nils_stan <- function() {
    if (is.null(cache$model))
      cache$model <- rstan::stan_model(file = model_file, auto_write = FALSE, save_dso = TRUE)
    cache$model
  }
  worker$fit_nils_gbi <- function(dat, predictor, model_type, seed = NULL, design = NULL) {
    prepared <- worker$prepare_nils_stan(dat, predictor, model_type, seed, design = design)
    sample_sim_nils_stan(worker, prepared, seed, chains, cores, control, worker$STAN_INIT_METHOD)
  }
  worker$run_nils_cell <- function(dgp, predictor, n, replicate, seed_base = 0,
                                  models = c("M0", "M2")) {
    if (!length(models) || anyDuplicated(models) || any(!models %in% c("M0", "M2")))
      stop("models must contain M0, M2, or both, without duplicates.")
    sim <- worker$generate_data(dgp, predictor, n, replicate, seed_base = seed_base)
    designs <- lapply(models, function(m) worker$make_design_M0(sim$dat, predictor, m == "M0"))
    fits <- lapply(models, function(m)
      worker$fit_nils_gbi(sim$dat, predictor, m, seed = sim$seed + if (m == "M2") 2000 else 0))
    names(designs) <- names(fits) <- tolower(models)
    list(scenario = sim, designs = designs, fits = fits)
  }
  worker
}
