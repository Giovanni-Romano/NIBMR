# Parallel-safe compilation for the existing simulation Stan adapter.
# From R/: source("src/sim-study/02-worker-stan-parallel.R")
#          enable_sim_study_stan_parallel(study)

enable_sim_study_stan_parallel <- function(study, project_dir = ".") {
  if (!is.environment(study) || !is.function(study$compile_nils_stan) ||
      !is.function(study$fit_nils_gbi))
    stop("study must be the environment returned by load_sim_study_stan().")
  model_file <- normalizePath(file.path(project_dir, "src/sim-study/stan/nils-gbi.stan"),
                              mustWork = TRUE)
  cache <- new.env(parent = emptyenv())
  cache$model <- NULL
  # Keep the explicit compiler entry point for earlier parallel experiments.
  # Both the production compiler and this one retain the DSO for PSOCK.
  study$compile_nils_stan <- function() {
    if (is.null(cache$model)) {
      cache$model <- rstan::stan_model(file = model_file, auto_write = FALSE,
                                       save_dso = TRUE)
    }
    cache$model
  }
  invisible(study)
}

# Optional entry point for a new study, with the original loader's arguments.
# The source adapter is loaded privately, leaving global worker functions alone.
load_sim_study_stan_parallel <- function(project_dir = ".", benchmarks = FALSE,
                                         chains = 1L, cores = chains,
                                         control = list(adapt_delta = 0.95,
                                                        max_treedepth = 12L),
                                         init_method = c("dispersed", "stan")) {
  init_method <- match.arg(init_method)
  adapter <- new.env(parent = parent.frame())
  sys.source(file.path(project_dir, "src/sim-study/02-worker-stan.R"), envir = adapter)
  study <- adapter$load_sim_study_stan(project_dir = project_dir,
                                      benchmarks = benchmarks, chains = chains,
                                      cores = cores, control = control,
                                      init_method = init_method)
  enable_sim_study_stan_parallel(study, project_dir = project_dir)
  study
}
