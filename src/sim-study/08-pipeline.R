# File-based study orchestration. Source from R/, or use the hpc entry scripts.
# Existing fitting/statistical routines remain the source of every estimate.
source("src/sim-study/00-progress.R")
source("src/sim-study/00-summaries.R")

study_config <- function() {
  # Each block of 24 rows covers every setting. Start with n=1000 / X3 so a
  # pilot sees the expensive cases early; seeds never depend on task order.
  out <- expand.grid(dgp = paste0("D", 1:4), predictor = c("X3", "X2", "X1"),
                     n = c(1000L, 250L), replicate = 1:100,
                     KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  out <- out[order(out$replicate, -out$n, match(out$predictor, c("X3", "X2", "X1")), out$dgp), ]
  rownames(out) <- NULL
  out[c("n", "dgp", "predictor", "replicate")]
}

study_key <- function(row) with(row, sprintf("n%d_%s_%s_r%03d", n, dgp, predictor, replicate))

# Save beside the final file, then rename: interrupted writes never look complete.
study_save <- function(object, path) {
  tmp <- paste0(path, ".", Sys.getpid(), ".tmp")
  saveRDS(object, tmp)
  if (!file.rename(tmp, path)) stop("Could not finish saving ", path)
  invisible(path)
}

# A broken checkpoint is unavailable and can be recomputed on retry.
study_read <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(readRDS(path), error = function(e) NULL)
}

study_sources <- function() {
  files <- c("src/01a-utils.R", "src/01b-mcmc.R", "sim-study/tanzania.gra",
    file.path("src/sim-study", c("00-progress.R", "00-summaries.R", "01-dgp.R", "02-worker.R", "02-worker-stan.R",
      "03-extract.R", "06-evaluate.R", "07-test-evaluate.R", "08-pipeline.R",
      "stan/nils-gbi.stan")))
  tools::md5sum(files)
}

# Reuse the live object when preparation and sampling share an R session.
# Re-reading its serialized external pointers while the same DLL is already
# loaded can confuse RStan. Fresh Slurm processes load the saved binary once.
study_model <- local({
  cache <- new.env(parent = emptyenv())
  function(path, model = NULL) {
    path <- normalizePath(path)
    if (!is.null(model)) cache[[path]] <- model
    if (is.null(cache[[path]])) cache[[path]] <- readRDS(path)
    cache[[path]]
  }
})

study_worker <- function(settings, model_file = NULL) {
  e <- new.env(parent = globalenv())
  # Routine package/graph setup is quiet; startup warnings still reach stderr.
  run_output(suppressPackageStartupMessages({
    source("src/sim-study/02-worker.R", local = e)
    for (file in c("03-extract.R", "06-evaluate.R", "07-test-evaluate.R"))
      sys.source(file.path("src/sim-study", file), e)
  }))
  e$MCMC_PROFILE <- "cluster"
  e$MCMC_TOTAL_ITER <- settings$total_iter; e$MCMC_BURNIN <- settings$burnin
  e$MCMC_THIN <- settings$thin
  e$MCMC_SETTINGS <- e$mcmc_schedule(settings$total_iter, settings$burnin, settings$thin)
  e$MCMC_KEEP <- e$MCMC_SETTINGS$keep; e$MCMC_NDRAWS <- e$MCMC_SETTINGS$ndraws
  fitting <- c("D_SPLINE", "SF", "M5_CHAINS", "M5_CONTROL", "M5_INIT_METHOD",
               names(settings)[grepl("^(M3_|OPTIM_|CV_)", names(settings))])
  for (name in unique(fitting))
    e[[name]] <- settings[[name]]
  e$M5_CORES <- 1L
  if (!is.null(settings$parameter_probs)) e$POSTERIOR_PROBS <- settings$parameter_probs
  options(mc.cores = 1L)
  if (!is.null(model_file)) {
    # The compiled model is built once on the cluster, not once per replicate.
    # save_dso=TRUE embeds its Linux binary. Never copy a compiled Mac model here.
    model <- NULL
    e$compile_nils_stan <- function() {
      if (is.null(model)) model <<- study_model(model_file)
      model
    }
  }
  e
}

prepare_study <- function(config, settings, directory) {
  stopifnot(all(c("n", "dgp", "predictor", "replicate") %in% names(config)),
    all(config$n %in% c(250L, 1000L)), all(config$dgp %in% paste0("D", 1:4)),
    all(config$predictor %in% paste0("X", 1:3)),
    all(config$replicate > 0 & config$replicate == floor(config$replicate)),
    !anyDuplicated(study_key(config)), !settings$test_replicate %in% config$replicate)
  e <- study_worker(settings)
  for (part in c("fits", "results", "extracted", "failures", "summary", "logs"))
    dir.create(file.path(directory, part), recursive = TRUE, showWarnings = FALSE)
  directory <- normalizePath(directory)
  manifest <- list(config = config, settings = settings, source_md5 = study_sources())
  manifest_file <- file.path(directory, "manifest.rds")
  if (file.exists(manifest_file)) {
    # A resumed study must use the same code and settings. Use a new directory
    # when changing the target, simulation design, chain length or evaluation.
    stopifnot(identical(manifest, readRDS(manifest_file)))
  } else {
    study_save(manifest, manifest_file)
    write.csv(config, file.path(directory, "config.csv"), row.names = FALSE)
    study_save(sessionInfo(), file.path(directory, "environment.rds"))
  }
  tests_file <- file.path(directory, "test-sets.rds")
  if (!file.exists(tests_file)) study_save(e$make_test_sets(settings$seed_base,
    settings$test_replicate, settings$n_test), tests_file)
  model_file <- file.path(directory, "stan-model.rds")
  if (!file.exists(model_file)) {
    model <- e$compile_nils_stan()
    study_save(model, model_file)
    study_model(model_file, model)
  }
  invisible(manifest)
}

# Console flags affect display only; they are not part of saved study settings.
run_study_task <- function(directory, task_id, worker = NULL,
                           verbose = TRUE, engine_output = FALSE) {
  directory <- normalizePath(directory)
  manifest <- readRDS(file.path(directory, "manifest.rds"))
  stopifnot(task_id >= 1L, task_id <= nrow(manifest$config), task_id == floor(task_id),
            identical(study_sources(), manifest$source_md5))
  settings <- manifest$settings; row <- manifest$config[task_id, ]
  key <- study_key(row)
  checkpoint <- file.path(directory, "fits", paste0(key, ".rds"))
  result_file <- file.path(directory, "results", paste0(key, ".rds"))
  completed <- study_read(result_file)
  if (!is.null(completed)) {
    # Recover even if a process stopped after saving its result but before
    # removing the temporary fit. Older results without summaries are untouched.
    if (isFALSE(settings$keep_draws) && !is.null(completed$model_summaries)) unlink(checkpoint)
    run_status(paste(key, "already complete; skipped."), verbose)
    return(invisible(completed))
  }
  if (is.null(worker)) worker <- study_worker(settings, file.path(directory, "stan-model.rds"))
  timing <- numeric(); stage <- "fitting"
  started <- proc.time()[["elapsed"]]
  timed <- function(label, expr) {
    stage <<- label; start <- proc.time()[["elapsed"]]
    value <- force(expr)
    timing[label] <<- proc.time()[["elapsed"]] - start
    value
  }
  run_status(sprintf("%s | %s | n=%d | replicate=%d | %d draws/chain",
    row$dgp, row$predictor, row$n, row$replicate, worker$MCMC_NDRAWS), verbose)
  report <- function(label, expr) run_stage(label, expr, verbose, engine_output)
  tryCatch({
    s <- if (file.exists(checkpoint)) timed("load_fit", report("Loading saved fits", study_read(checkpoint))) else NULL
    if (is.null(s)) {
      s <- timed("fitting", worker$run_cell(row$dgp, row$predictor, row$n,
                                          row$replicate, seed_base = settings$seed_base, verbose = verbose,
                                          engine_output = engine_output, show_header = FALSE))
      s$settings$study <- settings
      timed("save_fit", report("Saving fit checkpoint", study_save(s, checkpoint)))
    }
    extraction_file <- file.path(directory, "extracted", paste0(key, ".rds"))
    fitted <- if (file.exists(extraction_file)) timed("load_extraction", report("Loading saved extraction", study_read(extraction_file))) else NULL
    if (is.null(fitted)) {
      fitted <- timed("extraction", worker$results_on_train(s, xval_N = settings$grid_N,
        probs = settings$probs, B_M3 = settings$B_M3,
        boot_chunk_size = settings$grid_chunk_size, grid_chunk_size = settings$grid_chunk_size,
        verbose = verbose, engine_output = engine_output))
      if (settings$save_extracted) timed("save_extraction", report("Saving extraction", study_save(fitted, extraction_file)))
    }
    structural <- timed("structural", report("Evaluating structural effects",
      worker$evaluate_replicate(s, fitted, grid_N = settings$grid_N)))
    tests <- readRDS(file.path(directory, "test-sets.rds"))
    unseen <- timed("test", report("Evaluating unseen responses",
      worker$evaluate_test_sample(s, tests[[row$dgp]], chunk_size = settings$kde_chunk_size)))
    probs <- if (is.null(settings$parameter_probs)) POSTERIOR_PROBS else settings$parameter_probs
    summaries <- timed("summaries", report("Summarizing parameters and MCMC diagnostics",
      compact_models(s, probs)))
    timing["total"] <- proc.time()[["elapsed"]] - started
    out <- list(scenario = s$scenario[c("dgp", "predictor", "n", "replicate", "seed")],
      study_settings = settings, source_md5 = manifest$source_md5, status = "OK",
      model_summaries = summaries, fit_settings = s$settings,
      training_data = s$scenario[c("dat", "f0", "eps", "gamma", "region")],
      structural = structural, test = unseen, timing = timing,
      execution = list(job = Sys.getenv("SLURM_JOB_ID"), task = task_id,
                       host = Sys.info()[["nodename"]], session = sessionInfo()))
    stage <- "save_result"
    study_save(out, result_file)
    # Draws are needed through interval construction AND predictive scoring.
    # Delete only after the compact result is atomically saved. Failed runs
    # keep their checkpoints; completed-run skips repeat any interrupted cleanup.
    if (isFALSE(settings$keep_draws)) {
      unlink(checkpoint)
      if (file.exists(checkpoint)) warning("Result saved, but fit checkpoint could not be removed: ", checkpoint)
    }
    run_status(sprintf("Results saved: %s; elapsed %.1f s", key,
                       proc.time()[["elapsed"]] - started), verbose)
    invisible(out)
  }, error = function(error) {
    run_status(paste0("FAILED at ", stage, ": ", conditionMessage(error)), verbose)
    study_save(list(scenario = row, stage = stage, message = conditionMessage(error),
      timing = timing, time = Sys.time()), file.path(directory, "failures", paste0(key, ".rds")))
    stop(error)
  })
}

# Collect one (n, DGP, predictor) group at a time; loading all 2,400 spatial
# evaluations together would need much more RAM. Full fits are never loaded.
collect_study <- function(directory) {
  directory <- normalizePath(directory)
  manifest <- readRDS(file.path(directory, "manifest.rds"))
  config <- manifest$config; settings <- manifest$settings
  keys <- study_key(config)
  inventory <- config
  inventory$status <- "MISSING"
  inventory$message <- ""
  inventory$seconds <- NA_real_
  groups <- interaction(config$n, config$dgp, config$predictor, drop = TRUE, sep = "_")
  all_scores <- list(); timing <- list()
  for (group in levels(groups)) {
    indices <- which(groups == group)
    evaluations <- vector("list", length(indices))
    for (j in seq_along(indices)) {
      i <- indices[j]
      path <- file.path(directory, "results", paste0(keys[i], ".rds"))
      record <- study_read(path)
      # A missing/failed run stays in the denominator and in the retry list.
      if (is.null(record)) {
        if (file.exists(path)) inventory$status[i] <- "UNREADABLE"
        failure <- file.path(directory, "failures", paste0(keys[i], ".rds"))
        if (file.exists(failure)) {
          fail <- readRDS(failure)
          inventory$status[i] <- "FAILED"; inventory$message[i] <- paste(fail$stage, fail$message)
        }
      } else {
        stopifnot(identical(record$study_settings, settings),
          identical(record$source_md5, manifest$source_md5),
          identical(study_key(as.data.frame(record$scenario)), keys[i]))
        inventory$status[i] <- "OK"; inventory$seconds[i] <- record$timing[["total"]]
        evaluations[[j]] <- record$structural
        all_scores[[as.character(i)]] <- cbind(config[rep(i, nrow(record$test$scores)), ], record$test$scores)
        timing[[as.character(i)]] <- cbind(config[i, ], as.data.frame(as.list(record$timing)))
      }
    }
    # Empty evaluations preserve the requested replicate rows. Model-specific
    # failures (e.g. no M4 LOO winner) remain NULL within successful records.
    for (j in seq_along(indices)) if (is.null(evaluations[[j]])) {
      row <- config[indices[j], ]
      evaluations[[j]] <- list(scenario = as.list(row[c("dgp", "predictor", "n", "replicate")]),
        grid = seq(-.5, .5, length.out = settings$grid_N), interval_probs = settings$probs,
        models = setNames(vector("list", 6L), paste0("m", 0:5)))
    }
    summary <- aggregate_evaluations(evaluations)
    study_save(summary, file.path(directory, "summary", paste0(group, ".rds")))
    rm(evaluations, record, summary); gc(FALSE)
  }
  # Include missing runs explicitly in prediction-score tables, too.
  for (i in which(inventory$status != "OK")) all_scores[[as.character(i)]] <- cbind(config[rep(i, 6L), ],
    data.frame(model = paste0("m", 0:5), RMSE = NA_real_, MAE = NA_real_, log_score = NA_real_,
      score_type = NA_character_, bandwidth = NA_real_, status = inventory$status[i]))
  scores <- do.call(rbind, all_scores); rownames(scores) <- NULL
  scores <- scores[order(scores$n, scores$dgp, scores$predictor, scores$replicate, scores$model), ]
  score_groups <- interaction(scores$n, scores$dgp, scores$predictor, scores$model, drop = TRUE)
  score_summary <- do.call(rbind, lapply(split(scores, score_groups), function(x) {
    out <- x[1L, c("n", "dgp", "predictor", "model")]
    out$n_runs <- nrow(x); out$n_ok <- sum(x$status == "OK")
    for (metric in c("RMSE", "MAE", "log_score")) {
      valid <- !is.na(x[[metric]])
      out[[metric]] <- if (any(valid)) mean(x[[metric]][valid]) else NA_real_
      out[[paste0(metric, "_n")]] <- sum(valid)
    }
    out
  }))
  write.csv(inventory, file.path(directory, "summary", "inventory.csv"), row.names = FALSE)
  write.csv(config[inventory$status != "OK", ], file.path(directory, "summary", "retry.csv"), row.names = FALSE)
  write.csv(scores, file.path(directory, "summary", "test-scores.csv"), row.names = FALSE)
  write.csv(score_summary, file.path(directory, "summary", "test-summary.csv"), row.names = FALSE)
  if (length(timing)) {
    columns <- unique(unlist(lapply(timing, names)))
    timing <- lapply(timing, function(x) {x[setdiff(columns, names(x))] <- NA_real_; x[columns]})
    write.csv(do.call(rbind, timing), file.path(directory, "summary", "timings.csv"), row.names = FALSE)
  }
  cat(sum(inventory$status == "OK"), "of", nrow(config), "runs complete; summaries and retry.csv written.\n")
  invisible(inventory)
}
