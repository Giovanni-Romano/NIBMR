# Standalone diagnostic; run from R/. No simulation scripts or checkpoints change.
# First compare median solvers on the same bootstrap samples:
# Rscript src/sim-study/diagnose-bootstrap-median.R STUDY [KEY] [B]
# Add "refit" to also compare final M3 fits:
# Rscript src/sim-study/diagnose-bootstrap-median.R STUDY KEY B refit
# Defaults: KEY=n1000_D1_X3_r001, B=10. Reports go in STUDY/diagnostics/.
# In an R session, source this file and call diagnose_bootstrap_median(s, worker).
# relative_loss_gap is (FN - BR) / max(1, BR absolute-error loss).
# Relative differences divide by the larger magnitude of the two values;
# scales/k columns report the largest difference across their two components.
# RMSE columns compare the two fitted predictors, not either predictor to truth.
# Final M3 comparisons change BOTH the start and loss scales; they do not
# identify their separate contributions or establish a global optimum.

diagnose_bootstrap_median <- function(s, worker, B = 10L, refit = FALSE,
                                      seed = s$scenario$seed) {
  stopifnot(B >= 1L, B == as.integer(B))
  dat <- s$scenario$dat; predictor <- s$scenario$predictor
  method <- s$fits$m3$method
  n <- nrow(dat)
  strata <- if (predictor == "X3") split(seq_len(n), dat$region) else list(seq_len(n))
  grid <- seq(-.5, .5, length.out = 50L)
  cn <- worker$c_n(n); kn <- worker$k_n(n)

  # Both solvers minimize sum(abs(y - X %*% beta)). Capture warnings without
  # hiding solver errors. Different coefficients need not imply different fits.
  median_fit <- function(X, y, solver) {
    messages <- character()
    fit <- withCallingHandlers(quantreg::rq.fit(X, y, tau = .5, method = solver),
      warning = function(w) {
        messages <<- c(messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    fitted <- drop(X %*% fit$coefficients)
    residuals <- y - fitted
    scales <- worker$est_s1_s2(residuals)
    k <- kn / scales
    stopifnot(all(is.finite(c(fit$coefficients, fitted, scales, k))), all(scales > 0))
    list(start = fit$coefficients, fitted = fitted, s = scales, k = k,
         objective = sum(abs(residuals)), warnings = messages)
  }

  # This is the same optimization recipe as boot_CI(), using the method chosen
  # on the full dataset. Only the preliminary median solver differs.
  m3_fit <- function(X, y, setup) {
    loss <- function(beta) sum(worker$loss_asymm(y - drop(X %*% beta), k = setup$k, c = cn))
    gradient <- function(beta) -drop(worker$score_pop(seq_len(ncol(X)), beta, y, X,
      k = setup$k, c = cn, freq = 1, link = identity,
      link_prime = function(x) rep(1, length(x))))
    fit <- worker$optim_with_diagnostics(setup$start, fn = loss,
      gr = if (method %in% c("BFGS", "CG", "L-BFGS-B")) gradient else NULL,
      method = method, control = worker$m3_optim_control(method),
      retry = worker$OPTIM_RETRY, retry_factor = worker$OPTIM_RETRY_FACTOR,
      ntry = worker$OPTIM_NTRY_M3, warn = FALSE)
    stopifnot(fit$usable)
    fit
  }

  relative_difference <- function(a, b) max(abs(a - b) / pmax(abs(a), abs(b), 1e-12))
  cases <- vector("list", B)
  summary <- vector("list", B)
  # One seed for the whole procedure, as in boot_CI(). The sampling statements
  # below are copied from that procedure, including within-region resampling.
  set.seed(seed)
  for (b in seq_len(B)) {
    idx <- integer(n)
    for (group in strata)
      idx[group] <- group[sample.int(length(group), length(group), replace = TRUE)]
    datboot <- dat[idx, , drop = FALSE]
    dd <- worker$make_design_M3(datboot, predictor)
    X <- dd$X; y <- datboot$y
    br <- median_fit(X, y, "br")
    fn <- median_fit(X, y, "fn")

    # Check omitted training observations too: agreement only at duplicated
    # bootstrap rows can conceal different predictions elsewhere.
    T1 <- mgcv::PredictMat(dd$bs1, data.frame(x = dat$x1))
    T2 <- mgcv::PredictMat(dd$bs2, data.frame(x = dat$x2))
    R <- if (predictor == "X3")
      dd$region_contrasts[match(as.character(dat$region), dd$region_levels), , drop = FALSE] else
      matrix(numeric(), nrow = n, ncol = 0L)
    original_X <- cbind(1, T1, T2, R)
    delta <- fn$start - br$start
    summary[[b]] <- data.frame(
      bootstrap = b, unique_rows = length(unique(idx)),
      br_warning = paste(unique(br$warnings), collapse = " | "),
      fn_warning = paste(unique(fn$warnings), collapse = " | "),
      abs_loss_br = br$objective, abs_loss_fn = fn$objective,
      relative_loss_gap = (fn$objective - br$objective) / max(1, br$objective),
      fitted_rmse = sqrt(mean((fn$fitted - br$fitted)^2)),
      fitted_max_abs = max(abs(fn$fitted - br$fitted)),
      original_fitted_rmse = sqrt(mean(drop(original_X %*% delta)^2)),
      s1_br = br$s[1L], s1_fn = fn$s[1L], s2_br = br$s[2L], s2_fn = fn$s[2L],
      scales_relative_difference = relative_difference(br$s, fn$s),
      k1_br = br$k[1L], k1_fn = fn$k[1L], k2_br = br$k[2L], k2_fn = fn$k[2L],
      k_relative_difference = relative_difference(br$k, fn$k), row.names = NULL)

    # SANN consumes random numbers. Replay the original BR-based optimizer
    # even in median-only mode so subsequent resamples match boot_CI(). All FN
    # refits wait until AFTER the original sampling stream has been replayed.
    baseline <- if (refit || method == "SANN") m3_fit(X, y, br) else NULL
    cases[[b]] <- list(indices = idx, design = dd, y = y, original_X = original_X,
                       br = br, fn = fn, baseline = baseline)
  }

  comparison <- NULL
  if (refit) {
    comparison <- vector("list", B)
    for (b in seq_len(B)) {
      case <- cases[[b]]; dd <- case$design; X <- dd$X; y <- case$y
      alternative <- m3_fit(X, y, case$fn)
      beta_br <- case$baseline$par; beta_fn <- alternative$par
      delta <- beta_fn - beta_br
      # Compare centered component curves using the ORIGINAL training sample,
      # just as boot_CI() does when constructing the reported intervals.
      curve_difference <- function(basis, columns) {
        on_grid <- drop(mgcv::PredictMat(basis, data.frame(x = grid)) %*% delta[columns])
        on_train <- drop(case$original_X[, columns, drop = FALSE] %*% delta[columns])
        max(abs(on_grid - mean(on_train)))
      }
      loss_at <- function(beta, k)
        sum(worker$loss_asymm(y - drop(X %*% beta), k = k, c = cn))
      # If k changes, the two optimized objectives are different functions.
      # Cross-evaluate both fits under EACH k; do not compare only their own minima.
      comparison[[b]] <- data.frame(bootstrap = b,
        convergence_br = case$baseline$convergence, convergence_fn = alternative$convergence,
        original_fitted_rmse = sqrt(mean(drop(case$original_X %*% delta)^2)),
        f1_grid_max_abs = curve_difference(dd$bs1, dd$idx$x1),
        f2_grid_max_abs = curve_difference(dd$bs2, dd$idx$x2),
        br_loss_at_br = loss_at(beta_br, case$br$k), br_loss_at_fn = loss_at(beta_fn, case$br$k),
        fn_loss_at_br = loss_at(beta_br, case$fn$k), fn_loss_at_fn = loss_at(beta_fn, case$fn$k))
      cases[[b]]$alternative <- alternative
    }
    comparison <- do.call(rbind, comparison)
  }

  list(scenario = s$scenario[c("dgp", "predictor", "n", "replicate", "seed")],
       seed = seed, method = method, B = B, grid = grid,
       settings = list(D_SPLINE = worker$D_SPLINE, control = worker$m3_optim_control(method),
         retry = worker$OPTIM_RETRY, ntry = worker$OPTIM_NTRY_M3,
         retry_factor = worker$OPTIM_RETRY_FACTOR),
       summary = do.call(rbind, summary), m3_comparison = comparison,
       cases = cases, session = sessionInfo())
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args)) stop("Usage: Rscript src/sim-study/diagnose-bootstrap-median.R STUDY [KEY] [B] [refit]")
  directory <- args[1L]
  key <- if (length(args) >= 2L) args[2L] else "n1000_D1_X3_r001"
  B <- if (length(args) >= 3L) as.integer(args[3L]) else 10L
  refit <- length(args) >= 4L && args[4L] == "refit"
  source("src/sim-study/08-pipeline.R")
  manifest <- readRDS(file.path(directory, "manifest.rds"))
  worker <- study_worker(manifest$settings)
  s <- readRDS(file.path(directory, "fits", paste0(key, ".rds")))
  cat("Comparing median solvers for", key, "on", B, "bootstrap samples.\n")
  if (s$fits$m3$method == "SANN")
    cat("Replaying SANN to preserve resamples; refit differences also include optimizer randomness.\n")
  out <- diagnose_bootstrap_median(s, worker, B = B, refit = refit)
  print(out$summary[c("bootstrap", "br_warning", "fn_warning", "relative_loss_gap",
                      "fitted_rmse", "original_fitted_rmse",
                      "scales_relative_difference", "k_relative_difference")], row.names = FALSE)
  if (refit) print(out$m3_comparison, row.names = FALSE)
  reports <- file.path(directory, "diagnostics")
  dir.create(reports, showWarnings = FALSE)
  output <- file.path(reports, paste0("bootstrap-median-", key, if (refit) "-refit" else "", ".rds"))
  saveRDS(out, output)
  cat("Diagnostic saved to", output, "\n")
}
