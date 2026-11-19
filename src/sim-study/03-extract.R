# === === === === === === === === === === === === === === === === === === ===
# src/sim-study/03-extract.R
# Reconstruct fitted effects and the complete predictor from one run_cell()
# result, and construct comparable simulation truth on observations or a grid.
# This script defines extraction helpers; it does not fit models, aggregate
# simulation replicates, or write result files.
# === === === === === === === === === === === === === === === === === === ===

source("src/sim-study/00-progress.R", local = TRUE)
library(coda)

# Naming and centering conventions used below
# ------------------------------------------
# The DGP supplies raw effects g1(x1), g2(x2), and gamma[r] (X3 only).
# Its structural signal is F(x1, x2, r) = g1(x1) + g2(x2) + gamma[r].
# For one training sample, define
#   a1 = mean(g1(x1_train)), a2 = mean(g2(x2_train)),
#   aS = mean(gamma) over ALL regions, or 0 for X1/X2.
# The same signal can be written as
#   F = (a1 + a2 + aS) + (g1 - a1) + (g2 - a2) + (gamma - aS).
# Thus the centered truth has
#   int = a1 + a2 + aS, f1 = g1 - a1, f2 = g2 - a2,
#   spat = gamma - aS, fcont = f1 + f2, f = fcont + spat,
#   f_full = int + f = F.
# Omit spat from these sums for X1/X2. With unequal region counts, the
# observation mean of spat (and therefore of f) need not be zero.
#
# results_on_train() uses these field names for fitted effects. Truth helpers
# use them for known DGP values. Neither centering nor adding the intercept
# changes the fitted model or the structural signal. f_full in a truth object
# means F, not automatically the conditional mode/median: for example, the
# D4 error mode is nonzero, so its conditional mode is F + mode(epsilon).


# True effects and centering ----

#' Split a true effect into its projection on the space of linear functions and
#' the nonlinear remainder.
#' @param effect is a function;
#' @param x is the original training covariate vector
#' @param new_x is the vector of points at which to evaluate the projection and 
#'    nonlinear part. It only chooses evaluation points; it never determines 
#'    the projection.
#' @param method is the method for computing the projection: "training-data" or "numerical-integration".
#' @param lower and upper are the bounds of the integration interval.
#' @param density is the probability density function on [lower, upper].
#' 
#' The projection is intercept + slope*x; nonlinear is effect(x) minus that projection.
#' Slopes refer to the original x units, not standardized covariates.
#' Methods: 
#'  - "training-data" uses empirical expectations
#'  - "numerical-integration" uses the provided density on [lower, upper], 
#'     defaulting to the study's uniform covariates.
#' intercept is the affine intercept at x=0. mean_effect is the constant
#' moved into the global intercept when the full effect is centered.
split_true_effect <- function(effect, x = NULL, new_x = x,
                              method = c("training-data", "numerical-integration"),
                              lower = -.5, upper = .5,
                              density = function(x) rep(1 / (upper - lower), length(x))) {
  method <- match.arg(method)
  if (method == "training-data") {
    stopifnot(length(x) > 1L, var(x) > 0)
    mean_x <- mean(x); mean_effect <- mean(effect(x))
    slope <- sum((x - mean_x) * (effect(x) - mean_effect)) / sum((x - mean_x)^2)
  } else {
    mass <- integrate(density, lower, upper, rel.tol = 1e-10)$value
    stopifnot(mass > 0)
    expectation <- function(f) integrate(function(z) f(z) * density(z), lower, upper,
                                         rel.tol = 1e-10)$value / mass
    mean_x <- expectation(identity); mean_effect <- expectation(effect)
    variance_x <- expectation(function(z) (z - mean_x)^2)
    stopifnot(variance_x > 0)
    slope <- expectation(function(z) (z - mean_x) * (effect(z) - mean_effect)) /
      variance_x
  }
  intercept <- mean_effect - slope * mean_x
  linear <- intercept + slope * new_x
  list(intercept = intercept, slope = slope, mean_x = mean_x, mean_effect = mean_effect,
       linear = linear, nonlinear = effect(new_x) - linear, method = method)
}

#' Compute the three constants needed to center simulation truth.
#'
#' This function only calculates offsets; it does not alter any effects.
#' Continuous offsets use the true effects evaluated at the observations in
#' THIS training sample, not the covariate means or the grid-point means.
#' The spatial offset gives every region equal weight, regardless of how many
#' observations belong to it. No effect is divided by a standard deviation.
#'
#' @param scenario The scenario returned by generate_data().
#'   Uses dat, predictor, f0$f1, f0$f2, and gamma (for X3).
#' @return A named numeric vector c(f1 = a1, f2 = a2, spat = aS).
#'   For X1/X2, spat is 0. Offsets depend on the training sample; compute them
#'   once per scenario and reuse them for its training and grid truth.
#' @examples
#' offsets <- truth_centering(s$scenario)
truth_centering <- function(scenario) {
  f1 <- scenario$f0$f1
  f2 <- scenario$f0$f2
  stopifnot(length(f1) == nrow(scenario$dat), length(f2) == nrow(scenario$dat),
            length(f1) > 0L, all(is.finite(f1)), all(is.finite(f2)))
  
  # Linear offsets
  a1 <- mean(f1)
  a2 <- mean(f2)
  
  # Spatial offset
  aS <- 0
  if (scenario$predictor == "X3") {
    stopifnot(length(scenario$gamma) > 0L, all(is.finite(scenario$gamma)))
    aS <- mean(scenario$gamma)
  }
  
  c(f1 = a1, f2 = a2, spat = aS)
}

#' Apply supplied offsets to a raw truth object (shared internal helper).
#'
#' training_truth() and grid_truth() call this helper automatically. Usually
#' there is no need to call it directly. It does NOT estimate offsets from the
#' supplied truth values; this prevents accidental recentering on grid points.
#'
#' @param truth An uncentered truth list constructed by training_truth() or
#'   grid_truth() with centering = NULL. Its int is initially 0 and its f and
#'   f_full both contain the raw structural signal. This is not a fitted result.
#' @param centering Named numeric vector returned by truth_centering(), or NULL.
#' @return With NULL, returns truth unchanged. Otherwise returns a centered
#'   copy: subtracts a1 from f1, a2 from f2, aS from spat, a1+a2 from fcont,
#'   and a1+a2+aS from f; sets int = a1+a2+aS. A spatial grid's f is a list,
#'   so the last subtraction is applied to each regional surface separately.
#'   f_full is unchanged. The centering vector is stored as metadata, and a
#'   second non-NULL centering call on that returned object raises an error.
center_truth <- function(truth, centering) {
  # Return uncentered truth if no offsets are supplied
  if (is.null(centering)) return(truth)
  
  # Checks
  stopifnot(is.numeric(centering), length(centering) == 3L, 
            setequal(names(centering), c("f1", "f2", "spat")), 
            all(is.finite(centering)))
  if (!is.null(truth$centering)) stop("Truth has already been centered.")
  
  # Apply centering: see truth_centering() for the definitions of a1, a2, aS and
  # descriptions at the top of the script.
  truth$int <- truth$int + unname(sum(centering))
  truth$f1 <- truth$f1 - unname(centering["f1"])
  truth$f2 <- truth$f2 - unname(centering["f2"])
  truth$fcont <- truth$fcont - unname(sum(centering[c("f1", "f2")]))
  # We update spatial effect only if present
  if (!is.null(truth$spat)) truth$spat <- truth$spat - unname(centering["spat"])
  
  # truth$f is a list when computed on a grid for predictor X3, otherwise it is numeric
  truth$f <- if (is.list(truth$f)) {
    lapply(truth$f, function(z) z - truth$int)
  } else {
    truth$f - truth$int
  }
  
  # save centering as metadata
  truth$centering <- centering
  
  return(truth)
}

#' Construct true effects at the observed training covariates and regions.
#'
#' @param scenario The scenario returned by generate_data(), usually s$scenario.
#'   In addition to the fields used by truth_centering(), uses f0$f and region
#'   (numeric indices into gamma, as supplied by generate_data(), for X3).
#' @param centering Offsets for the training sample. Defaults to computing them
#'   with truth_centering(scenario). Supply NULL for uncentered truth instead.
#' @return A list with scalar int and length-n vectors f1, f2, fcont, f,
#'   f_full, and spat (NULL for X1/X2). spat[i] is the effect of observation i's
#'   region, so regions with more observations appear more often in this vector.
#'   f_full is always exactly scenario$f0$f. With non-NULL centering, the list
#'   also contains the offsets in centering. True f1L/f1NL decompositions and
#'   uncertainty intervals are not constructed by this function.
#' @examples
#' truth_train <- training_truth(s$scenario)
#' fitted <- results_on_train(s)
#' # Both vectors below refer to the same observed covariates and regions:
#' fitted$m0$point_est$f_full - truth_train$f_full
training_truth <- function(scenario, centering = truth_centering(scenario)) {
  f0 <- scenario$f0
  spat <- if (scenario$predictor == "X3") scenario$gamma[scenario$region] else NULL
  
  truth <- list(f1 = f0$f1, f2 = f0$f2, fcont = f0$f1 + f0$f2,
                f = f0$f, f_full = f0$f, int = 0, spat = spat)
  
  center_truth(truth, centering)
}

#' Construct true effects on the Cartesian grid new_x by new_x.
#'
#' @param pred Predictor specification: "X1", "X2", or "X3".
#' @param new_x Numeric covariate grid, shared by x1 and x2. A surface entry
#'   [i, j] corresponds to x1 = new_x[i] and x2 = new_x[j].
#' @param sp_eff For X3, the complete vector of true regional effects (normally
#'   s$scenario$gamma). Regional surfaces follow this vector's order; match it
#'   to the region order of the fitted summaries. Ignored for X1/X2.
#' @param centering Offsets returned by truth_centering(s$scenario). Unlike
#'   training_truth(), the default is NULL to preserve historical grid_truth()
#'   calls: the grid alone cannot determine training-sample offsets.
#' @return A truth list with length-G f1/f2, G-by-G fcont, scalar int, and
#'   spat (the length-R regional vector for X3; otherwise NULL). f and f_full
#'   are G-by-G matrices for X1/X2, or lists of R such matrices for X3.
#'   With centering = NULL, int is 0 and f1/f2/fcont/f have their original raw
#'   values. With offsets supplied, these components are centered using the
#'   TRAINING observations, int is the sum of offsets, and centering is stored
#'   in the list. f_full always contains the same uncentered structural signal.
#'
#' @details For component comparisons with results_on_train(), supply offsets;
#'   for full-function comparisons, use f_full on both sides without a shift.
#'   The fitted grid must use the same bounds and number of points as new_x.
#' @examples
#' new_x <- seq(-0.5, 0.5, length.out = 100)
#' offsets <- truth_centering(s$scenario)
#' truth_grid <- grid_truth(s$scenario$predictor, new_x,
#'                          sp_eff = s$scenario$gamma, centering = offsets)
#' # The default results_on_train(s) uses this same 100-by-100 grid.
#' # In X3, [[r]] selects region r's surface for f and f_full.
grid_truth <- function(pred, new_x, sp_eff = NULL, centering = NULL) {
  ft <- f_true(pred, new_x, new_x)
  f1 <- ft$f1
  f2 <- ft$f2
  fcont <- outer(f1, f2, "+")
  if (pred == "X3"){
    f <- lapply(sp_eff, function(S) S + fcont)
  } else {
    f <- fcont
  }
  
  truth <- list(f1 = ft$f1, f2 = ft$f2, fcont = fcont, f = f,
                f_full = f, int = 0, spat = if (pred == "X3") sp_eff else NULL)
  
  center_truth(truth, centering)
}


# Helpers for point estimates and CI ----
get_pe <- function(draws, X, eff_idx){
  if (length(eff_idx) == 1){
    mean(draws[ , eff_idx])
  } else {
    drop(X[ , eff_idx] %*% colMeans(draws[ , eff_idx]))
  }
}

get_CIquant <- function(draws, X = NULL, probs = c(0.025, 0.975)) {
  if (is.null(X)){
    out <- setNames(as.numeric(quantile(draws, probs)), c("LB", "UB"))
  } else {
    out <- t(apply(draws %*% t(X), 2, quantile, probs = probs))
    colnames(out) <- c("LB", "HB")
  }
  return(out)
}
get_CIhpd <- function(draws, X = NULL, prob = 0.95) {
  if (!is.matrix(draws)) {
    H <- HPDinterval(as.mcmc(draws), prob = prob)
  } else {
    H <- HPDinterval(as.mcmc(draws %*% t(X)), prob = prob)
  }
  
  names(H) <- c("LB", "UB")
  return(H)
}



#' Empty M3 intervals with the same shapes as the bootstrap output.
#' Also used when bootstrap is disabled with B_M3=0.
m3_na_intervals <- function(new_x, region_levels = NULL) {
  G <- length(new_x)
  bounds <- c("LB", "UP")
  curve <- matrix(NA_real_, G, 2L, dimnames = list(as.character(new_x), bounds))
  surface <- array(NA_real_, c(G, G, 2L),
    dimnames = list(x1 = as.character(new_x), x2 = as.character(new_x), bound = bounds))
  spatial <- !is.null(region_levels)
  surfaces <- if (spatial) setNames(rep(list(surface), length(region_levels)), region_levels) else surface
  list(int = setNames(rep(NA_real_, 2L), bounds),
       f1L = NULL, f1NL = curve, f1 = curve,
       f2L = NULL, f2NL = curve, f2 = curve,
       spat = if (spatial) matrix(NA_real_, length(region_levels), 2L,
         dimnames = list(region_levels, bounds)) else NULL,
       fcont = surface, fcont_full = surface, f = surfaces, f_full = surfaces,
       slope1 = setNames(rep(NA_real_, 2L), bounds),
       slope2 = setNames(rep(NA_real_, 2L), bounds))
}

# Pointwise percentile intervals for M3, using the method selected on the full data.
# Rebuild the basis and loss parameters in each resample. Center fitted curves
# over the ORIGINAL training covariates, moving their means to the intercept.
# X3 resamples within regions and uses equal-region spatial centering.
# Return only intervals; f_full includes the intercept before taking quantiles.
# seed sets the RNG once. SANN and resampling then use the same RNG stream.
boot_CI <- function(B = 100L, dat, predictor, new_x, method,
                    probs = c(.025, .975), seed = NULL, chunk_size = 512L) {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(dat); G <- length(new_x)
  regional <- predictor == "X3"
  region_levels <- if (regional) levels(factor(dat$region)) else NULL
  strata <- if (regional) split(seq_len(n), dat$region) else list(seq_len(n))
  intercepts <- numeric(B)
  slopes <- matrix(NA_real_, B, 2L)
  # The same original training reference defines true and bootstrap slopes.
  slope_weights <- lapply(dat[c("x1", "x2")], function(x) (x - mean(x)) / sum((x - mean(x))^2))
  f1 <- f2 <- matrix(NA_real_, B, G)
  spat <- if (regional) matrix(NA_real_, B, length(region_levels)) else NULL

  control <- m3_optim_control(method)

  for (b in seq_len(B)) {
    
    idx <- integer(n)
    for (group in strata) idx[group] <- group[sample.int(length(group), length(group), replace = TRUE)]
    datboot <- dat[idx, , drop = FALSE]
    dd <- make_design_M3(datboot, predictor)
    X <- dd$X; y <- datboot$y
    # Apply the same setup recipe to this resample's data and basis. Its start
    # and scales are re-estimated here; seeding remains once at the top.
    setup <- nils_setup(X, y)
    k <- setup$k; cn <- setup$c
    loss <- function(beta) sum(loss_asymm(y - drop(X %*% beta), k = k, c = cn))
    gradient <- function(beta) -drop(score_pop(seq_len(ncol(X)), beta, y, X,
      k = k, c = cn, freq = 1, link = identity, link_prime = function(x) rep(1, length(x))))
    fit <- optim_with_diagnostics(setup$start, fn = loss,
      gr = if (method %in% c("BFGS", "CG", "L-BFGS-B")) gradient else NULL,
      method = method, control = control, retry = OPTIM_RETRY,
      retry_factor = OPTIM_RETRY_FACTOR, ntry = OPTIM_NTRY_M3, warn = OPTIM_WARN)
    if (!fit$usable) stop("M3 bootstrap fit failed at replicate ", b, ".")
    beta <- fit$par

    # Use this resample's basis, but the original training sample for centering.
    train1 <- drop(mgcv::PredictMat(dd$bs1, data.frame(x = dat$x1)) %*% beta[dd$idx$x1])
    train2 <- drop(mgcv::PredictMat(dd$bs2, data.frame(x = dat$x2)) %*% beta[dd$idx$x2])
    mean1 <- mean(train1); mean2 <- mean(train2)
    slopes[b, ] <- c(sum(slope_weights$x1 * train1), sum(slope_weights$x2 * train2))
    f1[b, ] <- drop(mgcv::PredictMat(dd$bs1, data.frame(x = new_x)) %*% beta[dd$idx$x1]) - mean1
    f2[b, ] <- drop(mgcv::PredictMat(dd$bs2, data.frame(x = new_x)) %*% beta[dd$idx$x2]) - mean2
    intercepts[b] <- beta[1L] + mean1 + mean2
    if (regional) {
      gamma <- drop(dd$region_contrasts %*% beta[dd$idx$region])
      gamma <- gamma[match(region_levels, dd$region_levels)]
      spat[b, ] <- gamma - mean(gamma)
      intercepts[b] <- intercepts[b] + mean(gamma)
    }
  }

  intervals <- m3_na_intervals(new_x, region_levels)
  column_quantiles <- function(M) t(apply(M, 2L, quantile, probs = probs, names = FALSE))
  intervals$int[] <- quantile(intercepts, probs = probs, names = FALSE)
  intervals$slope1[] <- quantile(slopes[, 1L], probs = probs, names = FALSE)
  intervals$slope2[] <- quantile(slopes[, 2L], probs = probs, names = FALSE)
  intervals$f1[,] <- column_quantiles(f1); intervals$f1NL <- intervals$f1
  intervals$f2[,] <- column_quantiles(f2); intervals$f2NL <- intervals$f2
  if (regional) intervals$spat[,] <- column_quantiles(spat)

  # Summarize one grid chunk at a time, including dependence between effects.
  surface <- function(shift) {
    bounds <- matrix(NA_real_, G * G, 2L)
    for (first in seq.int(1L, G * G, by = chunk_size)) {
      points <- seq.int(first, min(G * G, first + chunk_size - 1L))
      i1 <- (points - 1L) %% G + 1L; i2 <- (points - 1L) %/% G + 1L
      combined <- sweep(f1[, i1, drop = FALSE] + f2[, i2, drop = FALSE], 1L, shift, "+")
      bounds[points, ] <- column_quantiles(combined)
    }
    array(as.numeric(bounds), c(G, G, 2L), dimnames = dimnames(intervals$fcont))
  }
  intervals$fcont <- surface(rep(0, B))
  intervals$fcont_full <- surface(intercepts)
  if (regional) {
    for (r in seq_along(region_levels)) {
      intervals$f[[r]] <- surface(spat[, r])
      intervals$f_full[[r]] <- surface(intercepts + spat[, r])
    }
  } else {
    intervals$f <- intervals$fcont
    intervals$f_full <- intervals$fcont_full
  }
  intervals
}


# Results on trainining observations ----
#' Summarize fitted effects and the complete predictor from one fitted cell.
#'
#' Reads stored coefficients/designs and never uses true DGP effects to adjust
#' predictions. With B_M3=0 it does not refit models or consume RNG draws.
#' By default, 100 bootstrap resamples refit M3 using its selected method. Despite the
#' historical name, point estimates are returned on training data AND a grid.
#'
#' @param s A run_cell() result with scenario, fits, and designs. The latter two
#'   lists must have matching model names/order. M0/M1/M2/M4/M5 use posterior draws;
#'   M3 uses the fitted coefficient vector in beta.
#' @param xval_l,xval_u Endpoints of the common x1/x2 grid.
#' @param xval_N Number G of grid points per covariate, giving G^2 combinations.
#' @param probs Lower/upper quantile probabilities. CI_quant uses these two
#'   quantiles; CI_hpd uses probability mass probs[2] - probs[1]. Intervals are
#'   pointwise intervals for fitted functions, not simultaneous bands or
#'   predictive intervals for a new response (no observation noise is added).
#' @param B_M3 Number of M3 bootstrap resamples (default 100); 0 disables refitting.
#' @param boot_seed Seed set once before the bootstrap; defaults to the saved
#'   scenario seed. Explicit NULL uses the current RNG stream. Bootstrap consumes
#'   RNG draws without restoring them.
#' @param boot_chunk_size Maximum grid combinations summarized at once for M3 bootstrap intervals.
#' @param grid_chunk_size Maximum grid combinations summarized at once for
#'   posterior surfaces (default 512). Every location still uses all draws.
#'   Smaller batches use less temporary memory; returned summaries are unchanged.
#' @param verbose Print concise extraction start/completion messages.
#' @param engine_output Show routine engine/bootstrap prints (default FALSE).
#'   Warnings and errors remain visible in either mode.
#' @return A list indexed by model (m0, ..., m5 as present in s). Each contains:
#'   - point_est: int (scalar), and observation-level effect vectors.
#'   - point_est_grid: int (scalar); f1/f2 and their available L/NL components
#'     (length G); spat (one value per region); fcont (G-by-G matrix); f/f_full
#'     (G-by-G matrices, or lists of regional matrices when spatial).
#'   - CI_quant and CI_hpd: intervals for the GRID summaries and the intercept.
#'     For f_full these have dimensions G-by-G-by-2, or are lists of such
#'     arrays by region; the third dimension contains lower/upper bounds.
#'     M3 has only CI_quant: bootstrap intervals, or NA when B_M3=0.
#'   Additional fields: point_est$slope1/slope2 are scalar slopes in original
#'   x units (M0/M2/M3/M4/M5), with matching interval entries. fcont_full means
#'   int+f1+f2, excluding any spatial effect, on training points and the grid.
#'   grid records the exact covariate grid used for all grid summaries.
#'   Absent effect components are NULL. Legacy nonspatial M4 fits in X3 retain
#'   a single common surface rather than acquiring regional coefficients.
#'   If M4 selection failed, its result is a named m4=NULL entry; other models
#'   are extracted normally. Raw family draws remain in s$fits$m4$losing_draws.
#'
#' @details The effect names mean:
#'   int = fitted intercept; f1/f2 = fitted centered covariate effects;
#'   f1L/f2L and f1NL/f2NL = linear/nonlinear components where available;
#'   spat = fitted regional effect (equal-region centering);
#'   fcont = f1 + f2; f = fcont + spat (omit spat when absent);
#'   f_full = int + f, the complete fitted predictor.
#'   Existing f remains intercept-free for compatibility. M1/M4 spatial
#'   recentering shifts their intercepts by the compensating amount, preserving
#'   complete predictions. With unequal region counts, f need not have zero
#'   observation mean even though its components obey their constraints.
#'
#'   For posterior models, point estimates are posterior means. f_full grid
#'   intervals are computed after adding each draw's intercept to that SAME
#'   draw's effect sum, retaining uncertainty and dependence. M3 bootstrap
#'   intervals use the same principle across resamples. Previous results do not acquire
#'   f_full automatically: rerun this function on their saved run_cell() object.
#' @examples
#' fitted <- results_on_train(s)
#' truth_train <- training_truth(s$scenario)
#' # Whole-function comparison; the structural signal is the chosen benchmark:
#' fitted$m0$point_est$f_full - truth_train$f_full
#' # For centered effect comparisons, use f1/f2/spat/f from aligned truth.
results_on_train <- function(s, xval_l = -0.5, xval_u = 0.5, xval_N = 100,
                             probs = c(0.025, 0.975), B_M3 = 100L, boot_seed = s$scenario$seed,
                             boot_chunk_size = 512L, grid_chunk_size = 512L,
                             verbose = TRUE, engine_output = FALSE) {
  label <- if ("m3" %in% names(s$fits) && !is.null(s$fits$m3) && B_M3 > 0L)
    sprintf("Extracting effects and intervals (M3 bootstrap B=%d)", B_M3) else
      "Extracting effects and intervals"
  run_stage(label, {
    stopifnot(length(grid_chunk_size) == 1L, is.finite(grid_chunk_size),
              grid_chunk_size >= 1, grid_chunk_size == floor(grid_chunk_size))

    stopifnot("Models in fits should match those in designs. The names of the two objects do not coincide." =
                all(names(s$fits) == names(s$designs)))

    # Extract main objects
    pred <- s$scenario$predictor
    dat <- s$scenario$dat
    des <- s$designs
    models <- names(s$fits)

    # Helpers for spatial effects
    region_levels <- levels(dat$region)
    region_id <- match(as.character(dat$region), region_levels)
    ## One observation representing each region, in factor-level order.
    first_row <- match(region_levels, as.character(dat$region))
    if (anyNA(first_row)) {stop("Regional extraction requires every region to be observed.")}

    # Define values for the fit on the grid (will be used both for x1 and x2)
    new_x <- seq(xval_l, xval_u, length = xval_N)
    ## Helpers to compute surface-type effects, since f(x1, x2) is defined on R^2
    ### Ordering agrees with as.vector(outer(f1, f2, "+")): x1 varies fastest, then x2.
    grid_i1 <- rep(seq_len(xval_N), times = xval_N)
    grid_i2 <- rep(seq_len(xval_N), each = xval_N)
    grid_names <- list(x1 = as.character(new_x), x2 = as.character(new_x))
    ### Convert xval_N^2 values into an x1-by-x2 surface.
    as_surface <- function(z) {
      stopifnot(length(z) == xval_N * xval_N)
      matrix(z, nrow = xval_N, ncol = xval_N, dimnames = grid_names)
    }
    ### Convert a (xval_N^2)-by-2 CI matrix into a surface of intervals.
    as_interval_surface <- function(z) {
      stopifnot(nrow(z) == xval_N * xval_N, ncol(z) == 2L)
      array(as.numeric(z), dim = c(xval_N, xval_N, 2L), dimnames = c(grid_names, list(bound = c("LB", "UP"))))
    }

    # Define nominal level from probs
    cred_mass = probs[2] - probs[1]

    # Keep summaries for the whole grid, but detailed draws for only one batch.
    # A batch contains grid locations, not a subset of posterior draws.
    summarize_surfaces <- function(f1, f2, intercept, spatial = NULL) {
      empty_summary <- function() list(point_est = numeric(xval_N^2),
        CI_quant = matrix(NA_real_, xval_N^2, 2L), CI_hpd = matrix(NA_real_, xval_N^2, 2L))
      save_batch <- function(out, D, points) {
        out$point_est[points] <- colMeans(D)
        out$CI_quant[points, ] <- t(apply(D, 2L, quantile, probs = probs))
        out$CI_hpd[points, ] <- coda::HPDinterval(coda::as.mcmc(D), prob = cred_mass)
        out
      }
      finish <- function(out) list(point_est = as_surface(out$point_est),
        CI_quant = as_interval_surface(out$CI_quant), CI_hpd = as_interval_surface(out$CI_hpd))
      finish_regions <- function(out) {
        out <- lapply(out, finish)
        list(point_est = lapply(out, `[[`, "point_est"),
             CI_quant = lapply(out, `[[`, "CI_quant"), CI_hpd = lapply(out, `[[`, "CI_hpd"))
      }
      continuous <- continuous_full <- empty_summary()
      n_regions <- if (is.null(spatial)) 1L else ncol(spatial)
      regional <- full <- replicate(n_regions, empty_summary(), simplify = FALSE)
      names(regional) <- names(full) <- colnames(spatial)

      for (first in seq.int(1L, xval_N^2, by = grid_chunk_size)) {
        points <- seq.int(first, min(xval_N^2, first + grid_chunk_size - 1L))
        # Compute the covariate sum once per batch and reuse it in every region.
        # Keep x1 varying fastest, as in the original complete surface matrix.
        D_cont <- f1[, grid_i1[points], drop = FALSE] + f2[, grid_i2[points], drop = FALSE]
        continuous <- save_batch(continuous, D_cont, points)
        if (!is.null(spatial)) {
          continuous_full <- save_batch(continuous_full, sweep(D_cont, 1L, intercept, "+"), points)
        }
        for (r in seq_len(n_regions)) {
          if (is.null(spatial)) {
            D_region <- D_cont
          } else {
            # Add region r's effect from the SAME draw to each grid location.
            D_region <- D_cont + spatial[, r]
            regional[[r]] <- save_batch(regional[[r]], D_region, points)
          }
          # Add each draw's intercept before intervals, retaining its dependence
          # with the other effects. Only the small batch is held in memory.
          D_full <- sweep(D_region, 1L, intercept, "+")
          full[[r]] <- save_batch(full[[r]], D_full, points)
        }
      }
      continuous <- finish(continuous)
      list(fcont = continuous,
           fcont_full = if (is.null(spatial)) finish(full[[1L]]) else finish(continuous_full),
           f = if (is.null(spatial)) continuous else finish_regions(regional),
           f_full = if (is.null(spatial)) finish(full[[1L]]) else finish_regions(full))
    }

    # Container for results
    res <- list()



    for (mod in models){
      if (mod == "m4" && is.null(s$fits[[mod]]$draws)) {
        res[mod] <- list(NULL)
        next
      }

      if (mod == "m3"){
        beta <- s$fits[[mod]]$beta
        draws <- matrix(beta, ncol = length(beta))
        colnames(draws) <- names(beta)
        draws
      } else {
        draws <- s$fits[[mod]]$draws
      }

      X <- s$designs[[mod]]$X

      new_xtilde1 <- (new_x - mean(dat[ , "x1"]))/sd(dat[ , "x1"])
      new_xtilde2 <- (new_x - mean(dat[ , "x2"]))/sd(dat[ , "x2"])

      spat_idx <- NULL
      if (pred == "X3") {
        spat_idx <- grep(if (mod == "m4") "^region" else "^gmrf", colnames(draws))
      }

      # M0/M2/M3/M5 use equal-region centering: sum_r gamma_r = 0.
      # M1 uses observation-weighted centering; M4 uses a reference region.
      # Recenter their regional effects and shift each draw's intercept so all
      # fitted predictions are preserved. Legacy M4 fits without regional
      # coefficients retain their previous, nonspatial extraction.
      if (pred == "X3" && (mod == "m1" || (mod == "m4" && length(spat_idx) > 0L))) {
        spatial_idx <- spat_idx
        intercept_idx <- match("Intercept", colnames(draws))
        stopifnot(length(spatial_idx) > 0L, !is.na(intercept_idx), !anyNA(first_row))

        # Mean regional basis row; each draw's shift is b %*% a.
        a <- colMeans(X[first_row, spatial_idx, drop = FALSE])
        shift <- drop(draws[, spatial_idx, drop = FALSE] %*% a)

        draws[, intercept_idx] <- draws[, intercept_idx] + shift
        X[, spatial_idx] <- sweep(X[, spatial_idx, drop = FALSE], 2L, a, "-")
      }

      if (mod %in% c("m0", "m5")){
        dr1 <- des[[mod]]$dr1
        dr2 <- des[[mod]]$dr2

        new_B1 <- Predict.matrix(dr1$smooth_object, data.frame(x = new_x))
        new_B2 <- Predict.matrix(dr2$smooth_object, data.frame(x = new_x))

        new_Z1 <- new_B1 %*% dr1$Trafo
        new_Z2 <- new_B2 %*% dr2$Trafo

        new_X <- cbind(1, new_xtilde1, new_xtilde2, new_Z1, new_Z2)
      } else if (mod == "m1"){
        sm <- des[[mod]]$full_info$mu$smooth.construct
        sm1 <- sm$`s(x1)`
        sm2 <- sm$`s(x2)`

        new_TP1 <- PredictMat(sm1, setNames(data.frame(new_x), sm1$term))
        new_TP2 <- PredictMat(sm2, setNames(data.frame(new_x), sm2$term))

        new_X <- cbind(1, new_TP1, new_TP2)
      } else if (mod == "m2" | mod == "m4") {
        new_X <- cbind(1, new_xtilde1, new_xtilde2)
      } else if (mod == "m3"){

        bs1 <- des[[mod]]$bs1
        bs2 <- des[[mod]]$bs2

        new_B1 <- PredictMat(bs1, data.frame(x = new_x))
        new_B2 <- PredictMat(bs2, data.frame(x = new_x))

        new_X <- cbind(1, new_B1, new_B2)
      }


      int_idx <- 1
      x1L_idx <- which(colnames(draws) == "x1")
      x1NL_idx <- grep("x1.s", colnames(draws))
      x2L_idx <- which(colnames(draws) == "x2")
      x2NL_idx <- grep("x2.s", colnames(draws))
      effects_idx <- list(int = int_idx,
                          f1L = x1L_idx, f1NL = x1NL_idx, f1 = c(x1L_idx, x1NL_idx),
                          f2L = x2L_idx, f2NL = x2NL_idx, f2 = c(x2L_idx, x2NL_idx),
                          spat = spat_idx)
      effects_idx[["fcont"]] <- c(effects_idx$f1, effects_idx$f2)
      effects_idx[["f"]] <- c(effects_idx$fcont, effects_idx$spat)

      effects_idx <- lapply(effects_idx, function(x) if (length(x) == 0) NULL else x)

      if (mod != "m3"){

        point_est <- point_est_grid <- CI_quant <- CI_hpd <- vector(mode = "list", length = length(effects_idx))
        names(point_est) <- names(point_est_grid) <- names(CI_quant) <- names(CI_hpd) <- names(effects_idx)

        D_eff <- D_eff_grid <- vector(mode = "list", length = length(effects_idx))
        names(D_eff) <- names(D_eff_grid) <- names(effects_idx)

        for (E in seq_along(effects_idx)){

          eff <- names(effects_idx)[E]
          eff_idx <- effects_idx[[E]]

          if (is.null(eff_idx)){
            NULL
          } else {
            if (eff == "int"){
              # Compute
              D <- draws[ , 1]
              D_grid <- draws[ , 1]
              pe <- mean(D)
              pe_grid <- mean(D_grid)
              quant <- quantile(D_grid, probs = probs)
              hpd <- HPDinterval(as.mcmc(D_grid), prob = cred_mass)
              names(quant) <- names(hpd) <- c("LB", "UP")

              # Save
              D_eff[[eff]] <- D; D_eff_grid[[eff]] <- D_grid
              point_est[[eff]] <- unname(pe); point_est_grid[[eff]] <- unname(pe_grid)
              CI_quant[[eff]] <- quant
              CI_hpd[[eff]] <- hpd
            } else if (eff %in% c("f1L", "f1NL", "f1", "f2L", "f2NL", "f2")){
              # Compute
              D <- draws[ , eff_idx] %*% t(X[ , eff_idx])
              D_grid <- draws[ , eff_idx] %*% t(new_X[ , eff_idx])
              pe <- colMeans(D)
              pe_grid <- colMeans(D_grid)
              quant <- t(apply(D_grid, 2, quantile, probs = probs))
              hpd <- HPDinterval(as.mcmc(D_grid), prob = cred_mass)
              colnames(quant) <- colnames(hpd) <- c("LB", "UP")

              # Save
              D_eff[[eff]] <- D; D_eff_grid[[eff]] <- D_grid
              point_est[[eff]] <- unname(pe); point_est_grid[[eff]] <- unname(pe_grid)
              CI_quant[[eff]] <- quant
              CI_hpd[[eff]] <- hpd
            } else if (eff == "fcont"){
              # Training draws stay as before. Grid surfaces are summarized in
              # batches after the component loop, once regional draws are ready.
              D <- D_eff$f1 + D_eff$f2
              pe <- point_est$f1 + point_est$f2
              D_eff[[eff]] <- D
              point_est[[eff]] <- unname(pe)
            } else if (eff == "spat"){
              # Compute
              V <- X[, eff_idx, drop = FALSE]
              ## Extract one row for each region, in factor-level (i.e. region) order
              V_region <- V[first_row, , drop = FALSE]
              rownames(V_region) <- region_levels
              D_region <- draws[, eff_idx, drop = FALSE]
              D <- D_region %*% t(V)
              D_grid <- D_region %*% t(V_region)
              colnames(D_grid) <- region_levels

              # Couple of checks
              ## Every observation must use its region's design row.
              stopifnot(isTRUE(all.equal(unname(V), unname(V_region[region_id, , drop = FALSE]))))
              ## Regional reconstruction must reproduce observation-level effects.
              stopifnot(isTRUE(all.equal(unname(D), unname(D_grid[, region_id, drop = FALSE]))))

              pe <- colMeans(D)
              pe_grid <- colMeans(D_grid)
              quant <- t(apply(D_grid, 2, quantile, probs = probs))
              hpd <- HPDinterval(as.mcmc(D_grid), prob = cred_mass)
              colnames(quant) <- colnames(hpd) <- c("LB", "UP")

              # Save
              D_eff[[eff]] <- D; D_eff_grid[[eff]] <- D_grid
              point_est[[eff]] <- unname(pe); point_est_grid[[eff]] <- unname(pe_grid)
              CI_quant[[eff]] <- quant
              CI_hpd[[eff]] <- hpd
            } else if (eff == "f"){
              if (is.null(effects_idx$spat)){
                D_eff[[eff]] <- D_eff[["fcont"]]
                point_est[[eff]] <- point_est[["fcont"]]
              } else {
                D <- D_eff$spat + D_eff$fcont
                pe <- point_est$fcont + point_est$spat
                D_eff[[eff]] <- D
                point_est[[eff]] <- unname(pe)
              }
            }
          }
        }

        # Retain only the compact surface summaries, not a draw matrix for every
        # region and grid location. Curve and regional component draws are small.
        surfaces <- summarize_surfaces(D_eff_grid$f1, D_eff_grid$f2, D_eff$int, D_eff_grid$spat)
        for (eff in names(surfaces)) {
          point_est_grid[[eff]] <- surfaces[[eff]]$point_est
          CI_quant[[eff]] <- surfaces[[eff]]$CI_quant
          CI_hpd[[eff]] <- surfaces[[eff]]$CI_hpd
        }
        res[[mod]] <- list(point_est = point_est, point_est_grid = point_est_grid,
                           CI_quant = CI_quant, CI_hpd = CI_hpd)

      } else if (mod == "m3") {
        point_est <- point_est_grid <- CI_quant <- setNames(
          vector("list", length(effects_idx)), names(effects_idx))
        boot_out <- m3_na_intervals(new_x, if (pred == "X3") region_levels else NULL)
        if (B_M3 > 0L) {
          boot_out <- boot_CI(B = B_M3, dat = dat, predictor = pred, new_x = new_x,
            method = s$fits$m3$method, probs = probs, seed = boot_seed, chunk_size = boot_chunk_size)
        }
        # Bootstrap intervals are calculated once, outside this point-estimate loop.
        for (eff in names(effects_idx)) {
          idx <- effects_idx[[eff]]
          if (is.null(idx)) next
          if (eff == "int") {
            point_est[[eff]] <- point_est_grid[[eff]] <- unname(beta[1L])
          } else if (eff %in% c("f1NL", "f1", "f2NL", "f2")) {
            point_est[[eff]] <- unname(drop(X[, idx, drop = FALSE] %*% beta[idx]))
            point_est_grid[[eff]] <- unname(drop(new_X[, idx, drop = FALSE] %*% beta[idx]))
          } else if (eff == "spat") {
            point_est[[eff]] <- unname(drop(X[, idx, drop = FALSE] %*% beta[idx]))
            point_est_grid[[eff]] <- unname(drop(des[[mod]]$region_contrasts %*% beta[idx]))
          } else if (eff == "fcont") {
            point_est[[eff]] <- point_est$f1 + point_est$f2
            point_est_grid[[eff]] <- outer(point_est_grid$f1, point_est_grid$f2, "+")
          } else if (eff == "f") {
            point_est[[eff]] <- point_est$fcont + if (is.null(effects_idx$spat)) 0 else point_est$spat
            point_est_grid[[eff]] <- if (is.null(effects_idx$spat)) point_est_grid$fcont else
              lapply(point_est_grid$spat, function(gamma_r) {
                surface <- point_est_grid$fcont + gamma_r
                dimnames(surface) <- grid_names
                surface
              })
          }
          CI_quant[eff] <- boot_out[eff]
        }
        res[[mod]] <- list(point_est = point_est, point_est_grid = point_est_grid,
                           CI_quant = CI_quant)
      }

      # f_full is an additional field: existing f remains intercept-free.
      # Posterior means add linearly, so the training mean can be reconstructed
      # as mean(intercept) + mean(f). This shortcut is NOT valid for intervals.
      res[[mod]]$point_est$f_full <- point_est$int + point_est$f
      res[[mod]]$point_est$fcont_full <- point_est$int + point_est$fcont
      if (mod == "m3") {
        res[[mod]]$point_est_grid$fcont_full <- point_est_grid$int + point_est_grid$fcont
        res[[mod]]$CI_quant$fcont_full <- boot_out$fcont_full
        res[[mod]]$point_est_grid$f_full <- if (is.list(point_est_grid$f)) {
          lapply(point_est_grid$f, function(z) point_est_grid$int + z)
        } else {
          point_est_grid$int + point_est_grid$f
        }
        # f_full intervals include each bootstrap intercept before quantiles.
        res[[mod]]$CI_quant$f_full <- boot_out$f_full
      }

      # Empirical affine projection gives a scalar slope in original x units.
      # For M0/M2/M4/M5 this equals coefficient / sd(x). For M3 it separates the
      # fitted spline without changing it. M1 has no requested slope summary.
      if (mod != "m1") for (j in 1:2) {
        x <- dat[[paste0("x", j)]]
        weights <- (x - mean(x)) / sum((x - mean(x))^2)
        key <- paste0("slope", j)
        effect <- paste0("f", j)
        if (mod == "m3") {
          res[[mod]]$point_est[[key]] <- sum(point_est[[effect]] * weights)
          res[[mod]]$CI_quant[[key]] <- boot_out[[key]]
        } else {
          slope <- drop(D_eff[[effect]] %*% weights)
          res[[mod]]$point_est[[key]] <- mean(slope)
          res[[mod]]$CI_quant[[key]] <- setNames(as.numeric(quantile(slope, probs)), c("LB", "UP"))
          res[[mod]]$CI_hpd[[key]] <- setNames(as.numeric(HPDinterval(as.mcmc(slope), prob = cred_mass)), c("LB", "UP"))
        }
      }

      res[[mod]]$grid <- new_x
      res[[mod]]$interval_probs <- unname(probs)
      pg <- res[[mod]]$point_est_grid
      stopifnot(isTRUE(all.equal(unname(pg$fcont), unname(outer(pg$f1, pg$f2, "+")))))

      if (!is.null(effects_idx$spat)) {
        stopifnot(length(pg$f) == length(pg$spat))

        for (r in seq_along(pg$spat)) {
          stopifnot(isTRUE(all.equal(unname(pg$f[[r]]), unname(pg$fcont + pg$spat[r]))))
        }
      }

    }

    res
  }, verbose, engine_output)
}
