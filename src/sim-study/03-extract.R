# === === === === === === === === === === === === === === === === === === ===
# src/sim-study/03-extract.R
# Reconstruct fitted effects and the complete predictor from one run_cell()
# result, and construct comparable simulation truth on observations or a grid.
# This script defines extraction helpers; it does not fit models, aggregate
# simulation replicates, or write result files.
# === === === === === === === === === === === === === === === === === === ===

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


boot_CI <- function(){matrix(c(NA, NA), ncol = 2)}


# Results on trainining observations ----
#' Summarize fitted effects and the complete predictor from one fitted cell.
#'
#' This reads stored coefficients and designs; it does not refit models or use
#' the true DGP effects to adjust predictions. Despite the historical function
#' name, it returns point estimates on BOTH training observations and a grid.
#'
#' @param s A run_cell() result with scenario, fits, and designs. The latter two
#'   lists must have matching model names/order. M0/M1/M2/M4 use posterior draws;
#'   M3 uses the fitted coefficient vector in beta.
#' @param xval_l,xval_u Endpoints of the common x1/x2 grid.
#' @param xval_N Number G of grid points per covariate, giving G^2 combinations.
#' @param probs Lower/upper quantile probabilities. CI_quant uses these two
#'   quantiles; CI_hpd uses probability mass probs[2] - probs[1]. Intervals are
#'   pointwise intervals for fitted functions, not simultaneous bands or
#'   predictive intervals for a new response (no observation noise is added).
#' @return A list indexed by model (m0, ..., m4 as present in s). Each contains:
#'   - point_est: int (scalar), and observation-level effect vectors.
#'   - point_est_grid: int (scalar); f1/f2 and their available L/NL components
#'     (length G); spat (one value per region); fcont (G-by-G matrix); f/f_full
#'     (G-by-G matrices, or lists of regional matrices when spatial).
#'   - CI_quant and CI_hpd: intervals for the GRID summaries and the intercept.
#'     For f_full these have dimensions G-by-G-by-2, or are lists of such
#'     arrays by region; the third dimension contains lower/upper bounds.
#'     M3 has only CI_quant, with NA placeholders pending bootstrap inference.
#'   Absent effect components are NULL. Legacy nonspatial M4 fits in X3 retain
#'   a single common surface rather than acquiring regional coefficients.
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
#'   draw's effect sum, retaining uncertainty and dependence. For M3 only the
#'   point estimate is available. Previous extracted results do not acquire
#'   f_full automatically: rerun this function on their saved run_cell() object.
#' @examples
#' fitted <- results_on_train(s)
#' truth_train <- training_truth(s$scenario)
#' # Whole-function comparison; the structural signal is the chosen benchmark:
#' fitted$m0$point_est$f_full - truth_train$f_full
#' # For centered effect comparisons, use f1/f2/spat/f from aligned truth.
results_on_train <- function(s, xval_l = -0.5, xval_u = 0.5, xval_N = 100, probs = c(0.025, 0.975)) {
  
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
  
  # Container for results
  res <- list()
  
  
  
  for (mod in models){
    
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

    # M0/M2/M3 use equal-region centering: sum_r gamma_r = 0.
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
    
    if (mod == "m0"){
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
            # Compute
            D <- D_eff$f1 + D_eff$f2
            D_grid <- D_eff_grid$f1[, grid_i1, drop = FALSE] + D_eff_grid$f2[, grid_i2, drop = FALSE]
            pe <- point_est$f1 + point_est$f2
            pe_grid <- as_surface(colMeans(D_grid))
            quant <- as_interval_surface(t(apply(D_grid, 2, quantile, probs = probs)))
            hpd <- as_interval_surface(HPDinterval(as.mcmc(D_grid), prob = cred_mass))
            
            # Save
            D_eff[[eff]] <- D; D_eff_grid[[eff]] <- D_grid
            point_est[[eff]] <- unname(pe); point_est_grid[[eff]] <- pe_grid
            CI_quant[[eff]] <- quant
            CI_hpd[[eff]] <- hpd
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
            
            counts <- tabulate(region_id, nbins = length(region_levels))
            
            centering_check <- c(
              equal_region = max(abs(colMeans(V_region))),
              observation_weighted = max(abs(
                drop(crossprod(counts / sum(counts), V_region))
              ))
            )
            
            print(list(model = mod, centering = centering_check))
            
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
              D_eff[[eff]] <- D_eff[["fcont"]]; D_eff_grid[[eff]] <- D_eff_grid[["fcont"]]
              point_est[[eff]] <- point_est[["fcont"]]; point_est_grid[[eff]] <- point_est_grid[["fcont"]]
              CI_quant[[eff]] <- CI_quant[["fcont"]]
              CI_hpd[[eff]] <- CI_hpd[["fcont"]]
            } else {
              # Compute
              D <- D_eff$spat + D_eff$fcont
              D_grid <- apply(D_eff_grid$spat, 2L, function(x) x + D_eff_grid$fcont, simplify = FALSE)
              pe <- point_est$fcont + point_est$spat
              pe_grid <- lapply(D_grid, function(z) {as_surface(colMeans(z))})
              quant <- lapply(D_grid, function(z) {as_interval_surface(t(apply(z, 2L, quantile, probs = probs)))})
              hpd <- lapply(D_grid, function(z) {as_interval_surface(coda::HPDinterval(coda::as.mcmc(z), prob = cred_mass))})
              
              # Save
              D_eff[[eff]] <- D; D_eff_grid[[eff]] <- D_grid
              point_est[[eff]] <- unname(pe); point_est_grid[[eff]] <- pe_grid
              CI_quant[[eff]] <- quant
              CI_hpd[[eff]] <- hpd
            }
          }
        }
      }
      
      res[[mod]] <- list(point_est = point_est, point_est_grid = point_est_grid,
                         CI_quant = CI_quant, CI_hpd = CI_hpd)
      
    } else if (mod == "m3") {
      
      point_est <- point_est_grid <- CI_quant <- vector(mode = "list", length = length(effects_idx))
      names(point_est) <- names(point_est_grid) <- names(CI_quant) <- names(CI_quant) <- names(effects_idx)
      
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
            D <- drop(draws[ , 1])
            D_grid <- drop(draws[ , 1])
            pe <- D
            pe_grid <- D_grid
            quant <- boot_CI() # still to be implemented
            names(quant) <- c("LB", "UP")
            
            # Save
            D_eff[[eff]] <- D; D_eff_grid[[eff]] <- D_grid
            point_est[[eff]] <- unname(pe); point_est_grid[[eff]] <- unname(pe_grid)
            CI_quant[[eff]] <- quant
          } else if (eff %in% c("f1NL", "f1", "f2NL", "f2")){ # there are no explicit linear components in this model
            # Compute
            D <- drop(draws[ , eff_idx] %*% t(X[ , eff_idx]))
            D_grid <- drop(draws[ , eff_idx] %*% t(new_X[ , eff_idx]))
            pe <- D
            pe_grid <- D_grid
            quant <- boot_CI() # still to be implemented
            colnames(quant) <- c("LB", "UP")
            
            # Save
            D_eff[[eff]] <- D; D_eff_grid[[eff]] <- D_grid
            point_est[[eff]] <- unname(pe); point_est_grid[[eff]] <- unname(pe_grid)
            CI_quant[[eff]] <- quant
          } else if (eff == "fcont"){
            # Compute
            D <- D_eff$f1 + D_eff$f2
            D_grid <- outer(D_eff_grid$f1, D_eff_grid$f2, "+")
            dimnames(D_grid) <- grid_names
            pe <- point_est$f1 + point_est$f2
            pe_grid <- D_grid
            # Placeholder until bootstrap intervals are implemented.
            quant <- array(NA_real_, dim = c(xval_N, xval_N, 2L), dimnames = c(grid_names, list(bound = c("LB", "UP"))))
            
            # Save
            D_eff[[eff]] <- D; D_eff_grid[[eff]] <- D_grid
            point_est[[eff]] <- unname(pe); point_est_grid[[eff]] <- unname(pe_grid)
            CI_quant[[eff]] <- quant
          } else if (eff == "spat"){
            # Compute
            D <- drop(draws[ , eff_idx] %*% t(X[ , eff_idx]))
            D_grid <- drop(des[[mod]]$region_contrasts %*% as.numeric(draws[, eff_idx]))
            pe <- D
            pe_grid <- D_grid
            quant <- boot_CI() # still to be implemented
            colnames(quant) <- c("LB", "UP")
            
            # Save
            D_eff[[eff]] <- D; D_eff_grid[[eff]] <- D_grid
            point_est[[eff]] <- unname(pe); point_est_grid[[eff]] <- unname(pe_grid)
            CI_quant[[eff]] <- quant
          } else if (eff == "f"){
            if (is.null(effects_idx$spat)){
              D_eff[[eff]] <- D_eff[["fcont"]]; D_eff_grid[[eff]] <- D_eff_grid[["fcont"]]
              point_est[[eff]] <- point_est[["fcont"]]; point_est_grid[[eff]] <- point_est_grid[["fcont"]]
              CI_quant[[eff]] <- CI_quant[["fcont"]]
            } else {
              # Compute
              D <- D_eff$spat + D_eff$fcont
              D_grid <- lapply(D_eff_grid$spat, function(gamma_r) {D_eff_grid$fcont + gamma_r})
              pe <- point_est$fcont + point_est$spat
              pe_grid <- D_grid
              # One placeholder interval surface per region.
              quant <- lapply(D_grid, function(z) {
                array(NA_real_, dim = c(xval_N, xval_N, 2L), 
                      dimnames = c(grid_names, list(bound = c("LB", "UP"))))})
              
              # Save
              D_eff[[eff]] <- D; D_eff_grid[[eff]] <- D_grid
              point_est[[eff]] <- unname(pe); point_est_grid[[eff]] <- unname(pe_grid)
              CI_quant[[eff]] <- quant
            }
          }
        }
      }
      
      res[[mod]] <- list(point_est = point_est, point_est_grid = point_est_grid,
                         CI_quant = CI_quant)
    }
    
    # f_full is an additional field: existing f remains intercept-free.
    # Posterior means add linearly, so the training mean can be reconstructed
    # as mean(intercept) + mean(f). This shortcut is NOT valid for intervals.
    res[[mod]]$point_est$f_full <- point_est$int + point_est$f
    if (mod == "m3") {
      res[[mod]]$point_est_grid$f_full <- if (is.list(point_est_grid$f)) {
        lapply(point_est_grid$f, function(z) point_est_grid$int + z)
      } else {
        point_est_grid$int + point_est_grid$f
      }
      # CI_quant$f consists of NA placeholders here. Reuse their shapes for
      # f_full; this does not treat the intercept as known for inference.
      res[[mod]]$CI_quant$f_full <- CI_quant$f
    } else {
      summarize_full <- function(centered_draws) {
        # Rows = posterior draws; columns = xval_N^2 grid combinations. sweep(...,
        # MARGIN = 1) adds intercept[b] to every grid value in draw b, preserving
        # its posterior dependence with the other effects. Summarize columns
        # only after this addition, then restore x1-by-x2 surface dimensions.
        full_draws <- sweep(centered_draws, 1L, D_eff$int, "+")
        list(point_est = as_surface(colMeans(full_draws)),
             CI_quant = as_interval_surface(t(apply(full_draws, 2L, quantile, probs = probs))),
             CI_hpd = as_interval_surface(coda::HPDinterval(coda::as.mcmc(full_draws), prob = cred_mass)))
      }
      if (is.list(D_eff_grid$f)) {
        # Summarize one region at a time instead of keeping another full set
        # of regional draw matrices in memory.
        full <- lapply(D_eff_grid$f, summarize_full)
        res[[mod]]$point_est_grid$f_full <- lapply(full, `[[`, "point_est")
        res[[mod]]$CI_quant$f_full <- lapply(full, `[[`, "CI_quant")
        res[[mod]]$CI_hpd$f_full <- lapply(full, `[[`, "CI_hpd")
      } else {
        full <- summarize_full(D_eff_grid$f)
        res[[mod]]$point_est_grid$f_full <- full$point_est
        res[[mod]]$CI_quant$f_full <- full$CI_quant
        res[[mod]]$CI_hpd$f_full <- full$CI_hpd
      }
    }

    pg <- res[[mod]]$point_est_grid
    stopifnot(isTRUE(all.equal(unname(pg$fcont), unname(outer(pg$f1, pg$f2, "+")))))
    
    if (!is.null(effects_idx$spat)) {
      stopifnot(length(pg$f) == length(pg$spat))
      
      for (r in seq_along(pg$spat)) {
        stopifnot(isTRUE(all.equal(unname(pg$f[[r]]), unname(pg$fcont + pg$spat[r]))))
      }
    }
    
  }
  
  return(res)
  
}
