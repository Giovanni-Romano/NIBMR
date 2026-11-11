# === === === === === === === === === === === === === === === === === === ===
# 04-simulation/03-extract.R
# Load saved per-replicate `s`, apply effect_intervals (B1) and test_eval (B2),
# aggregate over replicates, write result tables.
# === === === === === === === === === === === === === === === === === === ===

library(coda)

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

grid_truth <- function(pred, new_x, sp_eff = NULL) {
  ft <- f_true(pred, new_x, new_x)
  f1 <- ft$f1
  f2 <- ft$f2
  fcont <- outer(f1, f2, "+")
  if (pred == "X3"){
    f <- lapply(sp_eff, function(S) S + fcont)
  } else {
    f <- fcont
  }
  list(f1 = ft$f1, f2 = ft$f2, fcont = fcont, f = f)
}

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
    
    # The smooth constructor for GMRF in mgcv uses a constraint on the spatial
    # effect that differs from the one used in M0/M2. This mines the comparison
    # of the estimates across models.
    # Define "gamma_r" and "n_r" as the effect and the number of observations in region r.
    # M0/M2 use equal-region centering, meaning that sum_{r = 1}^s gamma_r = 0.
    # M1 uses weighted centering, meaning that       sum_{r = 1}^s n_r gamma_r = 0.
    # Here we express M1 spatial effects under equal-region centering.
    if (mod == "m1" && pred == "X3") {
      spatial_idx <- grep("^gmrf", colnames(draws))
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
    spat_idx <- NULL
    if (pred == "X3" & mod != "m4"){
      spat_idx <- grep("gmrf", colnames(draws))
    }
    
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
              D_eff[[eff]] <- D_eff[["fcont"]]; D_eff_grid[[eff]] <- D_eff[["fcont"]]
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
              D_eff[[eff]] <- D_eff[["fcont"]]; D_eff_grid[[eff]] <- D_eff[["fcont"]]
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