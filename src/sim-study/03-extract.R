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
  
  pred <- s$scenario$predictor
  dat <- s$scenario$dat
  des <- s$designs
  models <- names(s$fits)
  
  new_x <- seq(xval_l, xval_u, length = xval_N)
  
  res <- list()
  
  for (mod in models){
    
    if (mod != "m3"){
      draws <- s$fits[[mod]]$beta
    } else {
      draws <- s$fits[[mod]]$draws
    }
    
    X <- s$designs[[mod]]$X
    
    new_xtilde1 <- (new_x - mean(dat[ , "x1"]))/sd(dat[ , "x1"])
    new_xtilde2 <- (new_x - mean(dat[ , "x2"]))/sd(dat[ , "x2"])
    
    if (mod == "m0" | mod == "m3"){
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
    }
    
    
    int_idx <- 1
    x1L_idx <- which(colnames(draws) == "x1")
    x1NL_idx <- grep("x1.s", colnames(draws))
    x2L_idx <- which(colnames(draws) == "x2")
    x2NL_idx <- grep("x2.s", colnames(draws))
    spat_idx <- NULL
    if (pred == "X3" & mod != "m4"){
      spat_idx <- grep("gmrf.", colnames(draws))
    }
    
    effects_idx <- list(int = int_idx, 
                        f1L = x1L_idx, f1NL = x1NL_idx, f1 = c(x1L_idx, x1NL_idx),
                        f2L = x2L_idx, f2NL = x2NL_idx, f2 = c(x2L_idx, x2NL_idx),
                        spat = spat_idx)
    effects_idx[["fcont"]] <- c(effects_idx$int, 
                                effects_idx$x1, effects_idx$x2)
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
            hpd <- HPDinterval(as.mcmc(D_grid), prob = probs[1] + (1-probs[2]))
            names(quant) <- names(hpd) <- c("LB", "UP")
            
            # Save
            D_eff[[eff]] <- D; D_eff_grid[[eff]] <- D_grid
            point_est[[eff]] <- pe; point_est_grid[[eff]] <- pe_grid
            CI_quant[[eff]] <- quant
            CI_hpd[[eff]] <- hpd
          } else if (eff %in% c("f1L", "f1NL", "f1", "f2L", "f2NL", "f2")){
            # Compute
            D <- draws[ , eff_idx] %*% t(X[ , eff_idx])
            D_grid <- draws[ , eff_idx] %*% t(new_X[ , eff_idx])
            pe <- colMeans(D)
            pe_grid <- colMeans(D_grid)
            quant <- t(apply(D_grid, 2, quantile, probs = probs))
            hpd <- HPDinterval(as.mcmc(D_grid), prob = probs[1] + (1-probs[2]))
            colnames(quant) <- colnames(hpd) <- c("LB", "UP")
            
            # Save
            D_eff[[eff]] <- D; D_eff_grid[[eff]] <- D_grid
            point_est[[eff]] <- pe; point_est_grid[[eff]] <- pe_grid
            CI_quant[[eff]] <- quant
            CI_hpd[[eff]] <- hpd
          } else if (eff == "fcont"){
            # Compute
            D <- D_eff$f1 + D_eff$f2
            D_grid <- D_eff_grid$f1 + D_eff_grid$f2
            pe <- point_est$f1 + point_est$f2
            pe_grid <- point_est_grid$f1 + point_est_grid$f2
            quant <- t(apply(D_grid, 2, quantile, probs = probs))
            hpd <- HPDinterval(as.mcmc(D_grid), prob = probs[1] + (1-probs[2]))
            colnames(quant) <- colnames(hpd) <- c("LB", "UP")
            
            # Save
            D_eff[[eff]] <- D; D_eff_grid[[eff]] <- D_grid
            point_est[[eff]] <- pe; point_est_grid[[eff]] <- pe_grid
            CI_quant[[eff]] <- quant
            CI_hpd[[eff]] <- hpd
          } else if (eff == "spat"){
            # Compute
            reg <- dat$region
            V <- X[order(reg), eff_idx]
            uV <- unique(V)
            D <- draws[ , eff_idx] %*% t(uV)
            pe <- colMeans(D)
            quant <- t(apply(D, 2, quantile, probs = probs))
            hpd <- HPDinterval(as.mcmc(D), prob = probs[1] + (1-probs[2]))
            colnames(quant) <- colnames(hpd) <- c("LB", "UP")
            
            # Save
            D_eff[[eff]] <- D; D_eff_grid[[eff]] <- D
            point_est[[eff]] <- pe; point_est_grid[[eff]] <- pe
            CI_quant[[eff]] <- quant
            CI_hpd[[eff]] <- hpd
          } else if (eff == "f"){
            if (is.null(effects_idx$spat)){
              D_eff[[eff]] <- D_eff[["f"]]; D_eff_grid[[eff]] <- D_eff[["f"]]
              point_est[[eff]] <- pe; point_est_grid[[eff]] <- point_est[["f"]]
              CI_quant[[eff]] <- CI_quant[["f"]]
              CI_hpd[[eff]] <- CI_hpd[["f"]]
            } else {
              # Compute
              D <- apply(D_eff$spat, 2, function(x) x + D_eff$fcont)
              D_grid <- apply(D_eff_grid$spat, 2, function(x) x + D_eff_grid$fcont)
              pe <- point_est$fcont + point_est$fspat
              quant <- t(apply(D_grid, 2, quantile, probs = probs))
              hpd <- HPDinterval(as.mcmc(D_grid), prob = probs[1] + (1-probs[2]))
              colnames(quant) <- colnames(hpd) <- c("LB", "UP")
              
              # Save
              D_eff[[eff]] <- D; D_eff_grid[[eff]] <- D_grid
              point_est[[eff]] <- pe; point_est_grid[[eff]] <- pe_grid
              CI_quant[[eff]] <- quant
              CI_hpd[[eff]] <- hpd
            }
          }
        }
      }
      
      res[[mod]] <- list(point_est = point_est, point_est_grid = point_est_grid,
                         CI_quant = CI_quant, CI_hpd = CI_hpd)
    } else {
      
      
      
      res[[mod]] <- integer(0)
      cat("Computations of CI for m3 currently not implemented!!!")
    }
    
  }
  
  return(res)
  
}