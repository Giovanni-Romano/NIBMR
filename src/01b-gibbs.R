# MCMC SAMPLER ----
MCMC_IWLS <- function(niter,
                      X, y, 
                      k, k_deriv, c, c_deriv,
                      p, d, beta0,
                      a_tau, b_tau,
                      debug = FALSE, verbose = FALSE, print_step = 1000,
                      hess_eps = 1e-3,
                      K_fold = 50,
                      seed = NULL) {
  
  # Total number of regression coefficients:
  # intercept + p linear + sum(d) nonlinear basis coeffs
  P <- 1 + p + sum(d)
  
  if (P != length(beta0)) {
    stop("Length of beta0 is not equal to total parameter dimension P")
  }
  
  n <- nrow(X)
  loss_fun <- loss_asymm
  
  
  # --- Estimate w ---
  w_num <- 1/2
  if (!is.null(seed)) {set.seed(seed)}
  
  K <- K_fold
  if (verbose > 0) cat("Start ", K, "-fold CV. \t", sep = "")
  fold_id <- sample(rep(1:K, length.out = n))
  
  CV <- sapply(1:K, function(f) {
    test_idx  <- which(fold_id == f)
    train_idx <- setdiff(seq_len(n), test_idx)
    
    est <- optim(
      par = beta0,
      fn = function(b) {
        sum(loss_fun(y[train_idx] - X[train_idx, , drop = FALSE] %*% b, k, c))
      },
      method = "L-BFGS-B"
    )$par
    
    sum(loss_fun(y[test_idx] - X[test_idx, , drop = FALSE] %*% est, k, c))
  })
  
  mean_cv <- sum(CV) / n
  w <- w_num / mean_cv
  
  if (verbose > 0) cat("End CV.\n\n")
  
  
  # --- Storage ---
  beta_sample <- matrix(NA_real_, nrow = niter, ncol = P)
  colnames(beta_sample) <- c(
    "beta0",
    paste0("beta_lin", seq_len(p)),
    unlist(lapply(seq_len(p), function(i) {
      paste0("beta_nonlin", i, ".", seq_len(d[i]))
    }))
  )
  
  tau_sample <- matrix(NA_real_, nrow = niter, ncol = 2 * p)
  colnames(tau_sample) <- c(
    paste0("tau_lin", seq_len(p)),
    paste0("tau_nonlin", seq_len(p))
  )
  
  # alfa_sample <- matrix(NA_real_, nrow = niter, ncol = 1+p)
  logpost_sample <- matrix(NA_real_, nrow = niter, ncol = 1+p)
  
  
  # --- Initial values ---
  beta_curr <- drop(beta0)
  tau <- rep(10, 2 * p)
  
  # Maps beta[-1] to the corresponding tau component:
  # first p coefficients are linear effects,
  # next sum(d) coefficients are nonlinear groups
  tau_idx <- c(1:p, p + rep(seq_len(p), d))
  
  acc <- rep(0L, 1+p)
  acc_tmp <- rep(0L, 1+p)
  
  # --- Precompute outer X ---
  X2 <- X^2
  
  # --- MCMC loop ---
  for (m in seq_len(niter)) {
    
    if (verbose == 1 && m %% max(1, floor(niter / 10)) == 0) {
      cat(round(100 * m / niter), "%\t")
    }
    
    # --- IWLS update for beta ---
    
    ## --- Update intercept ---
    up <- iwls_beta_update(
      idx = 1,
      beta_curr = beta_curr,
      X = X, y = y, X2 = X2,
      k_post = k, k_deriv = k_deriv,
      c_post = c, c_deriv = c_deriv, 
      w = w,
      tau = tau, tau_idx = tau_idx,
      hess_eps = hess_eps
    )
    beta_curr <- up$beta
    # alfa_sample[m, 1] <- exp(min(up$log_alpha, 0))
    logpost_sample[m, 1] <- up$logpost
    acc[1] <- acc[1] + up$accepted
    acc_tmp[1] <- acc_tmp[1] + up$accepted
    
    ## --- Block update for each function/covariate ---
    for (j in 1:p){
      idx.b = 2+p+sum(d[1:j])-d[j]
      idx.e = idx.b+d[j]-1
      idx = c(1+j, idx.b:idx.e)
      up <- iwls_beta_update(
        idx = idx,
        beta_curr = beta_curr,
        X = X, y = y, X2 = X2,
        k_post = k, k_deriv = k_deriv,
        c_post = c, c_deriv = c_deriv, 
        w = w,
        tau = tau, tau_idx = tau_idx,
        hess_eps = hess_eps
      )
      beta_curr <- up$beta
      # alfa_sample[m, j+1] <- exp(min(up$log_alpha, 0))
      logpost_sample[m, j+1] <- up$logpost
      acc[j+1] <- acc[j+1] + up$accepted
      acc_tmp[j+1] <- acc_tmp[j+1] + up$accepted
    }
    
    # Save updated beta after full update
    beta_sample[m, ] <- beta_curr
    
    
    # --- Update tau ---
    for (j in seq_len(2 * p)) {
      beta_tmp <- beta_sample[m, 1 + which(tau_idx == j)]
      
      var_sample = 1 / rgamma(1,
                              shape = a_tau[j] + c(rep(1, p), d)[j] / 2,
                              rate  = b_tau[j] + sum(beta_tmp^2) / 2
      )
      
      tau_sample[m, j] <- sqrt(var_sample)
    }
    
    tau <- tau_sample[m, ]
    
    if (verbose == 2 && (m %% print_step == 0)) {
      # cat("iter =", m,
      #     "| acc rate =", round(acc / m, 3),
      #     "| alpha =", round(alfa_sample[m], 3),
      #     "| logpost =", round(logpost_sample[m], 3), "\n")
      cat("iter =", m,
          "| acc rate (last 1k) =", round(acc_tmp / print_step, 3),
          # "| alpha mean (last 1k) =", round(colMeans(alfa_sample[(m-print_step+1):m, ]), 3),
          "| acc rate =", round(acc / m, 3),
          # "| alpha mean =", round(mean(alfa_sample, na.rm = T), 3),
          "| logpost =", round(logpost_sample[m], 3),
          "\n")
      acc_tmp <- rep(0L, p+1)
    }
  }
  
  list(
    beta = beta_sample,
    tau = tau_sample,
    acc_prop = acc / niter,
    # alfa = alfa_sample,
    logpost = logpost_sample,
    w = w
  )
}



# MCMC SAMPLER BIKES APPLICATION ----
MCMC_IWLS_bikes <- function(niter,
                            X_dummy = NULL,
                            X_cont = NULL,
                            Z_cont = NULL,
                            Z_spat = NULL,
                            y,
                            freq = NULL,
                            link = function(x) x,
                            link_prime = function(x) rep(1, length(x)),
                            link_second = function(x) rep(0, length(x)),
                            k, k_deriv, c, c_deriv,
                            d, # splines dimension for every continuous covariate
                            beta0,
                            a_tau, b_tau,
                            debug = FALSE, verbose = 0, print_step = 1000,
                            hess_eps = 1e-3,
                            step = 1,
                            K_fold = 50,
                            seed = NULL) {
  
  cat("Fitting starts at: ", format(Sys.time(), "%H:%M:%S"), "\n")
  
  # Dimensions
  # number of binary covariates
  if (is.null(X_dummy)) {
    q = 0
  } else {
    q = ncol(X_dummy)
  }
  # number of continuous covariates
  p = ncol(X_cont)
  # number of spatial covariates
  if (is.null(Z_spat)) {
    s = 0
  } else {
    s = ncol(Z_spat)
  }
  n = length(y)
  
  
  if (!(length(d) == p)) {
    stop("Length of d must be equal to number of continuous covariates (i.e. # cols of X_cont)")
  }
  
  # Total number of regression coefficients:
  # intercept + q binary + p linear + sum(d) nonlinear basis coeffs
  P <- 1 + q + p + sum(d) + s
  
  if (P != length(beta0)) {
    stop("Length of beta0 is not equal to total parameter dimension P")
  }
  
  if (is.null(freq)) freq <- rep(1, length(y))
  
  loss_fun <- loss_asymm
  
  
  # Build design matrix X
  X = cbind(1, X_dummy, X_cont, Z_cont, Z_spat)
  # colnames(X) = c("Intercept", colnames(X_dummy), colnames(X_cont), 
  #                 c(sapply(1:p, function(j) paste0(colnames(X_cont)[j], "_DR", 1:d[j]))),
  #                 paste0("spatial_DR", 1:ncol(Z_spat)))
  
  # --- Estimate w ---
  w_num <- 1/2
  if (!is.null(seed)) {set.seed(seed)}
  
  start_Kfold <- Sys.time()
  K <- K_fold
  if (verbose > 0) cat("Start ", K, "-fold CV. ", sep = "")
  fold_id <- sample(rep(1:K, length.out = n))

  CV <- sapply(1:K, function(f) {
    test_idx  <- which(fold_id == f)
    train_idx <- setdiff(seq_len(n), test_idx)

    est <- optim(
      par = beta0,
      fn = function(b) {
        lin_pred = X[train_idx, , drop = FALSE] %*% b
        sum(freq[train_idx] * loss_fun(y[train_idx] - link(lin_pred), k, c))
      },
      method = "L-BFGS-B"
    )$par

    sum(freq[test_idx] * loss_fun(y[test_idx] - link(X[test_idx, , drop = FALSE] %*% est), k, c))
  })

  mean_cv <- sum(CV) / sum(freq)
  w <- w_num / mean_cv
  # w = 0.84
  end_Kfold <- Sys.time()
  diff_time_kfold = round(difftime(end_Kfold, start_Kfold, units = "mins"), 2)
  if (verbose > 0) cat("End CV. \t", "Exec. time: ", diff_time_kfold, " mins \n\n", sep = "")
  
  
  # --- Storage ---
  beta_sample <- matrix(NA_real_, nrow = niter, ncol = P)
  colnames(beta_sample) <- colnames(X)
  
  tau_sample <- matrix(NA_real_, nrow = niter, ncol = q + 2 * p + 1*(s>0))
  # colnames(tau_sample) <- c(
  #   paste0("tau_bin", seq_len(q)),
  #   paste0("tau_lin", seq_len(p)),
  #   paste0("tau_nonlin", seq_len(p)),
  #   "tau_spat")
  
  
  
  
  # --- Initial values ---
  beta_curr <- unname(drop(beta0))
  tau <- rep(10, q + 2 * p + 1*(s>0))
  
  # Maps beta[-1] to the corresponding tau component:
  # first q coefficients are binary covariates,
  # next p coefficients are linear effects,
  # next sum(d) coefficients are nonlinear groups
  # Corrected tau_idx mapping
  tau_idx <- c(
    if (q > 0) 1:q else integer(0),
    if (p > 0) (q + 1):(q + p) else integer(0),
    if (p > 0) q + p + rep(seq_len(p), d) else integer(0),
    if (s > 0) rep(q + 2 * p + 1, s) else integer(0)
  )
  
  # IG shape 
  IG_size_vec <- c(rep(1, q + p), d)
  if (s > 0) IG_size_vec <- c(IG_size_vec, s)
  
  # Storage for acceptance/rejection of MH
  n_sp_up = floor(s / 5) + 1*(s %% 5 > 0)
  acc <- rep(0L, 1+q+2*p+n_sp_up)
  # acc <- rep(0L, 1+q+2*p+1*(s>0))
  acc_tmp <- acc
  logpost_sample <- matrix(NA_real_, nrow = niter, ncol = length(acc))
  
  # --- Precompute outer X ---
  X2 <- X^2
  
  # Elapsed time initialization
  time_start = time_last = Sys.time()
  cat("Start Gibbs of ", niter, " iterations \n", sep = "")
  
  # --- MCMC loop ---
  for (m in seq_len(niter)) {
    
    if (verbose == 1 && m %% max(1, floor(niter / 10)) == 0) {
      cat(round(100 * m / niter), "%\t")
    }
    
    # --- IWLS update for beta ---
    
    ## --- Update intercept ---
    up <- iwls_beta_update(
      idx = 1,
      beta_curr = beta_curr,
      X = X, y = y, X2 = X2,
      freq = freq,
      link = link, link_prime = link_prime, link_second = link_second,
      k_post = k, k_deriv = k_deriv,
      c_post = c, c_deriv = c_deriv, 
      w = w,
      tau = tau, tau_idx = tau_idx,
      hess_eps = hess_eps,
      step = step
    )
    beta_curr <- up$beta
    # alfa_sample[m, 1] <- exp(min(up$log_alpha, 0))
    logpost_sample[m, 1] <- up$logpost
    acc[1] <- acc[1] + up$accepted
    acc_tmp[1] <- acc_tmp[1] + up$accepted
    
    ## --- Update binary covariates ---
    if (!is.null(X_dummy)){
      for (j_bin in 1:q){
        idx = 1+j_bin
        up <- iwls_beta_update(
          idx = idx,
          beta_curr = beta_curr,
          X = X, y = y, X2 = X2,
          freq = freq,
          link = link, link_prime = link_prime, link_second = link_second,
          k_post = k, k_deriv = k_deriv,
          c_post = c, c_deriv = c_deriv, 
          w = w,
          tau = tau, tau_idx = tau_idx,
          hess_eps = hess_eps,
          step = step
        )
        beta_curr <- up$beta
        logpost_sample[m, 1+j_bin] <- up$logpost
        acc[1+j_bin] <- acc[1+j_bin] + up$accepted
        acc_tmp[1+j_bin] <- acc_tmp[1+j_bin] + up$accepted
      }
    }
    
    ## --- Block update for each continuous covariate (linear effect) ---
    if (!is.null(X_cont)){
      for (j_cont in 1:p){
        idx = 1+q+j_cont
        up <- iwls_beta_update(
          idx = idx,
          beta_curr = beta_curr,
          X = X, y = y, X2 = X2,
          freq = freq,
          link = link, link_prime = link_prime, link_second = link_second,
          k_post = k, k_deriv = k_deriv,
          c_post = c, c_deriv = c_deriv, 
          w = w,
          tau = tau, tau_idx = tau_idx,
          hess_eps = hess_eps,
          step = step
        )
        beta_curr <- up$beta
        logpost_sample[m, 1+q+j_cont] <- up$logpost
        acc[1+q+j_cont] <- acc[1+q+j_cont] + up$accepted
        acc_tmp[1+q+j_cont] <- acc_tmp[1+q+j_cont] + up$accepted
      }
    }
    
    ## --- Block update for each continuous covariate (NON-linear effect) ---
    if (!is.null(X_cont)){
      for (j_DR in 1:p){
        idx.b = (1+q+p) + (sum(d[1:j_DR])-d[j_DR]) + 1
        idx.e = idx.b+d[j_DR]-1
        idx = idx.b:idx.e
        up <- iwls_beta_update(
          idx = idx,
          beta_curr = beta_curr,
          X = X, y = y, X2 = X2,
          freq = freq,
          link = link, link_prime = link_prime, link_second = link_second,
          k_post = k, k_deriv = k_deriv,
          c_post = c, c_deriv = c_deriv, 
          w = w,
          tau = tau, tau_idx = tau_idx,
          hess_eps = hess_eps,
          step = step
        )
        beta_curr <- up$beta
        logpost_sample[m, 1+q+p+j_DR] <- up$logpost
        acc[1+q+p+j_DR] <- acc[1+q+p+j_DR] + up$accepted
        acc_tmp[1+q+p+j_DR] <- acc_tmp[1+q+p+j_DR] + up$accepted
      }
    }
    
    
    ## --- Multi-blocks update for spatial effect ---
    if (!is.null(Z_spat)) {
      idx.e = 1 + q + p + sum(d)
      
      for (r in 1:n_sp_up){

        idx.b = idx.e + 1
        idx.e = min(idx.b + 4, P)
        idx = idx.b:idx.e

        up <- iwls_beta_update(
          idx = idx,
          beta_curr = beta_curr,
          X = X, y = y, X2 = X2,
          freq = freq,
          link = link, link_prime = link_prime, link_second = link_second,
          k_post = k, k_deriv = k_deriv,
          c_post = c, c_deriv = c_deriv,
          w = w,
          tau = tau, tau_idx = tau_idx,
          hess_eps = hess_eps,
          step = step/4
        )
        beta_curr <- up$beta
        logpost_sample[m, 1+q+2*p+r] <- up$logpost
        acc[1+q+2*p+r] <- acc[1+q+2*p+r] + up$accepted
        acc_tmp[1+q+2*p+r] <- acc_tmp[1+q+2*p+r] + up$accepted
      }

      
        # idx.b = idx.e + 1
        # idx.e = idx.e + s
        # idx = idx.b:idx.e
        # 
        # up <- iwls_beta_update(
        #   idx = idx,
        #   beta_curr = beta_curr,
        #   X = X, y = y, X2 = X2,
        #   freq = freq,
        #   link = link, link_prime = link_prime, link_second = link_second,
        #   k_post = k, k_deriv = k_deriv,
        #   c_post = c, c_deriv = c_deriv, 
        #   w = w,
        #   tau = tau, tau_idx = tau_idx,
        #   hess_eps = hess_eps,
        #   step = step/4
        # )
        # beta_curr <- up$beta
        # logpost_sample[m, 1+q+2*p+1] <- up$logpost
        # acc[1+q+2*p+1] <- acc[1+q+2*p+1] + up$accepted
        # acc_tmp[1+q+2*p+1] <- acc_tmp[1+q+2*p+1] + up$accepted
    }
    
    # Save updated beta after full update
    beta_sample[m, ] <- beta_curr
    
    
    # --- Update tau ---
    for (j in seq_len(q + 2*p + 1*(s>0))) {
      beta_tmp <- beta_sample[m, 1 + which(tau_idx == j)]
      
      var_sample = 1 / rgamma(1,
                              shape = a_tau[j] + IG_size_vec[j] / 2,
                              rate  = b_tau[j] + sum(beta_tmp^2) / 2
      )
      
      tau_sample[m, j] <- sqrt(var_sample)
    }
    
    tau <- tau_sample[m, ]
    
    if (verbose == 2 && (m %% print_step == 0)) {
      now = Sys.time()
      diff_time_start = round(difftime(now, time_start, units = "mins"), 2)
      diff_time_last = round(difftime(now, time_last, units = "mins"), 2)
      clock_time = format(Sys.time(), "%H:%M:%S")
      cat("iter = ", m, "\n",
          "print clock time: ", clock_time, " | ",
          "last ", print_step, " iter time: ", diff_time_last, " mins | ", 
          "total exe time: ", diff_time_start, " mins \n",
          "acc rate (last ", print_step, ") = ", paste(round(acc_tmp / print_step, 3), collapse = " "),
          # " | acc rate = ", paste(round(acc / m, 3), collapse = " "),
          ifelse(verbose > 2, paste0(" | logpost = ", round(logpost_sample[m], 3)), ""),
          "\n\n\n",
          sep = "")
      time_last = now
      acc_tmp <- rep(0L, length(acc))
    }
  }
  
  list(
    beta = beta_sample,
    tau = tau_sample,
    acc_prop = acc / niter,
    logpost = logpost_sample,
    w = w
  )
}