# REQUIRED LIBRARIES ----
require(mvnfast)
require(mvtnorm)
require(mgcv)

# UTILS ----
## Loss function ----
loss_asymm = function(eps, k, c, g = 1){
  k1 = k[1]
  k2 = k[2]
  h = 2*g
  (eps < 0) * (1 - exp(c^(1/h)-( (k1*eps)^h + c )^(1/h) )) + 
    (eps > 0) * (1 - exp(c^(1/h)-( (k2*eps)^h + c )^(1/h) ))
}

loss_asymm_pop = function(beta, y, X, k, c, g = 1){
  eps = y - X%*%beta
  k_vec = k[(eps > 0) + 1]
  h = 2*g
  sum(1 - exp(c^(1/h)-( (k_vec*eps)^h + c )^(1/h) ))
}


## Derivatives ----
H_pop = function(idx, beta, y, X, X2, k, c){
  eps = drop(y - drop(X %*% beta))
  k_vec = k[(eps>0)+1]
  u = k_vec*eps
  
  F1 = k_vec^2 * exp(sqrt(c) - sqrt(u^2+c)) * ((u^2)/(u^2+c) - c/((u^2+c)^(3/2)))
  out = colSums(X2[ , idx, drop = F] * F1, dims = 1)
  
  return(out)
}

score_pop = function(idx, beta, y, X, k, c){
  eps = drop((y - X %*% beta))
  k_vec = k[(eps>0)+1]
  u = k_vec*eps
  
  F1 = k_vec * exp(sqrt(c) - sqrt(u^2+c)) * u/sqrt(u^2+c)
  out = colSums(X[ , idx, drop = F]*F1)
  
  return(out)
}


## DR Decomposition ----
# This function returns the nonlinear Demmler-Reinsch spline basis functions
# Construct Demmler Reinsch basis according to Bach and Klein 2024
construct_DR_basis <- function(x, d = 10, rescale = TRUE) {
  
  n <- length(x)
  # set up P-spline design
  fit=smoothCon(s(x, bs="ps", k=d+2, m = c(2, 2)), data=data.frame(x),scale.penalty=FALSE)
  sm = fit[[1]]
  
  B=fit[[1]]$X
  K=fit[[1]]$S[[1]]
  
  # set up constraint matrix
  C <- t(cbind(1, scale(x))) %*% B
  V0 <- svd(C,nv=d+2)$v[, 3:(d+2), drop = FALSE]
  
  # transformed penalty and gram matrices
  K_tilde <- t(V0) %*% K %*% V0
  G_tilde <- t(V0) %*% (t(B) %*% B / n) %*% V0
  # solve generalized eigenvalue problem
  eig <- geigen::geigen(G_tilde, K_tilde)
  A <- eig$vectors
  Trf <- V0 %*% A
  
  # perform change of basis
  Z <- B %*% Trf
  # normalize to unit norm
  scaling_factor <- 1
  if (rescale) {
    scaling_factor <- sqrt(n / sum(Z^2))
    Z <- Z * scaling_factor
    Trf <- Trf * scaling_factor
  }
  # you don't really need the full list but I returned everything and kept what I needed
  list(Z = Z, Trafo = Trf, 
       smooth_object = sm,
       Gram = G_tilde, Penalty = K_tilde)
}


reconstruct_DR_for_grid <- function(x, Z_saved, x_grid, d = 20, rescale = TRUE) {
  
  dr_rec <- construct_DR_basis(x, d = d, rescale = rescale)
  
  # Reconstruct B_grid using the training smooth object
  B_grid <- mgcv::PredictMat(
    dr_rec$smooth_object,
    data.frame(x = x_grid)
  )
  
  Z_rec_train <- dr_rec$Z
  Z_rec_grid  <- B_grid %*% dr_rec$Trafo
  
  # Map reconstructed basis coordinates to the basis actually used in fitting
  M <- qr.solve(Z_rec_train, Z_saved)
  
  Z_grid_aligned <- Z_rec_grid %*% M
  
  # Diagnostic: should be very small
  err <- max(abs(Z_rec_train %*% M - Z_saved))
  
  x_center <- mean(x)#dr_rec$x_center
  x_scale  <- sd(x)#dr_rec$x_scale
  x_grid_scaled <- (x_grid - x_center) / x_scale
  
  list(
    x_scaled = as.numeric(x_grid_scaled),
    Z = Z_grid_aligned,
    reconstruction_error = err
  )
}


## IWLS utils ----
### IWLS single update ----
iwls_beta_update <- function(idx, beta_curr, 
                             X, y, X2,
                             k_post, k_deriv, 
                             c_post, c_deriv, 
                             w, tau, tau_idx,
                             hess_eps = 1e-3) {
  # Current state
  cur <- eval_logfullcond_beta(
    idx = idx, beta = beta_curr, 
    X = X, y = y, X2 = X2,
    k_post = k_post, k_deriv = k_deriv, 
    c_post = c_post, c_deriv = c_deriv,
    w = w,
    tau = tau, tau_idx = tau_idx
  )
  
  Hc <- make_negdef(cur$H, eps = hess_eps)
  Sigma_c <- -1/Hc
  mu_c <- drop(beta_curr[idx] + Sigma_c*cur$grad)
  
  # Proposal draw
  beta_prop <- beta_curr
  beta_prop[idx] <- rnorm(length(idx), mean = mu_c, sd = sqrt(Sigma_c))
  
  # Proposed state
  prop <- eval_logfullcond_beta(
    idx = idx, beta = beta_prop, 
    X = X, y = y, X2 = X2,
    k_post = k_post, k_deriv = k_deriv,
    c_post = c_post, c_deriv = c_deriv,
    w = w,
    tau = tau, tau_idx = tau_idx
  )
  
  Hp <- make_negdef(prop$H, eps = hess_eps)
  Sigma_p <- -1/Hp
  mu_p <- drop(beta_prop[idx] + Sigma_p*prop$grad)
  
  # MH acceptance probability
  log_q_prop_given_cur <- sum(dnorm(beta_prop[idx], mean = mu_c, sd = sqrt(Sigma_c), log = T))
  log_q_cur_given_prop <- sum(dnorm(beta_curr[idx], mean = mu_p, sd = sqrt(Sigma_p), log = T))
  
  log_alpha <- prop$logpost - cur$logpost +
    log_q_cur_given_prop - log_q_prop_given_cur
  
  if (log(runif(1)) < log_alpha) {
    list(
      beta = beta_prop,
      accepted = 1,
      log_alpha = log_alpha,
      logpost = prop$logpost
    )
  } else {
    list(
      beta = beta_curr,
      accepted = 0,
      log_alpha = log_alpha,
      logpost = cur$logpost
    )
  }
}

### Log-full-conditional ----
eval_logfullcond_beta <- function(idx, beta, 
                                  X, y, X2 = X2,
                                  k_post, k_deriv, 
                                  c_post, c_deriv,
                                  w, tau, tau_idx) {
  beta <- drop(beta)
  n <- nrow(X)
  P <- length(beta)
  
  loss_fun <- loss_asymm
  
  # --- log-likelihood ---
  r <- drop(y - X %*% beta)
  loglik <- -w * sum(loss_fun(r, k = k_post, c = c_post))
  
  grad_ll = w * score_pop(idx = idx, beta = beta, y = y, X = X, k = k_deriv, c = c_deriv) 
  H_ll = w * H_pop(idx = idx, beta = beta, y = y, X = X, X2 = X2, k = k_deriv, c = c_deriv)
  
  # --- Gaussian prior ---
  # No prior on intercept beta[1]
  # Prior on beta[-1] with group-specific variances tau[tau_idx]^2
  if (idx[1] == 1) {
    logprior <- 0
    grad_prior <- 0
    H_prior <- 0
  } else {
    tau_tmp <- tau[tau_idx[idx-1]]
    Omega <- 1 / (tau_tmp^2)
    logprior <- -0.5 * sum(beta[-1]^2 / (tau_tmp^2))
    grad_prior <- - (Omega * beta[idx])
    H_prior <- - Omega
  }
  
  # --- full conditional ---
  logpost <- loglik + logprior
  grad <- grad_ll + grad_prior
  H <- H_ll + H_prior
  
  list(
    logpost = as.numeric(logpost),
    grad = drop(grad),
    H = H
  )
}

### Enforce symmetry and negative definiteness in Hessian ----
make_negdef <- function(H, eps = 1e-3) {
  H[H>=0] = -eps
  return(H)
}

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
      method = "BFGS"
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
  # logpost_sample <- matrix(NA_real_, nrow = niter, ncol = 1+p)
  
  
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
    # logpost_sample[m, 1] <- up$logpost
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
          "| alpha mean (last 1k) =", round(colMeans(alfa_sample[(m-print_step+1):m, ]), 3),
          # "| acc rate =", round(acc / m, 3),
          # "| alpha mean =", round(mean(alfa_sample, na.rm = T), 3),
          # "| logpost =", round(logpost_sample[m], 3), 
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