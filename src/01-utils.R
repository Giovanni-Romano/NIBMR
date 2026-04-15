# REQUIRED LIBRARIES ----
require(mvnfast)
require(mvtnorm)
require(mgcv)

# UTILS ----
## Loss functions ----
loss_symm = function(eps, k, c, g = 1, theta = 0){
  h = 2*g
  1 - exp( c^(1/h)-(((k*eps)^h + c)^(1/h)) )
}

loss_asymm = function(eps, k, c, g = 1, theta = 0){
  h = 2*g
  m1 = tan(pi/4 - theta)
  1 - exp( c^(1/h)-(m1^sign(eps))*(((k*eps)^h + c)^(1/h)) )
}


loss_asymm2 = function(eps, k, c, g = 1){
  k1 = k[1]
  k2 = k[2]
  h = 2*g
  (eps < 0) * (1 - exp(c^(1/h)-( (k1*eps)^h + c )^(1/h) )) + 
    (eps > 0) * (1 - exp(c^(1/h)-( (k2*eps)^h + c )^(1/h) ))
}


loss_pop = function(beta, y, X, k, c= 1e-3, g = 1){
  n = length(y)
  sum(sapply(1:n, function(i) loss(y[i]-X[i, ]%*%beta, k, c, g)))
}

## Derivatives ----
H_i = function(beta, y_i, x_i, outerX_i, k, c){
  eps = drop((y_i - drop(crossprod(x_i, beta))))
  k_i = k[(eps>0)+1]
  u = k_i*eps
  
  F1 = exp(sqrt(c) - sqrt(u^2+c)) * ((u^2)/(u^2+c) - c/((u^2+c)^(3/2)))
  # F2 = outer(-k_i*x_i, -k_i*x_i)
  F2 = outerX_i*(k_i^2)
  
  return(F1 * F2)
}

H_pop = function(beta, y, X, outerX, k, c){
  eps = drop(y - drop(X %*% beta))
  k_vec = k[(eps>0)+1]
  u = k_vec*eps
  
  F1 = k_vec^2 * exp(sqrt(c) - sqrt(u^2+c)) * ((u^2)/(u^2+c) - c/((u^2+c)^(3/2)))
  out = colSums(outerX * F1, dims = 1)
  
  return(out)
}

score_i = function(beta, y_i, x_i, k, c){
  eps = drop((y_i - crossprod(x_i, beta)))
  k_i = k[(eps>0)+1]
  u = k_i*eps
  
  F1 = -exp(sqrt(c) - sqrt(u^2+c)) * u/sqrt(u^2+c)
  F2 = -k_i*x_i
  
  return(F1 * F2)
}

score_pop = function(beta, y, X, k, c){
  eps = drop((y - X %*% beta))
  k_vec = k[(eps>0)+1]
  u = k_vec*eps
  
  F1 = k_vec * exp(sqrt(c) - sqrt(u^2+c)) * u/sqrt(u^2+c)
  out = colSums(X * F1)
  
  return(out)
}


## DR Decomposition ----
# This function returns the nonlinear Demmler-Reinsch spline basis functions
# Construct Demmler Reinsch basis according to Bach and Klein 2024
construct_DR_basis <- function(x, d = 10, rescale = TRUE) {
  
  n <- length(x)
  # set up P-spline design
  fit=smoothCon(s(x, bs="ps", k=d+2, m = c(2, 2)), data=data.frame(x),scale.penalty=FALSE)
  
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
  list(Z = Z, Trafo = Trf, Gram = G_tilde, Penalty = K_tilde)
}

## IWLS utils ----
### IWLS single update ----
iwls_beta_update <- function(beta_curr, X, y, outerX,
                             k, c, w, tau, tau_idx,
                             asymm = TRUE, hess_eps = 1e-8,
                             score_i, H_i) {
  # Current state
  cur <- eval_logfullcond_beta(
    beta = beta_curr, X = X, y = y, outerX = outerX,
    k = k, c = c, w = w,
    tau = tau, tau_idx = tau_idx, asymm = asymm,
    score_i = score_i, H_i = H_i
  )
  
  Hc <- make_negdef(cur$H, eps = hess_eps)
  Sigma_c.tmp <- -solve(Hc)
  Sigma_c = (Sigma_c.tmp + t(Sigma_c.tmp))/2
  mu_c <- drop(beta_curr + drop(Sigma_c %*% cur$grad))
  
  # Proposal draw
  beta_prop <- drop(mvnfast::rmvn(1, mu = mu_c, sigma = Sigma_c, isChol = F))
  
  # Proposed state
  prop <- eval_logfullcond_beta(
    beta = beta_prop, X = X, y = y, outerX = outerX,
    k = k, c = c, w = w,
    tau = tau, tau_idx = tau_idx, asymm = asymm,
    score_i = score_i, H_i = H_i
  )
  
  Hp <- make_negdef(prop$H, eps = hess_eps)
  Sigma_p.tmp <- -solve(Hp)
  Sigma_p = (Sigma_p.tmp + t(Sigma_p.tmp))/2
  mu_p <- drop(beta_prop - solve(Hp, prop$grad))
  
  # MH acceptance probability
  log_q_prop_given_cur <- dmvnorm(beta_prop, mu_c, Sigma_c, log=TRUE,checkSymmetry=FALSE) 
  log_q_cur_given_prop <- dmvnorm(beta_curr, mu_p, Sigma_p, log=TRUE,checkSymmetry=FALSE)
  
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
eval_logfullcond_beta <- function(beta, X, y, outerX = outerX,
                                  k, c, w, tau, tau_idx,
                                  asymm = TRUE,
                                  score_i, H_i) {
  beta <- drop(beta)
  n <- nrow(X)
  P <- length(beta)
  
  loss_fun <- if (asymm) loss_asymm2 else loss_symm
  
  # --- log-likelihood ---
  r <- drop(y - X %*% beta)
  loglik <- -w * sum(loss_fun(r, k = k, c = c))
  
  grad_ll = w * score_pop(beta, y, X, k, c) 
  H_ll = w * H_pop(beta, y, X, outerX, k, c)
  
  # --- Gaussian prior ---
  # No prior on intercept beta[1]
  # Prior on beta[-1] with group-specific variances tau[tau_idx]^2
  tau_tmp <- tau[tau_idx]  # length P-1
  
  Omega_diag <- c(0, 1 / (tau_tmp^2))
  Omega <- diag(Omega_diag, nrow = P)
  
  logprior <- -0.5 * sum(beta[-1]^2 / (tau_tmp^2))
  grad_prior <- - drop(Omega %*% beta)
  H_prior <- - Omega
  
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
make_negdef <- function(H, eps = 1e-8) {
  H <- (H + t(H)) / 2
  ee <- eigen(H, symmetric = TRUE)
  
  # Force all eigenvalues to be <= -eps
  ee$values <- pmin(ee$values, -eps)
  H_new <- ee$vectors %*% diag(ee$values, nrow = nrow(H)) %*% t(ee$vectors)
  
  # # Modification in formula 3.44 of Nocedal&Wright2006
  # ee_val = ee$values
  # e_new = ee_val
  # e_new[ee_val >-eps] = -(eps + ee_val)[ee_val > -eps]
  # Lambda_new = diag(e_new)
  # H_new <- ee$vectors %*% Lambda_new %*% t(ee$vectors)
  # (H_new + t(H_new)) / 2
}

# MCMC SAMPLERS ----
MCMC_RWMH = function(niter, 
                     X, y, k, c, p, d, beta0,
                     a_tau = 1e-3, b_tau = 1e-3,
                     asymm = T, debug = F, verbose = F){
  
  P = 1+p+sum(d) # total size of the design matrix
  
  if (P != length(beta0)){stop("Length of beta0 is not equal to p")}
  
  n = nrow(X)
  
  loss = ifelse(asymm, loss_asymm2, loss_symm)
  
  
  if (verbose) {cat("Start 50-fold CV \t")}
  w_num = 1/2
  
  set.seed(123)  # for reproducibility
  
  K <- 50
  fold_id <- sample(rep(1:K, length.out = n))
  CV10 <- sapply(1:K, function(f) {
    test_idx  <- which(fold_id == f)
    train_idx <- setdiff(1:n, test_idx)
    
    est <- optim(
      par = beta0,
      fn = function(b) sum(loss(y[train_idx] - X[train_idx, ] %*% b, k, 1e-3)),
      method = "BFGS"
    )$par
    
    sum(loss(y[test_idx] - X[test_idx, ] %*% est, k, 1e-3))
  })
  mean_cv10 <- sum(CV10) / n
  w = w_num / mean_cv10
  if (verbose) {cat("End CV. \n\n")}
  
  # LOO <- sapply(1:n, function(i) {
  #   est <- optim(par = beta0,
  #                fn = function(b) sum(loss(y[-i] - X[-i, ] %*% b, k, 1e-3)),
  #                method = "BFGS")$par
  #   loss(y[i] - drop(X[i, ] %*% est), k, 1e-3)
  # })
  # w = w_num / (sum(LOO)/n)
  
  
  beta_sample = matrix(NA, nrow = niter, ncol = P)
  colnames(beta_sample) = c("beta0", 
                            paste0("beta_lin", 1:p), 
                            unlist(lapply(seq_len(p), function(i) {paste0("beta_nonlin", i, ".", seq_len(d[i]))}))
  )
  beta_init = beta0
  last = beta_init
  
  tau_sample = matrix(NA, nrow = niter, ncol = 2*p)
  tau = rep(1, 2*p)
  tau_idx = c(1:p, p+rep(1:p, d))
  
  acc = 0
  alfa_sample = rep(NA, niter)
  sd_prop_sample = c()
  
  Sigma0 = diag(1, nrow = P)
  cholSigma0 = chol(Sigma0)
  
  for (m in 1:niter){
    
    if (verbose & m %% (niter/10) == 0){
      cat(10 * m / (niter/10), "% \t")
    }
    
    if (m == 1){
      RM = log(2.38/sqrt(P))
    } else {
      RM = log(sd_prop) + (1/m)^(3/4)*(acc/(m-1) - 0.234)
    }
    
    sd_prop = exp(RM); sd_prop_sample = c(sd_prop_sample, sd_prop)
    prop = t(mvnfast::rmvn(n = 1, mu = last, sigma = sd_prop*cholSigma0, isChol = T))
    
    r_last = y - X%*%last
    r_prop = y - X%*%prop
    
    loss_last = loss(r_last, k = k, c = c)
    loss_prop = loss(r_prop, k = k, c = c)
    
    if (debug){
      plot(density(r_last)); lines(density(r_prop), col = 2)
      cat("\n\n Mean losses", colMeans(loss_sample[m, , ]), "\n\n")
    }
    
    loglik_last = -w*sum(loss_last)
    loglik_prop = -w*sum(loss_prop)
    
    if (debug){
      cat("loglik last:", loglik_last, "\t | \t loglik prop:", loglik_prop, "\n")
    }
    
    tau_tmp = tau[tau_idx]
    
    target_last = loglik_last - 0.5*sum(last[-1]^2/(tau_tmp^2))
    target_prop = loglik_prop - 0.5*sum(prop[-1]^2/(tau_tmp^2))
    
    alfa = alfa_sample[m] = exp(min(target_prop - target_last, 0))
    u = runif(1, 0, 1)
    
    if (u < alfa){
      beta_sample[m, ] = prop
      last = prop
      is_acc = 1
    } else {
      beta_sample[m, ] = last
      is_acc = 0
    }
    acc = acc + is_acc
    
    
    # Sample tau
    for (j in 1:(2*p)){
      beta_tmp = beta_sample[m, 1+which(tau_idx == j)]
      tau_sample[m, j] = 1/rgamma(1, 
                                  shape = a_tau + c(rep(1, p), d)[j]/2,
                                  rate = b_tau + sum(beta_tmp^2)/2)
    }
  }
  
  return(list(beta = beta_sample, tau = tau_sample,
              acc_prop = acc/niter, alfa = alfa_sample, sd_prop = sd_prop_sample,
              w = w))
}


MCMC_IWLS <- function(niter,
                      X, y, k, c, p, d, beta0,
                      a_tau = 1e-3, b_tau = 1e-3,
                      asymm = TRUE, debug = FALSE, verbose = FALSE,
                      hess_eps = 1e-8,
                      seed = 123) {
  
  # Total number of regression coefficients:
  # intercept + p linear + sum(d) nonlinear basis coeffs
  P <- 1 + p + sum(d)
  
  if (P != length(beta0)) {
    stop("Length of beta0 is not equal to total parameter dimension P")
  }
  
  n <- nrow(X)
  loss_fun <- if (asymm) loss_asymm2 else loss_symm
  
  
  # --- Estimate w ---
  if (verbose) cat("Start 50-fold CV\t")
  
  w_num <- 1/2
  set.seed(seed)
  
  K <- 10
  fold_id <- sample(rep(1:K, length.out = n))
  
  CV10 <- sapply(1:K, function(f) {
    test_idx  <- which(fold_id == f)
    train_idx <- setdiff(seq_len(n), test_idx)
    
    est <- optim(
      par = beta0,
      fn = function(b) {
        sum(loss_fun(y[train_idx] - X[train_idx, , drop = FALSE] %*% b, k, 1e-3))
      },
      method = "BFGS"
    )$par
    
    sum(loss_fun(y[test_idx] - X[test_idx, , drop = FALSE] %*% est, k, 1e-3))
  })
  
  mean_cv10 <- sum(CV10) / n
  w <- w_num / mean_cv10
  
  if (verbose) cat("End CV.\n\n")
  
  
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
  
  alfa_sample <- rep(NA_real_, niter)
  logpost_sample <- rep(NA_real_, niter)
  

  # --- Initial values ---
  beta_curr <- drop(beta0)
  tau <- rep(100, 2 * p)
  
  # Maps beta[-1] to the corresponding tau component:
  # first p coefficients are linear effects,
  # next sum(d) coefficients are nonlinear groups
  tau_idx <- c(1:p, p + rep(seq_len(p), d))
  
  acc <- 0L
  
  # --- Precompute outer X ----
  outerX <- aperm(simplify2array(lapply(1:nrow(X), function(i) tcrossprod(X[i, ]))), c(3, 1, 2))
  
  # --- MCMC loop ---
  for (m in seq_len(niter)) {
    
    if (verbose && m %% max(1, floor(niter / 10)) == 0) {
      cat(round(100 * m / niter), "%\t")
    }
    
    # --- IWLS update for beta ---
    up <- iwls_beta_update(
      beta_curr = beta_curr,
      X = X, y = y, outerX = outerX,
      k = c(1, 1), c = c, w = w,
      tau = tau, tau_idx = tau_idx,
      asymm = asymm, hess_eps = hess_eps,
      score_i = score_i, H_i = H_i
    )
    
    beta_curr <- up$beta
    beta_sample[m, ] <- beta_curr
    alfa_sample[m] <- exp(min(up$log_alpha, 0))
    logpost_sample[m] <- up$logpost
    acc <- acc + up$accepted
    
    # --- Update tau ---
    for (j in seq_len(2 * p)) {
      beta_tmp <- beta_sample[m, 1 + which(tau_idx == j)]
      
      var_sample = 1 / rgamma(1,
                              shape = a_tau + c(rep(1, p), d)[j] / 2,
                              rate  = b_tau + sum(beta_tmp^2) / 2
      )
      
      tau_sample[m, j] <- sqrt(var_sample)
    }
    
    tau <- tau_sample[m, ]
    
    if (debug && (m %% 100 == 0)) {
      cat("iter =", m,
          "| acc rate =", round(acc / m, 3),
          "| alpha =", round(alfa_sample[m], 3),
          "| logpost =", round(logpost_sample[m], 3), "\n")
    }
  }
  
  list(
    beta = beta_sample,
    tau = tau_sample,
    acc_prop = acc / niter,
    alfa = alfa_sample,
    logpost = logpost_sample,
    w = w
  )
}