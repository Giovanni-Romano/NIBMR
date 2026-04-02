# UTILS ----
# Loss function
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


MCMC_MH = function(niter, 
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
      tau_sample[m, ] = 1/rgamma(1, 
                                 shape = a_tau + c(rep(1:p), d)[j]/2,
                                 rate = b_tau + sum(beta_tmp^2)/2)
    }
  }
  
  return(list(beta = beta_sample, tau = tau_sample,
              acc_prop = acc/niter, alfa = alfa_sample, sd_prop = sd_prop_sample,
              w = w))
}


# Other objects ----
H_i = function(beta, y_i, x_i, k, c){
  eps = drop((y_i - drop(crossprod(x_i, beta))))
  k_i = k[(eps>0)+1]
  u = k_i*eps
  
  F1 = exp(sqrt(c) - sqrt(u^2+c)) * ((u^2)/(u^2+c) - c/((u^2+c)^(3/2)))
  F2 = outer(-k_i*x_i, -k_i*x_i)
  
  return(F1 * F2)
}

score_i = function(beta, y_i, x_i, k, c){
  eps = drop((y_i - crossprod(x_i, beta)))
  k_i = k[(eps>0)+1]
  u = k_i*eps
  
  F1 = -exp(sqrt(c) - sqrt(u^2+c)) * u/sqrt(u^2+c)
  F2 = -k_i*x_i
  
  return(F1 * F2)
}

# DR Decomposition ----
## this function returns the nonlinear Demmler-Reinsch spline basis functions
# Construct Demmler Reinsch basis according to Bach and Klein 2024
construct_DR_basis <- function(x, d = 10, rescale = TRUE) {
  
  require(mgcv)
  
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