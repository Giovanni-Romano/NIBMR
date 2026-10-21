# REQUIRED LIBRARIES ----
# require(mvnfast)
# require(mvtnorm)
require(mgcv)
# library(mgcv, lib.loc = "/Users/giovanni/Library/Caches/org.R-project.R/R/renv/sandbox/macos/R-4.5/aarch64-apple-darwin20/4cd76b74")

# UTILS ----
## Loss function ----
loss_asymm = function(eps, k, c, g = 1){
  k1 = k[1]
  k2 = k[2]
  h = 2*g
  (eps < 0) * (1 - exp(c^(1/h)-( (k1*eps)^h + c )^(1/h) )) + 
    (eps > 0) * (1 - exp(c^(1/h)-( (k2*eps)^h + c )^(1/h) ))
}

loss_asymm_pop = function(beta, y, X, k, c, g = 1, link = function(x) x){
  eps = y - link(X%*%beta)
  k_vec = k[(eps > 0) + 1]
  h = 2*g
  sum(1 - exp(c^(1/h)-( (k_vec*eps)^h + c )^(1/h) ))
}


## Derivatives ----
score_pop = function(idx, beta, y, X, k, c, freq = 1, link, link_prime){
  eta = drop(X %*% beta)
  eps = y - link(eta)
  k_vec = k[(eps>0)+1]
  u = k_vec*eps
  
  F1 = k_vec * exp(sqrt(c) - sqrt(u^2+c)) * u/sqrt(u^2+c) * link_prime(eta)
  # out = colSums(X[ , idx, drop = F]*F1)
  out = drop(crossprod(X[ , idx, drop = F], freq * F1))
  return(out)
}

H_pop = function(idx, beta, y, X, X2, k, c, freq = 1, 
                 link, link_prime, link_second){
  eta = drop(X %*% beta)
  eps = (y - link(eta))
  k_vec = k[(eps>0)+1]
  u = k_vec*eps
  
  F1 = k_vec^2 * exp(sqrt(c) - sqrt(u^2+c)) * (
    link_prime(eta)^2 * ((u^2)/(u^2+c) - c/((u^2+c)^(3/2))) +
      link_second(eta) * eps/sqrt(u^2+c)
  )
  # out = colSums(X2[ , idx, drop = F] * F1)
  out = drop(crossprod(X2[ , idx, drop = F], freq * F1))
  
  return(out)
}


## DR Decomposition ----
# This function returns the nonlinear Demmler-Reinsch spline basis functions
# Construct Demmler Reinsch basis according to Bach and Klein 2024
construct_DR_basis <- function(x, d = 10, rescale = TRUE, type_of_splines = "ps") {
  
  n <- length(x)
  # set up P-spline design
  if (type_of_splines == "ps"){
    fit=smoothCon(s(x, bs="ps", k=d+2, m = c(2, 2)), data=data.frame(x),scale.penalty=FALSE)
  } else if (type_of_splines == "gps"){
    require(gps.mgcv)
    fit=smoothCon(s(x, bs="gps", k=d+2, m = c(3, 2)), data=data.frame(x),scale.penalty=FALSE)
  } else {
    stop("Type of splines not supported")
  }
  sm = fit[[1]]
  
  B=fit[[1]]$X
  K=fit[[1]]$S[[1]]
  
  # set up constraint matrix
  # Orthogonality is enforced against both the intercept (column of 1s) and the
  # linear term. We use scale(x) (centred AND scaled), but this spans the same
  # 2D subspace as centred x: scale(x) = c0*1 + c1*(x - mean(x)), so the
  # nullspace of (1, scale(x))^T B is identical to that of (1, x-centred)^T B
  # and to that of original x, so (1, x)^T B.
  # The resulting DR basis Z is hence orthogonal to any linear/intercept effect.
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

# DR decomposition for continuous spatial components from Paul's code
construct_DR_Spatial<-function(x1, x2, k1 = 10, k2 = 10, normalize_cols = FALSE){
  
  require(geigen)
  
  n = length(x1)
  
  fit1=smoothCon(s(x1,bs="ps",k=k1,m=c(2,1)),data=data.frame(x1),scale.penalty=FALSE)
  fit2=smoothCon(s(x2,bs="ps",k=k2,m=c(2,1)),data=data.frame(x2),scale.penalty=FALSE)
  
  B1=fit1[[1]]$X
  B2=fit2[[1]]$X
  
  K1=fit1[[1]]$S[[1]]
  K2=fit2[[1]]$S[[1]]
  
  B=tensor.prod.model.matrix(list(B1,B2))
  K=Reduce("+",tensor.prod.penalties(list(K1,K2)))
  
  G=crossprod(B)/n
  
  C=t(rep(1,n))%*%B
  # NOTE: in the spatial extension we only enforce orthogonality against the
  # intercept (span of 1) and deliberately NOT against any linear term: the
  # spatial effect is defined as purely smooth/nonlinear and no linear spatial
  # covariate is included in the predictor, so orthogonality to a (meaningless)
  # linear term is not required.

  
  V0=svd(C,nv=k1*k2)$v[,2:(k1*k2)]
  
  G2=t(V0)%*%G%*%V0
  K2=t(V0)%*%K%*%V0
  
  fit=geigen(G2,K2)
  
  A=fit$vectors
  
  Z=B%*%V0%*%A
  
  frob=sqrt(sum(diag(t(Z)%*%Z/n)))
  
  Z=Z/frob ### unit Frobenius norm
  
  # sum(diag(t(Z)%*%Z)) ## check
  
  Trafo=V0%*%A/frob
  
  ret=list()
  ret[[1]]=Z
  ret[[2]]=Trafo
  ret[[3]]=B
  ret[[4]] = B1
  ret[[5]] = B2
  ret[[6]] = fit1[[1]]$knots
  ret[[7]] = fit2[[1]]$knots
  names(ret)=c("Z","Trafo","B", "B1", "B2", "knots1", "knots2")
  
  return(ret)
}

## IWLS utils ----
### IWLS single update ----
iwls_beta_update <- function(idx, beta_curr, u_curr,
                             y, X, X2, freq = 1,
                             link, link_prime, link_second,
                             k_post, k_deriv, 
                             c_post, c_deriv, 
                             w, tau, tau_idx,
                             hess_eps = 1e-3,
                             step = 1){
  # Current state
  cur <- eval_logfullcond_beta(
    idx = idx, beta = beta_curr, u = u_curr,
    X = X, y = y, X2 = X2, freq = freq,
    link = link, link_prime = link_prime, link_second = link_second,
    k_post = k_post, k_deriv = k_deriv, 
    c_post = c_post, c_deriv = c_deriv,
    w = w,
    tau = tau, tau_idx = tau_idx
  )
  
  Hc <- make_negdef(cur$H, eps = hess_eps)
  Sigma_c <- - 1/Hc
  mu_c <- drop(beta_curr[idx] + step*Sigma_c*cur$grad)
  
  # Proposal draw
  beta_prop <- beta_curr
  beta_prop[idx] <- rnorm(length(idx), mean = mu_c, sd = sqrt(Sigma_c))
  
  # Proposed state
  prop <- eval_logfullcond_beta(
    idx = idx, beta = beta_prop, u = u_curr,
    X = X, y = y, X2 = X2, freq = freq,
    link = link, link_prime = link_prime, link_second = link_second,
    k_post = k_post, k_deriv = k_deriv,
    c_post = c_post, c_deriv = c_deriv,
    w = w,
    tau = tau, tau_idx = tau_idx
  )
  
  Hp <- make_negdef(prop$H, eps = hess_eps)
  Sigma_p <- -1/Hp
  mu_p <- drop(beta_prop[idx] + step*Sigma_p*prop$grad)
  
  # MH acceptance probability
  log_q_prop_given_cur <- sum(dnorm(beta_prop[idx], mean = mu_c, sd = sqrt(Sigma_c), log = T))
  log_q_cur_given_prop <- sum(dnorm(beta_curr[idx], mean = mu_p, sd = sqrt(Sigma_p), log = T))
  
  # Check if it's the target ratio or the proposal ratio that brings acceptance to zero
  if (as.numeric(get("verbose", envir = parent.frame())) > 2){
    if (get("m", envir = parent.frame()) %% get("print_step", envir = parent.frame()) == 0){
      cat("idx:", idx, "| Target diff:", round(prop$logpost - cur$logpost, 2), 
          "| Proposal diff:", round(log_q_cur_given_prop - log_q_prop_given_cur, 2), "\n")
    } 
  }
  
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


iwls_re_update <- function(u_curr, y, fe, freq = 1,
                           link, link_prime, link_second,
                           k_post, k_deriv, 
                           c_post, c_deriv, 
                           w, sigma2_u,
                           hess_eps = 1e-3,
                           step = 1){
  # Current state
  cur <- eval_logfullcond_re(
    y = y, fe = fe, u = u_curr,
    freq = freq,
    link,  link_prime, link_second,
    k_post, k_deriv, 
    c_post, c_deriv,
    w, sigma2_u
  )
  
  Hc <- make_negdef(cur$H, eps = hess_eps)
  Sigma_c <- - 1/Hc
  mu_c <- drop(u_curr + step*Sigma_c*cur$grad)
  
  # Proposal draw
  u_prop <- rnorm(1, mean = mu_c, sd = sqrt(Sigma_c))
  
  # Proposed state
  prop <- eval_logfullcond_re(
    y = y, fe = fe, u = u_prop,
    freq = freq,
    link,  link_prime, link_second,
    k_post, k_deriv, 
    c_post, c_deriv,
    w, sigma2_u
  )
  
  Hp <- make_negdef(prop$H, eps = hess_eps)
  Sigma_p <- -1/Hp
  mu_p <- drop(u_prop + step*Sigma_p*prop$grad)
  
  # MH acceptance probability
  log_q_prop_given_cur <- sum(dnorm(u_prop, mean = mu_c, sd = sqrt(Sigma_c), log = T))
  log_q_cur_given_prop <- sum(dnorm(u_curr, mean = mu_p, sd = sqrt(Sigma_p), log = T))
  
  # Check if it's the target ratio or the proposal ratio that brings acceptance to zero
  if (as.numeric(get("verbose", envir = sys.frame(sys.parent(2)))) > 2){
    if (get("m", envir = parent.frame()) %% get("print_step", envir = parent.frame()) == 0){
      cat("RE update  |  Target diff:", round(prop$logpost - cur$logpost, 2), 
          "| Proposal diff:", round(log_q_cur_given_prop - log_q_prop_given_cur, 2), "\n")
    } 
  }
  
  log_alpha <- prop$logpost - cur$logpost +
    log_q_cur_given_prop - log_q_prop_given_cur
  
  if (log(runif(1)) < log_alpha) {
    list(
      u = u_prop,
      accepted = 1,
      log_alpha = log_alpha,
      logpost = prop$logpost
    )
  } else {
    list(
      u = u_curr,
      accepted = 0,
      log_alpha = log_alpha,
      logpost = cur$logpost
    )
  }
}


### Log-full-conditional ----
eval_logfullcond_beta <- function(idx, beta, 
                                  X, y, X2 = X2, freq = 1, u,
                                  link,  link_prime, link_second,
                                  k_post, k_deriv, 
                                  c_post, c_deriv,
                                  w, tau, tau_idx) {
  beta <- drop(beta)
  loss_fun <- loss_asymm
  
  # --- log-likelihood ---
  eta <- X %*% beta + u
  r <- drop(y - link(eta))
  loglik <- -w * sum(freq * loss_fun(r, k = k_post, c = c_post))
  
  grad_ll = w * score_pop(idx = idx, beta = beta, y = y, X = X, k = k_deriv, c = c_deriv, freq = freq,
                          link = link, link_prime = link_prime) 
  H_ll = w * H_pop(idx = idx, beta = beta, y = y, X = X, X2 = X2, k = k_deriv, c = c_deriv, freq = freq,
                   link = link, link_prime = link_prime, link_second = link_second)
  
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
    logprior <- -0.5 * sum(beta[idx]^2 / (tau_tmp^2))
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


eval_logfullcond_re <- function(y, fe, u,
                                freq = 1,
                                link,  link_prime, link_second,
                                k_post, k_deriv, 
                                c_post, c_deriv,
                                w, sigma2_u) {
  
  loss_fun <- loss_asymm
  
  # --- log-likelihood ---
  r <- drop(y - link(fe + u))
  loglik <- -w * sum(freq * loss_fun(r, k = k_post, c = c_post))
  
  grad_ll = w * score_pop(idx = 1, beta = u, 
                          y = y, X = matrix(rep(1, length(y)), ncol = 1), freq = freq,
                          k = k_deriv, c = c_deriv, 
                          link = link, link_prime = link_prime) 
  H_ll = w * H_pop(idx = 1, beta = u, 
                   y = y, X = matrix(rep(1, length(y)), ncol = 1), 
                   X2 = matrix(rep(1, length(y)), ncol = 1), freq = freq,
                   k = k_deriv, c = c_deriv, 
                   link = link, link_prime = link_prime, link_second = link_second)
  
  # --- Gaussian prior ---
  Omega <- 1 / (sigma2_u)
  logprior <- -0.5 * u^2 / (sigma2_u)
  grad_prior <- - (Omega * u)
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
make_negdef <- function(H, eps = 1e-3) {
  H[H>=0] = -eps
  return(H)
}


## Slice sampler for sigma_u (sd rand. eff.) DO NOT USE ----
update_sigma_u <- function(sigma_u_curr, u, 
                           A, nu = 1, 
                           slice_w = 1, m = 30) {
  
  
  # Work on z = log(sigma_u)
  z_curr <- log(sigma_u_curr)
  
  # Log posterior of sigma_u
  # p(sigma_u | u) propto p(u | sigma_u) p(sigma_u)
  # Transformation of density for z: p(z | u) = p(exp(z) | u) * exp(z)
  # - z = g(sigma_u) = log(sigma_u)
  # - d/dz g^{-1}(z) = d/dz exp(z) = exp(z)
  
  logpost_z = function(z) {
    
    sigma = exp(z)
    
    # Random-effects likelihood: u_i | sigma ~ N(0, sigma^2)
    loglik_u = -length(u) * log(sigma) - sum(u^2) / (2 * sigma^2)
    
    # Half-t prior.
    # Density of t with nu df, location 0 and scale A.
    # The t density should be multiplied by 2 because it's half-t, but
    # the factor 2 is constant and can be omitted.
    logprior_sigma = dt(sigma / A, df = nu, log = TRUE) - log(A)
    
    # Jacobian in logpost is log(exp(z)) = z
    log_jacobian = z
    
    loglik_u + logprior_sigma + log_jacobian
  }
  
  
  # Slice level: y ~ Uniform(0, logpost_z) --> log_y = logpost_z + log(U(0, 1)) = logpost_z - rexp(1)
  log_y = logpost_z(z_curr) - rexp(1)
  
  
  # Step out procedure
  u_rand <- runif(1, 0, slice_w)
  L = z_curr - u_rand
  R = L + slice_w
  J = floor(runif(1, 0, m))
  K = (m - 1) - J
  
  while (J > 0 && log_y < logpost_z(L)) {
    L = L - slice_w
    J = J - 1
  }
  
  while (K > 0 && logpost_z(R) > log_y) {
    R = R + slice_w
    K = K - 1
  }
  
  
  # Shrinkage
  repeat {
    
    # Sample z_prop from the interval [L, R]
    z_prop <- runif(1, L, R)
    
    # If it is in the slice, accept it and break the loop
    if (logpost_z(z_prop) >= log_y) {
      break
    }
    
    # If z_prop it's outside the slice, shrink the interval [L, R] accordingly
    if (z_prop < z_curr) {
      L <- z_prop
    } else {
      R <- z_prop
    }
  }
  
  
  # Return
  sigma_new <- exp(z_prop)
  
  return(sigma_new)
}


eval_logfullcond_SDP = function(log_tau2, K, rank_K, beta, theta){
  # u is the notation used in Klein and Kneib (2016)
  # I use log_tau as argument because for me u is the random effect
  u  = log_tau2
  aux = drop(beta%*%K%*%beta)/(2 * exp(u))
  logpost = u/2*(1 - rank_K) - aux - (exp(u)/theta)^(1/2)
  fd = (1 - rank_K)/2 + aux - (exp(u)/(4*theta))^(1/2)
  sd = -aux - (exp(u)/(16*theta))^(1/2)
  return(list(logpost = logpost, grad = fd, hess = sd))
}

iwls_SDP_update = function(log_tau2, K, rank_K, beta, theta, step = 1){
  
  # Current state
  cur <- eval_logfullcond_SDP(log_tau2 = log_tau2, K = K, rank_K = rank_K, beta = beta, theta = theta)
  
  Hc <- cur$hess
  Sigma_c <- - 1/Hc
  mu_c <- drop(log_tau2 + step*Sigma_c*cur$grad)
  
  # Proposal draw
  log_tau2_prop <- rnorm(1, mean = mu_c, sd = sqrt(Sigma_c))
  
  # Proposed state
  prop <- eval_logfullcond_SDP(log_tau2 = log_tau2_prop, K = K, rank_K = rank_K, 
                               beta = beta, theta = theta)
  
  Hp <- prop$hess
  Sigma_p <- -1/Hp
  mu_p <- drop(log_tau2_prop + step*Sigma_p*prop$grad)
  
  # MH acceptance probability
  log_q_prop_given_cur <- sum(dnorm(log_tau2_prop, mean = mu_c, sd = sqrt(Sigma_c), log = T))
  log_q_cur_given_prop <- sum(dnorm(log_tau2, mean = mu_p, sd = sqrt(Sigma_p), log = T))
  
  # Check if it's the target ratio or the proposal ratio that brings acceptance to zero
  if (as.numeric(get("verbose", envir = sys.frame(sys.parent(2)))) > 2){
    if (get("m", envir = sys.frame(sys.parent(2))) %% get("print_step", envir = sys.frame(sys.parent(2))) == 0){
      cat("SDP  |  Target diff:", round(prop$logpost - cur$logpost, 2), 
          "| Proposal diff:", round(log_q_cur_given_prop - log_q_prop_given_cur, 2), "\n")
    } 
  }
  
  log_alpha <- prop$logpost - cur$logpost +
    log_q_cur_given_prop - log_q_prop_given_cur
  
  if (log(runif(1)) < log_alpha) {
    list(
      log_tau2 = log_tau2_prop,
      accepted = 1,
      log_alpha = log_alpha,
      logpost = prop$logpost
    )
  } else {
    list(
      log_tau2 = log_tau2,
      accepted = 0,
      log_alpha = log_alpha,
      logpost = cur$logpost
    )
  }
}

## Random effects update ----
update_re <- function(u, sigma2_u, id, theta,
                      y, X, freq,
                      beta,
                      link, link_prime, link_second,
                      k_post, k_deriv, c_post, c_deriv, w,
                      hess_eps = 1e-3, step = 1
) {
  
  fe = X %*% beta
  id_unique = unique(id)
  n_id = length(id_unique)
  
  u_new = rep(NA, n_id)
  acc_u = rep(NA, n_id)
  
  # Update each random effect u_i separately
  for (i in 1:n_id){
    id_i = id_unique[i]
    y_i = y[id == id_i]
    fe_i = fe[id == id_i]
    freq_i = freq[id == id_i]
    u_i = u[i]
    up_re = iwls_re_update(u_curr = u_i, y = y_i, fe = fe_i, freq = freq_i,
                           link = link, link_prime = link_prime, link_second = link_second,
                           k_post = k_post, k_deriv = k_deriv, 
                           c_post = c_post, c_deriv = c_deriv, 
                           w = w, sigma2_u = sigma2_u,
                           hess_eps = hess_eps,
                           step = step)
    
    acc_u[i] = up_re$accepted
    u_new[i] = up_re$u
  }
  
  # Update sigma2_u using the Scale Dependent Prior (SDP) by Klein&Kneib (2016, BA)
  up_sigma2_u = iwls_SDP_update(log_tau2 = log(sigma2_u), 
                                K = diag(1, n_id), rank_K = n_id, 
                                beta = u_new, theta = theta, step = 1)
  
  sigma2_u_new = exp(up_sigma2_u$log_tau2)
  acc_sigma2_u = up_sigma2_u$accepted
  
  list(u = u_new, sigma2_u = sigma2_u_new, acc_u = acc_u, acc_sigma2_u = acc_sigma2_u)
}


## Link functions and derivatives ----
### Logit ----
logit_link = function(eta) {
  return(1 / (1 + exp(-eta)))
}

logit_prime = function(eta) {
  ll = logit_link(eta)
  return(ll * (1 - ll))
}

logit_second = function(eta) {
  ll = logit_link(eta)
  return(ll * (1 - ll) * (1 - 2 * ll))
}

## Plots ----
### Spatial splines ----
plot_spatial_splines = function(coord, grid, bases, idx_to_plot, normalize = FALSE,
                                type, title = NULL, boundaries = NULL){
  require(ggplot2); require(dplyr); require(tidyr)
  
  if (is.null(title)){title = type}
  
  # Prepare data in a df
  df_points <- data.frame(x1 = coord[, 1], x2 = coord[, 2])
  df_grid <- data.frame(x1 = grid[, 1], x2 = grid[, 2])
  
  # Extract bases to plot
  if (normalize){
    bases = apply(bases, 2, function(x) (x - min(x)) / (max(x) - min(x)))
  }
  df_bases <- as.data.frame(bases[, idx_to_plot])
  colnames(df_bases) <- paste0("Basis_", idx_to_plot)
  
  df_plot <- cbind(df_grid, df_bases) %>%
    pivot_longer(cols = starts_with("Basis_"), names_to = "Basis", values_to = "Value") |> 
    mutate(Basis = factor(Basis, levels = paste0("Basis_", idx_to_plot)))
  
  if (type == "DR"){
    M = max(abs(df_plot$Value))
    breaks = round(seq(-M-1e-4, M, length = 8), 3)
    m1 = breaks[4]; m2 = breaks[5]
    breaks = c(breaks[1:3], seq(m1, m2, length = 6), breaks[6:8])
    palette = RColorBrewer::brewer.pal(n = 11, name = "RdBu")[11:1]
  } else if (type == "B-Splines"){
    M = max(df_plot$Value)
    m = min(df_plot$Value)
    breaks = round(seq(m-1e-4, M, length = 12), 3)
    palette = colorRampPalette(c("white", "#a4161a"))(11)
  } else {
    stop("Type not supported")
  }
  
  
  ncol_wrap = ceiling(sqrt(length(idx_to_plot)) * 1.5)
  
  # Plot with contour lines + points
  plt = ggplot() +
    geom_contour_filled(data = df_plot, aes(x = x1, y = x2, z = Value),
                        alpha = 0.8,
                        breaks = breaks) +
    geom_contour(data = df_plot, aes(x = x1, y = x2, z = Value),
                 color = "black", linewidth = 0.2, 
                 breaks = breaks) +
    geom_point(data = df_points, aes(x = x1, y = x2), 
               color = "red", size = 1, shape = 16, alpha = 0.75) +
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_continuous(limits = c(0, 1)) +
    # scale_fill_brewer(palette = palette, direction = -1) +
    scale_fill_manual(values = palette, name = "") +
    facet_wrap(~ Basis, ncol = ncol_wrap) +
    coord_equal() +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(title = title, x = "x1", y = "x2", fill = "Value of Basis Function")
  
  if (!is.null(boundaries)){
    plt = plt +  geom_path(data = data.frame(x = boundaries[, 1], y = boundaries[, 2]),
                           mapping = aes(x = x, y = y))
  }
  
  plt
}


plot_effect_spat<-function(Xspat, beta_est, Trafo, xlim=c(-0.2,1.2), ylim=c(-0.2,1.2), 
                           convex_hull = FALSE, concave_hull = FALSE,
                           ncol = 30, ngrid = 100, k1 = 5, k2 = 5,
                           knots1, knots2,
                           name = "", cex = 1){
  
  require(plot3D)
  
  if (convex_hull & concave_hull){
    stop("Please choose either convex or concave hull, not both.")
  }
  
  x1=Xspat[,1]
  x2=Xspat[,2]
  
  t1=seq(min(x1),max(x1),length=ngrid)
  t2=seq(min(x2),max(x2),length=ngrid)
  
  TT=expand.grid(t1,t2)
  
  t1n=TT[,1]
  t2n=TT[,2]
  
  BT1=smoothCon(s(t1n,bs="ps",k=k1,m=c(2,1)), knots = data.frame(knots1),
                data=data.frame(t1n),scale.penalty=FALSE)[[1]]$X
  BT2=smoothCon(s(t2n,bs="ps",k=k2,m=c(2,1)), knots = data.frame(knots2),
                data=data.frame(t2n),scale.penalty=FALSE)[[1]]$X
  
  BT=tensor.prod.model.matrix(list(BT1,BT2))
  
  ##################################### overall effect
  
  XT=BT%*%Trafo
  f.hat=XT%*%beta_est
  fhat=matrix(f.hat,ngrid,ngrid)
  
  
  
  if (convex_hull){
    require(geometry)
    ch=convhulln(Xspat)
    inside=inhulln(ch,p=as.matrix(TT))
    fhat[!inside]=NA
  }
  
  if (concave_hull){
    require(concaveman)
    require(sp)
    ch=concaveman(Xspat, concavity = 2.3)
    inside=point.in.polygon(TT[,1],TT[,2],ch[,1],ch[,2])
    fhat[inside==0]=NA
  }
  
  colors = colorRampPalette(c("#4575b4", "#f7f7f7", "#d73027"))(ncol)
  M = max(abs(fhat), na.rm = TRUE)
  breaks = round(seq(-M, M, length.out = ncol + 1), 2)
  image2D(t1,t2,z=fhat,contour=FALSE,main="Spatial effect", 
          col=colors, breaks = breaks,
          xlab="Easting",ylab="Northing",xlim=xlim,ylim=ylim)
  points(x1,x2,pch=16,cex=cex)
}

################ this function plots the estimated overall effects for all continuous covariates

plot_continuous_effects<-function(fit,ylim=c(-1,1))
{
  
  X=fit$design$X
  indices=fit$design$indices
  dim=fit$design$dim
  beta=fit$samples$beta
  Xcont=fit$design$Xcont
  
  if(is.null(ncol(fit$design$Xdummy))){
    nd=0
  } else {
    nd=ncol(fit$design$Xdummy) ## number of dummy covariates  
  }
  
  nc=ncol(fit$design$Xcont) ## number of continuous covariates
  
  for(j in 1:nc){
    
    # xtemp=Xcont[,j]
    # betatemp=cbind(beta[,j+nd],beta[,indices[[j+(nd+nc)]]])
    # Xtemp=cbind(X[,j+nd],X[,indices[[j+(nd+nc)]]])
    # 
    # F=betatemp%*%t(Xtemp)
    # 
    # f.hat=apply(F,2,FUN=mean)
    # low=apply(F,2,FUN=quantile,prob=0.025)
    # up=apply(F,2,FUN=quantile,prob=1-0.025)
    # 
    # plot(xtemp,f.hat,ylim=ylim,type="p")
    # points(xtemp,low,lty="dashed")
    # points(xtemp,up,lty="dashed")
    
    x=Xcont[,j]
    t=seq(min(x),max(x),length=200)
    
    BT=smoothCon(s(t,bs="ps",k=dim[j+nd+nc]+2),data=data.frame(t))[[1]]$X
    
    m=mean(x)
    s=sd(x)
    
    XT=BT%*%fit$design$Trafo[[j]]
    XT=cbind((t-m)/s,XT)
    
    betatemp=cbind(beta[,j+nd],beta[,indices[[j+(nd+nc)]]])
    
    F=betatemp%*%t(XT)
    
    f.hat=apply(F,2,FUN=mean)
    low=apply(F,2,FUN=quantile,prob=0.025)
    up=apply(F,2,FUN=quantile,prob=1-0.025)
    
    plot(t,f.hat,ylim=ylim,type="l")
    lines(t,low,lty="dashed")
    lines(t,up,lty="dashed")
    
  }
}
