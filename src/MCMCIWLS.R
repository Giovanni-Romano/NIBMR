

# ============================================================
# One joint IWLS-MH beta update
# ============================================================

iwls_beta_update <- function(beta_curr, X, y, k, c, w, tau, tau_idx,
                             asymm = TRUE, hess_eps = 1e-8,
                             score_i, H_i) {
  # Current state
  cur <- eval_logfullcond_beta(
    beta = beta_curr, X = X, y = y, k = k, c = c, w = w,
    tau = tau, tau_idx = tau_idx, asymm = asymm,
    score_i = score_i, H_i = H_i
  )
  
  Hc <- make_negdef(cur$H, eps = hess_eps)
  Sigma_c.tmp <- -solve(Hc)
  Sigma_c = (Sigma_c.tmp + t(Sigma_c.tmp))/2
  mu_c <- drop(beta_curr + drop(Sigma_c %*% cur$grad))
  
  # Proposal draw
  beta_prop <- drop(rmvnorm(1, mean = mu_c, sigma = Sigma_c, checkSymmetry = F))
  
  # Proposed state
  prop <- eval_logfullcond_beta(
    beta = beta_prop, X = X, y = y, k = k, c = c, w = w,
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


# ============================================================
# Log full conditional, gradient, Hessian
# ============================================================

# IMPORTANT:
# score_i(beta, y_i, x_i, k, c) should be the gradient of the
# log-likelihood contribution for observation i
#
# H_i(beta, y_i, x_i, k, c) should be the Hessian of the
# log-likelihood contribution for observation i

eval_logfullcond_beta <- function(beta, X, y, k, c, w, tau, tau_idx,
                                  asymm = TRUE,
                                  score_i, H_i) {
  beta <- drop(beta)
  n <- nrow(X)
  P <- length(beta)
  
  loss_fun <- if (asymm) loss_asymm2 else loss_symm
  
  # ----- log-likelihood -----
  r <- drop(y - X %*% beta)
  loglik <- -w * sum(loss_fun(r, k = k, c = c))
  
  grad_ll <- rep(0, P)
  H_ll <- matrix(0, nrow = P, ncol = P)
  
  for (i in seq_len(n)) {
    grad_ll <- grad_ll + w * drop(score_i(beta, y[i], X[i, ], k, c))
    H_ll <- H_ll + w * H_i(beta, y[i], X[i, ], k, c)
  }
  
  # ----- Gaussian prior -----
  # No prior on intercept beta[1]
  # Prior on beta[-1] with group-specific variances tau[tau_idx]^2
  tau_tmp <- tau[tau_idx]  # length P-1
  
  Omega_diag <- c(0, 1 / (tau_tmp^2))
  Omega <- diag(Omega_diag, nrow = P)
  
  logprior <- -0.5 * sum(beta[-1]^2 / (tau_tmp^2))
  grad_prior <- - drop(Omega %*% beta)
  H_prior <- - Omega
  
  # ----- full conditional -----
  logpost <- loglik + logprior
  grad <- grad_ll + grad_prior
  H <- H_ll + H_prior
  
  list(
    logpost = as.numeric(logpost),
    grad = drop(grad),
    H = H
  )
}



# ============================================================
# Utilities
# ============================================================

# Enforce symmetry and negative definiteness
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

# Log-density of multivariate normal
log_dmvnorm_chol <- function(x, mean, Sigma) {
  x <- drop(x)
  mean <- drop(mean)
  
  R <- chol(Sigma)
  z <- backsolve(R, x - mean, transpose = TRUE)
  quad <- sum(z^2)
  logdet <- 2 * sum(log(diag(R)))
  
  -0.5 * (length(x) * log(2 * pi) + logdet + quad)
}

# ============================================================
# Main MCMC
# ============================================================

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
  
  # ----------------------------------------------------------
  # Calibrate w exactly as in your original code
  # ----------------------------------------------------------
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
  
  # ----------------------------------------------------------
  # Storage
  # ----------------------------------------------------------
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
  
  # ----------------------------------------------------------
  # Initial values
  # ----------------------------------------------------------
  beta_curr <- drop(beta0)
  tau <- rep(100, 2 * p)
  
  # Maps beta[-1] to the corresponding tau component:
  # first p coefficients are linear effects,
  # next sum(d) coefficients are nonlinear groups
  tau_idx <- c(1:p, p + rep(seq_len(p), d))
  
  acc <- 0L
  
  # ----------------------------------------------------------
  # MCMC loop
  # ----------------------------------------------------------
  for (m in seq_len(niter)) {
    
    if (verbose && m %% max(1, floor(niter / 10)) == 0) {
      cat(round(100 * m / niter), "%\t")
    }
    
    # ----------------------------
    # Joint IWLS update for beta
    # ----------------------------
    up <- iwls_beta_update(
      beta_curr = beta_curr,
      X = X, y = y, k = c(1, 1), c = c, w = w,
      tau = tau, tau_idx = tau_idx,
      asymm = asymm, hess_eps = hess_eps,
      score_i = score_i, H_i = H_i
    )
    
    beta_curr <- up$beta
    beta_sample[m, ] <- beta_curr
    alfa_sample[m] <- exp(min(up$log_alpha, 0))
    logpost_sample[m] <- up$logpost
    acc <- acc + up$accepted
    
    # ----------------------------
    # Update tau exactly as before
    # ----------------------------
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

out_MCMC = MCMC_IWLS(
  niter = 2000,
  X = X_scaled, y = y, k = k_multi, c = 1e-3, p = p, d = d, beta0 = rep(0, 1+p+sum(d)), #out_optim_scaled[, "SANN"],
  a_tau = 1e-3, b_tau = 1e-3,
  asymm = TRUE, debug = TRUE, verbose = TRUE,
  hess_eps = 1e-4,
  seed = 666
)
