# MCMC SAMPLER ----
#' MCMC sampler (IWLS proposals) for mode regression on a single additive
#' predictor, without random effects.
#'
#' Runs the Gibbs-within-Gibbs sampler of the notes for the model that has no
#' random effects (used e.g. by the simulation study scripts). At each iteration
#' it (i) updates the intercept and each coefficient block (one linear plus one
#' DR-nonlinear block per continuous covariate) with iwls_beta_update, (ii)
#' samples the block-specific variances tau_j^2 from their inverse-Gamma
#' full conditionals (conjugate inverse-gamma with shape a_tau[j]+size/2 and
#' rate b_tau[j]+sum(beta_tmp^2)/2), and (iii) stores the draws.
#'
#' The generalized-Bayes weight `w` is estimated by K-fold cross-validation of
#' the empirical loss: with the unit-information-loss criterion (UILC) the
#' numerator is the number of informative (prior-tasked) coefficients divided
#' by 2, i.e. w_num = (p + sum(d))/2 = (P - 1)/2, divided by the average
#' out-of-sample loss (see the BMR_NIL notes, UILC section).
#'
#' @param niter Integer, number of MCMC iterations.
#' @param X Numeric design matrix Phi = (1, x, Z), n x P with P = 1 + p + sum(d).
#' @param y Numeric response vector of length n.
#' @param k Flatness parameter k_ni of the nested interval loss for the posterior
#'   (target) distribution. Passed as the `k_post` argument of iwls_beta_update.
#' @param k_deriv Flatness parameter used only for the smoother proposal loss
#'   (gradient/Hessian), passed to iwls_beta_update.
#' @param c Differentiability parameter c_n of the nested interval loss for the
#'   posterior (target) distribution. Passed as the `c_post` argument.
#' @param c_deriv Differentiability parameter used only for the smoother proposal
#'   loss, passed to iwls_beta_update.
#' @param p Number of continuous covariates (linear + nonlinear DR blocks).
#' @param d Integer vector of length p with the DR dimension per covariate
#'   (number of nonlinear basis functions per smooth).
#' @param beta0 Numeric vector of length P with starting values for beta.
#' @param a_tau Numeric vector (length 2*p) of inverse-gamma shape parameters
#'   for the block variances.
#' @param b_tau Numeric vector (length 2*p) of inverse-gamma rate parameters.
#' @param debug Logical, reserved for extra diagnostics (currently unused).
#' @param verbose Integer verbosity level: 0 quiet, 1 progress by 10% steps,
#'   2 per-`print_step` summaries.
#' @param print_step Number of iterations between verbose progress reports.
#' @param hess_eps Small positive constant used by make_negdef (default 1e-3).
#' @param K_fold Number of folds for the CV estimate of `w` (default 50).
#' @param seed Optional RNG seed for reproducibility of the CV folds.
#'
#' @return A list with
#'   \item{beta}{niter x P matrix of sampled coefficients (named columns).}
#'   \item{tau}{niter x (2*p) matrix of sampled block standard deviations.}
#'   \item{acc_prop}{Vector of per-block acceptance rates.}
#'   \item{logpost}{niter x (1+p) matrix of log-full-conditionals per update.}
#'   \item{w}{The estimated generalized-Bayes weight.}
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
  # Generalized Bayes (Bissiri, Holmes & Walker 2016, JRSSB).
  # Unit-information-loss criterion (UILC): set w so ONE observation
  # carries the same "information" as the unit-information prior.
  #
  # Closed-form UILC numerator (see BMR_NIL notes, Section on w/ UILC):
  #   w_num = E_beta[ log pi(beta.hat^prior) / pi(beta) ]
  #         = (number of informative coefficients) / 2.
  # Each independent Gaussian N(0, tau_j^2) coefficient contributes
  # E[ beta_j^2 / (2 tau_j^2) ] = 1/2 to this expectation; the flat
  # intercept prior contribution cancels (intercept has no prior).
  #
  # For this sampler the informative coefficients are the p linear slope
  # coefficients plus the sum(d) nonlinear DR-basis coefficients, so
  #   w_num = (p + sum(d)) / 2.
  # Equivalently w_num = (P - 1) / 2, since P = 1 + p + sum(d) counts exactly
  # one flat intercept plus those informative coefficients.
  w_num <- (p + sum(d)) / 2
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
#' MCMC sampler (IWLS proposals) for the bikes application: additive mode
#' regression with binary, continuous (linear + DR-nonlinear) and spatial terms
#' plus random intercepts.
#'
#' Runs the full Gibbs-within-Gibbs sampler used by the application on the bikes
#' data. The linear predictor is
#'   eta = beta0 + X_dummy beta_d + X_cont beta_l + Z_cont beta_nl + Z_spat beta_s + u_id,
#' with u_i ~ N(0, sigma2_u) a random intercept per group and sigma2_u updated
#' under the Scale-Dependent Prior (see update_re / iwls_SDP_update). At each
#' iteration the sampler (i) IWLS-updates the intercept, each binary/linear block
#' and each DR-nonlinear block, plus the spatial coefficients in blocks of five
#' (with step/4, to lift the otherwise too-low acceptance), (ii) sweeps the random
#' intercepts and sigma2_u, (iii) samples the block variances from their
#' inverse-Gamma full conditionals, and (iv) stores the draws. The
#' generalized-Bayes weight `w` is again set by K-fold cross-validation with the
#' UILC numerator w_num = (P-1)/2 (P = 1 + q + p + sum(d) + s).
#'
#' @param niter Integer, number of MCMC iterations.
#' @param X_dummy Numeric design matrix of binary/indicator covariates (n x q);
#'   may be NULL if there are none.
#' @param X_cont Numeric design matrix of continuous covariates (n x p), one
#'   column per linear effect.
#' @param Z_cont List or matrix of DR bases for the nonlinear continuous effects
#'   (the columns are appended to X as one block per covariate).
#' @param Z_spat Numeric matrix of DR spatial basis functions (n x s); NULL if
#'   no spatial term.
#' @param id Group index vector of length n mapping each observation to a
#'   random-intercept group.
#' @param y Numeric response vector of length n.
#' @param freq Numeric vector of observation weights (default all 1).
#' @param link Link function mu = link(eta) (default identity).
#' @param link_prime First derivative of the link.
#' @param link_second Second derivative of the link.
#' @param k Flatness parameter k_ni of the nested interval loss for the posterior
#'   (target) distribution; passed as `k_post`.
#' @param k_deriv Flatness parameter of the smoother proposal loss (gradient/
#'   Hessian); passed as `k_deriv`.
#' @param c Differentiability parameter c_n of the nested interval loss for the
#'   posterior (target) distribution; passed as `c_post`.
#' @param c_deriv Differentiability parameter of the smoother proposal loss;
#'   passed as `c_deriv`.
#' @param d Integer vector of length p with the DR dimension per continuous
#'   covariate (nonlinear basis functions per smooth).
#' @param beta0 Numeric vector of length P = 1 + q + p + sum(d) + s of starting
#'   values for beta.
#' @param a_tau Numeric vector (length q + 2*p + s) of inverse-gamma shape
#'   parameters for the block variances.
#' @param b_tau Numeric vector (length q + 2*p + s) of inverse-gamma rate
#'   parameters for the block variances.
#' @param theta_re Score of the SDP prior placed on the random-effects variance
#'   sigma2_u (scale of the prior).
#' @param debug Logical, reserved for extra diagnostics (currently unused).
#' @param verbose Integer verbosity level: 0 quiet, 1 progress by 10% steps,
#'   2 (or more) detailed per-`print_step` summaries.
#' @param print_step Number of iterations between verbose progress reports.
#' @param hess_eps Small positive constant used by make_negdef (default 1e-3).
#' @param step Scale multiplying the gradient in the IWLS proposal mean.
#' @param K_fold Number of folds for the CV estimate of `w` (default 50).
#' @param seed Optional RNG seed for reproducibility of the CV folds.
#'
#' @return A list with
#'   \item{beta}{niter x P matrix of sampled coefficients (named as in X).}
#'   \item{tau}{niter x (q + 2*p + s) matrix of sampled block standard deviations.}
#'   \item{u}{niter x n_id matrix of sampled random intercepts.}
#'   \item{sigma2u}{Vector (length niter) of sampled random-effects variances.}
#'   \item{acc_prop}{List with per-block beta acceptance (`beta`), per-group
#'     random-intercept acceptance (`re_u`) and random-effects-var acceptance
#'     (`re_sigma2u`), each divided by niter.}
#'   \item{logpost}{niter x (1 + q + 2*p + n_sp_up) matrix of log-full-conditionals.}
#'   \item{w}{The estimated generalized-Bayes weight.}
MCMC_IWLS_bikes <- function(niter,
                            X_dummy = NULL,
                            X_cont = NULL,
                            Z_cont = NULL,
                            Z_spat = NULL,
                            id,
                            y,
                            freq = NULL,
                            link = function(x) x,
                            link_prime = function(x) rep(1, length(x)),
                            link_second = function(x) rep(0, length(x)),
                            k, k_deriv, c, c_deriv,
                            d, # splines dimension for every continuous covariate
                            beta0,
                            a_tau, b_tau,
                            theta_re,
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
  # Same UILC as in MCMC_IWLS (see comment there): the closed-form numerator
  # equals the number of informative (Gaussian-prior) coefficients divided by 2.
  # Here the informative coefficients are q binary + p linear + sum(d) nonlinear
  # DR-basis + s spatial coefficients (all have Gaussian priors; only the flat
  # intercept beta0 has no prior), so w_num = (q + p + sum(d) + s) / 2.
  w_num <- (P - 1) / 2
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
  
  # --- Storage for random effects ---
  # u_i ~ N(0, sigma_u^2)
  # sigma_u ~ half-t with nu_re df, location 0 and scale A_re
  id_fac <- match(id, unique(id))
  n_id <- max(id_fac)
  u_sample   <- matrix(NA_real_, nrow = niter, ncol = n_id)
  sigma2u_sample <- numeric(niter)
  
  
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
  
  
  # Initial values for parameter-expanded random effects
  sigma2u_curr <-rweibull(n = 1, shape = 1/2, scale = theta_re)
  u_curr <- rnorm(n_id, mean = 0, sd = sqrt(sigma2u_curr))
  
  
  
  # ---  Storage for acceptance/rejection of MH
  n_sp_up = floor(s / 5) + 1*(s %% 5 > 0)
  acc <- rep(0L, 1+q+2*p+n_sp_up)
  acc_tmp <- acc
  
  acc_re_u <- rep(0L, n_id)
  acc_re_sigma2u  <- 0L
  acc_re_u_tmp <- acc_re_u
  acc_re_sigma2u_tmp  <- acc_re_sigma2u
  
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
      u_curr = u_curr[id_fac],
      y = y, X = X, X2 = X2,
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
          u_curr = u_curr[id_fac],
          y = y, X = X, X2 = X2,
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
          u_curr = u_curr[id_fac],
          y = y, X = X, X2 = X2,
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
          u_curr = u_curr[id_fac],
          y = y, X = X, X2 = X2,
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
          u_curr = u_curr[id_fac],
          y = y, X = X, X2 = X2,
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
      
      # Single block update for spatial effects (acc. prob. too low)
      # idx.b = idx.e + 1
      # idx.e = idx.e + s
      # idx = idx.b:idx.e
      # 
      # up <- iwls_beta_update(
      #   idx = idx,
      #   beta_curr = beta_curr,
      #   u_curr = u_curr[id_fac],
      #   y = y, X = X, X2 = X2,
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
    
    
    # --- Random intercept update ---
    up_re <- update_re(
      u = u_curr, sigma2_u = sigma2u_curr, id = id,
      theta = theta_re,
      X = X, y = y, freq = freq,
      beta = beta_curr,
      link = link, link_prime = link_prime, link_second = link_second,
      k_post = k, k_deriv = k_deriv, c_post = c, c_deriv = c_deriv,
      w = w,
      hess_eps = hess_eps,
      step = step)
    
    # Update current values for rand. eff.
    u_curr <- up_re$u
    sigma2u_curr <- up_re$sigma2_u
    
    # Save sampled values for rand. eff.
    u_sample[m, ] <- u_curr
    sigma2u_sample[m] <- sigma2u_curr
    
    # Save acceptance rates for rand. eff.
    acc_re_u_tmp <- acc_re_u_tmp + up_re$acc_u
    acc_re_u <- acc_re_u + up_re$acc_u
    acc_re_sigma2u_tmp <- acc_re_sigma2u_tmp + up_re$acc_sigma
    acc_re_sigma2u <- acc_re_sigma2u + up_re$acc_sigma
    
    
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
          "acc rate (last ", print_step, ") = ", paste(round(acc_tmp / print_step, 3), collapse = " "), "\n",
          "acc rate re [quartiles] (last ", print_step, ") = ", paste(round(quantile(acc_re_u_tmp / print_step, probs = c(0, 0.25, 0.5, 0.75, 1)), 3), collapse = " "), "\n",
          "acc rate re var (last ", print_step, ") = ", paste(round(acc_re_sigma2u_tmp / print_step, 3), collapse = " "), "\n",
          ifelse(verbose > 2, paste0(" | logpost = ", round(logpost_sample[m], 3)), ""),
          "\n\n\n",
          sep = "")
      time_last = now
      acc_tmp <- rep(0L, length(acc))
      acc_re_u_tmp <- rep(0L, n_id)
      acc_re_sigma2u_tmp <- 0L
    }
  }
  
  list(
    beta = beta_sample,
    tau = tau_sample,
    u = u_sample,
    sigma2u = sigma2u_sample,
    acc_prop = list(beta = acc / niter, re_u = acc_re_u / niter, re_sigma2u = acc_re_sigma2u / niter),
    logpost = logpost_sample,
    w = w
  )
}