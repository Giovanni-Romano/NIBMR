# ==============================================================================
# 01a-utils.R --- Utility functions for Generalized Bayesian mode regression
#                 via the nested interval loss (NIL).
#
# This file implements the building blocks of the sampler in 01b-gibbs.R:
#   * the asymmetric nested interval loss and its derivatives,
#   * the (univariate and spatial bivariate) Demmler-Reinsch basis,
#   * the IWLS-type Metropolis-Hastings proposals and log-full-conditionals,
#   * link functions and derivative helpers (for the logit link),
#   * plotting helpers.
#
# The mathematical notation follows the accompanying notes
# Two naming conventions deserve attention:
#   * k/c          refer to the loss "flatness" (k_ni) and "differentiability" (c_n) parameters
#                  used to define the posterior distribution.
#   * k_deriv/c_deriv refer to the parameters of a smoother loss used only to
#                  build the IWLS proposal (smoother proposal -> hopefully better
#                  exploration), as described in the notes.
#
# The code implements the same mathematics as in the notes; the description of
# every function below states the correspondence between argument names and the
# TeX symbols whenever they do not match literally.
#
# Required libraries:
#   * mgcv     (P-spline design matrices, smoothCon / PredictMat)
#   * geigen   (generalized eigenvalue decomposition for the DR basis)
#   * gps.mgcv (only if type_of_splines = "gps")
#   * ggplot2/dplyr/tidyr, plot3D, ... (only for the plotting helpers)
# ==============================================================================

# REQUIRED LIBRARIES ----
require(mgcv)
# library(mgcv, lib.loc = "/Users/giovanni/Library/Caches/org.R-project.R/R/renv/sandbox/macos/R-4.5/aarch64-apple-darwin20/4cd76b74")

# UTILS ----
## Loss function ----
#' Asymmetric nested interval loss (approximate NIL).
#'
#' Computes, pointwise, the continuous approximation of the 0-1 loss used to
#' target the conditional mode. With g = 1 (the value assumed in the notes) the
#' loss equals
#'   loss_asymm(e; k, c) = 1 - exp(sqrt(c) - sqrt( k[2]^2 e^2 + c ))   for e >= 0,
#'   loss_asymm(e; k, c) = 1 - exp(sqrt(c) - sqrt( k[1]^2 e^2 + c ))   for e <  0,
#' i.e. an (asymmetric) version of the loss \ell(e_i; c_n, k_ni) in the notes,
#' with k[1] = k_n/s_1 (negative side) and k[2] = k_n/s_2 (positive side).
#' The exponent h = 2*g generalises the "2" appearing in the notes; the default
#' g = 1 reproduces the notes exactly. Larger k makes the loss steeper (more
#' peaked around 0), larger c makes it smoother (more differentiable in 0).
#'
#' @param eps Numeric vector of residuals e_i = y_i - f(x_i).
#' @param k Numeric vector of length 2, k = c(k_n/s_1, k_n/s_2): the loss-shape
#'   parameter k_ni on the negative and positive side of the mode, respectively
#'   (TeX: k_ni, asymmetric version).
#' @param c Positive scalar, the smoothness parameter c_n of the loss (TeX: c_n).
#' @param g Scalar >= 0 controlling the exponent h = 2*g; default 1 (notes case).
#'
#' @return Numeric vector of the same length as eps with the individual losses.
#' @examples
#' e <- seq(-3, 3, length.out = 100)
#' plot(e, loss_asymm(e, k = c(1, 1), c = 1e-3), type = "l")
loss_asymm = function(eps, k, c, g = 1){
  k1 = k[1]
  k2 = k[2]
  h = 2*g
  (eps < 0) * (1 - exp(c^(1/h)-( (k1*eps)^h + c )^(1/h) )) + 
    (eps > 0) * (1 - exp(c^(1/h)-( (k2*eps)^h + c )^(1/h) ))
}

#' DEPRECATED. Population (whole-vector) asymmetric nested interval loss.
#'
#' Sums loss_asymm over all residuals r_i = y_i - link(X_i %*% beta). It was
#' used only by the (now outdated) simulation-study scripts; the sampler itself
#' works with the vectorised loss_asymm directly. Kept for reference.
#'
#' @param beta Regression-coefficient vector at which to evaluate the loss.
#' @param y Numeric response vector.
#' @param X Numeric design matrix (a column of 1s plus covariates).
#' @param k Numeric vector of length 2 with the asymmetric loss-shape parameter
#'   (k_ni on the negative and positive side).
#' @param c Positive scalar smoothness parameter c_n of the loss.
#' @param g Scalar >= 0 controlling the exponent h = 2*g; default 1.
#' @param link Link function mapping the linear predictor to mu (default identity).
#'
#' @return Scalar: the summed empirical loss.

loss_asymm_pop = function(beta, y, X, k, c, g = 1, link = function(x) x){
  eps = y - link(X%*%beta)
  k_vec = k[(eps > 0) + 1]
  h = 2*g
  sum(1 - exp(c^(1/h)-( (k_vec*eps)^h + c )^(1/h) ))
}


## Derivatives ----
#' Score (gradient of the log-likelihood) for a block of regression coefficients.
#'
#' Computes the gradient of the log-likelihood part of the log-full-conditional
#' for the coefficients in `idx` (a single block). This is the "loss" contribution
#' from the notes' gradient equation (for the selected block),
#'     u_loss = w * sum_i k_ni^2 w_i (e_i / s_i) phi_i,
#' multiplied by the link-function chain-rule factor link_prime(eta_i) and by the
#' observation weights `freq`. The prior contribution -Omega*beta is added later
#' in eval_logfullcond_beta. The arguments `k`/`c` correspond to the loss shape
#' k_ni and smoothness c_n, and are meant to be the proposal parameters
#' (k_deriv / c_deriv) when called from the sampler.
#'
#' @param idx Integer vector selecting the columns (coefficients) of X whose
#'   gradient is returned (typically one covariate block).
#' @param beta Full coefficient vector at which to evaluate the gradient.
#' @param y Numeric response vector.
#' @param X Numeric design matrix Phi = (1, x, Z) with one column per coefficient.
#' @param k Numeric vector of length 2: the asymmetric loss-shape parameter
#'   k_ni used for the gradient.
#' @param c Positive scalar: the loss smoothness parameter c_n.
#' @param freq Numeric vector of observation weights (default all 1).
#' @param link Link function mapping the linear predictor to mu.
#' @param link_prime First derivative of the link function.
#'
#' @return Numeric vector of length length(idx): the gradient of the
#'   log-likelihood with respect to the selected coefficients.
score_pop = function(idx, beta, y, X, k, c, freq = 1, link, link_prime){
  eta = drop(X %*% beta)
  eps = y - link(eta)
  k_vec = k[(eps>0)+1]
  u = k_vec*eps
  
  F1 = k_vec * exp(sqrt(c) - sqrt(u^2+c)) * u/sqrt(u^2+c) * link_prime(eta)
  out = drop(crossprod(X[ , idx, drop = F], freq * F1))
  return(out)
}

#' Diagonal proposal-covariance / Hessian helper for a block of coefficients.
#'
#' Computes the diagonal of the Hessian of the log-likelihood part of the
#' log-full-conditional for the coefficients in `idx` (one block), i.e. the
#' "loss" contribution from the notes' Hessian equation. The full Hessian is
#'   sum_i k_ni^2 w_i [ (k_ni^2 e_i^2)/s_i^2 - c_n/s_i^3 ] phi_i phi_i^T ;
#' the code evaluates only the DIAGONAL of phi_i phi_i^T, which is why a second
#' design matrix X2 = X^2 (element-wise square) is passed. This matches the
#' diagonal-Hessian IWLS proposal described in the notes. The prior contribution
#' -Omega is added in eval_logfullcond_beta. Arguments `k`/`c` are the proposal
#' loss parameters (k_deriv / c_deriv) when called from the sampler.
#'
#' @param idx Integer vector selecting the columns (coefficients) of X whose
#'   Hessian diagonal is returned.
#' @param beta Full coefficient vector at which to evaluate the Hessian.
#' @param y Numeric response vector.
#' @param X Numeric design matrix Phi = (1, x, Z).
#' @param X2 Numeric matrix, the element-wise square of X (X^2), used for the
#'   diagonal of phi_i phi_i^T.
#' @param k Numeric vector of length 2: the asymmetric loss-shape parameter.
#' @param c Positive scalar: the loss smoothness parameter c_n.
#' @param freq Numeric vector of observation weights (default all 1).
#' @param link Link function.
#' @param link_prime First derivative of the link function.
#' @param link_second Second derivative of the link function.
#'
#' @return Numeric vector of length length(idx): the diagonal of the Hessian of
#'   the log-likelihood with respect to the selected coefficients.
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
  out = drop(crossprod(X2[ , idx, drop = F], freq * F1))
  
  return(out)
}


## DR Decomposition ----
#' Nonlinear Demmler-Reinsch (DR) basis (Bach & Klein, 2024).
#'
#' Constructs the nonlinear DR basis functions for a smooth effect of a single
#' covariate x. Starting from the P-/GPS-spline design B (with penalty K), the
#' basis is reparametrised so that (i) it is constrained to be orthogonal to the
#' intercept and the (scaled) linear term, i.e. Z lives in the span of
#' the nonlinear part of f, and (ii) in the DR basis the penalty becomes the
#' identity. Concretely one first discards the two directions of B aligned with
#' the intercept and the linear term (the first two principal directions after
#' the constraint (1, scale(x))), then solves a generalised eigen-decomposition
#' of the Gram matrix G_tilde against the penalty K_tilde, and finally rescales
#' the columns as requested.
#'
#' @param x Numeric vector of covariate values (length n).
#' @param d Number of nonlinear basis functions (the DR dimension), default 10.
#'   This corresponds to d_j in the notes (the dimension of the nonlinear
#'   coefficients of the j-th effect); the full spline uses k = d + 2 basis.
#' @param rescale Logical; if TRUE the basis columns are normalised to unit norm
#'   (scaled by sqrt(n / sum(Z^2))), default TRUE.
#' @param type_of_splines Character; "ps" (P-splines) or "gps" (for which the
#'   package gps.mgcv is required).
#'
#' @return A list with elements
#'   \item{Z}{n x d matrix of DR basis functions evaluated at x.}
#'   \item{Trafo}{Transformation matrix such that Z = B %*% Trafo.}
#'   \item{smooth_object}{The mgcv smooth object used to reconstruct B on a grid.}
#'   \item{Gram}{The Gram matrix G_tilde (B^T B / n) in the DR coordinates.}
#'   \item{Penalty}{The transformed penalty matrix K_tilde.}
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
  
  B = fit[[1]]$X
  K = fit[[1]]$S[[1]]
  
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


#' DEPRECATED. Reconstruct a DR basis on a new grid of points.
#'
#' Evaluates the DR basis on a grid x_grid aligned with the basis Z_saved that
#' was actually used in a fit (aligning via a least-squares map M between the
#' freshly built basis and the one used at fitting time). Also returns the
#' scaled grid values. Not referenced by the application pipeline; it was used
#' ad-hoc for plotting/prediction. Kept for reference.
#'
#' @param x Training covariate values used to build the DR basis.
#' @param Z_saved n x d matrix of DR basis used in the fit to be replicated.
#' @param x_grid Numeric vector of new covariate values at which to evaluate.
#' @param d Dimension of the DR basis to build on x_grid (default 20).
#' @param rescale Logical; passed to construct_DR_basis.
#'
#' @return A list with x_scaled (grid scaled as (x_grid - mean x)/sd x), Z (the
#'   aligned basis evaluated on the grid) and reconstruction_error (diagnostic).
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

#' Bivariate (continuous spatial) Demmler-Reinsch basis.
#'
#' Builds a tensor-product DR basis for a smooth effect of two spatial
#' coordinates (x1, x2). This is the spatial counterpart of construct_DR_basis
#' for a 2D P-spline smooth: the tensor-product design B = tensor(B1, B2) with
#' sum penalty K is constrained to be orthogonal to the constant only (the code
#' deliberately does NOT enforce orthogonality to a linear term, because no
#' linear spatial effect is included in the predictor). The Gram/penalty pair is
#' then diagonalised in a generalized eigen-problem and the columns normalised
#' to unit Frobenius norm. The returned basis/transformation are used to build
#' the spatial columns Z_spat of the design matrix in the application scripts.
#'
#' @param x1 Numeric vector of first spatial coordinate (length n).
#' @param x2 Numeric vector of second spatial coordinate (length n).
#' @param k1 Number of basis functions in the first marginal (TeX: ~k1 splines).
#' @param k2 Number of basis functions in the second marginal.
#' @param normalize_cols Logical; accepted for interface-consistency but currently
#'   NOT used inside the function (the columns are always normalised to unit
#'   Frobenius norm).
#'
#' @return A list with components
#'   \item{Z}{n x (k1*k2 - 1) DR basis matrix for the spatial effect.}
#'   \item{Trafo}{Transformation such that Z = B %*% Trafo.}
#'   \item{B}{Tensor-product B-spline design before the DR reparametrisation.}
#'   \item{B1, B2}{Marginal B-spline design matrices for x1 and x2.}
#'   \item{knots1, knots2}{Knots of the two marginal smooths.}
#'
#' @seealso construct_DR_basis
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
#' Single IWLS-type Metropolis-Hastings update for a block of regression
#' coefficients.
#'
#' Implements one MH step updating the coefficients in `idx` (intercept, a
#' binary/linear term, a whole nonlinear block, or a slice of the spatial
#' block) with the IWLS proposal described in the notes: conditionally on the
#' other coefficients and on the random effects `u_curr`, the proposal is
#'   beta[idx] ~ N( beta[idx] + step * Sigma * grad , Sigma ),  Sigma = -1/Hess,
#' where grad/Hess are the gradient and (diagonal) Hessian of the
#' log-full-conditional returned by eval_logfullcond_beta. `step` multiplies the
#' gradient (NOT the covariance), matching the notes' implementation comment.
#' The proposal is evaluated on both sides, so the needed densities
#' q(beta_prop|cur) and q(beta_curr|prop) are computed at their own gradients.
#'
#' @param idx Integer vector indexing the block of beta to update.
#' @param beta_curr Current full coefficient vector.
#' @param u_curr Current random-effects vector, added to eta_i (used to evaluate
#'   the full conditional at the correct linear predictor).
#' @param y Numeric response vector.
#' @param X Numeric design matrix Phi = (1, x, Z).
#' @param X2 Element-wise square of X (diagonal-Hessian helper).
#' @param freq Numeric vector of observation weights.
#' @param link Link function (HL: mu = link(eta)).
#' @param link_prime First derivative of the link.
#' @param link_second Second derivative of the link.
#' @param k_post, c_post Flatness parameter k_post and differentiability parameter
#'   c_post of the posterior (target) distribution (TeX: k_ni, c_n for the likelihood).
#' @param k_deriv, c_deriv Flatness parameter k_deriv and differentiability
#'   parameter c_deriv of the smoother loss used ONLY to build the proposal
#'   gradient/Hessian (larger c, smaller k for wider support).
#' @param w Generalized-Bayes weight (UILC-calibrated).
#' @param tau Vector of prior standard deviations (sd, so tau^2 is the variance).
#' @param tau_idx Maps each coefficient beta[k], k>1, to its tau component.
#' @param hess_eps Small positive constant used by make_negdef.
#' @param step Scale multiplying the gradient in the proposal mean (default 1).
#'
#' @return A list with `beta` (updated full coefficient vector, either the
#'   proposal or the previous one), `accepted` (0/1), `log_alpha` (log MH ratio)
#'   and `logpost` (log full conditional at the returned value).
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


#' Single IWLS-type Metropolis-Hastings update for one random-effect coefficient.
#'
#' Random-effect counterpart of iwls_beta_update: updates a single random
#' intercept u (scalar) associated with one group, using the same IWLS proposal
#' as for the fixed effects. The linear predictor is fe + u, where `fe` contains
#' the fixed-effects part for the observations of that group, and the log
#' posterior includes the N(0, sigma2_u) prior on u. The proposal mean is
#'   u_prop ~ N( u + step * Sigma * grad , Sigma ), Sigma = -1/Hess,
#' with grad/Hess from eval_logfullcond_re, and the MH ratio uses the proposal
#' densities evaluated at both u_curr and u_prop.
#'
#' @param u_curr Current value of the random effect for the group.
#' @param y Numeric response values for the observations in this group.
#' @param fe Fixed-effects linear predictor (X %*% beta) for this group's rows.
#' @param freq Numeric observation weights for this group's rows.
#' @param link Link function.
#' @param link_prime First derivative of the link.
#' @param link_second Second derivative of the link.
#' @param k_post, c_post Flatness parameter k_post and differentiability parameter
#'   c_post of the posterior loss.
#' @param k_deriv, c_deriv Flatness parameter k_deriv and differentiability
#'   parameter c_deriv of the smoother loss used for the proposal
#'   (gradient/Hessian).
#' @param w Generalized-Bayes weight.
#' @param sigma2_u Variance of the random effect prior N(0, sigma2_u).
#' @param hess_eps Small positive constant used by make_negdef.
#' @param step Scale multiplying the gradient in the proposal mean (default 1).
#'
#' @return A list with `u` (either the proposal or the previous value),
#'   `accepted` (0/1), `log_alpha` and `logpost`.
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
#' Log-full-conditional, gradient and (diagonal) Hessian for a block of beta.
#'
#' Evaluates the full conditional of the coefficients in `idx` (one MH block).
#' The log-posterior is
#'   logpost(beta[idx]) = -w * sum_i freq_i loss_asymm(r_i; k_post, c_post)
#'                        - 1/2 * sum(beta[idx]^2 / tau_tmp^2)   [+ const],
#' with residuals r = y - link(X %*% beta + u). The gradient/Hessian combine
#' the loss contribution (score_pop / H_pop, evaluated with the *smoother*
#' proposal parameters k_deriv, c_deriv) and the Gaussian prior contribution
#' -Omega*beta and -Omega, where Omega = 1/tau_tmp^2. The intercept (idx[1]==1)
#' has a flat prior (no prior contribution).
#'
#' @param idx Integer vector indexing the block of beta to evaluate.
#' @param beta Full coefficient vector.
#' @param X Numeric design matrix Phi = (1, x, Z).
#' @param y Numeric response vector.
#' @param X2 Element-wise square of X (defaults to own X, used for Hessian).
#' @param freq Numeric observation weights.
#' @param u Random-effects vector; eta = X %*% beta + u.
#' @param link Link function.
#' @param link_prime First derivative of the link.
#' @param link_second Second derivative of the link.
#' @param k_post, c_post Flatness parameter k_post and differentiability parameter
#'   c_post of the posterior (target) distribution.
#' @param k_deriv, c_deriv Flatness parameter k_deriv and differentiability
#'   parameter c_deriv of the smoother loss used for gradient/Hessian.
#' @param w Generalized-Bayes weight.
#' @param tau Vector of prior standard deviations (tau^2 = prior variance).
#' @param tau_idx Maps each coefficient beta[k], k>1, to its tau component.
#'
#' @return A list with `logpost` (scalar), `grad` (vector of length length(idx))
#'   and `H` (diagonal Hessian, vector of length length(idx)).
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


#' Log-full-conditional, gradient and Hessian for one random effect u.
#'
#' Random-effect counterpart of eval_logfullcond_beta: evaluates the full
#' conditional of a single random intercept `u` with fixed-effects predictor
#' `fe` and prior N(0, sigma2_u). The log-posterior is
#'   logpost(u) = -w * sum_i freq_i loss_asymm(y_i - link(fe_i + u); k_post, c_post)
#'                - 1/2 * u^2 / sigma2_u,
#' and the gradient/Hessian combine the loss contribution (score_pop / H_pop,
#' evaluated with the smoother proposal params k_deriv, c_deriv) with the prior
#' contribution -u/sigma2_u and -1/sigma2_u.
#'
#' @param y Numeric response values for this group.
#' @param fe Fixed-effects linear predictor (X %*% beta) for this group's rows.
#' @param u Current value of the random effect.
#' @param freq Numeric observation weights.
#' @param link Link function.
#' @param link_prime First derivative of the link.
#' @param link_second Second derivative of the link.
#' @param k_post, c_post Flatness parameter k_post and differentiability parameter
#'   c_post of the posterior (target) distribution.
#' @param k_deriv, c_deriv Flatness parameter k_deriv and differentiability
#'   parameter c_deriv of the smoother loss used for gradient/Hessian.
#' @param w Generalized-Bayes weight.
#' @param sigma2_u Variance of the random-effect prior N(0, sigma2_u).
#'
#' @return A list with `logpost` (scalar), `grad` (scalar) and `H` (scalar).
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
#' Enforce (diagonal) negative definiteness of the Hessian.
#'
#' The IWLS proposal uses covariance -1/H, so the diagonal Hessian H must be
#' (element-wise) negative. Any entry that is >= 0 (e.g. where the local
#' log-posterior is non-concave) is overwritten by the small negative constant
#' -eps, guaranteeing a valid proposal covariance. This implements the simple
#' correction described in the notes.
#'
#' @param H Numeric vector/matrix, the diagonal Hessian.
#' @param eps Positive constant used to replace non-negative entries.
#'
#' @return H with all entries strictly negative (non-negative ones set to -eps).
make_negdef <- function(H, eps = 1e-3) {
  H[H>=0] = -eps
  return(H)
}


## Slice sampler for sigma_u (sd rand. eff.) DO NOT USE ----
#' DEPRECATED. Slice sampler for the standard deviation sigma_u of the random
#' effects. DO NOT USE.
#'
#' Older sampler for sigma_u (on the log scale) under a half-t prior, using a
#' slice-with-step-out scheme. It has been superseded by the Scale-Dependent
#' Prior (SDP) update implemented in update_re / iwls_SDP_update, which is what
#' the application uses. Kept only for reference.
#'
#' @param sigma_u_curr Current value of sigma_u (> 0).
#' @param u Vector of current random effects (one per group).
#' @param A Scale of the half-t prior.
#' @param nu Degrees of freedom of the half-t prior (default 1).
#' @param slice_w Slice width for the stepping-out step.
#' @param m Maximum number of stepping-out iterations.
#'
#' @return A single draw for sigma_u (> 0).
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


#' Log-full-conditional, gradient and Hessian for the SDP variance parameter.
#'
#' Evaluates the full conditional of log_tau2 = log(tau^2) under the
#' Scale-Dependent Prior (SDP, Klein & Kneib 2016, Bayesian Analysis). With the
#' random coefficients beta and quadratic-form penalty beta' K beta,
#'   logpost = u/2*(1 - rank_K) - beta'K beta/(2 e^u) - sqrt(e^u/theta),
#' where u = log_tau2 and theta is the prior "shape" scale. Returns the log
#' full conditional together with its first (`fd`) and second (`sd`) derivatives
#' in u, used directly by iwls_SDP_update as its target.
#'
#' @param log_tau2 Scalar, current log-variance log(tau^2).
#' @param K Penalty matrix (p x p), e.g. identity or a roughness penalty.
#' @param rank_K Rank of K (number of prior "effective" parameters).
#' @param beta Current random/beta coefficient vector plugged into beta'K beta.
#' @param theta Positive scale of the SDP prior.
#'
#' @return A list with `logpost` (scalar), `grad` (scalar) and `hess` (scalar).
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

#' IWLS-type Metropolis-Hastings update for the SDP variance parameter.
#'
#' Updates log_tau2 = log(tau^2) with one MH step using the IWLS proposal built
#' from the full conditional returned by eval_logfullcond_SDP. The proposal is
#'   log_tau2_prop ~ N( log_tau2 + step * Sigma * grad , Sigma ),  Sigma = -1/hess,
#' the Hessian being the (already negative) second derivative `sd` of the full
#' conditional, so no make_negdef call is needed. Accepts with the usual MH
#' ratio computed at both current and proposed states.
#'
#' @param log_tau2 Scalar, current log-variance log(tau^2).
#' @param K Penalty matrix (p x p).
#' @param rank_K Rank of K.
#' @param beta Current coefficient vector plugged into the full conditional.
#' @param theta Positive scale of the SDP prior.
#' @param step Scale multiplying the gradient in the proposal mean (default 1).
#'
#' @return A list with `log_tau2` (either the proposal or the previous value),
#'   `accepted` (0/1), `log_alpha` (log MH ratio) and `logpost`.
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
#' Update all random effects and their common variance sigma2_u.
#'
#' Sweeps over every random intercept u_i (one per group, indexed by `id`),
#' updating each with iwls_re_update and recording per-group acceptance rates.
#' After all u_i are refreshed, the shared variance sigma2_u is updated on the
#' log scale with iwls_SDP_update (SDP prior for the random-effects variance,
#' using penalty K = I with full rank), and the corresponding acceptance is
#' recorded. Returns the new u and sigma2_u together with acceptance vectors.
#'
#' @param u Vector of current random effects, one per group.
#' @param sigma2_u Current value of the random-effects variance.
#' @param id Group index vector of length n mapping each observation to a group.
#' @param theta Positive scale of the SDP prior for sigma2_u.
#' @param y Numeric response vector.
#' @param X Numeric design matrix of fixed effects.
#' @param freq Numeric vector of observation weights.
#' @param beta Current fixed-effect coefficient vector.
#' @param link Link function.
#' @param link_prime First derivative of the link.
#' @param link_second Second derivative of the link.
#' @param k_post, c_post Flatness parameter k_post and differentiability parameter
#'   c_post of the posterior (target) distribution.
#' @param k_deriv, c_deriv Flatness parameter k_deriv and differentiability
#'   parameter c_deriv of the smoother loss used for gradient/Hessian.
#' @param w Generalized-Bayes weight.
#' @param hess_eps Small positive constant used by make_negdef (default 1e-3).
#' @param step Scale multiplying the gradient in the proposal mean (default 1).
#'
#' @return A list with `u` (updated random effects), `sigma2_u` (updated
#'   variance), `acc_u` (per-group acceptance 0/1 vector) and `acc_sigma2_u`.
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
#' Logit link function: mu = 1/(1 + exp(-eta)).
#'
#' @param eta Numeric vector of linear predictors.
#' @return The logit-link value in (0, 1).
logit_link = function(eta) {
  return(1 / (1 + exp(-eta)))
}

#' First derivative of the logit link: dmu/deta = mu*(1 - mu).
#'
#' @param eta Numeric vector of linear predictors.
#' @return The derivative of the logit link w.r.t. eta.
logit_prime = function(eta) {
  ll = logit_link(eta)
  return(ll * (1 - ll))
}

#' Second derivative of the logit link: d^2mu/deta^2 = mu*(1-mu)*(1-2mu).
#'
#' @param eta Numeric vector of linear predictors.
#' @return The second derivative of the logit link w.r.t. eta.
logit_second = function(eta) {
  ll = logit_link(eta)
  return(ll * (1 - ll) * (1 - 2 * ll))
}

## Plots ----
### Spatial splines ----
#' Plot individual basis functions of a spatial spline basis (DR or B-spline).
#'
#' Renders the chosen columns (functions) of a spatial basis as filled contours
#' over the plotting grid, with the observed locations overlaid as red points.
#' Used to visualise DR or B-spline basis functions over the x1-x2 plane. Works
#' for either "DR" (diverging RdBu palette, signed values) or "B-Splines"
#' (non-negative warm palette) bases, optionally restricting the plots to a
#' polygon `boundaries`. Any function-plots can also be normalised to [0, 1].
#'
#' @param coord n x 2 numeric matrix of observed coordinates (x1, x2).
#' @param grid Number of points or a grid of coordinates across which the basis
#'   is evaluated.
#' @param bases The basis matrix whose columns are plotted.
#' @param idx_to_plot Integer vector of column indices of `bases` to plot.
#' @param normalize Logical; if TRUE each plotted column is rescaled to [0, 1].
#' @param type Character, "DR" or "B-Splines", selecting the palette/legend style.
#' @param title Optional plot title (defaults to `type`).
#' @param boundaries Optional 2-column matrix of polygon vertices to overlay.
#'
#' @return A ggplot object (invisibly the plot).
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


#' Plot the fitted spatial effect surface via image2D.
#'
#' Reconstructs the bivariate spatial effect f(x1, x2) from the DR spatial
#' basis coordinates `Trafo`, the fitted coefficients `beta_est` and the grid
#' of marginal B-spline knots, then draws it as a filled contour map with the
#' observed locations overlaid. The reconstruction first evaluates the two
#' marginal P-spline designs on an `ngrid` x `ngrid` grid, forms their tensor
#' product, maps through `Trafo` and multiplies by `beta_est`. Optionally the
#' surface is masked to a (convex or concave) hull of the observed points.
#'
#' @param Xspat n x 2 numeric matrix of observed spatial coordinates.
#' @param beta_est Coefficient vector for the spatial effect (length = #cols of Trafo).
#' @param Trafo Transformation matrix mapping B-spline design to DR coordinates.
#' @param xlim,ylim Numeric vectors of length 2 giving the plotting limits
#'   (defaults c(-0.2, 1.2)).
#' @param convex_hull Logical; if TRUE mask the surface to the convex hull of Xspat.
#' @param concave_hull Logical; if TRUE mask the surface to a concave hull of Xspat.
#' @param ncol Number of colours in the diverging palette (default 30).
#' @param ngrid Points per axis of the square evaluation grid (default 100).
#' @param k1,k2 Numbers of basis functions in the two marginal B-splines
#'   (defaults 5; must match the pp-school mind fit used for the design).
#' @param knots1,knots2 Knot vectors of the two marginal smooths.
#' @param name Optional title text (reserved; currently not plotted).
#' @param cex Point expansion factor for the overlaid observations.
#'
#' @return NULL (draws the image plot as a side effect).
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

#' DEPRECATED. Plot the estimated overall effects of all continuous covariates.
#'
#' For each continuous covariate, reconstructs the fitted smooth effect on a
#' grid, then plots the posterior mean together with 2.5%/97.5% pointwise
#' credible bands. Reconstruction uses the stored DR transformation `Trafo`,
#' the appropriate design indices/blocks and the MCMC / posterior samples in
#' `fit$samples$beta`. Reads the design information directly from the `fit`
#' object; primarily used by the earlier simulation scripts and not by the
#' pipeline. Kept for reference.
#'
#' @param fit A fitted model object containing at least `fit$design` (with
#'   `X`, `indices`, `dim`, `Xcont`, `Xdummy`, `Trafo`) and `fit$samples$beta`
#'   (matrix of posterior draws, rows = iterations).
#' @param ylim Numeric vector of length 2 with the y-axis limits (default c(-1,1)).
#'
#' @return Invisible NULL (draws one plot per continuous covariate).
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
