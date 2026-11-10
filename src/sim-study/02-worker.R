# === === === === === === === === === === === === === === === === === === ===
# 04-simulation/02-worker.R
# Single-replicate worker: sample data, fit M0-M4 on ONE data set, write the
# per-replicate objects required to the result script.
#
# A worker is the unit of SLURM-array computation. It reads
#   {dgp, predictor, replicate, n}
# from command-line args (through SLURM_ARRAY_TASK_ID and config file), 
# samples the data for that cell, fits all five models (run_cell), evaluate
# goodness of fit on training (effect_intervals) and on test (test_eval) and 
# saves the results.
# === === === === === === === === === === === === === === === === === === ===



# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# 0. LIBRARIES / SOURCES / PATHS ####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
library(GUD)         # rTPSC (D2), densities, stanmodels (M4)
library(R2BayesX)    # read.gra (Tanzania graph)
library(bamlss)      # M1 (ALD median)
library(mgcv)        # smooth terms / DR basis
library(posterior)   # posterior::draws (M1, M4)
library(loo)         # ELPD selection (M4)
library(quantreg)    # preliminary median fit (s1/s2, beta0)
library(sdPrior)     # get_theta -> GMRF SDP scale
library(coda)        # hdi
library(stringr)     # str_replace

source("src/sim-study/01-dgp.R")    # generate_data, f_true, car_setup
source("src/01a-utils.R")           # loss_asymm, construct_DR_basis, ...
source("src/01b-gibbs.R")           # MCMC_IWLS

OUT_DIR  <- "sim-study/results/"
GRA_PATH <- "sim-study/tanzania.gra"



# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# 1. CONFIGURATION (loss policy, MCMC, splines) ####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

# Loss parameters as functions of n (fixed across error distributions for a given n).
c_n <- function(n) 1e-3
k_n <- function(n) {
  const <- n^(5/12) * log(n)^(7/12) / n^(7/12)
  (n / const)^(1/7)
}

# MCMC settings (every Bayesian model, every replicate).
MCMC_NITER  <- 1100   # total iterations
MCMC_BURNIN <- 100    # discard
MCMC_THIN   <- 1      # keep every 10th  => (11000-1000)/10 = 1000 draws

# Spline dimension per continuous covariate (M0/M3), 2nd-order penalty.
# To make bamlss and construct_DR_basis coherent, both use a P-spline with k = D_SPLINE + 2 
# and a 2nd-order difference penalty m = c(2, 2). D_SPLINE is the DR dimension (the
# number of nonlinear basis functions), so the B-spline has D_SPLINE + 2 bases.
D_SPLINE <- 5

# IWLS proposal smoothing. SF = 1 disables extra smoothing; k_deriv = k/SF.
# (For wider exploration one may raise SF, but SF = 1 is the default.)
SF <- 1

# Priors on block variances: weakly-informative inverse-Gamma. With
# a_tau = b_tau = 1e-3 the prior is vague (mode near 0, heavy tail); a common
# alternative is shape = rate = 1 (standard default).
# In src/09-find-hyper-IG.R I also implemented a selection obtained by putting
# a contstraint on the right tail.
a_tau_default <- 1
b_tau_default <- 0.005 # CHECK IF SCALE / RATE for IG



# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# 2. LOSS-HELPER: s1, s2 PER /docs ####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# s1, s2 are the two sides of the law-of-total-variance split of the
# error, from median-regression residuals.
est_s1_s2 <- function(resid) {
  m_all <- mean(resid)
  neg <- resid[resid < 0]
  pos <- resid[resid >= 0]
  if (length(neg) == 0 || length(pos) == 0)
    stop("s1/s2 need both positive and negative median residuals")
  c(s1 = sqrt(var(neg) + (mean(neg) - m_all)^2),
    s2 = sqrt(var(pos) + (mean(pos) - m_all)^2))
}



# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# 3. MMR GMRF DESIGN (fit side, for X3) ####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Mixed-model reparameterization of the 20-region graph
gmrf_design <- local({
  K  <- unclass(R2BayesX::read.gra(GRA_PATH, sorted = TRUE))
  
  ee <- eigen(K, symmetric = TRUE)
  tol <- 1e-10 * max(ee$values)
  
  keep_pen <- ee$values > tol
  keep_null <- !keep_pen
  
  lambda <- ee$values[keep_pen]
  
  Gamma_pen  <- ee$vectors[, keep_pen,  drop = FALSE]
  Gamma_null <- ee$vectors[, keep_null, drop = FALSE]
  
  k <- sum(keep_pen)
  q <- sum(keep_null)
  
  cat("Number of positive eigenvalues (k):", k, " | Number of zero eigenvalues (q):", q, "\n")
  
  U_tilde <- Gamma_null
  
  Omega = diag(sqrt(lambda))
  
  L = Gamma_pen %*% Omega
  V_tilde = L %*% solve(crossprod(L))
  
  if (!(all(round(crossprod(V_tilde, K %*% V_tilde), 2) == diag(1, k)))){
    warning("V_tilde' K V_tilde is not diag(1, k)")
  }
  
  if (!(all(round(crossprod(U_tilde, K %*% U_tilde), 2) == diag(0, q)))){
    warning("U_tilde' K U_tilde is not diag(0, q)")
  }
  
  Omega_minus = diag(0, k, k)
  dOmega = diag(Omega)
  diag(Omega_minus)[dOmega > 0 ] = 1 / dOmega[dOmega > 0 ]
  
  K_minus = Gamma_pen %*% Omega_minus %*% t(Gamma_pen)
  
  list(V_tilde = V_tilde, U_tilde = U_tilde, K = K,
       K_minus = K_minus, k = k, q = q)
})


# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# 4. DESIGN BUILDER ####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# splines=TRUE -> DR P-splines for x1,x2 (M0, M3); FALSE -> linear (M2).
# X3 always adds the GMRF; X1/X2 never do.
make_design <- function(dat, predictor, splines) {
  x1 <- dat$x1; x2 <- dat$x2
  X_cont <- scale(cbind(x1, x2))
  Z_cont <- NULL
  d <- c(0L, 0L)
  dr1 <- dr2 <- NULL
  if (splines) {
    dr1 <- construct_DR_basis(x1, d = D_SPLINE) # keeps $Z, $Trafo, $smooth_object
    dr2 <- construct_DR_basis(x2, d = D_SPLINE) # keeps $Z, $Trafo, $smooth_object
    Z_cont <- cbind(dr1$Z, dr2$Z)
    d <- c(D_SPLINE, D_SPLINE)
  }
  X_GMRF <- K_GMRF <- theta_GMRF <- NULL
  s <- 0L
  if (predictor == "X3") {
    ridx <- as.integer(dat$region)
    X_GMRF <- gmrf_design$V_tilde[ridx, , drop = FALSE] # equivalent to Z %*% V_tilde
    s <- gmfr_design$k
    K_GMRF <- diag(1, s)
    cat("Currently theta_GMRF is fixed to 1 for debugging purposes!!!! \n")
    theta_GMRF <- 1
    # theta_GMRF <- sdPrior::get_theta(Z = model.matrix(~-1 + dat$region), Kinv = gmfr_design$K_minus)
  }
  X <- cbind(1, X_cont, Z_cont, X_GMRF)
  list(X = X, X_cont = X_cont, Z_cont = Z_cont, X_GMRF = X_GMRF,
       K_GMRF = K_GMRF, theta_GMRF = theta_GMRF,
       dr1 = dr1, dr2 = dr2,               # for test_design
       d = d, p = ncol(X_cont), q = 0L, s = s)
}



extract_design_draws_bamlss <- function(fit, model = "mu") {
  
  # Model structure
  xm <- fit$x[[model]]
  
  # Parametric design matrix
  Xp <- xm$model.matrix
  
  # All smooth terms (including MRF if present)
  sm <- xm$smooth.construct
  if (is.null(sm)) sm <- list()
  
  # Smooth design matrices
  Xs <- lapply(sm, function(s) s$X)
  
  # Full design matrix
  X <- do.call(cbind, c(list(Xp), Xs))
  
  # Names of parametric coefficients
  names_p <- paste0(model, ".p.", colnames(Xp))
  
  # Names of smooth coefficients
  names_s <- unlist(
    Map(function(term, s) {
      paste0(
        model, ".s.", term, 
        ifelse(term == "s(region)", ".x", ".b"),
        if (term == "s(region)"){seq_len(ncol(s$X)) + 1} else {seq_len(ncol(s$X))}
      )
    }, names(sm), sm),
    use.names = FALSE
  )
  
  # All coefficient names, in design-matrix order
  coef_names <- c(names_p, names_s)
  
  # All posterior draws
  S <- as.matrix(fit$samples)
  
  # Check that coefficients are available
  missing <- setdiff(coef_names, colnames(S))
  if (length(missing) > 0) {
    stop("Missing coefficients: ",
         paste(missing, collapse = ", "))
  }
  
  # Select and order coefficient draws
  B <- S[, coef_names, drop = FALSE]
  
  stopifnot(ncol(X) == ncol(B))
  
  colnames(X) <- coef_names
  
  return(list(
    X = X,
    B = B,
    terms = names(sm)
  ))
}



# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# 5. MODEL FITS (M0-M4) ####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

## M0 / M2 : Bayesian NIL+GBI ----
# model_type "M0" = splines, "M2" = linear; X3 adds the GMRF in both.
fit_nils_gbi <- function(dat, predictor, model_type, seed = NULL) {
  dd <- make_design(dat, predictor, splines = (model_type == "M0"))
  n <- nrow(dat); y <- dat$y; X <- dd$X
  kn <- k_n(n); cn <- c_n(n)
  
  med <- quantreg::rq.fit(X, y, tau = 0.5)
  
  s12 <- est_s1_s2(as.vector(X %*% med$coefficients - y))
  k   <- c(kn / s12["s1"], kn / s12["s2"])
  k_deriv <- k / SF; c_deriv <- cn * SF
  
  P <- 1 + dd$q + dd$p + sum(dd$d) + dd$s
  
  # Initialization: the median fit on the FULL design X returns a coefficient
  # for every column (intercept, linear, DR-spline, and GMRF).
  beta0 <- unname(med$coefficients)
  
  ntau <- dd$q + dd$p + (model_type == "M0")*dd$p + 1 * (dd$s > 0)
  a_tau <- rep(a_tau_default, ntau); b_tau <- rep(b_tau_default, ntau)
  
  # RE-layer placeholder (option c): single neutral group. SWITCH HERE when
  # the optional random-intercept mode is available.
  id <- rep(1L, n); theta_re <- 1
  
  fit <- MCMC_IWLS(niter = MCMC_NITER + MCMC_BURNIN,
                   X_dummy = NULL, X_cont = dd$X_cont, Z_cont = dd$Z_cont,
                   Z_spat = NULL, X_GMRF = dd$X_GMRF,
                   id = id, y = y,
                   link = function(x) x,
                   link_prime = function(x) rep(1, length(x)),
                   link_second = function(x) rep(0, length(x)),
                   k = k, k_deriv = k_deriv, c = cn, c_deriv = c_deriv,
                   d = dd$d, beta0 = beta0,
                   a_tau = a_tau, b_tau = b_tau,
                   theta_re = theta_re, theta_GMRF = dd$theta_GMRF, K_GMRF = dd$K_GMRF,
                   verbose = 0, K_fold = 10, seed = seed)
  
  idx <- seq(MCMC_BURNIN + 1, MCMC_NITER + MCMC_BURNIN, by = MCMC_THIN)
  draws <- fit$beta[idx, , drop = FALSE]
  
  if (model_type == "M0"){
    variable_names <- switch(predictor,
                             X1 = c("Intercept", "x1", "x2", 
                                    paste0("x1.s", 1:D_SPLINE), 
                                    paste0("x2.s", 1:D_SPLINE)),
                             X2 = c("Intercept", "x1", "x2", 
                                    paste0("x1.s", 1:D_SPLINE), 
                                    paste0("x2.s", 1:D_SPLINE)),
                             X3 = c("Intercept", "x1", "x2", 
                                    paste0("x1.s", 1:D_SPLINE), 
                                    paste0("x2.s", 1:D_SPLINE),
                                    paste("gmrf", 1:dd$s)))
  } else if (model_type == "M2") {
    variable_names <- switch(predictor,
                             X1 = c("Intercept", "x1", "x2"),
                             X2 = c("Intercept", "x1", "x2"),
                             X3 = c("Intercept", "x1", "x2", paste("gmrf", 1:dd$s)))
  }
  
  
  dimnames(draws) <- list(iteration = 1:MCMC_NITER, 
                          variable = variable_names)
  
  list(draws = draws, fhat = as.vector(X %*% colMeans(draws)),
       w = fit$w, acc = fit$acc_prop, loss = list(k = k, c = cn, s = s12))
}

## M1 : Bayesian ALD median (bamlss) ----
fit_bamlss_median <- function(dat, predictor, seed = NULL) {
  form <- switch(predictor,
                 X1 = y ~ s(x1) + s(x2),
                 X2 = y ~ s(x1) + s(x2),
                 X3 = y ~ s(x1) + s(x2) +
                   s(region, bs = "mrf", xt = list(penalty = gmrf_design$K)))
  
  # bamlss self-initializes via its default posterior-mode optimizer (opt_bfit)
  # before the sampler; we do not inject our DR-basis coefficients as start.
  if (!is.null(seed)) set.seed(seed)
  fit <- bamlss(form, data = dat, family = ALD_bamlss(tau = 0.5),
                n.iter = MCMC_NITER + MCMC_BURNIN, burnin = MCMC_BURNIN+1, thin = MCMC_THIN)
  
  # Extract design matrix and draws
  out <- extract_design_draws_bamlss(fit, model = "mu")
  # Full design matrix: observations x coefficients
  X <- out$X
  # Posterior coefficients: draws x coefficients
  draws <- out$B
  
  variable_names <- colnames(out$B)
  variable_names <- stringr::str_replace(variable_names, 
                                         pattern = fixed("mu.p.(Intercept)"), replacement = "Intercept")
  variable_names <- stringr::str_replace(variable_names, 
                                         pattern = fixed("mu.s.s(x1)"), replacement = "x1")
  variable_names <- stringr::str_replace(variable_names, 
                                         pattern = fixed("mu.s.s(x2)"), replacement = "x2")
  variable_names <- stringr::str_replace(variable_names, 
                                         pattern = ".b", replacement = ".s")
  variable_names <- stringr::str_replace(variable_names, 
                                         pattern = "mu.s.s(region).x", replacement = "gmrf")
  
  dimnames(draws) <- list(iteration = 1:MCMC_NITER, 
                          variable = variable_names)
  
  list(draws = draws,
       fhat = as.vector(predict(fit, newdata = dat, type = "parameter")),
       design = list(X = X, full_info = fit$x),
       fit = fit)
}

## M3 : frequentist NIL (same loss & design as M0) ----
fit_nils_freq <- function(dat, predictor, seed = NULL) {
  dd <- make_design(dat, predictor, splines = TRUE)
  n <- nrow(dat); y <- dat$y; X <- dd$X
  kn <- k_n(n); cn <- c_n(n)
  med <- quantreg::rq.fit(X, y, tau = 0.5)
  s12 <- est_s1_s2(as.vector(X %*% med$coefficients - y))
  k   <- c(kn / s12["s1"], kn / s12["s2"])
  P   <- ncol(X)
  
  loss_sum  <- function(b) sum(loss_asymm(y - X %*% b, k = k, c = cn))
  
  # score_pop returns grad of -sum(loss); the loss gradient is its negative.
  grad_loss <- function(b)
    -score_pop(seq_len(P), b, y, X, k = k, c = cn, freq = 1,
               link = function(x) x,
               link_prime = function(x) rep(1, length(x)))
  
  if (!is.null(seed)) set.seed(seed)   # determinism (esp. SANN)
  start <- med$coefficients
  methods <- c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN")
  out <- lapply(methods, function(m) {
    if (m %in% c("BFGS", "CG", "L-BFGS-B")) {
      optim(start, fn = loss_sum, gr = grad_loss, method = m,
            control = list(maxit = 5000, reltol = 1e-10))
    } else {
      optim(start, fn = loss_sum, method = m,
            control = list(maxit = 5000, reltol = 1e-10))
    }
  })
  names(out) <- methods
  best <- which.min(vapply(out, function(o) o$value, 0))
  o <- out[[best]]
  
  estimates <- o$par
  variable_names <- switch(predictor,
                           X1 = c("Intercept", "x1", "x2", 
                                  paste0("x1.s", 1:D_SPLINE), 
                                  paste0("x2.s", 1:D_SPLINE)),
                           X2 = c("Intercept", "x1", "x2", 
                                  paste0("x1.s", 1:D_SPLINE), 
                                  paste0("x2.s", 1:D_SPLINE)),
                           X3 = c("Intercept", "x1", "x2", 
                                  paste0("x1.s", 1:D_SPLINE), 
                                  paste0("x2.s", 1:D_SPLINE),
                                  paste("gmrf", 1:dd$s)))
  names(estimates) <- variable_names
  
  list(beta = estimates, fhat = as.vector(X %*% o$par),
       convergence = o$convergence, message = o$message,
       value = o$value, counts = o$counts, method = methods[best],
       all = out, loss = list(k = k, c = cn, s = s12))
}

## M4 : GUD fully parametric modal regression (with loo) ----
# Design for GUD. For X3 the spatial effect enters through the original region
# factor (categorical), since GUD has no smooth GMRF term: 
# X = model.matrix(~ x1 + x2 + region). For X1/X2 it is linear only.
build_gud_design <- function(dat) {
  dat_sc <- dat
  dat_sc$x1 <- scale(dat$x1)
  dat_sc$x2 <- scale(dat$x2)
  list(X = model.matrix(~ x1 + x2, data = dat_sc))
}

modal_regression_loo <- function(X, N, y, model, seed, ...) {
  model <- match.arg(model, c("FG", "DTP", "TPSC"))
  dat <- list(N = N, P = ncol(X), y = y, X = X)
  stanfit <- switch(model,
                    TPSC = rstan::sampling(GUD:::stanmodels$TPSC, data = dat, seed = seed, ...),
                    DTP  = rstan::sampling(GUD:::stanmodels$DTP,  data = dat, seed = seed, ...),
                    FG   = rstan::sampling(GUD:::stanmodels$FG,   data = dat, seed = seed, ...))
  
  arr = posterior::as_draws(stanfit)
  draws = matrix(arr,
                 nrow = dim(arr)[1], ncol = dim(arr)[3])[ , substr(dimnames(arr)[3]$variable, 1, 4) == "beta"]
  
  list(draws = draws, stanfit = stanfit)
}

# Retains ALL three GUD members (draws + ELPD) + the selected one.
fit_gud <- function(dat, predictor, seed_base) {
  X <- build_gud_design(dat)$X
  N <- nrow(X); y <- dat$y
  mods <- c("TPSC", "FG", "DTP")
  res <- lapply(seq_along(mods), function(j) {
    m  <- mods[j]
    fit <- modal_regression_loo(X, N, y, model = m, seed = seed_base + j,
                                chains = 1, iter = MCMC_NITER+MCMC_BURNIN, warmup = MCMC_BURNIN,
                                thin = MCMC_THIN,  # default init, seeded -> reproducible
                                refresh = 0, verbose = FALSE, open_progress = FALSE, show_messages = FALSE)  
    elpd <- tryCatch(loo::loo(fit$stanfit)$estimates["elpd_loo", "Estimate"],
                     error = function(e) {
                       warning("loo() failed on ", m, ": ",
                               conditionMessage(e)); NA_real_ })
    list(model = m, fit = fit, elpd = elpd)
  })
  names(res) <- mods
  best <- mods[which.max(vapply(res, function(z) z$elpd, 0))]
  
  all_draws <- lapply(res, 
                      function(z) {
                        out <- z$fit$draws
                        dimnames(out) <- list(iteration = 1:MCMC_NITER,
                                              variable = c("Intercept", "x1", "x2"))
                        out})
  draws <- all_draws[[best]]
  losing_draws <- all_draws[names(all_draws) != best]
  
  list(draws = draws,
       losing_draws = losing_draws,
       selected = best,
       elpd = setNames(vapply(res, function(z) z$elpd, 0), mods))
}



# # --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# # 6. B1: EFFECT-WISE INTERVALS ON TRAINING ####
# # --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# effect_columns <- function(des, effect) {
#   p <- des$p; d <- des$d; s <- des$s
#   z_start <- 2 + p; z_idx <- seq_len(sum(d))
#   spat <- if (s > 0) (2 + p + sum(d)):(2 + p + sum(d) + s - 1) else integer(0)
#   switch(effect,
#          intercept = 1,
#          x1   = c(2, z_start + seq_len(d[1])),
#          x2   = c(3, z_start + d[1] + seq_len(d[2])),
#          spatial = spat,
#          whole = seq_len(1 + p + sum(d) + s))
# }
# 
# # Per-observation HPD + quantile interval. NOTE: HDInterval::hdi() applied to a
# # *density* object can return the union of multiple intervals; we deliberately
# # use the sample/vector form (single HPD) and do not use that option.
# interval_row <- function(v, cred = 0.95) {
#   a <- 1 - cred
#   c(lo_q = as.numeric(quantile(v, a / 2)),
#     hi_q = as.numeric(quantile(v, 1 - a / 2)),
#     lo_h = as.numeric(HDInterval::hdi(v, credMass = cred)[1]),
#     hi_h = as.numeric(HDInterval::hdi(v, credMass = cred)[2]))
# }
# effect_ci_from_draws <- function(eff_draws, ff, cred = 0.95) {
#   # eff_draws: n x ndraw matrix of effect-function posterior draws
#   ci <- t(apply(eff_draws, 1, interval_row, cred = cred))
#   list(point = rowMeans(eff_draws), 
#        lo_q = ci[, "lo_q"], hi_q = ci[, "hi_q"],
#        lo_h = ci[, "lo_h"], hi_h = ci[, "hi_h"])
# }
# 
# effect_intervals <- function(s, cred = 0.95, B_m3 = 50) {
#   scen <- s$scenario; dat <- scen$dat
#   effs <- c("intercept", "x1", "x2", if (scen$predictor == "X3") "spatial",
#             "whole")
#   one_model <- function(design, draws_or_beta, is_boot = FALSE) {
#     lapply(setNames(effs, effs), function(e) {
#       cc <- effect_columns(design, e)
#       fd <- drop(design$X[, cc, drop = FALSE] %*% t(draws_or_beta[cc, , drop = FALSE]))
#       effect_ci_from_draws(fd)
#     })
#   }
#   list(
#     M0 = one_model(s$designs$m0, s$fits$m0$draws),
#     M2 = one_model(s$designs$m2, s$fits$m2$draws),
#     M4 = one_model_from_modelmatrix(s$designs$m4, s$fits$m4, scen$predictor),
#     M3 = one_model_boot(s$designs$m0, s$fits$m3, dat, scen$predictor, B_m3),
#     M1 = one_model_bamlss(s$fits$m1, dat, scen$predictor))
# }
# 
# one_model <- function(design, draws, predictor, cred = 0.95) {
#   effs <- effects_for(predictor)
#   lapply(setNames(effs, effs), function(e) {
#     cc  <- effect_columns(design, e)
#     fd  <- drop(design$X[, cc, drop = FALSE] %*% t(draws[cc, , drop = FALSE]))
#     effect_ci_from_draws(fd, cred = cred)
#   })
# }

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# 7. MAIN WORKER ####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
run_cell <- function(dgp, predictor, n, replicate, seed_base = 0) {
  sim <- generate_data(dgp, predictor, n, replicate, seed_base = seed_base)
  dat <- sim$dat
  s   <- sim$seed
  
  # Designs retained and returned (resolves the previously-unused `design`).
  designs <- list(m0 = make_design(dat, predictor, splines = TRUE),   # M0 
                  m1 = NULL, # M1 design is directly created in fit_bamlss_median
                  m2 = make_design(dat, predictor, splines = FALSE),  # M2
                  m3 = make_design(dat, predictor, splines = TRUE),   # M3 (equal to M0)
                  m4 = build_gud_design(dat))              # M4
  
  fits <- list(m0 = fit_nils_gbi(dat, predictor, "M0", seed = s),
               m1 = fit_bamlss_median(dat, predictor, seed = s + 1000),
               m2 = fit_nils_gbi(dat, predictor, "M2", seed = s + 2000),
               m3 = fit_nils_freq(dat, predictor, seed = s + 2500),
               m4 = fit_gud(dat, predictor, seed_base = s + 3000))
  
  designs$m1 <- fits$m1$design  # retain the bamlss design for M1
  fits$m1$design <- NULL          # remove the bamlss design from fits (redundant)
  
  list(scenario = sim, designs = designs, fits = fits)
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# --- 8. ARG PARSING / SLURM ARRAY MAPPING ####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Direct invocation:  Rscript 02-worker.R <dgp> <predictor> <n> <replicate>
# SLURM array:        Rscript 02-worker.R --array <ARRAY_TASK_ID> [-n 250]
grid_from_index <- function(idx) {
  stopifnot(idx >= 1, idx <= 1200)
  list(dgp       = paste0("D", ((idx - 1) %/% 300) + 1),
       predictor = paste0("X", (((idx - 1) %/% 100) %% 3) + 1),
       replicate = ((idx - 1) %% 100) + 1)
}

main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  n <- 250
  if ("--array" %in% argv) {
    idx <- as.integer(argv[match("--array", argv) + 1])
    if ("-n" %in% argv) n <- as.integer(argv[match("-n", argv) + 1])
    g <- grid_from_index(idx)
    run_cell(g$dgp, g$predictor, n, g$replicate)
  } else {
    stopifnot(length(argv) >= 3)
    n <- if (length(argv) >= 4) as.integer(argv[4]) else n
    run_cell(argv[1], argv[2], n, as.integer(argv[3]))
  }
}



# When run as a script (Rscript) this executes; when source() into an
# interactive session it only defines the functions (so I can call
# run_cell() manually).
if (!interactive()) main()