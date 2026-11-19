# === === === === === === === === === === === === === === === === === === ===
# Single-replicate worker: generate one dataset and fit M0-M5.
# M5 samples M0's exact posterior with Stan/NUTS. run_cell returns the fits,
# designs and settings; evaluation and study-wide collection are separate.
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

source("src/sim-study/01-dgp.R", local = TRUE)    # generate_data, f_true, car_setup
source("src/01a-utils.R", local = TRUE)           # loss_asymm, construct_DR_basis, ...
source("src/01b-mcmc.R", local = TRUE)            # MCMC_IWLS
source("src/sim-study/02-worker-stan.R", local = TRUE) # M5 backend definitions; compile lazily

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

# Validate settings before constructing indices or starting any model fit.
mcmc_schedule <- function(total_iter, burnin, thin) {
  settings <- list(total_iter = total_iter, burnin = burnin, thin = thin)
  for (nm in names(settings)) {
    x <- settings[[nm]]
    if (!(is.numeric(x) && length(x) == 1L && is.finite(x) &&
          x > 0 && x == floor(x) && x <= .Machine$integer.max)) {
      stop(nm, " must be a positive integer within R's integer range.")
    }
  }
  if (total_iter <= burnin) stop("total_iter must exceed burnin.")
  if ((total_iter - burnin) %% thin != 0) {
    stop("The post-burn-in iteration count must be divisible by thin.")
  }
  list(total_iter = as.integer(total_iter), burnin = as.integer(burnin),
       thin = as.integer(thin),
       keep = seq.int(burnin + 1L, total_iter, by = thin),
       ndraws = as.integer((total_iter - burnin) / thin))
}

assert_mcmc_draws <- function(draws, expected, model) {
  if (!is.matrix(draws)) stop(model, " must return a matrix of draws.")
  if (nrow(draws) != expected) {
    stop(model, " returned ", nrow(draws), " draws; expected ", expected, ".")
  }
  invisible(TRUE)
}

# MCMC settings (every Bayesian model, every replicate).
# Select the profile before evaluating this configuration block.
MCMC_PROFILE <- "development" # "production" or "development"
MCMC_PROFILES <- list(production = c(total_iter = 11000L, burnin = 1000L, thin = 10L),
                      development = c(total_iter = 2020L, burnin = 10L, thin = 1L))
stopifnot(MCMC_PROFILE %in% names(MCMC_PROFILES))
MCMC_TOTAL_ITER <- MCMC_PROFILES[[MCMC_PROFILE]][["total_iter"]] # Includes burn-in
MCMC_BURNIN <- MCMC_PROFILES[[MCMC_PROFILE]][["burnin"]]         # Number discarded
MCMC_THIN <- MCMC_PROFILES[[MCMC_PROFILE]][["thin"]]
MCMC_SETTINGS <- mcmc_schedule(MCMC_TOTAL_ITER, MCMC_BURNIN, MCMC_THIN)
MCMC_KEEP <- MCMC_SETTINGS$keep
MCMC_NDRAWS <- MCMC_SETTINGS$ndraws

# M5 uses M0's exact posterior, sampled by Stan. One chain preserves the
# study's default retained-draw count; more chains can be selected explicitly.
M5_CHAINS <- 1L
M5_CORES <- 1L
M5_CONTROL <- list(adapt_delta = .95, max_treedepth = 12L)
M5_INIT_METHOD <- "dispersed"

# Spline dimension per continuous covariate for M0/M3/M5.
# All three use P-splines with k = D_SPLINE + 2 and a second-order difference
# penalty m = c(2, 2). D_SPLINE counts the nonlinear DR basis functions in M0/M5;
# M3 uses the corresponding original spline basis with centering absorbed.
# M1 intentionally uses bamlss/mgcv s() defaults, independently of D_SPLINE.
D_SPLINE <- 8

# IWLS proposal smoothing. SF = 1 disables extra smoothing; k_deriv = k/SF.
# (For wider exploration one may raise SF, but SF = 1 is the default.)
SF <- 1

# Optimization policy for K-fold calibration (M0/M2) and frequentist M3.
# Warnings can be disabled for large simulation runs; diagnostics are always
# retained. NTRY counts total attempts, separately per CV fold / M3 method.
# Stop early on a usable zero-code result. Retries start from the best usable
# coefficients so far (otherwise the original start), multiplying the preceding
# maxit by OPTIM_RETRY_FACTOR. OPTIM_RETRY=FALSE overrides both limits to one.
# Finite nonconverged results remain eligible for lowest-objective selection.
OPTIM_WARN <- TRUE
OPTIM_RETRY <- TRUE
OPTIM_NTRY_CV <- 5L
OPTIM_NTRY_M3 <- 5L
OPTIM_RETRY_FACTOR <- 5
CV_OPTIM_CONTROL <- list(maxit = 100L, reltol = sqrt(.Machine$double.eps)) # BFGS defaults
# Analytical training-loss gradient for M0/M2 calibration. FALSE restores
# numerical differentiation; the choice can change the calibrated w and draws.
CV_ANALYTIC_GRADIENT <- TRUE
M3_OPTIM_MAXIT <- 5000L
M3_OPTIM_RELTOL <- 1e-10
# L-BFGS-B previously ignored reltol and used factr=1e7. This explicitly sets
# factr*machine_epsilon=1e-10; it is not an exact equivalence of stopping rules
# across methods and may change M3 results. Set 1e7 to restore its former value.
M3_OPTIM_FACTR <- M3_OPTIM_RELTOL / .Machine$double.eps

# Priors on block variances: weakly-informative inverse-Gamma. With
# a_tau = b_tau = 1e-3 the prior is vague (mode near 0, heavy tail); a common
# alternative is shape = rate = 1 (standard default).
# In src/09-find-hyper-IG.R I also implemented a selection obtained by putting
# a contstraint on the right tail.
a_tau_default <- 1
b_tau_default <- 0.005 # IG scale; enters the Gamma full conditional as a rate.

# Fixed Gaussian intercept prior for the NIL models M0 and M2.
# The mean and standard deviation refer to the intercept of the centered
# design; the standard deviation is not a variance or a sampled tau component.
NIL_INTERCEPT_MEAN <- 0
NIL_INTERCEPT_SD <- 100

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# 2. SHARED FITTING HELPERS ####
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

# Shared recipe for the median start and NIL loss parameters. Call it afresh
# for each model or bootstrap resample, using that fit's own X and y.
# Keep coefficient names here; callers retain their existing naming convention.
nils_setup <- function(X, y) {
  med <- quantreg::rq.fit(X, y, tau = 0.5)
  s12 <- est_s1_s2(as.vector(y - X %*% med$coefficients))
  kn <- k_n(nrow(X))
  list(start = med$coefficients, s = s12,
       k = c(kn / s12["s1"], kn / s12["s2"]), c = c_n(nrow(X)))
}

# One prior recipe for both samplers. The last IG entries are unused when
# there is a GMRF: its variance has the SDP prior specified by theta_GMRF.
nils_prior <- function(dd, splines = TRUE) {
  ntau <- dd$q + dd$p + as.integer(splines) * dd$p + as.integer(dd$s > 0)
  a_tau <- rep(a_tau_default, ntau); b_tau <- rep(b_tau_default, ntau)
  if (splines) {
    nonlinear <- dd$q + dd$p + seq_len(dd$p)
    a_tau[nonlinear] <- 2
    b_tau[nonlinear] <- .10^2 * qgamma(.05, shape = 2, rate = 1)
  }
  list(a_tau = a_tau, b_tau = b_tau,
       intercept = list(mean = NIL_INTERCEPT_MEAN, sd = NIL_INTERCEPT_SD))
}

# M3 and its bootstrap use the same method-specific settings. L-BFGS-B uses
# factr, while SANN has no reltol stopping rule. Values come from the existing
# configuration; seeding and the choice of optimization methods stay in callers.
m3_optim_control <- function(method) {
  switch(method,
    "L-BFGS-B" = list(maxit = M3_OPTIM_MAXIT, factr = M3_OPTIM_FACTR),
    "SANN" = list(maxit = M3_OPTIM_MAXIT),
    list(maxit = M3_OPTIM_MAXIT, reltol = M3_OPTIM_RELTOL))
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
  
  Omega = diag(lambda)
  
  L = Gamma_pen %*% Omega^{1/2}
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
# splines=TRUE -> DR P-splines for x1,x2 (M0); FALSE -> linear (M2).
# X3 always adds the GMRF; X1/X2 never do.
make_design_M0 <- function(dat, predictor, splines) {
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
    X_GMRF <- gmrf_design$V_tilde[ridx, , drop = FALSE]
    s <- gmrf_design$k
    K_GMRF <- diag(1, s)
    
    # gt_GMRF <- sdPrior::get_theta(Z = diag(nrow(gmrf_design$K)),
    #                               Kinv = gmrf_design$K_minus)
    # theta_GMRF <- gt_GMRF$root
    cat("Currently theta_GMRF is fixed to 1 for debugging purposes!!!! \n")
    theta_GMRF <- 1
  }
  X <- cbind(1, X_cont, Z_cont, X_GMRF)
  list(X = X, X_cont = X_cont, Z_cont = Z_cont, X_GMRF = X_GMRF,
       K_GMRF = K_GMRF, theta_GMRF = theta_GMRF,
       dr1 = dr1, dr2 = dr2,               # for test_design
       d = d, p = ncol(X_cont), q = 0L, s = s)
}

make_design_M3 <- function(dat, predictor) {
  predictor <- match.arg(predictor, c("X1", "X2", "X3"))
  
  make_bs <- function(x) {
    mgcv::smoothCon(mgcv::s(x, bs = "ps", k = D_SPLINE + 2L, m = c(2, 2)),
                    data = data.frame(x = x),
                    absorb.cons = TRUE, scale.penalty = FALSE)[[1L]]
  }
  
  bs1 <- make_bs(dat$x1)
  bs2 <- make_bs(dat$x2)
  
  B1 <- bs1$X
  B2 <- bs2$X
  
  R <- matrix(numeric(0), nrow = nrow(dat), ncol = 0L)
  region_levels <- region_contrasts <- NULL
  
  if (predictor == "X3") {
    reg <- if (is.factor(dat$region)) {dat$region} else {factor(dat$region)}
    
    if (anyNA(reg) || any(table(reg) == 0L)) {stop("M3 requires observations in every declared region.")}
    if (nlevels(reg) < 2L) {stop("X3 requires at least two observed regions.")}
    
    region_levels <- levels(reg)
    
    # Sum-to-zero region effects: s regions, s - 1 coefficients.
    region_contrasts <- contr.sum(region_levels)
    colnames(region_contrasts) <- paste0("region.c", seq_len(ncol(region_contrasts)))
    
    R <- region_contrasts[as.integer(reg), , drop = FALSE]
  }
  
  X <- cbind(Intercept = 1, B1, B2, R)
  
  if (qr(X)$rank < ncol(X)) {stop("M3 design matrix is rank deficient.")}
  
  idx <- list(intercept = 1L, x1 = 1L + seq_len(ncol(B1)),
              x2 = 1L + ncol(B1) + seq_len(ncol(B2)),
              region = 1L + ncol(B1) + ncol(B2) + seq_len(ncol(R)))
  
  list(X = X, bs1 = bs1, bs2 = bs2, 
       region_levels = region_levels, region_contrasts = region_contrasts, s = ncol(R),
       idx = idx, predictor = predictor)
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
      if (term == "s(region)") {
        # bamlss stores regional columns in its own order (e.g. x10 before x2).
        # Match coefficient names to those actual columns, not numeric suffixes.
        basis_names <- colnames(s$X)
        if (is.null(basis_names) || anyNA(basis_names) || anyDuplicated(basis_names))
          stop("M1 regional design requires unique column names for coefficient matching.")
        paste0(model, ".s.", term, ".", basis_names)
      } else {
        paste0(model, ".s.", term, ".b", seq_len(ncol(s$X)))
      }
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
# Uses NIL_INTERCEPT_MEAN/SD for the fixed Gaussian global-intercept prior.
# Returns calibration$w, calibration$numerator, calibration$mean_cv and per-fold
# optimizer diagnostics (calibration$folds), including any retry attempts, along
# with coefficient/tau draws, fhat, acceptance summaries and loss parameters.
# tau contains STANDARD DEVIATIONS, with the same retained rows as draws.
# intercept_prior and tau_prior record the priors actually used by this fit.
# Optional design must correspond to the same dat/predictor/model settings.
# run_cell passes its stored design; standalone calls build their own as before.
fit_nils_gbi <- function(dat, predictor, model_type, seed = NULL, design = NULL) {
  mcmc <- mcmc_schedule(MCMC_TOTAL_ITER, MCMC_BURNIN, MCMC_THIN)
  dd <- if (is.null(design)) make_design_M0(dat, predictor, splines = (model_type == "M0")) else design
  y <- dat$y; X <- dd$X
  setup <- nils_setup(X, y)
  s12 <- setup$s; k <- setup$k; cn <- setup$c
  k_deriv <- k / SF; c_deriv <- cn * SF
  
  P <- 1 + dd$q + dd$p + sum(dd$d) + dd$s
  
  # Initialization: the median fit on the FULL design X returns a coefficient
  # for every column (intercept, linear, DR-spline, and GMRF).
  beta0 <- unname(setup$start)
  
  ntau <- dd$q + dd$p + (model_type == "M0")*dd$p + 1 * (dd$s > 0)
  prior <- nils_prior(dd, splines = model_type == "M0")
  a_tau <- prior$a_tau; b_tau <- prior$b_tau
  
  fit <- MCMC_IWLS(niter = mcmc$total_iter,
                   X_dummy = NULL, X_cont = dd$X_cont, Z_cont = dd$Z_cont,
                   Z_spat = NULL, X_GMRF = dd$X_GMRF,
                   id = NULL, y = y,
                   link = function(x) x,
                   link_prime = function(x) rep(1, length(x)),
                   link_second = function(x) rep(0, length(x)),
                   k = k, k_deriv = k_deriv, c = cn, c_deriv = c_deriv,
                   d = dd$d, beta0 = beta0,
                   a_tau = a_tau, b_tau = b_tau,
                   intercept_mean = NIL_INTERCEPT_MEAN, intercept_sd = NIL_INTERCEPT_SD,
                   theta_re = NULL, theta_GMRF = dd$theta_GMRF, K_GMRF = dd$K_GMRF,
                   verbose = 0, K_fold = 10, seed = seed,
                   optim_warn = OPTIM_WARN, optim_retry = OPTIM_RETRY,
                   optim_retry_factor = OPTIM_RETRY_FACTOR,
                   cv_optim_control = CV_OPTIM_CONTROL,
                   optim_ntry_cv = OPTIM_NTRY_CV,
                   cv_analytic_gradient = CV_ANALYTIC_GRADIENT)
  
  assert_mcmc_draws(fit$beta, mcmc$total_iter, paste(model_type, "unthinned chain"))
  draws <- fit$beta[mcmc$keep, , drop = FALSE]
  assert_mcmc_draws(draws, mcmc$ndraws, model_type)
  tau <- fit$tau[mcmc$keep, , drop = FALSE]
  assert_mcmc_draws(tau, mcmc$ndraws, paste(model_type, "tau"))
  tau_names <- c(paste0(colnames(dd$X_cont), ".linear"),
                 if (model_type == "M0") paste0(colnames(dd$X_cont), ".nonlinear"),
                 if (dd$s > 0) "GMRF")
  dimnames(tau) <- list(iteration = seq_len(nrow(tau)), parameter = tau_names)
  # IG priors apply to tau^2 for linear/nonlinear blocks. The GMRF block
  # instead uses an SDP prior, so omit its unused IG entries and retain theta.
  ig_idx <- seq_len(ntau - as.integer(dd$s > 0))
  tau_prior <- list(a_tau = setNames(a_tau[ig_idx], tau_names[ig_idx]),
                    b_tau = setNames(b_tau[ig_idx], tau_names[ig_idx]),
                    theta_GMRF = dd$theta_GMRF)
  
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
                                    paste0("gmrf", 1:dd$s)))
  } else if (model_type == "M2") {
    variable_names <- switch(predictor,
                             X1 = c("Intercept", "x1", "x2"),
                             X2 = c("Intercept", "x1", "x2"),
                             X3 = c("Intercept", "x1", "x2", paste0("gmrf", 1:dd$s)))
  }
  
  
  dimnames(draws) <- list(iteration = seq_len(nrow(draws)),
                          variable = variable_names)
  
  list(draws = draws, tau = tau, fhat = as.vector(X %*% colMeans(draws)),
       intercept_prior = fit$intercept_prior, tau_prior = tau_prior,
       calibration = fit$calibration, acc = fit$acc_prop, loss = list(k = k, c = cn, s = s12))
}

## M1 : Bayesian ALD median (bamlss) ----
# Intentional benchmark choice: retain the package defaults for continuous
# terms s(x1) and s(x2), including the basis family and dimension. With the
# current package versions these are thin-plate regression splines with k = 10
# (9 coefficients per effect after centering). They are not matched to the
# D_SPLINE-controlled P-spline bases of M0/M3; their priors and smoothing
# treatment are also not constrained to match. In X3, the regional MRF term
# remains explicitly specified using the shared graph penalty.
fit_bamlss_median <- function(dat, predictor, seed = NULL) {
  mcmc <- mcmc_schedule(MCMC_TOTAL_ITER, MCMC_BURNIN, MCMC_THIN)
  form <- switch(predictor,
                 X1 = y ~ s(x1) + s(x2),
                 X2 = y ~ s(x1) + s(x2),
                 X3 = y ~ s(x1) + s(x2) +
                   s(region, bs = "mrf", xt = list(penalty = gmrf_design$K)))
  
  # bamlss self-initializes via its default posterior-mode optimizer (opt_bfit)
  # before the sampler; we do not inject our DR-basis coefficients as start.
  if (!is.null(seed)) set.seed(seed)
  fit <- bamlss(form, data = dat, family = ALD_bamlss(tau = 0.5),
                n.iter = mcmc$total_iter, burnin = mcmc$burnin + 1L, thin = mcmc$thin)
  
  # Extract design matrix and draws
  out <- extract_design_draws_bamlss(fit, model = "mu")
  # Full design matrix: observations x coefficients
  X <- out$X
  # Posterior coefficients: draws x coefficients
  draws <- out$B
  assert_mcmc_draws(draws, mcmc$ndraws, "M1")
  
  variable_names <- colnames(out$B)
  variable_names <- str_replace(variable_names, 
                                pattern = fixed("mu.p.(Intercept)"), replacement = "Intercept")
  variable_names <- str_replace(variable_names, 
                                pattern = fixed("mu.s.s(x1)"), replacement = "x1")
  variable_names <- str_replace(variable_names, 
                                pattern = fixed("mu.s.s(x2)"), replacement = "x2")
  variable_names <- str_replace(variable_names, 
                                pattern = fixed(".b"), replacement = ".s")
  variable_names <- str_replace(variable_names, 
                                pattern = fixed("mu.s.s(region).x"), replacement = "gmrf")
  
  dimnames(draws) <- list(iteration = seq_len(nrow(draws)),
                          variable = variable_names)
  
  list(draws = draws,
       fhat = as.vector(predict(fit, newdata = dat, model = "mu", type = "parameter")),
       design = list(X = X, full_info = fit$x),
       fit = fit)
}

## M3 : frequentist NIL ----
# Each method retains all optimizer attempts in all[[method]]$attempts.
# The minimum objective across numerically usable results remains the selection
# rule, including finite results with a nonzero termination code. SANN's zero
# code means budget completion, not a convergence certificate.
# design optionally reuses run_cell's basis for the same data and settings.
# return_design=TRUE returns the design actually used for fitting.
# Bootstrap resamples still need their own designs, built from each resample.
fit_nils_freq <- function(dat, predictor, seed = NULL, return_design = FALSE, design = NULL) {
  dd <- if (is.null(design)) make_design_M3(dat, predictor) else design
  
  y <- dat$y; X <- dd$X
  setup <- nils_setup(X, y)
  s12 <- setup$s; k <- setup$k; cn <- setup$c
  P   <- ncol(X)
  
  loss_sum <- function(b) {
    sum(loss_asymm(y - drop(X %*% b), k = k, c = cn))
  }
  
  grad_loss <- function(b) {
    -drop(score_pop(seq_len(P), b, y, X, k = k, c = cn, freq = 1,
                    link = function(x) x, link_prime = function(x) rep(1, length(x))))
  }
  
  if (!is.null(seed)) set.seed(seed)   # determinism (esp. SANN)
  start <- setup$start
  methods <- c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN")
  out <- lapply(methods, function(m) {
    optim_with_diagnostics(start, fn = loss_sum,
      gr = if (m %in% c("BFGS", "CG", "L-BFGS-B")) grad_loss else NULL,
      method = m, control = m3_optim_control(m),
      retry = OPTIM_RETRY, retry_factor = OPTIM_RETRY_FACTOR,
      warn = OPTIM_WARN, label = paste0("M3 ", m), ntry = OPTIM_NTRY_M3)
  })
  names(out) <- methods
  if (!any(vapply(out, function(o) o$usable, logical(1)))) {
    stop(structure(list(
      message = "M3 has no finite coefficients/objective after optimization attempts for any method.",
      call = NULL, optimizations = out),
      class = c("mode_optim_error", "error", "condition")))
  }
  best <- which.min(vapply(out, function(o) o$value, 0))
  o <- out[[best]]
  
  estimates <- o$par
  variable_names <- switch(predictor,
                           X1 = c("Intercept",
                                  paste0("x1.s", 1:(D_SPLINE+1)), 
                                  paste0("x2.s", 1:(D_SPLINE+1))),
                           X2 = c("Intercept",
                                  paste0("x1.s", 1:(D_SPLINE+1)), 
                                  paste0("x2.s", 1:(D_SPLINE+1))),
                           X3 = c("Intercept",
                                  paste0("x1.s", 1:(D_SPLINE+1)), 
                                  paste0("x2.s", 1:(D_SPLINE+1)),
                                  paste0("gmrf", 1:dd$s)))
  names(estimates) <- variable_names
  
  result <- list(beta = estimates, fhat = as.vector(X %*% o$par),
       convergence = o$convergence, message = o$message,
       value = o$value, counts = o$counts, method = methods[best],
       optim_settings = list(warn = OPTIM_WARN, retry = OPTIM_RETRY,
         retry_factor = OPTIM_RETRY_FACTOR, ntry = OPTIM_NTRY_M3, maxit = M3_OPTIM_MAXIT,
         reltol = M3_OPTIM_RELTOL, factr = M3_OPTIM_FACTR),
       all = out, loss = list(k = k, c = cn, s = s12))
  if (isTRUE(return_design)) result$design <- dd
  result
}

## M4 : GUD fully parametric modal regression (with loo) ----
# Design for GUD. For X3 the spatial effect enters through the original region
# factor (categorical), since GUD has no smooth GMRF term: 
# X = model.matrix(~ x1 + x2 + region). For X1/X2 it is linear only.
build_gud_design <- function(dat, predictor) {
  predictor <- match.arg(predictor, c("X1", "X2", "X3"))
  dat_sc <- dat
  dat_sc$x1 <- scale(dat$x1)
  dat_sc$x2 <- scale(dat$x2)
  if (predictor == "X3") {
    if (is.null(dat$region)) stop("M4 requires region for X3.")
    reg <- if (is.factor(dat$region)) dat$region else factor(dat$region)
    if (anyNA(reg) || any(table(reg) == 0L)) {
      stop("M4 requires observations in every declared region.")
    }
    if (nlevels(reg) < 2L) stop("X3 requires at least two observed regions.")
    dat_sc$region <- reg
    # Keep reference-region contrasts in the fitted model. Equal-region
    # centering is performed only during result extraction.
    X <- model.matrix(~ x1 + x2 + region, data = dat_sc,
                      contrasts.arg = list(region = contr.treatment(levels(reg))))
  } else {
    X <- model.matrix(~ x1 + x2, data = dat_sc)
  }
  list(X = X)
}

modal_regression_loo <- function(X, N, y, model, seed, ...) {
  model <- match.arg(model, c("FG", "DTP", "TPSC"))
  dat <- list(N = N, P = ncol(X), y = y, X = X)
  stanfit <- switch(model,
                    TPSC = rstan::sampling(GUD:::stanmodels$TPSC, data = dat, seed = seed, ...),
                    DTP  = rstan::sampling(GUD:::stanmodels$DTP,  data = dat, seed = seed, ...),
                    FG   = rstan::sampling(GUD:::stanmodels$FG,   data = dat, seed = seed, ...))
  
  # Select coefficients explicitly in design-column order, including beta[10]
  # and later regional coefficients. Flatten iterations/chains as draw rows.
  samples <- posterior::as_draws_matrix(stanfit)
  beta_names <- paste0("beta[", seq_len(ncol(X)), "]")
  if (!all(beta_names %in% colnames(samples))) {
    stop("GUD draws do not contain all design coefficients.")
  }
  draws <- unclass(samples[, beta_names, drop = FALSE])
  # Keep nuisance parameters from the SAME rows as beta for predictive scoring.
  distribution_names <- switch(model, TPSC = c("w", "sigma", "delta"),
    FG = c("w", "sigma1", "sigma2"), DTP = c("sigma1", "sigma2", "delta1", "delta2"))
  distribution_draws <- unclass(samples[, distribution_names, drop = FALSE])
  list(draws = draws, distribution_draws = distribution_draws, stanfit = stanfit)
}

# Retains coefficient draws for all three GUD members and compact distribution
# draws for the selected member. If LOO provides no finite ELPD, return no
# winner and retain every family's draws in losing_draws for inspection.
# Reuse the stored design when supplied; standalone calls retain their builder.
fit_gud <- function(dat, predictor, seed_base, design = NULL) {
  mcmc <- mcmc_schedule(MCMC_TOTAL_ITER, MCMC_BURNIN, MCMC_THIN)
  dd <- if (is.null(design)) build_gud_design(dat, predictor) else design
  X <- dd$X
  N <- nrow(X); y <- dat$y
  mods <- c("TPSC", "FG", "DTP")
  # Separate candidate streams by 1000, rather than +1/+2/+3, which reused
  # seeds between adjacent replicates. run_cell supplies data seed + 4000.
  seeds <- seed_base + (seq_along(mods) - 1L) * 1000
  res <- lapply(seq_along(mods), function(j) {
    m  <- mods[j]
    fit <- modal_regression_loo(X, N, y, model = m, seed = seeds[j],
                                chains = 1, iter = mcmc$total_iter, 
                                warmup = mcmc$burnin, thin = mcmc$thin,  
                                # default init, seeded -> reproducible
                                refresh = 0, verbose = FALSE, open_progress = FALSE, show_messages = FALSE)  
    assert_mcmc_draws(fit$draws, mcmc$ndraws, paste("M4", m))
    elpd <- tryCatch(loo::loo(fit$stanfit)$estimates["elpd_loo", "Estimate"],
                     error = function(e) {
                       warning("loo() failed on ", m, ": ",
                               conditionMessage(e)); NA_real_ })
    list(model = m, fit = fit, elpd = elpd)
  })
  names(res) <- mods
  elpd <- setNames(vapply(res, function(z) z$elpd, numeric(1)), mods)
  
  all_draws <- lapply(res, 
                      function(z) {
                        out <- z$fit$draws
                        stopifnot(ncol(out) == ncol(X))
                        variable_names <- colnames(X)
                        variable_names[variable_names == "(Intercept)"] <- "Intercept"
                        dimnames(out) <- list(iteration = seq_len(nrow(out)),
                                              variable = variable_names)
                        out})
  valid <- is.finite(elpd)
  if (!any(valid)) return(list(draws = NULL, losing_draws = all_draws,
                               selected = "FAILURE", elpd = elpd, distribution_draws = NULL))
  best <- mods[valid][which.max(elpd[valid])]
  draws <- all_draws[[best]]
  losing_draws <- all_draws[names(all_draws) != best]
  
  list(draws = draws,
       losing_draws = losing_draws,
       selected = best,
       elpd = elpd,
       distribution_draws = res[[best]]$fit$distribution_draws)
}


# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# 7. MAIN WORKER ####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
run_cell <- function(dgp, predictor, n, replicate, seed_base = 0, include_m5 = TRUE) {
  sim <- generate_data(dgp, predictor, n, replicate, seed_base = seed_base)
  dat <- sim$dat
  s   <- sim$seed
  # Save the current values, rather than the cached schedule or later globals.
  # Loss/calibration/optimizer details already live in the individual fits.
  seeds <- s + c(data = 0, m0 = 0, m1 = 1000, m2 = 2000, m3 = 3000,
                 m4_TPSC = 4000, m4_FG = 5000, m4_DTP = 6000, bootstrap_default = 0)
  settings <- list(mcmc = mcmc_schedule(MCMC_TOTAL_ITER, MCMC_BURNIN, MCMC_THIN)[
                     c("total_iter", "burnin", "thin", "ndraws")],
                   D_SPLINE = D_SPLINE, SF = SF, seed_base = seed_base,
                   seeds = seeds, rng_kind = RNGkind())
  if (include_m5) {
    seeds <- c(seeds, m5 = s + 7000)
    settings$seeds <- seeds
    settings$m5 <- list(chains = M5_CHAINS, cores = M5_CORES,
                        control = M5_CONTROL, init_method = M5_INIT_METHOD)
  }
  
  # Build each design once and pass the same object to fitting and extraction.
  # Keep the existing construction order and each model's own parameterization.
  designs <- list(m0 = make_design_M0(dat, predictor, splines = TRUE),   # M0 
                  m1 = NULL, # M1 design is directly created in fit_bamlss_median
                  m2 = make_design_M0(dat, predictor, splines = FALSE),  # M2
                  m3 = make_design_M3(dat, predictor),     # M3
                  m4 = build_gud_design(dat, predictor))   # M4
  
  # Regular 1000-step offsets; M4 uses 4000/5000/6000 and M5 uses 7000.
  # The data seed block reserves space for all of these offsets.
  fits <- list(m0 = fit_nils_gbi(dat, predictor, "M0", seed = seeds[["m0"]], design = designs$m0),
               m1 = fit_bamlss_median(dat, predictor, seed = seeds[["m1"]]),
               m2 = fit_nils_gbi(dat, predictor, "M2", seed = seeds[["m2"]], design = designs$m2),
               m3 = fit_nils_freq(dat, predictor, seed = seeds[["m3"]], design = designs$m3),
               m4 = fit_gud(dat, predictor, seed_base = seeds[["m4_TPSC"]], design = designs$m4))
  if (include_m5) {
    designs$m5 <- designs$m0
    fits$m5 <- fit_nils_stan(dat, predictor, seed = seeds[["m5"]], design = designs$m5,
                            reference_fit = fits$m0)
  }
  
  designs$m1 <- fits$m1$design  # retain the bamlss design for M1
  fits$m1$design <- NULL          # remove the bamlss design from fits (redundant)
  
  list(scenario = sim, designs = designs, fits = fits, settings = settings)
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
    stopifnot(length(argv) == 4L)
    run_cell(argv[1], argv[2],
             n = as.integer(argv[3]),
             replicate = as.integer(argv[4]))
  }
}



# Execute only when this file itself is the Rscript entry point. Sourcing
# from another batch script must define functions without starting another run.
if (sys.nframe() == 0L) main()
