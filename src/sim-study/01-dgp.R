# === === === === === === === === === === === === === === === === === === ===
# 04-simulation/01-dgp.R
# Data-generating utilities for the Bayesian mode-regression simulation study.
# === === === === === === === === === === === === === === === === === === ===

# The study has 12 settings = 4 error distributions (D1-D4) x 3 predictor
# specifications (X1-X3), and sample sizes n in {250, 1000}.
#
# Design decisions (locked with the author):
#   * Covariates x1, x2 ~ iid U(-0.5, 0.5), i.e. already centered.
#   * Within a DGP, the covariate/error draws are identical across the three
#     predictors X1/X2/X3; only the functional form f changes. The seed is
#     therefore a function of (dgp, replicate) only, not of the predictor.
#   * For each (dgp, replicate) the full sample of size n_max = 1000 is drawn
#     once; the n = 250 setting uses the first 250 rows. This nests the smaller
#     sample in the larger one ("information coherence") and avoids sampling
#     twice.
#   * X3 adds a discrete spatial effect over the 20 Tanzania regions
#     (from tanzania.gra). The gamma vector is generated with seed 1
#     and held fixed across DGPs, predictors, and replicates.
#   * Region assignments are generated after x1, x2, eps, preserving
#     the covariate/error stream shared across predictors.
#   * The complete region vector is constructed independently of n:
#     the first 250 assignments and the remaining assignments are
#     randomly permuted separately. Calls with the same scenario seed
#     and n_max therefore have identical first-250 assignments.
#   * With n_max = 1000, every region has exactly 50 observations.
#     Among the first 250 observations, regions 1-10 have 13 each
#     and regions 11-20 have 12 each.
#
# Errors (Var(eps) = 1 in all cases):
#   D1  N(0,1)                                    (median-fair)
#   D2  TPSC-Student-t(theta=0, w=0.7, sigma=0.48, delta=3)   (Liu-Huang-Bai)
#   D3  Gamma(1.5, sqrt(1.5)) - (1.5-1)/sqrt(1.5) (mode 0, right-skewed)
#   D4  Gamma(5, sqrt(5))     - 0.5               (mode ~1.29, median ~1.59)
# === === === === === === === === === === === === === === === === === === ===

# --- REQUIRED LIBRARIES ---
library(GUD)      # rTPSC for D2
library(R2BayesX) # read.gra / read.bnd / plotmap (spatial)
# --- --- --- --- --- --- ---

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# 1. Error samplers ####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

#' Draw n i.i.d. errors for a given DGP.
error_draw <- function(dgp, n) {
  switch(dgp,
         D1 = rnorm(n, 0, 1),
         D2 = GUD::rTPSC(n = n, w = 0.7, theta = 0, sigma = 0.48, delta = 3),
         D3 = rgamma(n, shape = 1.5, rate = sqrt(1.5)) - (1.5 - 1) / sqrt(1.5),
         D4 = rgamma(n, shape = 5,    rate = sqrt(5))    - 0.5
  )
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# 2. Predictor (true regression function) ####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

#' Exact true regression function f0 for a given predictor specification.
#'   X1: -4 x1 + sqrt(32) x2
#'   X2: (2.75 x1)^3 + [2 sin(2 pi x2) + 3.25 exp(-(10 x2)^2)]
#'   X3: sqrt(0.8) * f_X2 + gamma[region]      (variance scaled by 0.8)
#' gamma  : length-20 vector of region effects (only used by X3).
#' region : length-n vector of region memberships in 1:20 (only used by X3).
f_true <- function(predictor, x1, x2, gamma = NULL, region = NULL) {
  f_X1.x1 <- function(x) -4 * x
  f_X1.x2 <- function(x) sqrt(32) * x
  f_X1 <- function(x1, x2) f_X1.x1(x1) + f_X1.x2(x2)
  f_X2.x1 <- function(x) (2.75 * x1)^3
  f_X2.x2 <- function(x) (2 * sin(2 * pi * x) + 3.25 * exp(-(10 * x)^2))
  f_X2 <- function(x1, x2) f_X2.x1(x1) + f_X2.x2(x2)
  
  switch(predictor,
         X1 = list(f1 = f_X1.x1(x1), f2 = f_X1.x2(x2), f = f_X1(x1, x2)),
         X2 = list(f1 = f_X2.x1(x1), f2 = f_X2.x2(x2), f = f_X2(x1, x2)),
         X3 = list(f1 = sqrt(0.8) * f_X2.x1(x1), f2 = sqrt(0.8) * f_X2.x2(x2), 
                   f = sqrt(0.8) * f_X2(x1, x2) + gamma[region]),
  )
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# 3. Spatial configuration (X3 only) ####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

# CAR / spatial setup (X3 only)
GRA_PATH   <- "sim-study/tanzania.gra"  # path to the graphical structure file
VAR_REGION <- 0.78            # requested marginal variance of gamma
RHO        <- 0.9             # CAR dependence; 0 <= rho < 1 (it does 
                              # not coincide with correlation, which is 
                              # usually much lower than this)


#' Build the proper-CAR sampler, computed once for a fixed graph and rho
#' read.gra returns the graph Laplacian D - W; sorted so gamma[i] is region i.
car_setup <- local({
  L <- unclass(R2BayesX::read.gra(GRA_PATH, sorted = TRUE))
  N_REGIONS <- nrow(L)
  D <- diag(diag(L))
  W <- D - L
  stopifnot(
    isSymmetric(W),
    all(diag(W) == 0),
    all(rowSums(W) > 0)
  )
  # Proper-CAR precision, then standardize to the target marginal variance
  Q0     <- D - RHO * W
  C0     <- chol2inv(chol(Q0))
  Sigma  <- VAR_REGION * cov2cor(C0)   # each region marginal var = VAR_REGION
  U      <- chol(Sigma)                # t(U) %*% U = Sigma
  draw_gamma <- function(seed) {
    set.seed(seed)
    setNames(drop(t(U) %*% rnorm(N_REGIONS)), rownames(L))
  }
  list(
    N_REGIONS = N_REGIONS,
    W         = W,
    D         = D,
    Sigma     = Sigma,
    U         = U,
    draw_gamma = draw_gamma
  )
})

N_REGIONS <- car_setup$N_REGIONS

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# 4. Main data-generator ####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

#' Generate one (predictor-specific) simulated data set.
#'
#' @param dgp       One of "D1", "D2", "D3", "D4".
#' @param predictor One of "X1", "X2", "X3".
#' @param n         Sample size; must be in {250, 1000}.
#' @param replicate Positive integer replicate index (1..R).
#' @param n_max     Size of the single full draw (default 1000) of which the
#'                  first `n` rows are kept (nesting / information coherence).
#' @param seed_base Base seed; the per-cell seed is derived from (dgp, replicate)
#'                  only, so all three predictors share the x1/x2/eps stream.
#'
#' @return A list with
#'   dat    : data.frame with columns y, x1, x2 (and region, factor, for X3),
#'   f0     : length-n true regression function,
#'   eps    : length-n errors,
#'   gamma  : (X3 only) length-20 region effects,
#'   region : (X3 only) length-n region membership in 1:20,
#'   dgp, predictor, n, replicate, seed : bookkeeping.
generate_data <- function(dgp, predictor, n, replicate,
                          n_max = 1000, seed_base = 0) {
  stopifnot(dgp %in% c("D1", "D2", "D3", "D4"),
            predictor %in% c("X1", "X2", "X3"),
            n %in% c(250, 1000),
            n <= n_max)
  
  # --- Spatial pieces (X3 only) ---
  gamma  <- NULL
  region <- NULL
  if (predictor == "X3") {
    gamma  <- car_setup$draw_gamma(1)
  }
  
  # --- Seed depends only on (dgp, replicate) ---
  dgp_idx  <- match(dgp, c("D1", "D2", "D3", "D4"))
  seed     <- seed_base + (dgp_idx - 1L) * 1e2 + replicate
  set.seed(seed)
  
  # --- Draw covariates and errors for the full sample ---
  x1  <- runif(n_max, -0.5, 0.5)
  x2  <- runif(n_max, -0.5, 0.5)
  eps <- error_draw(dgp, n_max)
  # I sample the region *after* x1, x2, eps so that the covariate/error stream is identical across X1/X2/X3.
  if (predictor == "X3") {
    # Construct the full allocation independently of the requested n.
    n_first <- 250L
    reg_val <- rep(seq_len(N_REGIONS), length.out = n_max)
    region <- integer(n_max)
    
    # First block: balanced allocation for the smaller sample.
    first_idx <- seq_len(n_first)
    region[first_idx] <- reg_val[first_idx][sample.int(length(first_idx))]
    
    # Remaining observations: complementary allocation.
    if (n_max > n_first) {
      rest_idx <- seq.int(n_first + 1L, n_max)
      region[rest_idx] <- reg_val[rest_idx][sample.int(length(rest_idx))]
    }
  }
  
  # --- true function ---
  f0 <- f_true(predictor, x1, x2, gamma, region)
  
  # --- response (full sample) ---
  y  <- f0$f + eps
  
  # --- Keep the first `n` rows (nesting: n=250 subset of n=1000) ---
  idx <- seq_len(n)
  dat <- data.frame(y = y[idx], x1 = x1[idx], x2 = x2[idx])
  if (predictor == "X3") dat$region <- factor(region[idx], levels = 1:20)
  
  list(
    dat       = dat,
    f0        = lapply(f0, function(x) x[idx]),
    eps       = eps[idx],
    gamma     = gamma,                    # NULL unless X3
    region    = region[idx],              # NULL unless X3
    dgp       = dgp,
    predictor = predictor,
    n         = n,
    replicate = replicate,
    seed      = seed
  )
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# 5. Optional sanity checks (run once, not part of the pipeline) ####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
if (FALSE) {
  # Empirical unit-variance / mode checks for the four error distributions.
  set.seed(11)
  chk <- function(dgp) {
    e <- error_draw(dgp, 2e6)
    c(var = var(e), mode_hat = density(e, n = 2^12)$x[which.max(
      density(e, n = 2^12)$y)])
  }
  sapply(c("D1", "D2", "D3", "D4"), chk)   # var ~ 1 in all cases
  
  # Same stream across predictors within (d, r): x1/x2/eps identical.
  g1 <- generate_data("D2", "X1", 250, 1)
  g2 <- generate_data("D2", "X2", 250, 1)
  stopifnot(identical(g1$dat$x1, g2$dat$x1),
            identical(g1$dat$x2, g2$dat$x2),
            identical(g1$eps,     g2$eps))
  
  # Nesting: n=250 is the first 250 rows of n=1000.
  g_full <- generate_data("D3", "X2", 1000, 2)
  g_sub  <- generate_data("D3", "X2", 250, 2)
  stopifnot(identical(unname(g_sub$dat$x1), unname(g_full$dat$x1[1:250])),
            identical(g_sub$eps, g_full$eps[1:250]))
  
  # Region nesting for X3: n=250 region is the first 250 of n=1000.
  rg_full <- generate_data("D1", "X3", 1000, 2)
  rg_sub  <- generate_data("D1", "X3", 250, 2)
  stopifnot(identical(unname(rg_sub$region), unname(rg_full$region[1:250])))
}