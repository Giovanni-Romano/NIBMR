# === === === === === === === === === === === === === === === === === === === ===
# 03c-run-2025-logit.R
#
# Initial values (quantile regression -> asymmetric-loss optimisation), tuning
# of the loss asymmetry k, and the IWLS MCMC run.
#
# Input : application/bike-sharing/rds_prep/prep_<MODEL_TAG>.rds
# Output: application/bike-sharing/rds_mcmc/out_MCMC_<MODEL_TAG>.rds
# === === === === === === === === === === === === === === === === === === === ===

suppressPackageStartupMessages({
  library(tidyverse)
  library(lqmm)
  library(sdPrior)
})

source("src/bike-sharing/03a-config-run.R")
source("src/01a-utils.R")
source("src/01b-gibbs.R")

prep <- readRDS(PATH_PREP)

n <- prep$dims$n; p <- prep$dims$p; q <- prep$dims$q
s <- prep$dims$s; d <- prep$dims$d

df_lqmm            <- prep$df_lqmm
X_complete_grouped <- prep$X_complete_grouped
X_complete_unique  <- prep$X_complete_unique
X_design_grouped   <- drop_id_col(X_complete_grouped)
X_design_unique    <- drop_id_col(X_complete_unique)

y_grouped    <- prep$grouped$y
freq_grouped <- prep$grouped$freq
labs_unique  <- prep$labels$unique

# ---- Median regression starting point ----------------------------------------
# Linear-quantile mixed model on the logit scale, with a station random intercept.
covariate_names <- setdiff(colnames(df_lqmm), c("y", "id_re"))
qr_formula <- as.formula(
  paste("qlogis(y) ~", paste(covariate_names, collapse = " + "))
)

qr_reg <- lqmm(
  fixed   = qr_formula,
  random  = ~ 1,
  group   = id_re,
  tau     = 0.5,
  data    = df_lqmm,
  weights = freq_grouped
)

res_qr  <- residuals(qr_reg)
beta_qr <- coef(qr_reg)
names(beta_qr) <- colnames(X_design_grouped)

unique_fit_qr <- plogis(
  drop(X_design_unique[, names(beta_qr)] %*% beta_qr) +
    ranef(qr_reg)[["(Intercept)"]][X_complete_unique[, "id_re"]]
)
names(unique_fit_qr) <- labs_unique

# ---- Asymmetry/flatness parameter k -------------------------------------------
# k_init follows the n^(5/12) log(n)^(7/12) / n^(7/12) rate; the two-sided
# multipliers rescale it by the second moment of the positive / negative
# quantile-regression residuals.
const  <- (n^(5 / 12)) * (log(n)^(7 / 12)) / (n^(7 / 12))
k_init <- (n / const)^(1 / 7)
k      <- k_init
# k_final <- k_init * n^(1/7)   # asymptotic target, currently unused

res_pos <- res_qr[res_qr > 0]
res_neg <- res_qr[res_qr < 0]

factor_pos <- sqrt(var(res_pos) + (mean(res_pos) - mean(res_qr))^2)
factor_neg <- sqrt(var(res_neg) + (mean(res_neg) - mean(res_qr))^2)

k_multi <- c(sqrt(k) / factor_neg, sqrt(k) / factor_pos)

# ---- Refine the starting point under the asymmetric loss ---------------------
asymm_objective <- function(b) {
  predictor <- logit_link(X_design_grouped %*% b)
  sum(freq_grouped * loss_asymm(y_grouped - predictor, k_multi, C_POST))
}

grad_fn <- function(b) -score_pop(1:length(b), b, y_grouped, X_design_grouped,
                                  k_multi, C_POST, freq_grouped, logit_link, logit_prime)
fit_optim <- optim(par = beta_qr, fn = asymm_objective, gr = grad_fn,
                   method = "BFGS",
                   control = list(maxit = 2000, reltol = 1e-10))
b_init <- fit_optim$par

unique_fit_init <- logit_link(drop(X_design_unique %*% b_init))
names(unique_fit_init) <- labs_unique

# ---- Scale-dependent prior for the random-effect variance --------------------
n_re      <- max(prep$grouped$id_re)
theta_SDP <- get_theta(Z = diag(1, n_re), Kinv = diag(1, n_re))

# ---- MCMC --------------------------------------------------------------------
# One variance per binary effect (q) and per linear effect (p), then one per
# non-linear spline block (p) and one for the spatial block (1).
a_tau <- c(rep(A_TAU_LINEAR, q + p), rep(A_TAU_SMOOTH, p + 1))
b_tau <- c(rep(B_TAU_LINEAR, q + p), rep(B_TAU_SMOOTH, p + 1))

out_MCMC <- MCMC_IWLS_bikes(
  niter       = N_ITER,
  X_dummy     = prep$grouped$X_dummy,
  X_cont      = prep$grouped$X_cont,
  Z_cont      = prep$grouped$Z_cont,
  Z_spat      = prep$grouped$Z_spat,
  id          = prep$grouped$id_re,
  y           = y_grouped,
  freq        = freq_grouped,
  link        = logit_link,
  link_prime  = logit_prime,
  link_second = logit_second,
  k           = k_multi,
  k_deriv     = k_multi^K_PROP_EXP,
  c           = C_POST,
  c_deriv     = C_PROP,
  d           = d,
  a_tau       = a_tau,
  b_tau       = b_tau,
  theta_re    = theta_SDP$root,
  verbose     = MCMC_VERBOSE,
  K_fold      = MCMC_K_FOLD,
  print_step  = MCMC_PRINT_EVERY,
  step        = MCMC_STEP,
  beta0       = b_init,
  seed        = SEED
)

# Carry the initial fits along so 03c needs neither lqmm nor a re-optimisation.
out_MCMC$init <- list(
  beta_qr         = beta_qr,
  unique_fit_qr   = unique_fit_qr,
  b_init          = b_init,
  unique_fit_init = unique_fit_init,
  k_multi         = k_multi,
  theta_SDP       = theta_SDP$root
)

saveRDS(out_MCMC, PATH_MCMC)
message("Saved MCMC output to ", PATH_MCMC)
