# ============================================================
# 03-simulation/03-run-simulation.R
# Simulation driver: 4 distributions x 3 predictors x 2 sample sizes
# ============================================================

library(quantreg)
library(qgam)     # additive median regression
library(mgcv)
library(loo)      # ELPD for M4

source("03-simulation/01-dgp.R")
source("03-simulation/02-models.R")

## One replicate: generate data, fit all models, compute RMSEs.
run_replicate <- function(distr, pred, n, seed) {
  sim <- generate_data(distr = distr, pred = pred, n = n, seed = seed)
  dat <- sim$dat
  parts_true <- f0_parts(dat, pred)

  t0 <- Sys.time()
  fits <- list(
    M1 = tryCatch(fit_M1(dat, pred), error = function(e) e),
    M2 = tryCatch(fit_M2(dat, pred, seed = seed), error = function(e) e),
    M3 = tryCatch(fit_M3(dat, pred), error = function(e) e),
    M4 = tryCatch(select_gud_family(dat), error = function(e) e)
  )
  runtime_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  get_f <- function(f) if (inherits(f, "error")) NA else
    if (!is.null(f$fhat)) f$fhat else
    if (!is.null(f$fits)) f$fits$TPSC$fhat else NA  # M4: TPSC fit's fhat

  rmse <- vapply(fits, function(f) {
    fh <- get_f(f)
    if (length(fh) == 1 && is.na(fh)) return(NA_real_)
    sqrt(mean((fh - sim$f0)^2))
  }, 0)

  ## Effectwise RMSE for M2 (linear / nonlinear parts of each effect)
  eff_rmse <- if (inherits(fits$M2, "error")) c(NA, NA) else
    c(sqrt(mean((fits$M2$eff[[1]]$overall - parts_true$lin1 -
                 parts_true$non1)^2)),
      sqrt(mean((fits$M2$eff[[2]]$overall - parts_true$lin2 -
                 parts_true$non2)^2)))

  ## Coverage/width at nominal 90%.
  ## M2: credible interval of the mode function. M4: credible interval of
  ## X beta (ci) AND posterior predictive ystar interval.
  covwid <- function(ci_low, ci_high, truth)
    c(cov = mean(truth >= ci_low & truth <= ci_high),
      wid = mean(ci_high - ci_low))
  cw <- list(
    M2 = if (inherits(fits$M2, "error")) NA else
      covwid(fits$M2$ci[1, ], fits$M2$ci[2, ], sim$f0),
    M4 = if (inherits(fits$M4, "error")) NA else
      covwid(fits$M4$fits$TPSC$ci[1, ], fits$M4$fits$TPSC$ci[2, ], sim$f0),
    M4_ystar = if (inherits(fits$M4, "error") ||
                   is.null(fits$M4$fits$TPSC$ystar_ci)) NA else
      covwid(fits$M4$fits$TPSC$ystar_ci[1, ],
             fits$M4$fits$TPSC$ystar_ci[2, ], dat$y))

  sel <- if (inherits(fits$M4, "error")) NA else fits$M4$selected
  elpd <- if (inherits(fits$M4, "error")) NA else fits$M4$elpd

  data.frame(distr = distr, pred = pred, n = n, seed = seed,
             rmse_M1 = rmse[["M1"]], rmse_M2 = rmse[["M2"]],
             rmse_M3 = rmse[["M3"]], rmse_M4 = rmse[["M4"]],
             eff1_rmse = eff_rmse[1], eff2_rmse = eff_rmse[2],
             cov_M2 = cw$M2["cov"], wid_M2 = cw$M2["wid"],
             cov_M4 = cw$M4["cov"], wid_M4 = cw$M4["wid"],
             cov_M4_ystar = cw$M4_ystar["cov"],
             wid_M4_ystar = cw$M4_ystar["wid"],
             gud_selected = sel,
             elpd_TPSC = unname(elpd["TPSC"]),
             elpd_FG = unname(elpd["FG"]),
             elpd_DTP = unname(elpd["DTP"]),
             runtime_s = runtime_s)
}

## Full grid: 4 x 3 x 2 = 24 configurations, R replicates each.
run_grid <- function(distrs = c("D1", "D2", "D3", "D4"),
                     preds   = c("P1", "P2", "P3"),
                     ns      = c(250, 1000),
                     R       = 100, seed0 = 2024) {
  grid <- expand.grid(distr = distrs, pred = preds, n = ns,
                      rep = seq_len(R), stringsAsFactors = FALSE)
  grid$seed <- seed0 + seq_len(nrow(grid))
  results <- do.call(rbind, Map(function(d, p, nn, sd)
    run_replicate(d, p, nn, sd),
    grid$distr, grid$pred, grid$n, grid$seed))
  attr(results, "grid") <- grid
  results
}

## Example: one setting first, then the full grid.
res_test <- run_replicate("D3", "P2", n = 250, seed = 1)
# res <- run_grid()
# saveRDS(res, "results/sim_grid.rds")

## Parallelization (once everything runs): replace Map() with
## future.apply::future_mapply() and plan(multisession, workers = 4).