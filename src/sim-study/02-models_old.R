# ============================================================
# 03-simulation/02-models.R
# Fitting wrappers for M1-M4
# ============================================================

# ---- NIL loss (shared by M2/M3)

## Oelker et al. (2015): k_initial = (n/const)^(1/7),
## const = n^(5/12) log(n)^(7/12) / n^(7/12).
## Algebraically this is n^(1/6) / log(n)^(1/12); both forms below cross-check.
k_initial <- function(n) {
  const  <- n^(5/12) * log(n)^(7/12) / n^(7/12)
  ki     <- (n / const)^(1/7)
  stopifnot(abs(ki - n^(1/6) / log(n)^(1/12)) < 1e-10)  # sanity check
  ki
}

## s1, s2: spread of the negative / positive side of the error distribution,
## estimated from the residuals of a preliminary median regression.
## TODO: align the exact definition with `loss_asymm` in 01a-utils.R.
estimate_s1_s2 <- function(resid) {
  neg <- resid[resid < 0]
  pos <- resid[resid >= 0]
  c(s1 = median(abs(neg)), s2 = median(pos))
}

## Individual asymmetric NIL loss: k_i = k_n/s1 if e < 0, k_n/s2 if e >= 0.
## TODO: reuse `loss_asymm` from 01a-utils.R instead of this copy.
nil_loss_sum <- function(beta, X, y, k_n, c_n, s1, s2) {
  e   <- y - as.vector(X %*% beta)
  k_i <- ifelse(e < 0, k_n / s1, k_n / s2)
  sum(1 - exp(sqrt(c_n) - sqrt(k_i^2 * e^2 + c_n)))
}

# ---- Design matrix for M3 (frequentist NIL)
## Unpenalized B-splines to keep `optim` manageable.
## TODO: if we want penalized estimation, switch to the IRLS machinery
## of Oelker et al. (2015) (NILS) instead.
make_design_M3 <- function(dat, pred) {
  if (pred == "P1") return(cbind(1, dat$x1, dat$x2))
  B1 <- splines::bs(dat$x1, df = 8)
  B2 <- splines::bs(dat$x2, df = 8)
  if (pred == "P2") return(cbind(1, B1, dat$x2))
  # P3: tensor-product smooth for the spatial coordinates via mgcv
  sm <- mgcv::smoothCon(mgcv::s(spat1, spat2, k = 25), data = dat)[[1]]
  cbind(1, B1, B2, sm$X)
}

# ---- M1: median regression
fit_M1 <- function(dat, pred) {
  if (pred == "P1") {
    return(quantreg::rq(y ~ x1 + x2, tau = 0.5, data = dat))
  }
  ## Additive median regression via qgam (must be installed).
  ## TODO: confirm predict() semantics for qgam objects.
  form <- switch(pred,
                 "P2" = y ~ s(x1) + s(x2),
                 "P3" = y ~ s(x1) + s(x2) + te(spat1, spat2))
  qgam::qgam(form, qu = 0.5, data = dat)
}

fhat_M1 <- function(fit, dat, pred) {
  if (pred == "P1") as.vector(predict(fit, newdata = dat))
  else as.vector(predict(fit, newdata = dat))
}

# ---- M2: our Bayesian model (linear and full)
# ============================================================
# 03-simulation/02-models.R  — FINAL fit_M2 / fit_M4 (replaces stubs)
# ============================================================

## ---- Design for M2 ---------------------------------------------------
## X = (1, xbar_1..xbar_p, Z_1..Z_p) as expected by MCMC_IWLS.
## d = c(0,0) for the linear predictor (P1): sum(d) = 0, purely linear.
make_design_M2 <- function(dat, pred, dj = 9) {
  x1_bar <- dat$x1 - mean(dat$x1)
  x2_bar <- dat$x2 - mean(dat$x2)
  switch(pred,
    "P1" = list(X = cbind(1, x1_bar, x2_bar), p = 2, d = c(0, 0),
                Z = list(NULL, NULL)),
    "P2" = {
      Z1 <- construct_DR_basis(dat$x1, d = dj)$Z
      list(X = cbind(1, x1_bar, x2_bar, Z1), p = 2, d = c(dj, 0),
           Z = list(Z1, NULL))
    },
    "P3" = {
      Z1 <- construct_DR_basis(dat$x1, d = dj)$Z
      Z2 <- construct_DR_basis(dat$x2, d = dj)$Z
      Zs <- construct_DR_Spatial(dat$spat1, dat$spat2)$Z   # check args vs docs
      ## Spatial columns folded into the last nonlinear block:
      ## d[2] = dj + ncol(Zs). The spatial coefficients then share the
      ## tau of the x2-block. TODO (cleaner): extend MCMC_IWLS with a
      ## separate spatial tau; MCMC_IWLS_bikes has Z_spat but drags in
      ## random-effect machinery (id, theta_re) we do not need here.
      list(X = cbind(1, x1_bar, x2_bar, Z1, Z2, Zs), p = 2,
           d = c(dj, dj + ncol(Zs)), Z = list(Z1, cbind(Z2, Zs)))
    })
}

## True linear/nonlinear parts of each effect, for effectwise RMSE
## (Bach & Klein 2024, Def. 1; both parts are centered by construction).
f0_parts <- function(dat, pred) {
  x1 <- dat$x1; x2 <- dat$x2
  lin1 <- -6 / pi * (x1 - 0.5)
  non1 <- sin(2 * pi * x1) + 6 / pi * (x1 - 0.5)
  switch(pred,
    "P1" = list(lin1 = 2 * x1 - 2 * 0.5, non1 = rep(0, length(x1)),
                lin2 = -1.5 * x2 + 1.5 * 0.5, non2 = rep(0, length(x2)),
                spat = rep(0, length(x1))),
    "P2" = list(lin1 = lin1, non1 = non1,
                lin2 = -1.5 * x2 + 1.5 * 0.5, non2 = rep(0, length(x2)),
                spat = rep(0, length(x1))),
    "P3" = list(lin1 = lin1, non1 = non1,
                lin2 = rep(0, length(x2)),           # x2 effect is purely
                non2 = 4 * (x2 - 0.5)^2 - 1/3,       # nonlinear in (P3)
                spat = f_spat(cbind(dat$spat1, dat$spat2))))
}

## ---- M2: generalized Bayesian NIL mode regression --------------------
fit_M2 <- function(dat, pred, niter = 5000, burnin = 1000, dj = 9,
                   c_n = 1e-5, sf = 10, level = 0.90,
                   seed = NULL, verbose = 0) {
  des <- make_design_M2(dat, pred, dj)
  X <- des$X; p <- des$p; d <- des$d
  n <- nrow(X)

  kn <- k_initial(n)
  med <- quantreg::rq.fit(X, dat$y, tau = 0.5)
  sp  <- estimate_s1_s2(as.vector(X %*% med$coefficients - dat$y))
  k   <- c(kn / sp[["s1"]], kn / sp[["s2"]])
  ## Smoother proposal loss for IWLS (notes: increase c, decrease k).
  ## TODO: align sf with the convention used in the project scripts.
  k_deriv <- k / sf
  c_deriv <- c_n * sf

  a_tau <- rep(1, 2 * p)   # TODO: confirm IG hypervalues with the notes
  b_tau <- rep(1, 2 * p)

  fit <- MCMC_IWLS(niter, X, dat$y, k, k_deriv, c_n, c_deriv,
                   p, d, med$coefficients, a_tau, b_tau,
                   seed = seed, verbose = verbose)

  beta <- fit$beta[(burnin + 1):niter, , drop = FALSE]
  fhat <- as.vector(X %*% colMeans(beta))
  ff   <- X %*% t(beta)                       # n x ndraw function draws
  ci   <- apply(ff, 1, quantile,
                probs = c((1 - level) / 2, 1 - (1 - level) / 2))

  ## Effectwise estimates: overall effect j = beta_lin_j * xbar_j + Z_j b,
  ## split into linear part (beta_lin_j * xbar_j) and nonlinear part (Z_j b).
  eff <- lapply(seq_len(p), function(j) {
    Zj <- des$Z[[j]]
    lin <- beta[, 1 + j] * (X[, 1 + j])          # beta_lin_j * xbar_j
    non <- if (is.null(Zj)) rep(0, n) else
      as.vector(Zj %*% colMeans(beta[, grep(
        paste0("beta_nonlin", j, "\\."), colnames(beta)), drop = FALSE]))
    list(lin = lin, non = non, overall = lin + non)
  })

  list(fhat = fhat, ci = ci, eff = eff, beta = beta,
       k = k, c_n = c_n, s = sp, w_gbi = fit$w, acc = fit$acc_prop)
}

# ---- M3: frequentist NIL via numerical optimization
fit_M3 <- function(dat, pred, c_n = 1e-5) {
  X  <- make_design_M3(dat, pred)
  y  <- dat$y
  n  <- nrow(dat)
  kn <- k_initial(n)

  # s1, s2 from median residuals on the SAME design
  med_fit  <- quantreg::rq.fit(X, y, tau = 0.5)
  sp       <- estimate_s1_s2(as.vector(X %*% med_fit$coefficients - y))
  s1 <- sp[["s1"]]; s2 <- sp[["s2"]]

  # start from the median fit (the loss is spiky for large k_n,
  # so a good starting point matters)
  start <- med_fit$coefficients
  opt <- optim(par     = start,
               fn      = nil_loss_sum,
               X = X, y = y, k_n = kn, c_n = c_n, s1 = s1, s2 = s2,
               method  = "BFGS",
               control = list(maxit = 5000))
  # TODO: multi-start / optimr fallback if BFGS stalls (see TeX, model M3)
  beta <- opt$par
  list(beta = beta, fhat = as.vector(X %*% beta),
       k_n = kn, c_n = c_n, s1 = s1, s2 = s2)
}

## ---- M4: GUD package (Liu, Huang & Bai 2024)
## API per GitHub README / CRAN: modal_regression(formula, data, model,
## chains, iter) -> rstan object with w, delta, sigma, coefficients,
## and posterior predictive draws ystar[1]..ystar[n].
fit_M4 <- function(dat, model = c("TPSC", "FG", "DTP"),
                   chains = 2, iter = 2000, level = 0.90, ...) {
  model <- match.arg(model)
  fit <- GUD::modal_regression(y ~ x1 + x2, data = dat, model = model,
                               chains = chains, iter = iter, ...)
  draws <- rstan::as.matrix(fit)               # post-warmup, all chains
  cf    <- c("(Intercept)", "x1", "x2")
  fhat  <- as.vector(draws[, cf[1]] +
                       draws[, cf[2]] * dat$x1 +
                       draws[, cf[3]] * dat$x2)
  ff    <- draws[, cf[1]] + outer(draws[, cf[2]], dat$x1) +
           outer(draws[, cf[3]], dat$x2)       # ndraw x n
  ci    <- apply(ff, 2, quantile,
                 probs = c((1 - level) / 2, 1 - (1 - level) / 2))
  ## Posterior predictive intervals (ystar) — only exists for M4;
  ## see note below on the asymmetry with M2's credible intervals.
  ystar_cols <- grep("^ystar\\[", colnames(draws), value = TRUE)
  ystar_ci <- if (length(ystar_cols) == nrow(dat))
    apply(draws[, ystar_cols], 2, quantile,
          probs = c((1 - level) / 2, 1 - (1 - level) / 2)) else NULL

  ## ELPD via PSIS-LOO, as proposed by Liu et al.
  elpd <- tryCatch(
    loo::loo(fit)$estimates["elpd_loo", "Estimate"],
    error = function(e) {
      warning("loo() failed; is log_lik exposed in the Stan model? ", e)
      NA_real_
    })
  list(fhat = fhat, ci = ci, ystar_ci = ystar_ci, elpd = elpd,
       model = model, draws = draws)
}

select_gud_family <- function(dat, ...) {
  mods <- c("TPSC", "FG", "DTP")
  fits <- lapply(mods, function(m) fit_M4(dat, model = m, ...))
  elpd <- vapply(fits, function(f) f$elpd, 0)
  list(selected = mods[which.max(elpd)], elpd = setNames(elpd, mods),
       fits = fits)
}