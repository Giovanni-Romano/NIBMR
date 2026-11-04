# ============================================================
# 03-simulation/01-dgp.R
# Data generating processes: modular predictors + error distributions
# ============================================================

# ---- (P) Predictors ----
## Each function takes (x1, x2, s) and returns the true mode function f0.
## s is a two-column matrix of spatial coordinates, only used by (P3).

# Effect of x1 decomposed into linear + nonlinear parts,
# following the L2 projection convention of Bach & Klein (2024, Def. 1).
f1_lin    <- function(x1) -6 / pi * (x1 - 0.5)
f1_nonlin <- function(x1) sin(2 * pi * x1) + 6 / pi * (x1 - 0.5)

f0_linear    <- function(x1, x2, s = NULL) {
  1 + 2 * x1 - 1.5 * x2
}

f0_nonlinear <- function(x1, x2, s = NULL) {
  f1_lin(x1) + f1_nonlin(x1) - 1.5 * x2
}

f0_spatial   <- function(x1, x2, s) {
  f1_lin(x1) + f1_nonlin(x1) + 4 * (x2 - 0.5)^2 - 1/3 + f_spat(s)
}

f_spat <- function(s) {
  center <- matrix(c(0.4, 0.6), nrow = nrow(s), ncol = 2, byrow = TRUE)
  1.5 * exp(-rowSums((s - center)^2) / 0.05)
}

# ---- (D) Error distributions ----
## Each function takes (n, x1) and returns n errors with mode 0.
## x1 is only needed by the heteroskedastic distribution (D4).

error_normal <- function(n, x1 = NULL) {
  rnorm(n)
}

error_t3 <- function(n, x1 = NULL) {
  rt(n, df = 3) / sqrt(3)          # Var(t_3) = 3, rescaled to unit variance
}

## ---- TPSC-Student-t (Liu, Huang & Bai 2024, eq. 7) ------------------
## f_TPSC(y | w, theta, sigma, delta) =
##   w * f_LT(y | theta, sigma*sqrt(w/(1-w)), delta) +
##   (1-w) * f_RT(y | theta, sigma*sqrt((1-w)/w), delta)
## with f_LT / f_RT the left/right theta-truncated Student-t densities
## sharing df = delta. Note: w > 0.5  =>  LEFT-skewed (Fernandez-Steel
## convention); w < 0.5 => right-skewed. At w = 0.5, delta -> Inf the
## density is N(theta, sigma^2).

## Exact inverse-CDF sampler for a one-sided truncated Student-t.
## For side = "left":  y = theta + scale * qt(u * 0.5,          df)
## For side = "right": y = theta + scale * qt(0.5 + u * 0.5,    df)
## (valid because pt(0, df) = 0.5 for the centered t).
r_t_trunc <- function(n, theta = 0, scale = 1, df = 3,
                      side = c("left", "right")) {
  side <- match.arg(side)
  u <- runif(n)
  p <- if (side == "left") 0.5 * u else 0.5 + 0.5 * u
  theta + scale * qt(p, df = df)
}

## TPSC-Student-t with mode theta, following the paper's sampling
## algorithm: Z ~ Bernoulli(w) selects the component, then a draw
## from the corresponding truncated t.
r_tpsc_t <- function(n, theta = 0, sigma = 1, delta = 3, w = 0.7) {
  stopifnot(w > 0, w < 1, delta > 0, sigma > 0)
  z     <- rbinom(n, 1, w)
  left  <- z == 1                       # left-skewed component, prob w
  scale_l <- sigma * sqrt(w / (1 - w))
  scale_r <- sigma * sqrt((1 - w) / w)
  out <- numeric(n)
  out[left]  <- r_t_trunc(sum(left),  theta, scale_l, delta, "left")
  out[!left] <- r_t_trunc(sum(!left), theta, scale_r, delta, "right")
  out
}

## D3 error distribution: TPSC-Student-t, mode 0.
## w = 0.7 => left-skewed (median < 0 = mode), delta = 3 => fat tails.
error_tpsc <- function(n, x1 = NULL, sigma = 1, delta = 3, w = 0.7) {
  r_tpsc_t(n, theta = 0, sigma = sigma, delta = delta, w = w)
}

## Density for the mode/sanity checks (and for plotting against
## the estimated densities, as in Fig. 5 of Liu et al.).
d_tpsc_t <- function(y, theta = 0, sigma = 1, delta = 3, w = 0.7) {
  scale_l <- sigma * sqrt(w / (1 - w))
  scale_r <- sigma * sqrt((1 - w) / w)
  w * dt((y - theta) / scale_l, df = delta) / scale_l * (y <= theta) +
    (1 - w) * dt((y - theta) / scale_r, df = delta) / scale_r * (y > theta)
}

## Shifted Gamma (mode 0), homoskedastic.
## mode(Gamma(shape, rate)) = (shape - 1)/rate, so we shift by that amount.
error_gamma_homosk <- function(n, x1 = NULL, shape = 1.25, rate = 0.5) {
  rgamma(n, shape = shape, rate = rate) - (shape - 1) / rate
}

## Shifted Gamma with multiplicative heteroskedasticity (D4).
## mode stays 0 (positive scaling preserves the mode);
## heteroskedasticity violates the s1, s2 assumption of our loss.
error_gamma_hetero <- function(n, x1 = stop("x1 required"),
                               shape = 1.25, rate = 0.5) {
  h <- 0.5 + 2 * x1
  (rgamma(n, shape = shape, rate = rate) - (shape - 1) / rate) * h
}

## ---- Dispatchers ----------------------------------------------------

ERROR_DISTS <- list(
  D1 = list(fun = error_normal,           needs_x = FALSE),
  D2 = list(fun = error_t3,               needs_x = FALSE),
  D3 = list(fun = error_tpsc,             needs_x = FALSE),
  D4 = list(fun = error_gamma_hetero,     needs_x = TRUE)
)

PREDICTORS <- list(
  P1 = list(fun = f0_linear,    needs_s = FALSE),
  P2 = list(fun = f0_nonlinear, needs_s = FALSE),
  P3 = list(fun = f0_spatial,   needs_s = TRUE)
)

generate_data <- function(distr = "D1", pred = "P1", n = 250, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  x1 <- runif(n)
  x2 <- runif(n)
  s  <- if (PREDICTORS[[pred]]$needs_s) cbind(runif(n), runif(n)) else NULL

  f0  <- do.call(PREDICTORS[[pred]]$fun, list(x1 = x1, x2 = x2, s = s))
  eps <- if (ERROR_DISTS[[distr]]$needs_x) {
    ERROR_DISTS[[distr]]$fun(n, x1 = x1)
  } else {
    ERROR_DISTS[[distr]]$fun(n)
  }

  dat <- data.frame(y = f0 + eps, x1 = x1, x2 = x2)
  if (!is.null(s)) { dat$spat1 <- s[, 1]; dat$spat2 <- s[, 2] }

  list(dat = dat, f0 = f0, eps = eps, x1 = x1, x2 = x2, s = s,
       distr = distr, pred = pred, n = n)
}

## Sanity check: empirical mode of the generated errors should be ~0.
distr_mode_check <- function(distr, n = 2e6, seed = 1) {
  set.seed(seed)
  x1  <- runif(n)
  eps <- if (ERROR_DISTS[[distr]]$needs_x) {
    ERROR_DISTS[[distr]]$fun(n, x1 = x1)
  } else ERROR_DISTS[[distr]]$fun(n)
  d <- density(eps, n = 2^12)
  d$x[which.max(d$y)]
}
# sapply(c("D1","D2","D3","D4"), distr_mode_check)   # all should be close to 0