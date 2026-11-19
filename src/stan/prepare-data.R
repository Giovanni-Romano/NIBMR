# Standalone data preparation for MCMC_IWLS.stan; sourcing has no side effects.
# Uses the design/block order and CV objective in ../01b-mcmc.R.
prepare_mcmc_iwls_stan_data <- function(
    y, k, c, d = integer(0), a_tau = numeric(0), b_tau = numeric(0),
    X_dummy = NULL, X_cont = NULL, Z_cont = NULL, Z_spat = NULL,
    X_GMRF = NULL, K_GMRF = NULL, theta_GMRF = NULL,
    id = NULL, theta_re = NULL, freq = NULL, link = c("identity", "logit"),
    intercept_mean = 0, intercept_sd = 100,
    w = NULL, beta0 = NULL, K_fold = 50L, seed = NULL) {
  link <- match.arg(link)
  positive_scalar <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0)
      stop(name, " must be a finite positive scalar.")
  }
  if (!is.numeric(y) || !length(y) || any(!is.finite(y)))
    stop("y must be a nonempty finite numeric vector.")
  y <- as.numeric(y)
  n <- length(y)
  positive_scalar(c, "c")
  positive_scalar(intercept_sd, "intercept_sd")
  if (!is.numeric(intercept_mean) || length(intercept_mean) != 1L ||
      !is.finite(intercept_mean)) stop("intercept_mean must be finite and scalar.")
  if (!is.numeric(k) || length(k) != 2L || any(!is.finite(k)) || any(k < 0))
    stop("k must contain two finite nonnegative values: negative then positive side.")
  if (!is.null(Z_spat) && !is.null(X_GMRF))
    stop("Supply only one of Z_spat and X_GMRF.")
  if (is.null(X_cont) && !is.null(Z_cont))
    stop("Z_cont requires X_cont.")

  design <- function(x, name) {
    if (is.null(x)) return(matrix(numeric(0), nrow = n, ncol = 0L))
    x <- as.matrix(x)
    if (!is.numeric(x) || nrow(x) != n || !ncol(x) || any(!is.finite(x)))
      stop(name, " must be a finite numeric matrix with length(y) rows and at least one column.")
    x
  }
  X_dummy <- design(X_dummy, "X_dummy")
  X_cont <- design(X_cont, "X_cont")
  q <- ncol(X_dummy)
  p <- ncol(X_cont)
  if (!is.numeric(d) || length(d) != p || any(!is.finite(d)) ||
      any(d < 0 | d != floor(d))) stop("d must contain p nonnegative integers.")
  d <- as.integer(d)
  use_nonlinear <- !is.null(Z_cont)
  if (use_nonlinear && any(d == 0L))
    stop("Each active nonlinear block must have positive dimension d[j].")
  if (!use_nonlinear && sum(d) != 0L)
    stop("Set d to rep(0L, ncol(X_cont)) when Z_cont is NULL.")
  if (is.list(Z_cont) && !is.data.frame(Z_cont)) {
    if (length(Z_cont) != p) stop("A Z_cont list must contain p matrices.")
    blocks <- lapply(seq_len(p), function(j) {
      z <- design(Z_cont[[j]], paste0("Z_cont[[", j, "]]"))
      if (ncol(z) != d[j]) stop("Z_cont block dimensions must match d.")
      z
    })
    Z_cont <- do.call(cbind, blocks)
  }
  Z_cont <- design(Z_cont, "Z_cont")
  if (ncol(Z_cont) != sum(d)) stop("ncol(Z_cont) must equal sum(d).")
  use_gmrf <- !is.null(X_GMRF)
  spatial <- design(if (use_gmrf) X_GMRF else Z_spat, "spatial design")
  s <- ncol(spatial)
  p_nonlin <- if (use_nonlinear) p else 0L
  X <- cbind(Intercept = 1, X_dummy, X_cont, Z_cont, spatial)
  P <- ncol(X)
  n_tau <- q + p + p_nonlin + as.integer(s > 0L)
  tau_idx <- as.integer(c(
    seq_len(q), q + seq_len(p),
    if (use_nonlinear) q + p + rep(seq_len(p), d) else integer(0),
    if (s > 0L) rep(n_tau, s) else integer(0)))
  # The final GMRF a_tau/b_tau are unused, just as in MCMC_IWLS.
  for (nm in c("a_tau", "b_tau")) {
    value <- get(nm)
    if (!is.numeric(value) || length(value) != n_tau || any(!is.finite(value)) ||
        any(value < 0)) stop(nm, " must contain n_tau finite nonnegative values.")
    ig_idx <- seq_len(n_tau - as.integer(use_gmrf))
    if (any(value[ig_idx] <= 0)) stop(nm, " must be positive for inverse-gamma blocks.")
  }
  if (is.null(freq)) freq <- rep(1, n)
  if (length(freq) == 1L) freq <- rep(freq, n)
  if (!is.numeric(freq) || length(freq) != n || any(!is.finite(freq)) ||
      any(freq < 0) || sum(freq) <= 0) stop("freq must have n nonnegative weights and positive total.")

  use_re <- !is.null(id)
  id_levels <- if (use_re) unique(id) else NULL
  if (use_re) {
    if (length(id) != n || anyNA(id)) stop("id must have length(y) entries without missing values.")
    positive_scalar(theta_re, "theta_re")
    id_index <- match(id, id_levels)
  } else {
    id_index <- rep(0L, n)
    theta_re <- 1
  }
  if (use_gmrf) {
    positive_scalar(theta_GMRF, "theta_GMRF")
    K_GMRF <- as.matrix(K_GMRF)
    if (!is.numeric(K_GMRF) || !identical(dim(K_GMRF), c(s, s)) ||
        any(!is.finite(K_GMRF)) || !isSymmetric(K_GMRF, tol = 1e-10))
      stop("K_GMRF must be a finite symmetric s by s matrix.")
    if (inherits(try(chol(K_GMRF), silent = TRUE), "try-error"))
      stop("K_GMRF must be positive definite: the R sampler assumes rank_K = s. See README.md for intrinsic GMRFs.")
  } else {
    K_GMRF <- matrix(numeric(0), 0L, 0L)
    theta_GMRF <- 1
  }

  mean_cv <- NA_real_
  folds <- NULL
  if (is.null(w)) {
    if (!is.numeric(beta0) || length(beta0) != P || any(!is.finite(beta0)))
      stop("Supply finite beta0 of length P for calibration, or supply w directly.")
    if (length(K_fold) != 1L || !is.finite(K_fold) || K_fold != floor(K_fold) ||
        K_fold < 2L || K_fold > n) stop("K_fold must be an integer between 2 and n.")
    if (!is.null(seed)) set.seed(seed)
    folds <- sample(rep(seq_len(K_fold), length.out = n))
    response_mode <- if (link == "identity") function(x) x else plogis
    # Same g=1 loss and unpenalized, fixed-effect-only CV as MCMC_IWLS.
    loss <- function(residual) {
      kr <- k[1L + as.integer(residual > 0)] * residual
      1 - exp(sqrt(c) - sqrt(kr^2 + c))
    }
    cv <- vapply(seq_len(K_fold), function(f) {
      test <- which(folds == f)
      train <- setdiff(seq_len(n), test)
      est <- optim(par = beta0, fn = function(b) {
        eta <- X[train, , drop = FALSE] %*% b
        sum(freq[train] * loss(y[train] - response_mode(eta)))
      }, method = "L-BFGS-B")$par
      sum(freq[test] * loss(y[test] - response_mode(X[test, , drop = FALSE] %*% est)))
    }, numeric(1))
    mean_cv <- sum(cv) / sum(freq)
    positive_scalar(mean_cv, "mean held-out loss")
    w <- (P / 2) / mean_cv
  }
  positive_scalar(w, "w")
  # Preserve vector dimensions even at length one (not a scalar in Stan).
  stan_vector <- function(x) array(as.numeric(x), dim = length(x))
  stan_index <- function(x) array(as.integer(x), dim = length(x))

  list(
    data = list(n = n, P = P, X = X, y = stan_vector(y), freq = stan_vector(freq),
                k = stan_vector(k), c = c, w = w,
                link_id = if (link == "identity") 1L else 2L,
                n_tau = n_tau, tau_idx = stan_index(tau_idx),
                a_tau = stan_vector(a_tau), b_tau = stan_vector(b_tau),
                intercept_mean = intercept_mean, intercept_sd = intercept_sd,
                use_re = as.integer(use_re), n_id = length(id_levels),
                id = stan_index(id_index), theta_re = theta_re,
                use_gmrf = as.integer(use_gmrf), s_gmrf = if (use_gmrf) s else 0L,
                K_GMRF = K_GMRF, theta_GMRF = theta_GMRF),
    calibration = list(w = w, numerator = P / 2, mean_cv = mean_cv),
    fold_id = folds, id_levels = id_levels, coefficient_names = colnames(X))
}
