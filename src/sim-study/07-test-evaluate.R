# B2: independent response prediction on a fixed test realization per DGP.
# Source after the worker definitions. This file defines helpers only.

# Replicate 101 is reserved for testing in the R=100 study. Each DGP has ONE
# common covariate/error/region realization; y changes with predictor X1-X3.
make_test_sets <- function(seed_base = 0, replicate = 101L, n_test = 500L) {
  stopifnot(n_test > 0, n_test == floor(n_test))
  gamma <- car_setup$draw_gamma(1)
  dgps <- paste0("D", 1:4)
  setNames(lapply(dgps, function(dgp) {
    seed <- data_seed(dgp, replicate, seed_base)
    set.seed(seed)
    x1 <- runif(n_test, -.5, .5); x2 <- runif(n_test, -.5, .5)
    eps <- error_draw(dgp, n_test)
    # Assign regions after covariates/errors; at n_test=500 each has 25 rows.
    region <- rep(seq_len(N_REGIONS), length.out = n_test)[sample.int(n_test)]
    dat <- data.frame(x1 = x1, x2 = x2, region = factor(region, levels = seq_len(N_REGIONS)))
    y <- setNames(lapply(paste0("X", 1:3), function(pred)
      unname(f_true(pred, x1, x2, gamma, region)$f + eps)), paste0("X", 1:3))
    list(dgp = dgp, replicate = replicate, seed = seed, n = n_test,
         dat = dat, eps = eps, gamma = gamma, y = y)
  }), dgps)
}

# Reuse the fitted coordinate system. Scaling, spline constraints/knots and
# regional coding come from TRAINING, never from the new covariate sample.
newdata_design <- function(s, model, newdata) {
  dat <- s$scenario$dat; dd <- s$designs[[model]]
  spatial <- s$scenario$predictor == "X3"
  if (spatial) {
    newdata$region <- factor(as.character(newdata$region), levels = levels(dat$region))
    stopifnot(!anyNA(newdata$region))
  }
  if (model == "m1") {
    sm <- dd$full_info$mu$smooth.construct
    return(cbind(1, do.call(cbind, lapply(sm, mgcv::PredictMat, data = newdata))))
  }
  if (model == "m3") {
    X <- cbind(1, mgcv::PredictMat(dd$bs1, data.frame(x = newdata$x1)),
                  mgcv::PredictMat(dd$bs2, data.frame(x = newdata$x2)))
    if (spatial) X <- cbind(X, dd$region_contrasts[match(as.character(newdata$region), dd$region_levels), , drop = FALSE])
    return(X)
  }
  scaled <- newdata
  for (key in c("x1", "x2")) scaled[[key]] <- (newdata[[key]] - mean(dat[[key]])) / sd(dat[[key]])
  if (model == "m4") {
    form <- if (spatial) ~ x1 + x2 + region else ~ x1 + x2
    return(model.matrix(form, scaled, contrasts.arg = attr(dd$X, "contrasts")))
  }
  X <- cbind(1, scaled$x1, scaled$x2)
  if (model %in% c("m0", "m5")) X <- cbind(X,
    mgcv::PredictMat(dd$dr1$smooth_object, data.frame(x = newdata$x1)) %*% dd$dr1$Trafo,
    mgcv::PredictMat(dd$dr2$smooth_object, data.frame(x = newdata$x2)) %*% dd$dr2$Trafo)
  if (spatial) X <- cbind(X, dd$X_GMRF[match(as.character(newdata$region), as.character(dat$region)), , drop = FALSE])
  X
}

# Rows are joint posterior draws, columns are new observations. M3 has one row.
# Use raw stored coefficients with the raw design, including its intercept.
predict_unseen <- function(s, model, newdata) {
  fit <- s$fits[[model]]
  B <- if (model == "m3") matrix(fit$beta, nrow = 1L) else fit$draws
  if (is.null(B)) return(NULL)
  X <- newdata_design(s, model, newdata)
  stopifnot(ncol(X) == ncol(B))
  unname(B %*% t(X))
}

# Stable log averages: log(mean(exp(x))) differs from mean(x).
# All -Inf rows retain -Inf, avoiding the undefined subtraction -Inf - -Inf.
row_log_mean_exp <- function(x) {
  m <- apply(x, 1L, max)
  finite <- is.finite(m)
  out <- m
  out[finite] <- m[finite] + log(rowMeans(exp(x[finite, , drop = FALSE] - m[finite])))
  out
}

log_add_exp <- function(a, b) {
  m <- pmax(a, b); finite <- is.finite(m)
  out <- m
  out[finite] <- m[finite] + log(exp(a[finite] - m[finite]) + exp(b[finite] - m[finite]))
  out
}

# Full normalized GUD log densities, matching the fitted Stan likelihoods.
# residual is draw-by-observation; each parameter vector shares those draw rows.
gud_log_density <- function(residual, family, parameters) {
  p <- as.data.frame(parameters)
  left <- residual < 0
  if (family == "FG") {
    z1 <- residual / p$sigma1; z2 <- residual / p$sigma2
    a <- log(p$w) - log(p$sigma1) + z1 - exp(z1)
    b <- log1p(-p$w) - log(p$sigma2) - z2 - exp(-z2)
    return(log_add_exp(a, b))
  }
  if (family == "TPSC") {
    scale <- p$sigma * sqrt(ifelse(left, p$w / (1 - p$w), (1 - p$w) / p$w))
    return(log(2) + ifelse(left, log(p$w), log1p(-p$w)) +
             dt(residual / scale, df = p$delta, log = TRUE) - log(scale))
  }
  # DTP's mixture weight depends on both scales and degrees of freedom.
  a <- log(p$sigma1) + dt(0, df = p$delta2, log = TRUE)
  b <- log(p$sigma2) + dt(0, df = p$delta1, log = TRUE)
  normalizer <- log_add_exp(a, b)
  log(2) + ifelse(left,
    a - normalizer + dt(residual / p$sigma1, df = p$delta1, log = TRUE) - log(p$sigma1),
    b - normalizer + dt(residual / p$sigma2, df = p$delta2, log = TRUE) - log(p$sigma2))
}

# Normalize the native smoothed median-ALD kernel after fitting. For a=sqrt(eps)
# and z=a/(2*sigma^2), its missing integral is Z=z*K_1(z). Using exp(z)*K_1(z)
# avoids underflow; rationalizing sqrt(r^2+eps)-a avoids cancellation near r=0.
# The normalized log density is -log(2*a)-log(exp(z)*K_1(z))
#                              -(sqrt(r^2+eps)-a)/(2*sigma^2).
smoothed_ald_log_density <- function(residual, sigma, eps = .01) {
  if (eps == 0) return(log(.25) - 2*log(sigma) - .5*abs(residual)/sigma^2)
  a <- sqrt(eps)
  log_z <- log(a) - log(2) - 2*log(sigma)
  # Outside the representable Bessel argument range, its limiting forms are
  # accurate to floating-point precision. This is a numerical range safeguard.
  cutoff <- log(.Machine$double.xmax) - 2
  log_k <- -log_z
  large <- log_z > cutoff
  log_k[large] <- .5 * (log(pi/2) - log_z[large])
  regular <- abs(log_z) <= cutoff
  log_k[regular] <- log(besselK(exp(log_z[regular]), nu = 1, expon.scaled = TRUE))
  root <- sqrt(residual^2 + eps)
  delta <- abs(residual) * (abs(residual) / (root + a))
  delta[is.infinite(root)] <- abs(residual)[is.infinite(root)]
  loss <- exp(log(delta) - log(2) - 2*log(sigma))
  -log(2*a) - log_k - loss
}

# bamlss uses sigma^2 as the ALD scale. Predict all scale draws in native row
# order, then transpose to the coefficient-draw convention. Use the saved eps;
# normalization affects scoring only, not the native fit or posterior draws.
m1_log_density <- function(fit, newdata, eta, y) {
  sigma <- as.matrix(predict(fit$fit, newdata = newdata, model = "sigma", type = "parameter", FUN = identity))
  sigma <- t(matrix(as.numeric(sigma), nrow = nrow(newdata), ncol = nrow(eta)))
  eps <- get("eps", environment(fit$fit$family$d))
  smoothed_ald_log_density(sweep(-eta, 2L, y, "+"), sigma, eps)
}

# One Gaussian residual KDE per fit. Leave residuals uncentered, preserving
# location shifts in misspecified scenarios. bw=NULL selects bw.nrd0.
fit_residual_kde <- function(residuals, bw = NULL) {
  if (is.null(bw)) bw <- stats::bw.nrd0(residuals)
  stopifnot(length(bw) == 1L, is.finite(bw), bw > 0)
  list(residuals = residuals, bandwidth = bw)
}

# Evaluate the exact KDE, not a truncated/interpolated density grid. Gaussian
# kernels are combined on the log scale so distant test residuals stay usable.
# Chunks limit temporary memory; they contain queries, not subsampled kernels.
kde_log_density <- function(x, kde, chunk_size = 256L) {
  stopifnot(chunk_size >= 1, chunk_size == floor(chunk_size))
  out <- numeric(length(x))
  for (first in seq.int(1L, length(x), by = chunk_size)) {
    i <- seq.int(first, min(length(x), first + chunk_size - 1L))
    kernels <- dnorm(outer(x[i], kde$residuals, "-") / kde$bandwidth, log = TRUE) - log(kde$bandwidth)
    out[i] <- row_log_mean_exp(kernels)
  }
  out
}

# B2 scores use observed test RESPONSES. No refitting, bootstrap, test-based
# family selection or test-based bandwidth tuning is performed here.
evaluate_test_sample <- function(s, test_set, bw = NULL, chunk_size = 256L) {
  stopifnot(identical(s$scenario$dgp, test_set$dgp))
  if (s$scenario$replicate == test_set$replicate)
    stop("The reserved test replicate cannot also be used as training data.")
  dat <- test_set$dat
  y <- test_set$y[[s$scenario$predictor]]
  scores <- do.call(rbind, lapply(names(s$fits), function(model) {
    fit <- s$fits[[model]]
    eta <- predict_unseen(s, model, dat)
    type <- if (model %in% c("m0", "m2", "m5")) "KDE_posterior_predictive" else
      if (model == "m3") "KDE_plugin" else "posterior_predictive"
    row <- data.frame(model = model, RMSE = NA_real_, MAE = NA_real_, log_score = NA_real_,
                      score_type = type, bandwidth = NA_real_, status = "OK")
    if (is.null(eta)) {row$status <- "FAILURE"; return(row)}
    error <- y - colMeans(eta)
    row$RMSE <- sqrt(mean(error^2)); row$MAE <- mean(abs(error))
    residual <- sweep(-eta, 2L, y, "+")
    if (model %in% c("m0", "m2", "m3", "m5")) {
      B <- if (model == "m3") matrix(fit$beta, nrow = 1L) else fit$draws
      train_hat <- drop(s$designs[[model]]$X %*% colMeans(B))
      kde <- fit_residual_kde(s$scenario$dat$y - train_hat, bw)
      row$bandwidth <- kde$bandwidth
      log_density <- matrix(kde_log_density(as.numeric(residual), kde, chunk_size), nrow = nrow(eta))
    } else if (model == "m1" && !is.null(fit$fit)) {
      log_density <- m1_log_density(fit, dat, eta, y)
    } else if (model == "m4" && !is.null(fit$distribution_draws)) {
      stopifnot(nrow(fit$distribution_draws) == nrow(eta))
      log_density <- gud_log_density(residual, fit$selected, fit$distribution_draws)
    } else {
      # Older saved M4 objects may have beta draws but no distribution draws.
      row$status <- "MISSING_DISTRIBUTION_DRAWS"
      return(row)
    }
    row$log_score <- mean(row_log_mean_exp(t(log_density)))
    row
  }))
  rownames(scores) <- NULL
  list(scenario = s$scenario[c("dgp", "predictor", "n", "replicate", "seed")],
       test = test_set[c("dgp", "replicate", "seed", "n")],
       settings = s$settings, selected_m4 = s$fits$m4$selected,
       kde_bandwidth_rule = if (is.null(bw)) "nrd0" else "fixed", scores = scores)
}
