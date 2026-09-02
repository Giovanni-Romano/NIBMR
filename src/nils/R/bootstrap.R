#' Function performing the bootstrap
#'
#' This function performs the bootstrap for the bootstrap samples generated via \code{\link{bootstrap_samples}}.
#' It can be called by \code{\link{m}} and by \code{\link{m_gam}} and allows for unpenalized as well as penalized models.
#' If \code{lambda} is chosen by cross-validation, the optimal value is estimated for each bootstrap sample.
#' @param B number of bootstrap samples
#' @param B_fit logical; if \code{FALSE}, the estimated coefficients for each bootstrap sample are returned in \code{ciraw}; if \code{TRUE}, the fitted values per sample are returned
#' @param x for \code{\link{m}}: design matrix; for \code{\link{m_gam}}: NULL
#' @param fitted fitted values
#' @param residuals residuals
#' @param d derivative of chosen loss function divided by its argument
#' @param tuning set of tuning parameters
#' @param lambda for \code{\link{m}}: initial lambda; for \code{\link{m_gam}}: smoothing parameter estimation method
#' as required by \code{\link[mgcv]{gam}}
#' @param K for \code{\link{m}}: quadratic penalty matrix; for \code{\link{m_gam}}: preceeding gam model with \code{fit=FALSE}
#' @param ... further arguments which can be passed to \code{\link[mgcv]{gam}}
#' @inheritParams m
#' @return \describe{
#' \item{\code{ciquantiles}}{bootstrap confidence intervals}
#' \item{\code{ciraw}}{estimated coefficients or fitted values for each bootstrap sample}
#' }
#' @import stats
#' @importFrom pbmcapply pbmclapply pbmcmapply
#' @importFrom purrr map
#' @importFrom mgcv gam gam.control
#' @export

bootstrap <- function(
  B,
  WB = FALSE,
  fu = function(x) x,
  type = "gumbel",
  B_fit = FALSE,
  alpha = NULL,
  x,
  fitted,
  residuals,
  start,
  loss,
  d = deriv_loss(loss),
  tuning,
  maxi,
  tau,
  lambda,
  K,
  ncores = parallel::detectCores() - 1,
  ...
) {
  # definitions
  n <- length(fitted)
  ciquantiles <- if (is.null(alpha)) c(.025, .975) else c(alpha / 2, 1 - alpha / 2)

  # generate bootstrap samples
  samples <- bootstrap_samples(B = B, WB = WB, fu = fu, type = type, x = x, fitted = fitted, residuals = residuals)

  if (length(lambda) == 1 && !is.character(lambda)) {
    # model unpenalized or penalized with fixed value of lambda, m will be applied
    print(paste("Fitting of the according models for the", B, "bootstrap samples using", ncores, "cores:"), quote = FALSE)
    fits <- pbmclapply(samples, m_fit, x = x, start = start, loss = loss, d = d, tuning = tuning,
                       maxi = maxi, tau = tau, K = K, lambda = lambda, mc.cores = ncores)
    error <- error_fits <- rep(FALSE, B)
  } else if (!is.null(x)) {
    # model penalized, lambda chosen by cross-validation via m_cv, m will be applied
    print(paste("Cross-validation for the", B, "bootstrap samples using", ncores, "cores:"), quote = FALSE)
    cv <- pbmclapply(samples, try_parallel(m_cv), x = x, start = start, loss = loss, d = d, tuning = tuning,
                     maxi = maxi, tau = tau, K = K, lambda = lambda, mc.cores = ncores)
    if ("warning" %in% names(fits)) cv <- cv$value
    error_cv <- unlist(lapply(cv, inherits, 'error'))
    print(paste("Fitting of the according models (with lambda chosen by cross-validation) for the", B, "bootstrap samples using", ncores, "cores:"), quote = FALSE)
    fits <- pbmcmapply(try_parallel(m_fit), y = samples[!error_cv], lambda = map(cv, "lambda_opt")[!error_cv],
                       MoreArgs = list(x = x, start = start, loss = loss, d = d, tuning = tuning,
                                       maxi = maxi, tau = tau, K = K), mc.cores = ncores, SIMPLIFY = FALSE)
    if ("warning" %in% names(fits)) fits <- fits$value
    error_fits <- unlist(lapply(fits, inherits, 'error'))
    error <- error_cv
  } else {
    # model penalized, lambda chosen inside m_fit, m_gam will be applied
    m_gam_fit_y <- function(z, G, ...) {G$y <- z; return(m_fit(K = G, ...))}
    print(paste("Fitting of the according models for the", B, "bootstrap samples using", ncores, "cores:"), quote = FALSE)
    fits <- pbmclapply(samples, m_gam_fit_y, G = K, x = NULL, y = NULL, start = start, loss = loss, d = d,
                       tuning = tuning, maxi = maxi, tau = tau, lambda = lambda, mc.cores = ncores)
    error <- error_fits <- rep(FALSE, B)
  }

  conv <- unlist(map(fits, "converged"))
  nconv <- sum(as.numeric(conv))

  error[which(error == FALSE)] <- error_fits_conv <- error_fits | !conv

  # save estimated coefficients and fitted values for every bootstrap iteration
  coefs_bootstrap <- matrix(NA, nrow = length(start), ncol = B)
  fitted_bootstrap <- matrix(NA, nrow = n, ncol = B)
  coefs_bootstrap[, !error] <- matrix(unlist(map(fits, "coefficients")), nrow = length(start), byrow = FALSE)[, !error_fits_conv]
  fitted_bootstrap[, !error] <- matrix(unlist(map(fits, "fitted")), nrow = n, byrow = FALSE)[, !error_fits_conv]

  # quantiles of the empirical distribution of the bootstrap estimates
  ciquantiles <- t(apply(coefs_bootstrap, 1, quantile, probs = ciquantiles, na.rm = TRUE, type = 7))
  rownames(ciquantiles) <- colnames(x)

  ciraw <- if (B_fit == FALSE) coefs_bootstrap else fitted_bootstrap

  if (nconv != B) warning(paste("The IRLS algorithm did not converge for", B - nconv, "of", B, "bootstrap iterations."))

  return(list(ci = ciquantiles, ciraw = ciraw))
}
