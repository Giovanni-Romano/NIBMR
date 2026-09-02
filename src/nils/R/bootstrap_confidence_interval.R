#' Function calculating the bootstrap confidence intervals
#'
#' This function calculates the bootstrap confidence intervals. Different construction methods can be chosen:
#' bootstrap reverse percentile, equal-tailed bootstrap percentile-t and symmetric bootstrap percentile-t.
#' The percentile confidence interval is already calculated in the function \code{\link{bootstrap}}.
#'
#' @param n number of observations
#' @param est_coefs estimated coefficients
#' @param bs_coefs estimated coefficients for each bootstrap sample
#' @param bs_ci_method bootstrap confidence interval construction method:
#' \code{"rev_percentile"} = bootstrap reverse percentile, \code{"et_percentile-t"} = equal-tailed bootstrap percentile-t,
#' \code{"sym_percentile-t"} = symmetric bootstrap percentile-t
#' @param B_nested for \code{bs_ci_method = "et_percentile-t"}: number of nested bootstrap iterations per original bootstrap iteration
#' @return \describe{
#' \item{\code{bootstrap_ci}}{matrix containing the bootstrap confidence intervals for the different coefficients as rows}
#' }
#' @inheritParams m
# #' @importFrom coxed bca
#' @export

bootstrap_confidence_interval <- function (
  n,
  alpha = NULL,
  est_coefs,
  bs_coefs,
  bs_ci_method = c("rev_percentile", "bca", "et_percentile-t", "sym_percentile-t"),
  B_nested = B,
  WB = FALSE,
  fu = function(x) x,
  type = "gumbel",
  x = NULL,
  y = NULL,
  start = NULL,
  loss = "nils",
  tuning = NULL,
  maxi = 5500,
  tau = 10^(-5),
  lambda = 0,
  K = NULL
) {
  # bs_ci_method <- match.arg(bs_ci_method)

  ciquantiles <- if (is.null(alpha)) c(.025, .975) else c(alpha / 2, 1 - alpha / 2)

  if ("rev_percentile" %in% bs_ci_method) {
    # ciquantiles <- t(apply(sqrt(n)*(bs_coefs - as.vector(est_coefs)), 1, quantile, probs = ciquantiles, na.rm = TRUE, type = 7))
    # bootstrap_ci <- cbind(est_coefs - n^(-0.5)*ciquantiles[, 2], est_coefs - n^(-0.5)*ciquantiles[, 1])
    ciquantiles <- t(apply(bs_coefs, 1, quantile, probs = ciquantiles, na.rm = TRUE, type = 7))
    bootstrap_ci <- cbind(2*est_coefs - ciquantiles[, 2], 2*est_coefs - ciquantiles[, 1])
  } # else if (bs_ci_method == "bca") {
    # bootstrap_ci <- t(apply(bs_coefs, 1, bca, conf.level = 1 - alpha))
  # }
  else if (sum(bs_ci_method %in% c("et_percentile-t", "sym_percentile-t")) >= 1) {
    B <- ncol(bs_coefs)
    bs_t_stat <- matrix(NA, length(est_coefs), B)

    for (i in 1:B) {
      fitted <- x %*% bs_coefs[, i]
      residuals <- y - fitted
      control <- m_control(x = x, y = y, loss = loss, tuning = tuning, maxi = maxi, start = start, lambda = lambda, K = K)
      bs <- bootstrap(B = B_nested, WB = WB, fu = fu, type = type, B_fit = FALSE, alpha = alpha, x = x, fitted = fitted, residuals = residuals,
                      start = control$start, loss = loss, d = deriv_loss(loss), tuning = control$tuning, maxi = control$maxi,
                      tau = tau, lambda = lambda, K = control$K, ncores = 1)
      se <- apply(bs$ciraw, 1, sd)
      bs_t_stat[, i] <- (bs_coefs[, i] - est_coefs)/se
    }

    se <- apply(bs_coefs, 1, sd)

    bootstrap_ci <- as.list(rep(NA, 2))
    names(bootstrap_ci) <- c("et_percentile-t", "sym_percentile-t")

    if ("et_percentile-t" %in% bs_ci_method) {
      ciquantiles_et <- t(apply(bs_t_stat, 1, quantile, probs = ciquantiles, na.rm = TRUE, type = 7))
      bootstrap_ci[[1]] <- cbind(est_coefs - se*ciquantiles_et[, 2], est_coefs - se*ciquantiles_et[, 1])
    }

    if ("sym_percentile-t" %in% bs_ci_method) {
      ciquantiles_sym <- t(apply(abs(bs_t_stat), 1, quantile, probs = 1 - alpha, na.rm = TRUE, type = 7))
      ciquantiles_sym <- as.vector(ciquantiles_sym)
      bootstrap_ci[[2]] <- cbind(est_coefs - se*ciquantiles_sym, est_coefs + se*ciquantiles_sym)
    }

    if (sum(is.na(bootstrap_ci)) == 1) bootstrap_ci <- bootstrap_ci[!is.na(bootstrap_ci)][[1]]
  } else {
    stop(strwrap("You need to choose one of the following
         bootstrap confidence interval construction methods:
         rev_percentile (for bootstrap reverse percentile), et_percentile-t (for equal-tailed bootstrap percentile-t),
         sym_percentile-t (for symmetric bootstrap percentile-t).", prefix = " ", initial = ""))
  }

  return(bootstrap_ci)
}
