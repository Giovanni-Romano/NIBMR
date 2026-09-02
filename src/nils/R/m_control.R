#' Function for defining the tuning parameters and the starting values
#'
#' This function defines the tuning parameters for the chosen loss function.
#' For the loss functions \code{"nils"} and \code{"med"} the input tuning parameters are preserved if the input parameter \code{tuning} is not \code{NULL}.
#' For \code{"kemp"} the tuning parameters are chosen as described in Kemp et al. (2012),
#' for \code{"mean"} they are chosen in a way that the IRLS algorithm yields the least-squares estimate.
#' Hence, for all four loss functions the same algorithm can be used.
#'
#' @param lambda penalty parameter
#' @inheritParams m
#' @return \describe{
#' \item{\code{tuning}}{set of tuning parameters}
#' \item{\code{maxi}}{maximal number of iterations in the IRLS algorithm}
#' \item{\code{start}}{initial values of coefficients in IRLS algorithm}
#' \item{\code{K}}{penalty matrix}
#' }
#' @import stats
#' @export

m_control <- function(
  x,
  y,
  loss,
  tuning,
  maxi,
  start,
  lambda,
  K
) {
  n <- nrow(x)
  if (any(is.null(start))) start <- rep(1, ncol(x))
  if (length(lambda) > 1) lambda <- lambda[1]
  if (any(is.null(K))) K <- diag(1, ncol(x))

  if (is.null(tuning) && loss == "nils") {
    # tuning parameters for mode regression via NILS
    tuning <- list(
      k = NULL,
      g = c(10:2, rep(1, maxi - 9)),
      smallc = 10^(-5),
      steps = rep(10, maxi),
      nu = rep(.25, maxi)
    )
    # adaptive tuning of parameters
    # fitted median regression
    med_tuning <- m_control(x = x, y = y, loss = "med", tuning = NULL, maxi = maxi, start = start, lambda = lambda, K = K)
    med_fit <- m_fit(x = x, y = y, start = start, loss = "med", d = deriv_loss("med"), tuning = med_tuning$tuning,
                     maxi = 500, tau = 10^(-5), K = K, lambda = lambda)
    # initial k
    q <- n^(-1/6)*log(n)^(7/12)
    k_initial <- (n/q)^(1/7)
    # final k
    k_final <- (k_initial*n^(1/7)) / sd(med_fit$residuals)
    # initial k is increased up to the final k in ten steps
    tuning$k <- c(seq(k_initial, k_final, length.out = 10), rep(k_final, maxi - 10))
    start <- med_fit$coefficients # median regression coefficients are employed as starting values
  } else if (loss == "kemp") {
    # tuning parameters for mode regression as described in Kemp et al. (2012)
    ols <- lsfit(x = x, y = y, intercept = FALSE)
    # if the factor k for the bandwith is not passed by the input variable tuning, choose 1.2
    if (is.null(tuning)) {k <- 1.2} else {k <- tuning$k}
    tuning <- list(
      k = rep(k, maxi),
      g = NULL,
      smallc = median(abs(ols$residuals - median(ols$residuals)))*n^(-0.143),
      steps = rep(1, maxi),
      nu = rep(1, maxi)
    )
    start <- ols$coef # ols coefficients are employed as starting values
  } else if (loss == "mean") {
    # tuning parameters for mean regression
    # performing only one iteration in the NILS algorithm with nu=1 and identity matrix as weights
    # yields the least squares estimator
    maxi <- 1
    tuning <- list(
      k = NULL,
      g = NULL,
      smallc = NULL,
      steps = 1,
      nu = 1
    )
  } else if (is.null(tuning) && loss == "med") {
    # tuning parameters for median regression via IRLS algorithm
    tuning <- list(
      k = NULL,
      g = NULL,
      smallc = 10^(-5),
      steps = rep(1, maxi),
      nu = rep(1, maxi)
    )
  } else {
    if (is.null(tuning)) stop("You need to choose one of the following losses: nils, kemp, mean, med.")
  }
  return(list(tuning = tuning, maxi = maxi, start = start, K = K))
}
