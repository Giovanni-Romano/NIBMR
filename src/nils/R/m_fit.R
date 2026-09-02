#' Function fitting the model
#'
#' This function fits the model by performing the IRLS algorithm proposed by Oelker et al. (2015).
#' @param x for \code{\link{m}}: design matrix, for \code{m_gam()}: \code{NULL}
#' @param y for \code{\link{m}}: response vector, for \code{m_gam()}: \code{NULL}
#' @param d derivative of chosen loss function divided by its argument
#' @param tuning set of tuning parameters
#' @param lambda for \code{\link{m}}: smoothing parameter, for \code{\link{m_gam}}: smoothing parameter estimation method
#' as required by \code{\link[mgcv]{gam}}
#' @param K for \code{\link{m}}: quadratic penalty matrix; for \code{\link{m_gam}}: preceeding gam model with \code{fit=FALSE}
#' @inheritParams m
#' @return \describe{
#' \item{\code{coefficients}}{estimated coefficients}
#' \item{\code{residuals}}{residuals}
#' \item{\code{fitted}}{fitted values}
#' \item{\code{weights}}{weights of final iteration in IRLS algorithm}
#' \item{\code{iter}}{number of iterations}
#' \item{\code{j}}{internal index of tuning parameters used in last IRLS iteration}
#' \item{\code{converged}}{logical; if \code{TRUE}, the IRLS algorithm converged}
#' \item{\code{gam_object}}{for \code{\link{m}}: \code{NULL}, for \code{\link{m_gam}}: gam object in final iteration}
#' }
#' @importFrom mgcv gam gam.control
#' @export

m_fit <- function(
  x = NULL,
  y = NULL,
  start,
  loss,
  d,
  tuning,
  maxi,
  tau,
  lambda,
  K
) {
  conv <- FALSE
  coefs <- start
  gam_object <- NULL

  if (is.null(x) && is.null(y)) {
    gam_boolean <- TRUE
    x <- K$X
    y <- K$y
  } else {
    gam_boolean <- FALSE
  }

  res <- y - x %*% coefs

  for (i in 1:maxi) {
    j <- which(cumsum(tuning$steps) >= i)[1]
    coefs_old <- coefs
    D <- d(xi = res, k = tuning$k[j], g = tuning$g[j], smallc = tuning$smallc)
    if (gam_boolean) {
      K$mf["(weights)"] <- K$w <- D
      gam_object <- gam(G = K, gamma = 1, scale = 0, control = gam.control(epsilon = tau),
                        method = lambda, optimizer = c("outer","newton"), fit = TRUE)
      update <- gam_object$coefficients
    } else {
      a <- crossprod(x, D*y)
      A <- crossprod(sqrt(D)*x)
      chol_matrix <- chol(A + lambda*K)
      update <- drop(chol2inv(chol_matrix)) %*% a
    }
    coefs <- (1 - tuning$nu[j])*coefs_old + tuning$nu[j]*update
    if (any(is.na(coefs))) {
      fitted <- x %*% coefs_old
      coefs <- coefs_old
      break
    }
    fitted <- x %*% coefs
    res <- y - fitted
    conv_crit <- sum(abs(coefs - coefs_old)) / sum(abs(coefs_old))
    if (conv_crit <= tau) {
      conv <- TRUE
      D <- d(xi = res, k = tuning$k[j], g = tuning$g[j], smallc = tuning$smallc)
      break
    }
  }

  if (loss == "mean") conv <- TRUE
  if (conv == FALSE) warning("The algorithm did not converge.")
  return(list(coefficients = coefs, residuals = res, fitted = fitted, weights = D, iter = i, j = j,
              converged = conv, gam_object = gam_object))
}
