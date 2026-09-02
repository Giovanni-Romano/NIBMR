#' Function returning the transformation function for wild bootstrap
#'
#' This function returns the transformation function for wild bootstrap depending on the chosen type of function.
#' The following types are available: \code{"identity"} (identity function), \code{"absolute"} (absolute value function),
#' \code{"sqrt_hat"} (function deviding the residuals by the square root of 1 minus the respective diagonal entry of the hat matrix),
#' \code{"hat"} (function deviding the residuals by 1 minus the respective diagonal entry of the hat matrix) and
#' \code{"df"} (function multiplying the residuals with \code{sqrt(n/(n-p))}.
#'
#' @param trafo_type type of transformation function to be returned, currently six options are implemented:
#' \code{"identity"}, \code{"absolute"}, \code{"sqrt_hat"}, \code{"hat"}, \code{"df"}, \code{"feng"}
#' @param n number of observations
#' @param p number of coefficients (including the intercept)
#' @param x design matrix
#' @return \describe{
#' \item{\code{f}}{transformation function for wild bootstrap}
#' }
#' @importFrom quantreg akj
#' @export

transformation_function <- function (
  trafo_type = c("identity", "absolute", "sqrt_hat", "hat", "df", "feng"),
  n = NULL,
  p = NULL,
  x = NULL
) {
  trafo_type <- match.arg(trafo_type)

  if (trafo_type == "identity") {
    f <- function (u) u
  } else if (trafo_type == "absolute") {
    f <- function (u) abs(u)
  } else if (trafo_type == "sqrt_hat") {
    chol_matrix <- chol(crossprod(x))
    hat_matrix <- x %*% chol2inv(chol_matrix) %*% t(x)
    f <- function (u) u / sqrt(1 - diag(hat_matrix))
  } else if (trafo_type == "hat") {
    chol_matrix <- chol(crossprod(x))
    hat_matrix <- x %*% chol2inv(chol_matrix) %*% t(x)
    f <- function (u) u / 1 - diag(hat_matrix)
  } else if (trafo_type == "df") {
    factor <- sqrt(n/(n - p))
    f <- function (u) factor*u
  } else if (trafo_type == "feng") {
    chol_matrix <- chol(crossprod(x))
    hat_matrix <- x %*% chol2inv(chol_matrix) %*% t(x)
    f <- function (u) {
      # fhat0 <- approxfun(density(u))(0)
      fhat0 <- akj(u, z = 0)$dens
      abs(u + (fhat0)^(-1)*diag(hat_matrix)*(0.5 - 1*(u < 0)))
    }
  } else {
    stop("You need to choose one of the following transformation function types:
         identity, absolute, sqrt_hat, hat, df, feng.")
  }

  return(f)
}
