#' Function needed for calculating the estimates in the IRLS algorithm
#'
#' The function calculates the derivative of the chosen loss function divided by its argument, constant factors are omitted.
#' This is needed to calculate the estimates of the coefficients in every iteration of the IRLS algorithm.
#'
#' @inheritParams m
#' @return \describe{
#' \item{\code{d}}{derivative of chosen loss function divided by its argument as function}
#' }
#' @export

deriv_loss <- function(
  loss
) {
  # derivatives of loss functions divided by its argument, constant factors are omitted
  if (loss == "nils") {
    d <- function(xi, k, g, smallc) {
      g <- 2 * g
      as.vector(k^g * ((k*xi)^g + smallc)^(1/g-1) * xi^(g-2) * exp(smallc^(1/g) - ((k*xi)^g + smallc)^(1/g)) + 1e-15)
    }
  } else if (loss == "kemp") {
    d <- function(xi, k, g, smallc) {as.vector(exp(-0.5 * (xi / (smallc*k))^2) + 1e-15)}
  } else if (loss == "mean") {
    d <- function(xi, k, g, smallc) {as.vector(rep(1, length(xi)))}
  } else if (loss == "med") {
    d <- function(xi, k, g, smallc) {as.vector((xi^2 + smallc)^(-1/2))}
  } else {
    stop("You need to choose one of the following losses: nils, kemp, mean, med.")
  }
  return(d)
}
