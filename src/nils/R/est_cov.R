#' Function for estimating the asymptotic covariance matrix
#'
#' This function calculates the asymptotic covariance matrix but returns the matrix
#' multiplied with the factor (k_n^3/n) to simplify the calculation of the asymptotic confidence intervals.
#'
#' @param x design matrix
#' @param residuals residuals
#' @param n number of observations
#' @param loss type of loss function: \code{"nils"} = approach of Oelker et al. (2015),
#' \code{"kemp"} = approach of Kemp et al. (2012)
#' @param tuning set of tuning parameters
#' @param j internal index for tuning the IRLS algorithm
#' @param correction logical; if TRUE, the finite sample corrected residuals are returned
#' @return \describe{
#' \item{\code{Omegahat}}{estimated asymptotic covariance matrix multiplied by \code{k^3/n}}
#' \item{\code{k}}{for \code{loss="nils"}: tuning parameter in last IRLS iteration, for \code{loss="kemp"}: inverse bandwith in last IRLS iteration}
#' }
#' @export

est_cov <- function(
  x,
  residuals,
  n,
  loss,
  tuning,
  j,
  correction = FALSE
) {
  if (loss == "nils") {
    k <- tuning$k[j]
    dK_nils <- function(u, smallc = tuning$smallc) {-0.5 * exp(-sqrt(u^2 + smallc)) * (u/sqrt(u^2 + smallc))}
    d2K_nils <- function(u, smallc = tuning$smallc) {0.5 * exp(-sqrt(u^2 + smallc)) * (u^2/(u^2 + smallc) - smallc / (u^2 + smallc)^1.5)}
  } else if (loss == "kemp") {
    delta <- tuning$smallc * tuning$k[j]
    k <- 1 / delta
    dK_kemp <- function(u) {-1/sqrt(2*pi)*u*exp(-u^2/2)}
    d2K_kemp <- function(u) {1/sqrt(2*pi)*(u^2-1)*exp(-u^2/2)}
  } else {
    stop(paste("The asymptotic covariance matrix can not be computed for the loss", loss, "."))
  }

  dK <- get(paste("dK_", loss, sep = ""))
  d2K <- get(paste("d2K_", loss, sep = ""))
  Bhat <- k / n * t(x) %*% diag(dK(k * as.vector(residuals))^2) %*% x
  Chat <- k^3 / n * t(x) %*% diag(d2K(k * as.vector(residuals))) %*% x
  Chatinv <- solve(Chat)
  Omegahat <- Chatinv %*% Bhat %*% Chatinv * (k^3 / n)
  if (correction) Omegahat <- residuals - k^2/n * x %*% Chatinv %*% colSums(dK(k * as.vector(residuals)) * x)

  return(list(Omegahat = Omegahat, k = k))
}
