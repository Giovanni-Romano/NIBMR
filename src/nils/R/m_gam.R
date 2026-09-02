#' Function performing the penalized mode, mean and median regression with lambda chosen inside the IRLS algorithm
#'
#' This function performs the penalized mode, median and mean regression where the penalization parameter lambda
#' is chosen inside the IRLS algorithm and the estimation of the coefficent vector and of the penalty parameter is interlaced.
#'
#' @param formula model formula as required by \code{\link[mgcv]{gam}}
#' @param data data set
#' @param method smoothing parameter estimation method as required by \code{\link[mgcv]{gam}}
#' @param ... further arguments which can be passed to \code{\link[mgcv]{gam}}
#' @inheritParams m
#' @return \describe{
#' \item{\code{coefficients}}{estimated coefficients}
#' \item{\code{residuals}}{residuals}
#' \item{\code{fitted}}{fitted values}
#' \item{\code{weights}}{weights of final iteration in IRLS algorithm}
#' \item{\code{iter}}{number of iterations in IRLS algorithm}
#' \item{\code{converged}}{logical; if \code{TRUE}, the IRLS algorithm converged}
#' \item{\code{tuning}}{set of tuning parameters}
#' \item{\code{lambda}}{estimated smoothing parameters for the gam model in the last IRLS iteration}
#' \item{\code{gam_object}}{gam object in the last IRLS iteration}
#' \item{\code{bs_ci}}{bootstrap confidence intervals}
#' \item{\code{bs_ciraw}}{estimated coefficients or fitted values for each bootstrap sample}
#' \item{\code{k}}{for \code{loss="nils"}: tuning parameter in last IRLS iteration, for \code{loss="kemp"}: inverse bandwith in last IRLS iteration}
#' \item{\code{x}}{design matrix (generated internally in \code{\link[mgcv]{gam}})}
#' }
#' @importFrom mgcv gam gam.control
#' @export

m_gam <- function(
  formula,
  data = list(),
  start = NULL,
  loss = "nils",
  tuning = NULL,
  maxi = 5500,
  tau = 10^(-5),
  B = 0,
  WB = FALSE,
  fu = function(x) x,
  type = "gumbel",
  B_fit = FALSE,
  alpha = NULL,
  method = "REML",
  ncores = parallel::detectCores() - 1,
  ...
) {
  bs <- list(ci = NULL, ciraw = NULL)
  k <- NA
  G <- gam(formula = formula, family = gaussian(), data = data, na.action = na.omit, offset = NULL, method = method,
           optimizer = c("outer","newton"), control = gam.control(epsilon = tau), scale = 0, select = FALSE,
           fit = FALSE, G = NULL, ...)

  # get derivative of chosen loss function divided by its argument as function
  d <- deriv_loss(loss)

  # define the tuning parameters and the starting values
  n <- nrow(G$X)
  control <- m_control(x = G$X, y = G$y, loss = loss, tuning = tuning, maxi = maxi, start = start, lambda = .5, K = diag(ncol(G$X)))

  # fit the according model
  fit <- m_fit(x = NULL, y = NULL, start = control$start, loss = loss, d = d, tuning = control$tuning, maxi = control$maxi,
               tau = tau, lambda = method, K = G)

  # estimate the bootstrap confidence intervals if the number B of bootstrap samples is larger than 0
  if (B > 0) bs <- bootstrap(B = B, WB = WB, fu = fu, type = type, B_fit = B_fit, alpha = alpha, x = NULL, fitted = fit$fitted,
                             residuals = fit$residuals, start = control$start, loss = loss, d = d, tuning = control$tuning,
                             maxi = control$maxi, tau = tau, lambda = method, K = G, ncores = ncores, ...)

  # calculate the tuning parameter / inverse bandwith in the last IRLS iteration
  if (loss == "nils") {k <- control$tuning$k[fit$j]}
  else if (loss == "kemp") {k <- (control$tuning$smallc*control$tuning$k[fit$j])^(-1)}

  return(list(
    coefficients = fit$coefficients,
    residuals = fit$residuals,
    fitted = fit$fitted,
    weights = fit$weights,
    iter = fit$iter,
    converged = fit$converged,
    tuning = control$tuning,
    lambda = fit$gam_object$sp,
    gam_object = fit$gam_object,
    bs_ci = bs$ci,
    bs_ciraw = bs$ciraw,
    k = k,
    x = G$X
  ))
}
