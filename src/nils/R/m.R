#' Function performing the (penalized) mode, mean and median regression with lambda chosen outside the IRLS algorithm
#'
#' This function performs the (possibly penalized) mode, median and mean regression where the penalization parameter lambda
#' can be a given value or can be chosen outside the IRLS algorithm via cross-validation.
#'
#' @param x design matrix
#' @param y response vector
#' @param start initial values of coefficients in IRLS algorithm
#' @param loss type of loss function: \code{"nils"} = approach of Oelker et al. (2015),
#' \code{"kemp"} = approach of Kemp et al. (2012), \code{"mean"} = mean regression, \code{"med"} = approximation of median
#' regression
#' @param tuning set of tuning parameters, if \code{NULL}, the tuning parameters are not set manually
#' @param maxi maximal number of iterations in the IRLS algorithm
#' @param tau criterion for convergence of IRLS algorithm
#' @param lambda vector of smoothing parameters; if \code{length(lambda>1)} they will be compared via cross-validation;
#' if \code{length(lambda)==2}, lambda is cross-validated on a log-scale between those values;
#' if \code{length(lambda>2)}, the values in \code{lambda} are cross-validated
#' @param K quadratic penalty matrix
#' @param B number of bootstrap samples; if \code{B=0} no bootstrap will be performed
#' @param WB logical; if \code{TRUE}, wild bootstrap instead of residual bootstrap is performed
#' @param fu transformation function for the residuals
#' @param type type of weight distribution, currently seven options are implemented:
#' \code{"gumbel"}, \code{"alaplace"}, \code{"mammen"}, \code{"normal"}, \code{"rademacher"}, \code{"mammen_two"}, \code{"webb_six"}
#' @param B_fit logical; if \code{FALSE}, the estimated coefficients for each bootstrap sample are returned in \code{bs_ciraw}; if \code{TRUE}, the fitted values per sample are returned
#' @param alpha confidence level for bootstrap confidence intervals; if \code{alpha=NULL} the confidence level is set to 0.05
#' @param cov logical; if \code{TRUE} and the type of loss function is either \code{"nils"} or \code{"kemp"},
#' the estimated asymptotic covariance matrix multiplied by k^3/n will be calculated
#' @param ncores number of cores to be used for parallelization
#' @return \describe{
#' \item{\code{coefficients}}{estimated coefficients}
#' \item{\code{residuals}}{residuals}
#' \item{\code{fitted}}{fitted values}
#' \item{\code{weights}}{weights of final iteration in IRLS algorithm}
#' \item{\code{iter}}{number of iterations in IRLS algorithm}
#' \item{\code{converged}}{logical; if \code{TRUE}, the IRLS algorithm converged}
#' \item{\code{tuning}}{set of tuning parameters}
#' \item{\code{K}}{penalty matrix}
#' \item{\code{lambda}}{penalization parameter used for fitting (might be determined via cross-validation)}
#' \item{\code{cv}}{list containing the GCV scores and the estimated coefficients for the different cross-validated smoothing parameter values}
#' \item{\code{bs_ci}}{bootstrap confidence intervals}
#' \item{\code{bs_ciraw}}{estimated coefficients or fitted values for each bootstrap sample}
#' \item{\code{Omegahat}}{estimated asymptotic covariance matrix multiplied by \code{k^3/n}}
#' \item{\code{k}}{for \code{loss="nils"}: tuning parameter in last IRLS iteration, for \code{loss="kemp"}: inverse bandwith in last IRLS iteration}
#' }
#' @export

m <- function(
  x,
  y,
  start = NULL,
  loss = "nils",
  tuning = NULL,
  maxi = 5500,
  tau = 10^(-5),
  lambda = 0,
  K = NULL,
  B = 0,
  WB = FALSE,
  fu = function(x) x,
  type = "gumbel",
  B_fit = FALSE,
  alpha = NULL,
  cov = FALSE,
  ncores = parallel::detectCores() - 1
) {
  cv <- list(lambda_opt = lambda, cv = NULL)
  bs <- list(ci = NULL, ciraw = NULL)
  cov_hat <- list(Omegahat = NULL, k = NULL)
  n <- nrow(x)

  # get derivative of chosen loss function divided by its argument as function
  d <- deriv_loss(loss)

  # define the tuning parameters and the starting values
  control <- m_control(x = x, y = y, loss = loss, tuning = tuning, maxi = maxi, start = start, lambda = lambda, K = K)

  # perform the cross-validation if length(lambda)>1
  if (length(lambda) > 1) cv <- m_cv(x = x, y = y, start = control$start, loss = loss, d = d, tuning = control$tuning, maxi = control$maxi,
                                     tau = tau, lambda = lambda, K = control$K, ncores = ncores)

  # fit the according model with lambda=cv$lambda_opt
  fit <- m_fit(x = x, y = y, start = control$start, loss = loss, d = d, tuning = control$tuning, maxi = control$maxi,
               tau = tau, lambda = cv$lambda_opt, K = control$K)

  # estimate the bootstrap confidence intervals if the number B of bootstrap samples is larger than 0
  if (B > 0) bs <- bootstrap(B = B, WB = WB, fu = fu, type = type, B_fit = B_fit, alpha = alpha, x = x, fitted = fit$fitted, residuals = fit$residuals,
                             start = control$start, loss = loss, d = d, tuning = control$tuning, maxi = control$maxi,
                             tau = tau, lambda = lambda, K = control$K, ncores = ncores)

  # if cov==TRUE, calculate the estimated asymptotic covariance matrix multiplied by k^3/n and
  # the tuning parameter / inverse bandwith in the last IRLS iteration
  # (needed for calculating the estimated asymptotic confidence intervals)
  if (cov == TRUE) {cov_hat <- est_cov(x = x, residuals = fit$residuals, n = n, loss = loss, tuning = control$tuning, j = fit$j)}
  # if cov==FALSE, calculate the tuning parameter / inverse bandwith in the last IRLS iteration
  else if (loss == "nils") {cov_hat$k <- control$tuning$k[fit$j]}
  else if (loss == "kemp") {cov_hat$k <- (control$tuning$smallc*control$tuning$k[fit$j])^(-1)}

  return(list(
    coefficients = fit$coefficients,
    residuals = fit$residuals,
    fitted = fit$fitted,
    weights = fit$weights,
    iter = fit$iter,
    converged = fit$converged,
    tuning = control$tuning,
    K = control$K,
    lambda = cv$lambda_opt,
    cv = cv$cv,
    bs_ci = bs$ci,
    bs_ciraw = bs$ciraw,
    Omegahat = cov_hat$Omegahat,
    k = cov_hat$k
  ))
}
