#' Function for calculating the GCV scores and finding the optimal lambda
#'
#' This function takes a vector \code{lambda} of possible values for lambda as input and fits the according model for each lambda.
#' Finally, it calculates the GCV scores for each lambda from the fits and detects the optimal lambda
#' (meaning the lambda which yields the minimal GCV score). If the minimum is not unique, the maximal optimal lambda is chosen.
#' @param lambda vector of possible smoothing parameters
#' @inheritParams m_cv
#' @return \describe{
#' \item{\code{lambda_opt}}{optimal lambda chosen by cross-validation via GCV criterion}
#' \item{\code{lambda}}{vector of smoothing parameters which will be compared via cross-validation}
#' \item{\code{scores}}{GCV scores for the different values in \code{lambda}}
#' \item{\code{coefs}}{estimated coefficients for the different values in \code{lambda}}
#' }
#' @importFrom pbmcapply pbmclapply
#' @importFrom purrr map
#' @export

m_cv_calc <- function(
  x,
  y,
  start,
  loss,
  d,
  tuning,
  maxi,
  tau,
  lambda,
  K,
  ncores = parallel::detectCores() - 1
) {
  gcv_score <- function(x, y, weights, lambda, K, coefs) {
    weighted_x <- sqrt(weights)*x
    chol_matrix <- chol(crossprod(weighted_x) + lambda*K)
    hat_matrix <- weighted_x %*% chol2inv(chol_matrix) %*% t(weighted_x)
    return(nrow(x) * t(weights) %*% as.vector(y - x %*% coefs)^2 / (nrow(x) - sum(diag(hat_matrix)))^2)
  }

  print(paste("Fitting of the according models for the current smoothing parameter values lambda =", paste(lambda, collapse = ", "), "using", ncores, "cores:"), quote = FALSE)
  fits <- pbmclapply(lambda, m_fit, x = x, y = y, start = start, loss = loss, d = d, tuning = tuning,
                     maxi = maxi, tau = tau, K = K, mc.cores = ncores)
  coefs <- matrix(unlist(map(fits, "coefficients")), nrow = length(start), byrow = FALSE)
  weights <- map(fits, "weights")
  scores <- mapply(weights = weights, lambda = lambda, coefs = coefs, gcv_score, MoreArgs = list(x = x, y = y, K = K))
  names(scores) <- colnames(coefs) <- as.character(lambda)
  lambda_opt <- max(lambda[which(scores == min(scores))])
  return(list(lambda_opt = lambda_opt, lambda = lambda, scores = scores, coefs = coefs))
}
