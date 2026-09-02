#' Function performing the cross-validation
#'
#' This function performs the cross-validation for models fitted with \code{\link{m}} and \code{length(lambda)>0}.
#' The cross-validation is done outside of the IRLS algorithm.
#' If \code{length(lambda)==2}, lambda is cross-validated on a log-scale between those values,
#' otherwise the values in \code{lambda} are cross-validated.
#' The GCV scores and optimal lambdas are calculated with \code{\link{m_cv_calc}}.
#'
#' @param lambda vector of smoothing parameters which will be compared via cross-validation;
#' if \code{length(lambda)==2}, lambda is cross-validated on a log-scale between those values;
#' if \code{length(lambda>2)}, the values in \code{lambda} are cross-validated
#' @param ncores number of cores to be used for parallelization
#' @inheritParams m_fit
#' @return \describe{
#' \item{\code{lambda_opt}}{optimal lambda chosen by cross-validation via GCV criterion}
#' \item{\code{cv}}{list containing the GCV scores and the estimated coefficients for the different cross-validated smoothing parameter values}
#' }
#' @export

m_cv <- function(
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
  if (length(lambda) > 2) {
    cv <- m_cv_calc(x = x, y = y, start = start, loss = loss, d = d, tuning = tuning, maxi = maxi, tau = tau,
                    lambda = lambda, K = K, ncores = ncores)
    lambda_opt <- cv$lambda_opt
    scores <- cv$scores
    coefs <- cv$coefs
  } else if (length(lambda) == 2) {
    lambda <- sort(lambda)
    basis <- (lambda[2] - lambda[1])^(1 / 4)
    acc <- .01
    lambda_temp <- if (basis > 1) {round(lambda[1] + c(.1, basis^(0:4)), digits = 2)} else {seq(lambda[1], lambda[2], by = acc)}
    coefs <- list()
    scores <- list()
    lambda_list <- list()
    i <- 1
    while (max(lambda_temp) - min(lambda_temp) >= acc && i <= 10) {
      lambda_temp <- unique(lambda_temp[which(lambda_temp <= lambda[2])])
      if (any(lambda_temp %in% unlist(lambda_list))) lambda_temp <- lambda_temp[-which(lambda_temp %in% unlist(lambda_list))]
      if (length(lambda_temp) > 0) {
        lambda_list[[i]] <- lambda_temp
        cv <- m_cv_calc(x = x, y = y, start = start, loss = loss, d = d, tuning = tuning, maxi = maxi, tau = tau,
                        lambda = lambda_temp, K = K, ncores = ncores)
        scores[[i]] <- cv$scores
        coefs[[i]] <- cv$coefs
        lambda_opt <- max(unlist(lambda_list)[which(unlist(scores) == min(unlist(scores)))])
      }
      j <- log(lambda_opt - lambda[1]) / log(basis)
      lambda_temp <- round(lambda[1] + basis^(c(j - 2^(-i), j + 2^(-i))), digits = 2)
      i <- i + 1
    }
    if (lambda[1] == 0) {
      lambda_temp <- round(seq(0, .5, length.out = 10), digits = 2)
      if (any(lambda_temp %in% unlist(lambda_list))) lambda_temp <- lambda_temp[-which(lambda_temp %in% unlist(lambda_list))]
      if (length(lambda_temp) > 0) {
        lambda_list[[length(lambda_list) + 1]] <- lambda_temp
        cv <- m_cv_calc(x = x, y = y, start = start, loss = loss, d = d, tuning = tuning, maxi = maxi, tau = tau,
                        lambda = lambda_temp, K = K, ncores = ncores)
        scores[[length(lambda_list) + 1]] <- cv$scores
        coefs[[length(lambda_list) + 1]] <- cv$coefs
        lambda_opt <- max(unlist(lambda_list)[which(unlist(scores) == min(unlist(scores)))])
      }
    }
  } else {
    stop("The length of the vector lambda is smaller or equal to 1. Hence, the values can not be cross-validated.")
  }
  coefs <- do.call(cbind, as.list(coefs))
  scores <- unlist(scores)
  coefs <- coefs[, order(as.numeric(colnames(coefs)))]
  scores <- scores[order(as.numeric(names(scores)))]
  return(list(lambda_opt = lambda_opt, cv = list(scores = scores, coefs = coefs)))
}
