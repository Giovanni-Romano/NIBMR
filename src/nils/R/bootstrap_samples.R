#' Function for generating the bootstrap samples
#'
#' This function generates the bootstrap samples via residual bootstrap.
#' If \code{WB = TRUE}, wild bootstrap instead of residual bootstrap is performed.
#'
#' @param x design matrix
#' @inheritParams bootstrap
#' @return \describe{
#' \item{\code{samples}}{bootstrap samples}
#' }
#' @export

bootstrap_samples <- function(
  B,
  WB = FALSE,
  fu = function(x) x,
  type = "gumbel",
  x,
  fitted,
  residuals
) {
  n <- length(fitted)
  samples <- list()
  if (WB == TRUE) {
    # wild bootstrap
    for (i in 1:B) {
      samples[[i]] <- fitted + trafo_WB(residuals, fu = fu, type = type)
    }
  } else {
    # residual bootstrap
    for (i in 1:B) {
      indices <- sample(1:n, n, replace = TRUE, prob = NULL)
      samples[[i]] <- fitted + residuals[indices]
    }
  }
  return(samples)
}
