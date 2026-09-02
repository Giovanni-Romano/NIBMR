#' Function calculating the coverage rates of the asymptotic and bootstrap confidence intervals
#'
#' This function calculates the coverage rates of the asymptotic and bootstrap confidence intervals.
#'
#' @param coefs vector of true coefficients
#' @param ci list of matrices containing the asymptotic confidence intervals as rows for the different repetitions
#' @param bs_ci list of matrices containing the bootstrap confidence intervals as rows for the different repetitions
#' @return \describe{
#' \item{\code{rates}}{matrix containing the coverage rates of the asymptotic confidence intervals and
#' of the bootstrap confidence intervals as columns}
#' }
#' @export

coverage_rates <- function(
  coefs,
  ci,
  bs_ci
) {
  ci <- ci[!is.na(ci)]
  bs_ci <- bs_ci[!is.na(bs_ci)]

  cover_fct <- function(coefs, ci) 1*(coefs <= ci[, 2] & coefs >= ci[, 1])

  cover <- lapply(ci, cover_fct, coefs = coefs)
  cover_sum <- Reduce(`+`, cover)
  cover_rate <- cover_sum / length(cover)

  bs_cover <- lapply(bs_ci, cover_fct, coefs = coefs)
  bs_cover_sum <- Reduce(`+`, bs_cover)
  bs_cover_rate <- bs_cover_sum / length(bs_cover)

  rates <- cbind(cover_rate, bs_cover_rate)

  return(rates)
}
