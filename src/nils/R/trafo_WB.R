#' Function transforming the residuals for the wild bootstrap
#'
#' This function transforms the residuals if \code{WB="TRUE"} is chosen.
#' If no transformation function for the residuals is specified, the identity function is used.
#' Currently, the weight distributions \code{"gumbel"} (Gumbel distribution with mode = 0, variance = 1, skewness = 1),
#' \code{"alaplace"} (asymmetric Laplace distribution with approx. mode = 0, variance = 1, skewness = 1),
#' \code{"mammen"} (continuous distribution proposed by Mammen (1993)),
#' \code{"normal"} (standard normal distribution), \code{"rademacher"} (Rademacher distribution),
#' \code{"mammen_two"} (Mammen's two-point distribution) and \code{"webb_six"} (Webb's six-point distribution) are implemented.
#' The parameter values are chosen such that the weight distribution has approximately
#' zero mode, variance one and skewness one.
#'
#' @param residuals original residuals
#' @param fu transformation function for the residuals
#' @param type type of weight distribution, currently eight options are implemented:
#' \code{"gumbel"}, \code{"alaplace"}, \code{"mammen"}, \code{"normal"}, \code{"rademacher"},
#' \code{"mammen_two"}, \code{"webb_six"}, \code{"feng"}
#' @return \describe{
#' \item{\code{residuals_WB}}{transformed residuals}
#' }
#' @importFrom gamlss.dist rGU
#' @importFrom LaplacesDemon ralaplace
# #' @importFrom extraDistr rsign
#' @importFrom stats rbinom
#' @export

trafo_WB <- function(
  residuals,
  fu = function(x) x,
  type = c("gumbel", "alaplace", "mammen", "normal", "rademacher", "mammen_two", "webb_six", "feng", "test")
)
{
  type <- match.arg(type)

  n <- length(residuals)

  if (type == "gumbel") {
    v <- rGU(n, mu = 0, sigma = sqrt(6/pi^2))
  } else if (type == "alaplace") {
    lambda <- 1.51
    kappa <- 0.78
    v <- ralaplace(n, location = 0, scale = sqrt(2)/lambda, kappa = kappa)
  } else if (type == "mammen"){
    u <- rnorm(n)
    w <- rnorm(n)
    v <- 1/sqrt(2)*u + 0.5*(w^2 - 1)
  } else if (type == "normal") {
    v <- rnorm(n)
  } else if (type == "rademacher") {
    # v <- rsign(n)
    v <- sample(c(1, -1), prob = c(0.5, 0.5), size = n, replace = TRUE)
  } else if (type == "mammen_two") {
    u <- rbinom(n, 1, (sqrt(5) + 1)/(2*sqrt(5)))
    v <- (-(sqrt(5) - 1)/2)*u + ((sqrt(5) + 1)/2)*(1 - u)
  } else if (type == "webb_six") {
    mass_points <- c(-sqrt(1.5), -1, -sqrt(0.5), sqrt(0.5), 1, sqrt(1.5))
    v <- replicate(n, sample(x = mass_points, size = 1, prob = rep(1/6, 6)))
  } else if (type == "feng") {
    u <- runif(n)
    v <- -sqrt((25/16 - 2*u)*1*(u < 0.5)) + sqrt((2*u - 7/16)*1*(u > 0.5))
  } else if (type == "test") {
    v <- sample(c(1.5, -1.5), prob = c(0.5, 0.5), size = n, replace = TRUE)
  } else {
    stop("You need to choose one of the following weight distributions:
         gumbel, alaplace, mammen, normal, rademacher, mammen_two, webb_six, feng.")
  }

  residuals_WB <- fu(residuals) * v
  return(residuals_WB)
}
