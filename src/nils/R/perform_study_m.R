#' Function performing a simulation study using \code{\link{m}}
#'
#' This function performs a simulation study for given data, true coefficients and a specific loss function.
#' The results containing the true coefficients, the estimated coefficients, the fitted values, the residuals,
#' the asymptotic and bootstrap confidence intervals, the lengths and coverage rates of the confidence intervals,
#' the bandwiths used in the final IRLS iterations and the numbers of IRLS iterations are saved to an .RData file.
#'
#' @param data list containing the data for the different repetitions of the study
#' @param coefs true coefficients
#' @param bs_ci_method vector of bootstrap confidence interval construction methods to be used: \code{"percentile"} = percentile bootstrap,
#' \code{"rev_percentile"} = bootstrap reverse percentile, \code{"et_percentile-t"} = equal-tailed bootstrap percentile-t,
#' \code{"sym_percentile-t"} = symmetric bootstrap percentile-t
#' @param B_nested for \code{bs_ci_method = "et_percentile-t"}: number of nested bootstrap iterations per original bootstrap iteration
#' @param trafo_type type of transformation function to be used, currently six options are implemented:
#' \code{"identity"}, \code{"absolute"}, \code{"sqrt_hat"}, \code{"hat"}, \code{"df"}, \code{"feng"}
#' @param exact logical; for \code{loss="mean"}: \code{exact="TRUE"} can be chosen to calculate the exact
#' confidence interval under the normal distribution assumption
#' @param kernel for \code{loss="med"}: kernel for estimating the asymptotic covariance matrix (\code{gaussian} or \code{uniform})
#' @param path path were the .Rdata file containing the results should be saved
#' @param filename filename of results file
#' @inheritParams m
#' @importFrom pbmcapply pbmclapply
#' @importFrom purrr map
#' @export

perform_study_m <- function (
  data,
  coefs,
  start = NULL,
  loss = "nils",
  tuning = NULL,
  maxi = 5500,
  tau = 10^(-5),
  lambda = 0,
  K = NULL,
  B = 0,
  B_nested = B,
  bs_ci_method = c("percentile", "rev_percentile", "et_percentile-t", "sym_percentile-t"),
  WB = FALSE,
  trafo_type = "identity",
  type = "gumbel",
  B_fit = FALSE,
  alpha = 0.05,
  cov = TRUE,
  exact = FALSE,
  kernel = "gaussian",
  ncores = parallel::detectCores() - 1,
  path = './results',
  filename
) {
  nrep <- length(data)

  fu <- est_coefs  <- fitted <- residuals <- ci <- ci_length <- bs_ci <- bs_ci_length <- bs_coefs <- k <- iter <- as.list(rep(NA, nrep))

  x <- lapply(data, function(mat) {mat[, -which(colnames(mat)=="y")]})
  y <- map(data, "y")
  x <- lapply(x, function (x, intercept) as.matrix(cbind(intercept, x)), intercept = rep(1, length(y[[1]])))

  print(paste("Fitting of the according models for the", nrep, "iterations using", ncores, "cores:"), quote = FALSE)

  if (loss == "mean" || loss == "med") cov = FALSE

  if (trafo_type == "nk" && loss %in% c("nils", "kemp")) {
    controls <- mapply(m_control, x = x, y = y,
                       MoreArgs = list(loss = loss, tuning = tuning, maxi = maxi, start = start, lambda = lambda, K = K),
                       SIMPLIFY = FALSE)
    initial_fits <- pbmcmapply(m_fit, x = x, y = y, start = map(controls, "start"),
                               tuning = map(controls, "tuning"), maxi = map(controls, "maxi"), K = map(controls, "K"),
                               MoreArgs = list(loss = loss, d = deriv_loss(loss), tau = tau, lambda = lambda),
                               SIMPLIFY = FALSE)
    corrected_residuals <- mapply(est_cov, x = x, residuals = map(initial_fits, "residuals"), tuning = map(controls, "tuning"),
                   j = map(initial_fits, "j"),
                   MoreArgs = list(n = dim(x[[1]])[1], loss = loss, correction = TRUE), SIMPLIFY = FALSE)
    corrected_residuals <- map(corrected_residuals, "Omegahat")
    fu <- lapply(corrected_residuals, function(cr) {f <- function(u) abs(cr)})
  } else if (trafo_type == "nk" && !(loss %in% c("nils", "kemp"))) {
    stop("trafo_type = nk can only be chosen for the losses nils and kemp.")
  } else {
    fu <- lapply(x, transformation_function, trafo_type = trafo_type, n = length(y[[1]]), p = ncol(x[[1]]))
  }

  RNGkind("L'Ecuyer-CMRG")
  fits <- pbmcmapply(try_parallel(m), x = x, y = y, fu = fu,
                     MoreArgs = list(start = start, loss = loss, tuning = tuning,
                                     maxi = maxi, tau = tau, lambda = lambda, K = K, B = B,
                                     WB = WB, type = type, B_fit = B_fit, alpha = alpha,
                                     cov = cov, ncores = 1), mc.cores = ncores, mc.set.seed = TRUE, SIMPLIFY = FALSE)

  if ("warning" %in% names(fits)) fits <- fits$value
  error_fits <- unlist(lapply(fits, inherits, 'error'))
  conv <- unlist(map(fits, "converged"))
  fits_conv <- !error_fits & conv
  est_coefs[fits_conv] <- map(fits, "coefficients")[fits_conv]
  fitted[fits_conv] <- map(fits, "fitted")[fits_conv]
  residuals[fits_conv] <- map(fits, "residuals")[fits_conv]

  ci[fits_conv] <- mapply(confidence_interval, x = x[fits_conv], obj = fits[fits_conv],
                               MoreArgs = list(loss = loss, alpha = alpha, exact = exact, kernel = kernel),
                               SIMPLIFY = FALSE)
  ci_length[fits_conv] <- map(ci[fits_conv], "length")
  ci[fits_conv] <- map(ci[fits_conv], "ci")

  bs_coefs[fits_conv] <- map(fits[fits_conv], "bs_ciraw")[fits_conv]

  # if (bs_ci_method == "percentile") {
  #   bs_ci[fits_conv] <- map(fits[fits_conv], "bs_ci")
  #   bs_ci_length[fits_conv] <- lapply(bs_ci[fits_conv], function (v) abs(v[, 2] - v[, 1]))
  # } else if (B_fit == FALSE) {
  #   # bs_coefs[fits_conv] <- map(fits[fits_conv], "bs_ciraw")[fits_conv]
  #   bs_ci <- mapply(bootstrap_confidence_interval, bs_coefs[fits_conv], est_coefs = est_coefs, fu = fu, x = x, y = y,
  #                   MoreArgs = list(n = length(y[[1]]), alpha = alpha,
  #                                   bs_ci_method = bs_ci_method, B_nested = B_nested, WB = WB,
  #                                   type = type, start = start, loss = loss, tuning =  tuning,
  #                                   maxi = maxi, tau = tau, lambda = lambda, K = K),
  #                   SIMPLIFY = FALSE)
  #   bs_ci_length[fits_conv] <- lapply(bs_ci[fits_conv], function (v) abs(v[, 2] - v[, 1]))
  # }

  bs_ci_all <- bs_ci_length_all <- coverage_rates_all <- as.list(rep(NA, length(bs_ci_method)))
  names(bs_ci_all) <- names(bs_ci_length_all) <- names(coverage_rates_all) <- bs_ci_method

  bs_ci_method_all <- bs_ci_method
  if (sum(bs_ci_method %in% c("et_percentile-t", "sym_percentile-t")) == 2) {
    bs_ci <- mapply(bootstrap_confidence_interval, bs_coefs[fits_conv], est_coefs = est_coefs, fu = fu, x = x, y = y,
                    MoreArgs = list(n = length(y[[1]]), alpha = alpha,
                                    bs_ci_method = c("et_percentile-t", "sym_percentile-t"), B_nested = B_nested, WB = WB,
                                    type = type, start = start, loss = loss, tuning =  tuning,
                                    maxi = maxi, tau = tau, lambda = lambda, K = K),
                    SIMPLIFY = FALSE)
    bs_ci_all[["et_percentile-t"]] <- map(bs_ci, "et_percentile-t")
    bs_ci_all[["sym_percentile-t"]] <- map(bs_ci, "sym_percentile-t")
    bs_ci_length[fits_conv] <- lapply(bs_ci_all[["et_percentile-t"]][fits_conv], function (v) abs(v[, 2] - v[, 1]))
    bs_ci_length_all[["et_percentile-t"]] <- bs_ci_length
    bs_ci_length[fits_conv] <- lapply(bs_ci_all[["sym_percentile-t"]][fits_conv], function (v) abs(v[, 2] - v[, 1]))
    bs_ci_length_all[["sym_percentile-t"]] <- bs_ci_length
    coverage_rates_all[["et_percentile-t"]] <- coverage_rates(coefs, ci, bs_ci_all[["et_percentile-t"]])
    coverage_rates_all[["sym_percentile-t"]] <- coverage_rates(coefs, ci, bs_ci_all[["sym_percentile-t"]])

    bs_ci_method <- bs_ci_method[!bs_ci_method %in% c("et_percentile-t", "sym_percentile-t")]
  }

  if (length(bs_ci_method[!bs_ci_method %in% c("et_percentile-t", "sym_percentile-t")]) > 0) {
    for (method in bs_ci_method) {
      if (method == "percentile") {
        bs_ci[fits_conv] <- map(fits[fits_conv], "bs_ci")
        bs_ci_length[fits_conv] <- lapply(bs_ci[fits_conv], function (v) abs(v[, 2] - v[, 1]))
      } else if (B_fit == FALSE) {
        bs_ci <- mapply(bootstrap_confidence_interval, bs_coefs[fits_conv], est_coefs = est_coefs, fu = fu, x = x, y = y,
                        MoreArgs = list(n = length(y[[1]]), alpha = alpha,
                                        bs_ci_method = method, B_nested = B_nested, WB = WB,
                                        type = type, start = start, loss = loss, tuning =  tuning,
                                        maxi = maxi, tau = tau, lambda = lambda, K = K),
                        SIMPLIFY = FALSE)
        bs_ci_length[fits_conv] <- lapply(bs_ci[fits_conv], function (v) abs(v[, 2] - v[, 1]))
      }
      bs_ci_all[[method]] <- bs_ci
      bs_ci_length_all[[method]] <- bs_ci_length
      coverage_rates_all[[method]] <- coverage_rates(coefs, ci, bs_ci)
    }
  }

  k[fits_conv] <- map(fits[fits_conv], "k")

  iter[fits_conv] <- map(fits[fits_conv], "iter")

  # coverage_rates <- coverage_rates(coefs, ci, bs_ci)

  # save(coefs, est_coefs, bs_coefs, fitted, residuals, ci, bs_ci, ci_length, bs_ci_length, coverage_rates, k, iter,
  #      file = paste(path, "/", filename, "_", loss, ".Rdata", sep = ""))

  for (method in bs_ci_method_all) {
    filename_split <- strsplit(filename, '_')[[1]]
    filename_split[grepl("method", filename_split, fixed = TRUE)] <- paste("method=", method, sep = "")
    filename_split <- filename_split[c(1:grep("method", filename_split, fixed = TRUE, value=FALSE),
                                      grep("WB", filename_split, fixed = TRUE, value=FALSE):length(filename_split))]
    filename <- paste(filename_split, collapse = "_")
    bs_ci <- bs_ci_all[[method]]
    bs_ci_length <- bs_ci_length_all[[method]]
    coverage_rates <- coverage_rates_all[[method]]
    save(coefs, est_coefs, bs_coefs, fitted, residuals, ci, bs_ci, ci_length,
         bs_ci_length, coverage_rates, k, iter,
         file = paste(path, "/", filename, "_", loss, ".Rdata", sep = ""))
  }

}
