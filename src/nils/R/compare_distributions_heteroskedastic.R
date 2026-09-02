#' Function comparing the kdes of the bootstrapped residuals and the true densities for the different losses
#' for the heteroskedasticity cases in the Kemp study
#'
#' This function plots the kdes of the bootstrapped residuals and the true densities for the different losses
#' for the heteroskedasticity cases in the Kemp study. For those cases, the true density depends on the regressor values.
#' Hence, there is for every observation another true density which needs to be calculated inside this function.
# #' Furthermore, it calculates the Kullback-Leibler divergences and plots QQ plots.
#' As result, a table of 4 plots is saved.
#'
# #' @param support support used for calculating the Kullback-Leibler divergences
#' @param n_plot number of observations to plot
#' @param bins parameter bins for \code{\link[ggplot2]{geom_histogram}}
#' @param data_path path to the folder which contains the data used for the preceeding simulation studies
#' @param data_filename filename of the data used for the preceeding simulation studies
#' @param results_path path to the folder which contains the results of the preceeding simulation studies
#' for the four different losses
#' @param results_filename filename (without the "_loss" ending) of the results files
#' @param plots_path path to the folder where the plot should be saved
#' @param plots_filename filename of the resulting plot
#' @importFrom purrr map
# #' @importFrom LaplacesDemon KLD
#' @importFrom ggplot2 ggplot aes aes_ ggtitle theme element_text geom_histogram geom_density stat_function stat_qq stat_qq_line
#' @importFrom gridExtra grid.arrange
#' @export

compare_distributions_heteroskedastic <- function (
  # support,
  n_plot = NULL,
  bins = bins,
  data_path = './data',
  data_filename,
  results_path = './results',
  results_filename,
  plots_path = './plots',
  plots_filename = results_filename
) {
  load(paste(data_path, "/", data_filename, ".Rdata", sep = ""))
  x <- lapply(data, function(mat) {mat[, -which(colnames(mat)=="y")]})
  y <- map(data, "y")
  x <- lapply(x, function (x, intercept) as.matrix(cbind(intercept, x)), intercept = rep(1, length(y[[1]])))
  x <- x1 <- x[[1]]
  y1 <- y[[1]]

  if (is.null(n_plot)) n_plot <- length(y[[1]])

  x <- x1 <- x1[1:n_plot, ]
  y1 <- y1[1:n_plot]

  losses <- c("nils", "kemp", "mean", "med")
  bs_coefs_all <- kld_all <- plots <- sapply(losses, function(x) NULL)

  for (loss in c("nils", "kemp", "mean", "med")) {
    load(paste(results_path, "/", results_filename, "_", loss, ".Rdata", sep = ""))
    bs_coefs <- bs_coefs[[1]]
    bs_coefs_all[[loss]] <- bs_coefs
  }

  residuals_all <- lapply(bs_coefs_all, function (matrix) {apply(matrix, 2, function(beta) y1 - x1 %*% beta)})
  residuals_df <- lapply(residuals_all, function(matrix) {data.frame(obs = as.factor(rep(1:length(y1), each = ncol(bs_coefs_all[[1]]))),
                                                                      residuals = as.vector(t(matrix)))})

  # kld_fct <- function (vec) {
  #   vec <- as.vector(vec)
  #   kde <- approxfun(density(vec, bw = "nrd0"))
  #   kde_values <- kde(seq(support[1], support[2], 0.01))
  #   density_values <- true_density(seq(support[1], support[2], 0.01))
  #   kld <- KLD(kde_values[!is.na(kde_values) & !is.na(density_values)], density_values[!is.na(kde_values) & !is.na(density_values)])
  #   return(c(round(kld$sum.KLD.px.py, digits = 3), round(kld$sum.KLD.py.px, digits = 3)))
  # }
  #
  # kld_all <- lapply(residuals_all, kld_fct)

  results_info <- strsplit(results_filename, '_')[[1]]
  alpha <- results_info[grepl("alpha=" , results_info, fixed = TRUE)]
  alpha <- as.numeric(strsplit(alpha, '=')[[1]][2])
  kappa <- alpha
  v <- results_info[grepl("v=" , results_info, fixed = TRUE)]
  v <- as.numeric(strsplit(v, '=')[[1]][2])
  m <- mean(x)
  m_sq <- mean(x^2)
  lambda <- ((log(kappa) - digamma(alpha))^2 * v^2 * (m_sq - m^2) +
               trigamma(alpha) * (1 + 2*v*m + v^2*m_sq))^(-0.5)

  true_density <- function(eps, x) {
    1/(lambda*(1 + v*x)) * (kappa^alpha)/gamma(alpha) * exp(-(alpha*eps)/(lambda*(1 + v*x))
                                                            - kappa*exp(-eps/(lambda*(1+v*x))))
  }

  plot_fct <- function (df, loss) {
    distr_plot <- ggplot(df, aes(x = residuals)) +
      geom_histogram(alpha = 0.5, aes(y = ..density.., fill = obs), position = 'identity', bins = bins) +
      Map(function(params, obs){stat_function(mapping = aes_(color = obs),
                                              fun = true_density, args = params)},
          params = as.list(x[, 2]),
          obs = levels(df$obs)) +
      # geom_density(colour = "darkgreen", lwd = 0.75) +
      # stat_function(fun = true_density, n = 1000, colour = "darkblue", lwd = 0.75) +
      ggtitle(loss) +
      theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8))
    # qq_plot <- ggplot(df, aes(sample = residuals)) +
    #   stat_qq(distribution = true_quantile, size = 0.75) +
    #   stat_qq_line(distribution = true_quantile, colour = "darkred", lwd = 0.75) +
    #   ggtitle(loss) +
    #   theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8))
    # return(list(distr_plot = distr_plot, qq_plot = qq_plot))
    return(list(distr_plot = distr_plot))
  }

  plots <- mapply(plot_fct, residuals_df, as.list(losses), SIMPLIFY = FALSE)

  results_info <- results_info[grepl("=" , results_info, fixed = TRUE)]
  results_info[grepl("dist" , results_info, fixed = TRUE)] <- paste("error_", results_info[grepl("dist" , results_info, fixed = TRUE)], sep = "")
  results_info[grepl("method" , results_info, fixed = TRUE)] <- paste("bs_ci_", results_info[grepl("method" , results_info, fixed = TRUE)], sep = "")
  results_info[grepl("type" , results_info, fixed = TRUE)][1] <- paste("trafo_", results_info[grepl("type" , results_info, fixed = TRUE)][1], sep = "")
  results_info[grepl("webb" , results_info, fixed = TRUE)] <- "type=webb_six"
  if (grepl("mammen_two", results_filename, fixed = TRUE)) results_info[grepl("mammen", results_info, fixed = TRUE)] <- paste(results_info[grepl("mammen", results_info, fixed = TRUE)], "_two", sep = "")
  if (grepl("sqrt_hat", results_filename, fixed = TRUE)) results_info[grepl("sqrt", results_info, fixed = TRUE)] <- "trafo_type=sqrt_hat"
  if (grepl("WB=FALSE", results_filename, fixed = TRUE)) results_info <- results_info[1:which(results_info == "WB=FALSE")]
  results_info <- paste(results_info, collapse = ', ')
  plot <- grid.arrange(plots$nils$distr_plot, # plots$nils$qq_plot,
                       plots$kemp$distr_plot, # plots$kemp$qq_plot,
                       plots$mean$distr_plot, # plots$mean$qq_plot,
                       plots$med$distr_plot, # plots$med$qq_plot,
                       ncol = 1,
                       top = paste("Histograms of bootstrapped residuals\nand true density of error terms for different losses\n(",
                                   paste(strwrap(results_info, width = 60), collapse = "\n"), ")", sep = ""))
  ggsave(filename = paste(plots_path, "/", plots_filename, "_comparision_distributions.png", sep = ""),
         plot = plot, width = 15, height = 40, units = "cm")
}
