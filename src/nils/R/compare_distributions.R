#' Function comparing the kdes of the bootstrapped residuals and the true density for the different losses
#'
#' This function plots the kdes of the bootstrapped residuals and the true density for the different losses.
#' Furthermore, it calculates the Kullback-Leibler divergences and plots QQ plots.
#' As result, a table of 4x2 plots is saved.
#'
#' @param n_plot number of observations to plot
#' @param true_density true density function
#' @param true_quantile true quantile function
#' @param support support used for calculating the Kullback-Leibler divergences
#' @param data_path path to the folder which contains the data used for the preceeding simulation studies
#' @param data_filename filename of the data used for the preceeding simulation studies
#' @param results_path path to the folder which contains the results of the preceeding simulation studies
#' for the four different losses
#' @param results_filename filename (without the "_loss" ending) of the results files
#' @param plots_path path to the folder where the plot should be saved
#' @param plots_filename filename of the resulting plot
#' @importFrom purrr map
#' @importFrom LaplacesDemon KLD
#' @importFrom ggplot2 ggplot aes ggtitle theme element_text geom_histogram geom_density stat_function stat_qq stat_qq_line
#' @importFrom gridExtra grid.arrange
#' @export

compare_distributions <- function (
  n_plot = NULL,
  true_density,
  true_quantile,
  support,
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
  x1 <- x[[1]]
  y1 <- y[[1]]

  if (is.null(n_plot)) n_plot <- length(y[[1]])

  x1 <- x1[1:n_plot, ]
  y1 <- y1[1:n_plot]

  losses <- c("nils", "kemp", "mean", "med")
  bs_coefs_all <- kld_all <- plots <- sapply(losses, function(x) NULL)

  for (loss in c("nils", "kemp", "mean", "med")) {
    load(paste(results_path, "/", results_filename, "_", loss, ".Rdata", sep = ""))
    bs_coefs <- bs_coefs[[1]]
    bs_coefs_all[[loss]] <- bs_coefs
  }

  residuals_all <- lapply(bs_coefs_all, function (matrix) {apply(matrix, 2, function(beta) y1 - x1 %*% beta)})
  residuals_df <- lapply(residuals_all, function(matrix) {data.frame(residuals = as.vector(matrix))})

  kld_fct <- function (vec) {
    vec <- as.vector(vec)
    kde <- approxfun(density(vec, bw = "nrd0"))
    kde_values <- kde(seq(support[1], support[2], 0.01))
    density_values <- true_density(seq(support[1], support[2], 0.01))
    kld <- KLD(kde_values[!is.na(kde_values) & !is.na(density_values)], density_values[!is.na(kde_values) & !is.na(density_values)])
    return(c(round(kld$sum.KLD.px.py, digits = 3), round(kld$sum.KLD.py.px, digits = 3)))
  }

  kld_all <- lapply(residuals_all, kld_fct)

  plot_fct <- function (df, loss) {
    distr_plot <- ggplot(df, aes(x = residuals)) +
      geom_histogram(aes(y = ..density..), colour = "black", fill = "white") +
      geom_density(colour = "darkgreen", lwd = 0.75) +
      stat_function(fun = true_density, n = 1000, colour = "darkblue", lwd = 0.75) +
      ggtitle(paste(loss, ", KL divergence: ", kld_all[loss], sep = "")) +
      theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8))
    qq_plot <- ggplot(df, aes(sample = residuals)) +
      stat_qq(distribution = true_quantile, size = 0.75) +
      stat_qq_line(distribution = true_quantile, colour = "darkred", lwd = 0.75) +
      ggtitle(loss) +
      theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8))
    return(list(distr_plot = distr_plot, qq_plot = qq_plot))
  }

  plots <- mapply(plot_fct, residuals_df, as.list(losses), SIMPLIFY = FALSE)
  results_info <- strsplit(results_filename, '_')[[1]]
  results_info <- results_info[grepl("=" , results_info, fixed = TRUE)]
  results_info[grepl("dist" , results_info, fixed = TRUE)] <- paste("error_", results_info[grepl("dist" , results_info, fixed = TRUE)], sep = "")
  results_info[grepl("method" , results_info, fixed = TRUE)] <- paste("bs_ci_", results_info[grepl("method" , results_info, fixed = TRUE)], sep = "")
  results_info[grepl("type" , results_info, fixed = TRUE)][1] <- paste("trafo_", results_info[grepl("type" , results_info, fixed = TRUE)][1], sep = "")
  results_info[grepl("webb" , results_info, fixed = TRUE)] <- "type=webb_six"
  if (grepl("mammen_two", results_filename, fixed = TRUE)) results_info[grepl("mammen", results_info, fixed = TRUE)] <- paste(results_info[grepl("mammen", results_info, fixed = TRUE)], "_two", sep = "")
  results_info <- paste(results_info, collapse = ', ')
  plot <- grid.arrange(plots$nils$distr_plot, plots$nils$qq_plot,
                       plots$kemp$distr_plot, plots$kemp$qq_plot,
                       plots$mean$distr_plot, plots$mean$qq_plot,
                       plots$med$distr_plot, plots$med$qq_plot,
                       ncol = 2,
                       top = paste("Left: Histograms and kdes (in green) of bootstrapped residuals\nand true density of error terms (in blue) for different losses\nRight: QQ plots of the bootstrapped residuals for different losses\n(n_plot=",
                                   n_plot, ", ",
                                   paste(strwrap(results_info, width = 60), collapse = "\n"), ")", sep = ""))
  ggsave(filename = paste(plots_path, "/", plots_filename, "_comparision_distributions.png", sep = ""),
         plot = plot, width = 15, height = 30, units = "cm")
}
