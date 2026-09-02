#' Function plotting histograms, kdes and bootstrap confidence intervals for given bootstrapped coefficients
#'
#' This function plots the histograms, kdes and bootstrap confidence intervals for given bootstrapped coefficients.
#' The true coefficients are indicated by blue vertical dashed lines. As result, a table containing the plots
#' for the different coefficients is saved.
#'
# #' @param all logical; if TRUE, all repitions are used, otherwise only the first repitition is used
#' @param nrep_plot number of repetitions to be plotted
#' @param results_path path to the folder which contains the results of the preceeding simulation study
#' @param results_filename filename of the result file
#' @param plots_path path to the folder where the plot should be saved
#' @param plots_filename filename of the resulting plot
#' @importFrom ggplot2 ggplot aes ggtitle theme element_text geom_histogram geom_density geom_vline
#' @importFrom gridExtra grid.arrange
#' @export

bs_coefs_distributions <- function (
  nrep_plot = NULL,
  results_path = './results',
  results_filename,
  plots_path = './plots',
  plots_filename = results_filename
) {
  load(paste(results_path, "/", results_filename, ".Rdata", sep = ""))

  if (is.null(nrep_plot)) {
    bs_coefs <- bs_coefs[[1]]
    bs_ci <- bs_ci[[1]]
  } else {
    bs_coefs <- do.call(cbind, bs_coefs)[, 1:nrep_plot]
    bs_ci <- do.call(cbind, bs_ci)[, 1:nrep_plot]
  }

  bs_coefs_df <- as.list(data.frame(t(bs_coefs)))
  bs_coefs_df <- lapply(bs_coefs_df, function (vec) {data.frame(bs_coefs = vec)})

  names(coefs) <- paste("beta_", 0:(length(coefs)-1), sep = "")

  plot_fct <- function (df, ciquantiles, coef, name) {
    distr_plot <- ggplot(df, aes(x = bs_coefs)) +
      geom_histogram(aes(y = ..density..), colour = "black", fill = "white") +
      geom_density(colour = "darkgreen", lwd = 0.75) +
      geom_vline(xintercept = ciquantiles, color = "darkred", linetype = "dashed", size = 0.75) +
      geom_vline(xintercept = coef, color="darkblue", linetype = "dashed", size = 0.75) +
      ggtitle(paste("Coefficient", name)) +
      theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8))
  }

  plots <- mapply(plot_fct, bs_coefs_df, as.list(data.frame(t(bs_ci))),
                  as.list(coefs), as.list(names(coefs)), SIMPLIFY = FALSE)
  names(plots) <- names(coefs)
  plot <- do.call("grid.arrange", c(plots, ncol = 2,
                                    top = "Histograms, kdes and bootstrap confidence intervals for the different coefficients\n(blue vertical line indicates the true coefficient)"))
  ggsave(filename = paste(plots_path, "/", plots_filename, "_bs_coefs_distributions.png", sep = ""),
         plot = plot, width = 15, height = 15, units = "cm")
}
