#' Function collecting and plotting the results of a simulation study performed with \code{\link{perform_study_m}}
#'
#' This function plots and saves boxplots of the estimated coefficients, of the lengths of the
#' asymptotic confidence intervals as well as of the bootstrap confidence intervals
#' and of the root mean squared errors. Furthermore, it collects and saves the coverage rates and the means
#' and standard deviations of the estimated coefficients.
#'
#' @param results_path path to the folder which contains the results of the preceeding simulation studies
#' for the four different losses
#' @param results_filename filename (without the "_loss" ending) of the results files
#' @param plots_path path to the folder where the plots should be saved
#' @param plots_filename filename of the resulting plots
#' @param tables_path path to the folder where the tables should be saved
#' @param tables_filename filename of generated tables
#' @importFrom ggplot2 ggplot ggtitle geom_boxplot geom_segment geom_text ggsave aes aes_string position_dodge
#' @importFrom dplyr '%>%' summarise group_by
# #' @importFrom xtable xtable print.xtable
#' @export

results_m <- function (
  results_path = './results',
  results_filename,
  plots_path = './plots',
  plots_filename = results_filename,
  tables_path = './tables',
  tables_filename = results_filename
) {
  results_info <- strsplit(results_filename, '_')[[1]]
  results_info <- results_info[grepl("=" , results_info, fixed = TRUE)]
  results_info[grepl("dist" , results_info, fixed = TRUE)] <- paste("error_", results_info[grepl("dist" , results_info, fixed = TRUE)], sep = "")
  results_info[grepl("method" , results_info, fixed = TRUE)] <- paste("bs_ci_", results_info[grepl("method" , results_info, fixed = TRUE)], sep = "")
  results_info[grepl("=rev" , results_info, fixed = TRUE)] <- "bs_ci_method=rev_percentile"
  results_info[grepl("=et" , results_info, fixed = TRUE)] <- "bs_ci_method=et_percentile-t"
  results_info[grepl("=sym" , results_info, fixed = TRUE)] <- "bs_ci_method=sym_percentile-t"
  results_info[grepl("type" , results_info, fixed = TRUE)][1] <- paste("trafo_", results_info[grepl("type" , results_info, fixed = TRUE)][1], sep = "")
  results_info[grepl("webb" , results_info, fixed = TRUE)] <- "type=webb_six"
  if (grepl("mammen_two", results_filename, fixed = TRUE)) results_info[grepl("mammen", results_info, fixed = TRUE)] <- paste(results_info[grepl("mammen", results_info, fixed = TRUE)], "_two", sep = "")
  if (grepl("sqrt_hat", results_filename, fixed = TRUE)) results_info[grepl("sqrt", results_info, fixed = TRUE)] <- "trafo_type=sqrt_hat"
  if (grepl("WB=FALSE", results_filename, fixed = TRUE)) results_info <- results_info[1:which(results_info == "WB=FALSE")]
  results_info <- paste(results_info, collapse = ', ')

  losses <- c("nils", "kemp", "mean", "med")
  est_coefs_all <- residuals_all <- ci_length_all <- bs_ci_length_all <- coverage_rates_all <- sapply(losses, function(x) NULL)

  for (loss in c("nils", "kemp", "mean", "med")) {
    load(paste(results_path, "/", results_filename, "_", loss, ".Rdata", sep = ""))
    est_coefs_all[[loss]] <- matrix(unlist(est_coefs), nrow = length(est_coefs), byrow = TRUE)
    ci_length_all[[loss]] <- matrix(unlist(ci_length), nrow = length(ci_length), byrow = TRUE)
    bs_ci_length_all[[loss]] <- matrix(unlist(bs_ci_length), nrow = length(bs_ci_length), byrow = TRUE)
    residuals_all[[loss]] <- matrix(unlist(residuals), nrow = length(residuals[[1]]), byrow = FALSE)
    coverage_rates_all[[loss]] <- coverage_rates
  }

  nrep <- dim(est_coefs_all$nils)[1]

  loss <- rep(losses, each = length(coefs)*nrep)
  coefficients <- rep(rep(paste("beta_", 0:(length(coefs) - 1), sep = ""), each = nrep), 4)

  estimates <- c(as.vector(est_coefs_all$nils), as.vector(est_coefs_all$kemp),
                 as.vector(est_coefs_all$mean), as.vector(est_coefs_all$med))
  ci_lengths <- c(as.vector(ci_length_all$nils), as.vector(ci_length_all$kemp),
                  as.vector(ci_length_all$mean), as.vector(ci_length_all$med))
  bs_ci_lengths <- c(as.vector(bs_ci_length_all$nils), as.vector(bs_ci_length_all$kemp),
                     as.vector(bs_ci_length_all$mean), as.vector(bs_ci_length_all$med))

  df_coverage_rates <- data.frame(loss = rep(losses, each = length(coefs)),
                                  coefficient = rep(rep(paste("beta_", 0:(length(coefs) - 1), sep = "")), 4),
                                  do.call(rbind, coverage_rates_all), row.names = NULL)

  df_est_coefs <- data.frame(loss, coefficients, estimates)
  df_est_coefs$loss <- factor(df_est_coefs$loss, levels = losses)
  summary_est_coefs <- df_est_coefs %>% group_by(loss, coefficients) %>%
    summarise(mean = mean(estimates), sd = sd(estimates))

  df_ci_lengths <- data.frame(loss, coefficients, ci_lengths)
  df_ci_lengths$yloc <- rep(min(df_ci_lengths$ci_lengths)
                            -(max(df_ci_lengths$ci_lengths)-min(df_ci_lengths$ci_lengths))/20, each = 4*length(coefs)*nrep)
  df_ci_lengths$label <- rep(NA, each = 4*length(coefs)*nrep)
  df_ci_lengths$label[1+nrep*(0:(4*length(coefs)-1))] <- df_coverage_rates$cover_rate
  df_ci_lengths$loss <- factor(df_ci_lengths$loss, levels = losses)

  df_bs_ci_lengths <- data.frame(loss, coefficients, bs_ci_lengths)
  df_bs_ci_lengths$yloc <- rep(min(df_bs_ci_lengths$bs_ci_lengths)
                               -(max(df_bs_ci_lengths$bs_ci_lengths)-min(df_bs_ci_lengths$bs_ci_lengths))/20, each = 4*length(coefs)*nrep)
  df_bs_ci_lengths$label <- rep(NA, each = 4*length(coefs)*nrep)
  df_bs_ci_lengths$label[1+nrep*(0:(4*length(coefs)-1))] <- df_coverage_rates$bs_cover_rate
  df_bs_ci_lengths$loss <- factor(df_bs_ci_lengths$loss, levels = losses)

  rmse_list <- lapply(residuals_all, function (res) {sqrt(colMeans(res^2))})
  df_rmse <- data.frame(loss = rep(losses, each = nrep), rmse = unlist(rmse_list))
  df_rmse$loss <- factor(df_rmse$loss, levels = losses)

  # save summary data frame for estimated coefficients

  save(summary_est_coefs, file = paste(tables_path, "/", tables_filename, "_summary_est_coefs.Rdata", sep = ""))
  # print(xtable(summary_est_coefs, digits = 3, type = "latex"),
  #       file = paste(tables_path, "/", tables_filename, "_summary_est_coefs.tex", sep = ""))

  # save data frame containing coverage rates

  save(df_coverage_rates, file = paste(tables_path, "/", tables_filename, "_coverage_rates.Rdata", sep = ""))
  # print(xtable(df_coverage_rates, digits = 3, type = "latex"),
  #       file = paste(tables_path, "/", tables_filename, "_coverage_rates.tex", sep = ""))

  # boxplots of estimated coefficients

  boxplots_est_coefs <- ggplot(df_est_coefs, aes(x = coefficients, y = estimates, fill = loss)) +
    ggtitle("Boxplots of estimated coefficients", subtitle = paste("(nrep=", nrep, ", ", strsplit(results_info, ', B=')[[1]][1], ")", sep="")) +
    geom_boxplot(outlier.shape = 1)

  for (i in 1:length(coefs)) {
    boxplots_est_coefs <- boxplots_est_coefs +
      geom_segment(aes_string(x = i-1+0.6, xend = i-1+1.4, y = coefs[i], yend = coefs[i]))
  }

  # boxplots of asymptotic confidence interval lengths

  boxplots_ci_lengths <- ggplot(df_ci_lengths, aes(x = coefficients, y = ci_lengths, fill = loss)) +
    ggtitle("Boxplots of asymptotic confidence interval lengths", subtitle = paste("(nrep=", nrep, ", ", strsplit(results_info, ', B=')[[1]][1], ")", sep="")) +
    geom_boxplot(outlier.shape = 1) +
    geom_text(aes(y = yloc, label = label), position = position_dodge(width = .75), size = 3, angle = 90)

  boxplots_ci_lengths_nils_kemp <- ggplot(df_ci_lengths[df_ci_lengths$loss == "nils" | df_ci_lengths$loss == "kemp", ], aes(x = coefficients, y = ci_lengths, fill = loss)) +
    ggtitle("Boxplots of asymptotic confidence interval lengths", subtitle = paste("for", nrep, "repetitions")) +
    geom_boxplot(outlier.shape = 1)

  # boxplots of bootstrap confidence interval lengths

  boxplots_bs_ci_lengths <- ggplot(df_bs_ci_lengths, aes(x = coefficients, y = bs_ci_lengths, fill = loss)) +
    ggtitle("Boxplots of bootstrap confidence interval lengths", subtitle =
              paste("(nrep=", nrep, ", ", paste(strwrap(results_info, width = 75), collapse = "\n"), ")", sep = "")) +
    geom_boxplot(outlier.shape = 1) +
    geom_text(aes(y = yloc, label = label), position = position_dodge(width = .75), size = 3, angle = 90)

  boxplots_bs_ci_lengths_nils_kemp <- ggplot(df_bs_ci_lengths[df_bs_ci_lengths$loss == "nils" | df_bs_ci_lengths$loss == "kemp", ], aes(x = coefficients, y = bs_ci_lengths, fill = loss)) +
    ggtitle("Boxplots of bootstrap confidence interval lengths", subtitle = paste("for", nrep, "repetitions")) +
    geom_boxplot(outlier.shape = 1)

  # boxplots of RMSE

  boxplots_rmse <- ggplot(df_rmse, aes(x = loss, y = rmse)) +
    ggtitle("Boxplots of root mean squared errors", subtitle = paste("(nrep=", nrep, ", ", strsplit(results_info, ', B=')[[1]][1], ")", sep="")) +
    geom_boxplot(outlier.shape = 1)

  # save boxplots

  ggsave(filename = paste(plots_path, "/", plots_filename, "_boxplots_est_coefs.png", sep = ""), plot = boxplots_est_coefs)
  ggsave(filename = paste(plots_path, "/", plots_filename, "_boxplots_ci_lengths.png", sep = ""), plot = boxplots_ci_lengths)
  ggsave(filename = paste(plots_path, "/", plots_filename, "_boxplots_bs_ci_lengths.png", sep = ""), plot = boxplots_bs_ci_lengths)
  ggsave(filename = paste(plots_path, "/", plots_filename, "_boxplots_rmse.png", sep = ""), plot = boxplots_rmse)
}
