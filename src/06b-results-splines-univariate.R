# rm(list = ls())
# 
# suppressPackageStartupMessages(library(tidyverse))
# library(reshape2)
# library(gridExtra)
# 
# input_final = readRDS("sim_study_nonC-GBI/splines_univ/sim_study_nonCalibrated_asymm_TL_final.RDS")
# input_init = readRDS("sim_study_nonC-GBI/splines_univ/sim_study_nonCalibrated_asymm_TL_init.RDS")
# 
# summ_final = lapply(input_final,
#                     function(I)
#                       abind::abind(lapply(I, function(x) x$summ), along = 3))
# summ_init = lapply(input_init,
#                    function(I)
#                      abind::abind(lapply(I, function(x) x$summ), along = 3))
# 
# PATH_IMG = "sim_study_nonC-GBI/splines_univ/"
# 
# # TRUE BETA AND SAMPLE SIZE ----
# n = 1e3
# 
# 
# # SIMULATION SETTINGS ----
# error_distr = c("Gaussian", "Gamma")
# cov_distr = c("Unif_rad2", "Beta_2_4") #"Gamma_1.5_1.5") #, "BVN")
# transform_x = c("parabola", "cubic", "trigonometric")
# SNR = c(2, 1)
# sim_settings = as.matrix(expand.grid(error_distr, cov_distr, transform_x, SNR))
# colnames(sim_settings) = c("err_distr", "cov_distr", "transf_x", "SNR")
# 
# # COLORS ----
# col_final = "#023e8a"
# col_init = "#588157"
# 
# 
# # FIT vs TRUE FUNCTION ----
# ## Fit ----
# fit_final = 
#   fit_init = array(NA, c(nrow(sim_settings), n, 100),
#                    dimnames = list("setting" = 1:24, "unit" = 1:n, "replicate" = 1:100))
# 
# for (s in 1:nrow(sim_settings)){
#   
#   cat(s, "\t")
#   
#   I_f = input_final[[s]]
#   I_i = input_init[[s]]
#   
#   fit_f = sapply(I_f, function(I) rowMeans(I$X_scaled %*% t(I$draws)))
#   fit_i = sapply(I_i, function(I) rowMeans(I$X_scaled %*% t(I$draws)))
#   
#   fit_final[s, , ] = fit_f
#   fit_init[s, , ] = fit_i
# }
# 
# ## X ----
# X_final = 
#   X_init = array(NA, c(nrow(sim_settings), n, 100),
#                  dimnames = list("setting" = 1:24, "unit" = 1:n, "replicate" = 1:100))
# 
# for (s in 1:nrow(sim_settings)){
#   
#   cat(s, "\t")
#   
#   I_f = input_final[[s]]
#   I_i = input_init[[s]]
#   
#   X_f = sapply(I_f, function(I) I$x1, simplify = "array")
#   
#   X_i = sapply(I_i, function(I) I$x1, simplify = "array")
#   
#   X_final[s, , ] = X_f
#   X_init[s, , ] = X_i
# }
# 
# 
# ## True functions ----
# true_f = data.frame(setting = integer(0L), x = numeric(0L), y = numeric(0L))
# for (s in 1:nrow(sim_settings)){
#   xmin = min(c(X_final[s, , ], X_init[s, , ]), na.rm = T)
#   xmax = max(c(X_final[s, , ], X_init[s, , ]), na.rm = T)
#   
#   transf = switch(sim_settings[s, "transf_x"],
#                   "parabola" = function(x) x^2,
#                   "cubic" = function(x) x^3,
#                   "trigonometric" = function(x) sin(2.5*x) + 2*exp(-(5^2)*((x)^2)))
#   
#   x = seq(xmin, xmax, length = 100)
#   y = transf(x)
#   
#   true_f = bind_rows(true_f, data.frame(setting = s, x = x, y = y))
# }
# true_f = left_join(true_f, cbind.data.frame(setting = 1:24, sim_settings), by = "setting")
# 
# 
# df_fit_final = left_join(fit_final %>% melt(value.name = "fit"),
#                          X_final %>% melt(value.name = "x"),
#                          by = c("setting", "unit", "replicate")) %>% 
#   left_join(cbind.data.frame(setting = 1:24, sim_settings), by = "setting")
# 
# df_fit_init = left_join(fit_init %>% melt(value.name = "fit"),
#                         X_init %>% melt(value.name = "x"),
#                         by = c("setting", "unit", "replicate")) %>% 
#   left_join(cbind.data.frame(setting = 1:24, sim_settings), by = "setting")
# 
# plt_fit_final2 = df_fit_final %>% 
#   filter(SNR == 2) %>% 
#   ggplot() +
#   geom_line(aes(x = x, y = fit, group = replicate), alpha = 0.15) +
#   geom_line(data = true_f %>% filter(SNR == 2),
#             mapping = aes(x = x, y = y), col = "red") +
#   facet_wrap(~setting, scales = "free", labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
#   labs(title = "Final - SNR = 2")
# 
# plt_fit_init2 = df_fit_init %>% 
#   filter(SNR == 2) %>% 
#   ggplot() +
#   geom_line(aes(x = x, y = fit, group = replicate), alpha = 0.15) +
#   geom_line(data = true_f %>% filter(SNR == 2),
#             mapping = aes(x = x, y = y), col = "red") +
#   facet_wrap(~setting, scales = "free", labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
#   labs(title = "Initial - SNR = 2")
# 
# plt_fit_final1 = df_fit_final %>% 
#   filter(SNR == 1) %>% 
#   ggplot() +
#   geom_line(aes(x = x, y = fit, group = replicate), alpha = 0.15) +
#   geom_line(data = true_f %>% filter(SNR == 1),
#             mapping = aes(x = x, y = y), col = "red") +
#   facet_wrap(~setting, scales = "free", labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
#   labs(title = "Final - SNR = 1")
# 
# plt_fit_init1 = df_fit_init %>% 
#   filter(SNR == 1) %>% 
#   ggplot() +
#   geom_line(aes(x = x, y = fit, group = replicate), alpha = 0.15) +
#   geom_line(data = true_f %>% filter(SNR == 1),
#             mapping = aes(x = x, y = y), col = "red") +
#   facet_wrap(~setting, scales = "free", labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
#   labs(title = "Initial - SNR = 1")
# 
# ggsave(filename = "fit.pdf",
#        plot = marrangeGrob(list(plt_fit_final2, plt_fit_init2, plt_fit_final1, plt_fit_init1), nrow = 1, ncol = 1, top = NULL),
#        path = PATH_IMG, width = 12, height = 9)
# 
# # DENSITY RESIDUALS ----
# y_init = 
#   y_final = array(NA, c(nrow(sim_settings), n, 100),
#                   dimnames = list("setting" = 1:24, "unit" = 1:n, "replicate" = 1:100))
# 
# for (s in 1:nrow(sim_settings)){
#   
#   cat(s, "\t")
#   
#   I_f = input_final[[s]]
#   I_i = input_init[[s]]
#   
#   y_f = sapply(I_f, function(I) I$y, simplify = "array")
#   
#   y_i = sapply(I_i, function(I) I$y, simplify = "array")
#   
#   y_final[s, , ] = y_f
#   y_init[s, , ] = y_i
# }
# 
# df_y_init = left_join(y_init %>% melt(value.name = "y"),
#                       fit_init %>% melt(value.name = "fit"),
#                       by = c("setting", "unit", "replicate")) %>% 
#   left_join(cbind.data.frame(setting = 1:24, sim_settings), by = "setting")
# 
# df_y_final = left_join(y_final %>% melt(value.name = "y"),
#                        fit_final %>% melt(value.name = "fit"),
#                        by = c("setting", "unit", "replicate")) %>% 
#   left_join(cbind.data.frame(setting = 1:24, sim_settings), by = "setting")
# 
# 
# plt_res_final2 = df_y_final %>%
#   filter(SNR == 2) %>% 
#   ggplot() +
#   geom_density(aes(x = y - fit, group = replicate), color = alpha("black", alpha = 0.2)) +
#   facet_wrap(~setting, scales = "free", labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
#   labs(title = "Final - SNR = 2")
# 
# plt_res_init2 = df_y_init %>%
#   filter(SNR == 2) %>% 
#   ggplot() +
#   geom_density(aes(x = y - fit, group = replicate), color = alpha("black", alpha = 0.2)) +
#   facet_wrap(~setting, scales = "free", labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
#   labs(title = "Initial - SNR = 2")
# 
# plt_res_final1 = df_y_final %>%
#   filter(SNR == 1) %>% 
#   ggplot() +
#   geom_density(aes(x = y - fit, group = replicate), color = alpha("black", alpha = 0.2)) +
#   facet_wrap(~setting, scales = "free", labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
#   labs(title = "Final - SNR = 1")
# 
# plt_res_init1 = df_y_init %>%
#   filter(SNR == 1) %>%
#   ggplot() +
#   geom_density(aes(x = y - fit, group = replicate), color = alpha("black", alpha = 0.2)) +
#   facet_wrap(~setting, scales = "free", labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
#   labs(title = "Initial - SNR = 1")
# 
# ggsave(filename = "residuals.pdf",
#        plot = marrangeGrob(list(plt_res_final2, plt_res_init2, plt_res_final1, plt_res_init1), nrow = 1, ncol = 1),
#        path = PATH_IMG, width = 12, height = 9)
# 
# 
# # DIAGNOSTICS ----
# ## Parabola ----
# I = input_final[[1]][[1]]
# 
# # Plot DR basis ----
# pdf(paste0(PATH_IMG, "DR_basis_parabola_example.pdf"), width = 12, height = 8)
# names_X = c("Intercept", "Linear", paste0("DR basis ", 1:20))
# par(mfrow = c(4, 6), mar = c(3, 2, 2, 1))
# for (j in 1:22){
#   plot(I$x1[order(I$x1)], I$X_scaled[order(I$x1), j], type = "l",
#        main = names_X[j],
#        xlab = "", ylab = "0")
# }
# dev.off()
# 
# 
# # Traceplots coefficients
# pdf(paste0(PATH_IMG, "traceplot_parabola_example.pdf"), width = 12, height = 8)
# par(mfrow = c(5, 6), mar = c(3, 2, 2, 1))
# for (j in 1:22){
#   plot(I$draws[ , j], type = "l", 
#        main = colnames(I$draws)[j], xlab = "MCMC Iteration", ylab = "")
# }
# plot.new(); plot.new()
# err = colMeans(I$X_scaled %*% t(I$draws) - matrix(I$y, nrow = 1000, ncol = 1000, byrow = F))
# plot(err, type = "l", main = "y - fit", xlab = "MCMC Iteration", ylab = "")
# err_noInt = colMeans(I$X_scaled[ , -1] %*% t(I$draws[ , -1]) - matrix(I$y, nrow = 1000, ncol = 1000, byrow = F))
# plot(err_noInt, type = "l", main = "y - fit (no beta0)", xlab = "MCMC Iteration", ylab = "")
# err_Square = colMeans(I$X_scaled[ , 22] %*% t(I$draws[ , 22]) - matrix(I$y, nrow = 1000, ncol = 1000, byrow = F))
# plot(err_Square, type = "l", main = "y - fit (only square)", xlab = "MCMC Iteration", ylab = "")
# dev.off()
# 
# # "Traceplot" reconstruction function
# df_complete = cbind(1:1e3, scale(I$x1), I$X_scaled %*% t(I$draws)) %>% 
#   as.data.frame() %>% 
#   pivot_longer(cols = -(V1:V2), names_to = "Iter", values_to = "x^2") %>%
#   rename(Unit = V1, x = V2) %>% 
#   mutate(Iter = as.integer(str_remove(Iter, "V"))-2)
# 
# df_noInt = cbind(1:1e3, scale(I$x1), I$X_scaled[ , -1] %*% t(I$draws[, -1])) %>% 
#   as.data.frame() %>% 
#   pivot_longer(cols = -(V1:V2), names_to = "Iter", values_to = "x^2") %>%
#   rename(Unit = V1, x = V2) %>% 
#   mutate(Iter = as.integer(str_remove(Iter, "V"))-2)
# 
# df_Z = cbind(1:1e3, scale(I$x1), I$Z %*% t(I$draws[, -(1:2)])) %>% 
#   as.data.frame() %>% 
#   pivot_longer(cols = -(V1:V2), names_to = "Iter", values_to = "x^2") %>%
#   rename(Unit = V1, x = V2) %>% 
#   mutate(Iter = as.integer(str_remove(Iter, "V"))-2)
# 
# df_Square = cbind(1:1e3, scale(I$x1), I$X_scaled[ , c(1, 22)] %*% t(I$draws[, c(1, 22)])) %>% 
#   as.data.frame() %>% 
#   pivot_longer(cols = -(V1:V2), names_to = "Iter", values_to = "x^2") %>%
#   rename(Unit = V1, x = V2) %>% 
#   mutate(Iter = as.integer(str_remove(Iter, "V"))-2)
# 
# traceplot_fit = bind_rows(
#   # df_Z %>% mutate(grp = "No Int + No Linear"),
#   df_noInt %>% mutate(grp = "No Int"),
#   df_complete %>% mutate(grp = "Complete"),
#   df_Square %>% mutate(grp = "Square + Int")) %>% 
#   mutate(grp = factor(grp, levels = c("Complete", "Square + Int", "No Int"))) %>% 
#   ggplot() +
#   geom_line(aes(x = x, y = `x^2`, group = Iter), alpha = 0.25) + 
#   geom_line(data = cbind.data.frame(x = scale(I$x1), transf = I$x1^2),
#             mapping = aes(x = x, y = transf), col = "red") + 
#   facet_wrap(~grp) +
#   labs(title = "Predictor during MCMC iterations")
# 
# ggsave(path = PATH_IMG,
#        filename = "traceplot_fit_parabola_example.pdf",
#        plot = traceplot_fit,
#        width = 12, height = 6)
# 
# 
# 
# ## Cubic ----
# I = input_final[[5]][[1]]
# 
# # Plot DR basis ----
# pdf(paste0(PATH_IMG, "DR_basis_cubic_example.pdf"), width = 12, height = 8)
# names_X = c("Intercept", "Linear", paste0("DR basis ", 1:20))
# par(mfrow = c(4, 6), mar = c(3, 2, 2, 1))
# for (j in 1:22){
#   plot(I$x1[order(I$x1)], I$X_scaled[order(I$x1), j], type = "l",
#        main = names_X[j],
#        xlab = "", ylab = "0")
# }
# dev.off()
# 
# 
# # Traceplots coefficients
# pdf(paste0(PATH_IMG, "traceplot_cubic_example.pdf"), width = 12, height = 8)
# par(mfrow = c(5, 6), mar = c(3, 2, 2, 1))
# for (j in 1:22){
#   plot(I$draws[ , j], type = "l", 
#        main = colnames(I$draws)[j], xlab = "MCMC Iteration", ylab = "")
# }
# plot.new(); plot.new()
# err = colMeans(I$X_scaled %*% t(I$draws) - matrix(I$y, nrow = 1000, ncol = 1000, byrow = F))
# plot(err, type = "l", main = "y - fit", xlab = "MCMC Iteration", ylab = "")
# err_noInt = colMeans(I$X_scaled[ , -1] %*% t(I$draws[ , -1]) - matrix(I$y, nrow = 1000, ncol = 1000, byrow = F))
# plot(err_noInt, type = "l", main = "y - fit (no beta0)", xlab = "MCMC Iteration", ylab = "")
# err_Square = colMeans(I$X_scaled[ , 22] %*% t(I$draws[ , 22]) - matrix(I$y, nrow = 1000, ncol = 1000, byrow = F))
# plot(err_Square, type = "l", main = "y - fit (only square)", xlab = "MCMC Iteration", ylab = "")
# dev.off()
# 
# # "Traceplot" reconstruction function
# df_complete = cbind(1:1e3, scale(I$x1), I$X_scaled %*% t(I$draws)) %>% 
#   as.data.frame() %>% 
#   pivot_longer(cols = -(V1:V2), names_to = "Iter", values_to = "x^2") %>%
#   rename(Unit = V1, x = V2) %>% 
#   mutate(Iter = as.integer(str_remove(Iter, "V"))-2)
# 
# df_noInt = cbind(1:1e3, scale(I$x1), I$X_scaled[ , -1] %*% t(I$draws[, -1])) %>% 
#   as.data.frame() %>% 
#   pivot_longer(cols = -(V1:V2), names_to = "Iter", values_to = "x^2") %>%
#   rename(Unit = V1, x = V2) %>% 
#   mutate(Iter = as.integer(str_remove(Iter, "V"))-2)
# 
# df_Partial = cbind(1:1e3, scale(I$x1), I$X_scaled[ , c(1:2, 21:22)] %*% t(I$draws[, c(1:2, 21:22)])) %>% 
#   as.data.frame() %>% 
#   pivot_longer(cols = -(V1:V2), names_to = "Iter", values_to = "x^2") %>%
#   rename(Unit = V1, x = V2) %>% 
#   mutate(Iter = as.integer(str_remove(Iter, "V"))-2)
# 
# traceplot_fit = bind_rows(
#   df_noInt %>% mutate(grp = "No Int"),
#   df_complete %>% mutate(grp = "Complete"),
#   df_Partial %>% mutate(grp = "Lin + DR19 + DR20 + Int")) %>% 
#   mutate(grp = factor(grp, levels = c("Complete", "Lin + DR19 + DR20 + Int", "No Int"))) %>% 
#   ggplot() +
#   geom_line(aes(x = x, y = `x^2`, group = Iter), alpha = 0.25) + 
#   geom_line(data = cbind.data.frame(x = scale(I$x1), transf = I$x1^3),
#             mapping = aes(x = x, y = transf), col = "red") + 
#   facet_wrap(~grp) +
#   labs(title = "Predictor during MCMC iterations")
# 
# ggsave(path = PATH_IMG,
#        filename = "traceplot_fit_cubic_example.pdf",
#        plot = traceplot_fit,
#        width = 12, height = 6)
# 
# 
# 
# ## Trigonometric ----
# I = input_final[[9]][[1]]
# 
# # Plot DR basis ----
# pdf(paste0(PATH_IMG, "DR_basis_trigo_example.pdf"), width = 12, height = 8)
# names_X = c("Intercept", "Linear", paste0("DR basis ", 1:20))
# par(mfrow = c(4, 6), mar = c(3, 2, 2, 1))
# for (j in 1:22){
#   plot(I$x1[order(I$x1)], I$X_scaled[order(I$x1), j], type = "l",
#        main = names_X[j],
#        xlab = "", ylab = "0")
# }
# dev.off()
# 
# 
# # Traceplots coefficients
# pdf(paste0(PATH_IMG, "traceplot_trigo_example.pdf"), width = 12, height = 8)
# par(mfrow = c(5, 6), mar = c(3, 2, 2, 1))
# for (j in 1:22){
#   plot(I$draws[ , j], type = "l", 
#        main = colnames(I$draws)[j], xlab = "MCMC Iteration", ylab = "")
# }
# plot.new(); plot.new()
# err = colMeans(I$X_scaled %*% t(I$draws) - matrix(I$y, nrow = 1000, ncol = 1000, byrow = F))
# plot(err, type = "l", main = "y - fit", xlab = "MCMC Iteration", ylab = "")
# err_noInt = colMeans(I$X_scaled[ , -1] %*% t(I$draws[ , -1]) - matrix(I$y, nrow = 1000, ncol = 1000, byrow = F))
# plot(err_noInt, type = "l", main = "y - fit (no beta0)", xlab = "MCMC Iteration", ylab = "")
# err_Square = colMeans(I$X_scaled[ , 22] %*% t(I$draws[ , 22]) - matrix(I$y, nrow = 1000, ncol = 1000, byrow = F))
# plot(err_Square, type = "l", main = "y - fit (only square)", xlab = "MCMC Iteration", ylab = "")
# dev.off()
# 
# # "Traceplot" reconstruction function
# df_complete = cbind(1:1e3, scale(I$x1), I$X_scaled %*% t(I$draws)) %>% 
#   as.data.frame() %>% 
#   pivot_longer(cols = -(V1:V2), names_to = "Iter", values_to = "x^2") %>%
#   rename(Unit = V1, x = V2) %>% 
#   mutate(Iter = as.integer(str_remove(Iter, "V"))-2)
# 
# df_noInt = cbind(1:1e3, scale(I$x1), I$X_scaled[ , -1] %*% t(I$draws[, -1])) %>% 
#   as.data.frame() %>% 
#   pivot_longer(cols = -(V1:V2), names_to = "Iter", values_to = "x^2") %>%
#   rename(Unit = V1, x = V2) %>% 
#   mutate(Iter = as.integer(str_remove(Iter, "V"))-2)
# 
# trigo = function(x) sin(2.5*x) + 3*exp(-(5^2)*((x)^2))
# 
# df_Partial = cbind(1:1e3, scale(I$x1), I$X_scaled[ , c(1:2, 21:22)] %*% t(I$draws[, c(1:2, 21:22)])) %>% 
#   as.data.frame() %>% 
#   pivot_longer(cols = -(V1:V2), names_to = "Iter", values_to = "x^2") %>%
#   rename(Unit = V1, x = V2) %>% 
#   mutate(Iter = as.integer(str_remove(Iter, "V"))-2)
# 
# 
# traceplot_fit = bind_rows(
#   df_noInt %>% mutate(grp = "No Int"),
#   df_complete %>% mutate(grp = "Complete"),
#   df_Partial %>% mutate(grp = "Lin + DR19 + DR20 + Int")) %>% 
#   mutate(grp = factor(grp, levels = c("Complete", "Lin + DR19 + DR20 + Int", "No Int"))) %>% 
#   ggplot() +
#   geom_line(aes(x = x, y = `x^2`, group = Iter), alpha = 0.25) + 
#   geom_line(data = cbind.data.frame(x = scale(I$x1), transf = trigo(I$x1)),
#             mapping = aes(x = x, y = transf), col = "red") + 
#   facet_wrap(~grp) +
#   labs(title = "Predictor during MCMC iterations")
# 
# ggsave(path = PATH_IMG,
#        filename = "traceplot_fit_trigo_example.pdf",
#        plot = traceplot_fit,
#        width = 12, height = 6)



rm(list = ls())

suppressPackageStartupMessages(library(tidyverse))
library(reshape2)
library(gridExtra)

input_init = readRDS("sim_study_nonC-GBI/splines_univ/RWMH/sim_study_nonCalibrated_asymm_TL_init.RDS")

summ_init = lapply(input_init,
                   function(I)
                     abind::abind(lapply(I, function(x) x$summ), along = 3))

PATH_IMG = "sim_study_nonC-GBI/splines_univ/RWMH/"

# TRUE BETA AND SAMPLE SIZE ----
n = 1e3


# SIMULATION SETTINGS ----
error_distr = c("Gaussian", "Gamma")
cov_distr = c("Unif_rad2", "Beta_2_4") #"Gamma_1.5_1.5") #, "BVN")
transform_x = c("parabola", "cubic", "trigonometric")
SNR = c(2, 1)
sim_settings = as.matrix(expand.grid(error_distr, cov_distr, transform_x, SNR))
colnames(sim_settings) = c("err_distr", "cov_distr", "transf_x", "SNR")

# COLORS ----
col_init = "#588157"


# FIT vs TRUE FUNCTION ----
## Fit ----
fit_init = array(NA, c(nrow(sim_settings), n, 100),
                 dimnames = list("setting" = 1:24, "unit" = 1:n, "replicate" = 1:100))

for (s in 1:nrow(sim_settings)){
  
  cat(s, "\t")
  
  I_i = input_init[[s]]
  
  fit_i = sapply(I_i, function(I) rowMeans(I$X_scaled %*% t(I$draws)))
  
  fit_init[s, , ] = fit_i
}

## X ----
X_init = array(NA, c(nrow(sim_settings), n, 100),
               dimnames = list("setting" = 1:24, "unit" = 1:n, "replicate" = 1:100))

for (s in 1:nrow(sim_settings)){
  
  cat(s, "\t")
  
  I_i = input_init[[s]]
  
  X_i = sapply(I_i, function(I) I$x1, simplify = "array")
  
  X_init[s, , ] = X_i
}


## True functions ----
true_f = data.frame(setting = integer(0L), x = numeric(0L), y = numeric(0L))
for (s in 1:nrow(sim_settings)){
  xmin = min(X_init[s, , ], na.rm = T)
  xmax = max(X_init[s, , ], na.rm = T)
  
  transf = switch(sim_settings[s, "transf_x"],
                  "parabola" = function(x) x^2,
                  "cubic" = function(x) x^3,
                  "trigonometric" = function(x) sin(2.5*x) + 2*exp(-(5^2)*((x)^2)))
  
  x = seq(xmin, xmax, length = 100)
  y = transf(x)
  
  true_f = bind_rows(true_f, data.frame(setting = s, x = x, y = y))
}
true_f = left_join(true_f, cbind.data.frame(setting = 1:24, sim_settings), by = "setting")


df_fit_init = left_join(fit_init %>% melt(value.name = "fit"),
                        X_init %>% melt(value.name = "x"),
                        by = c("setting", "unit", "replicate")) %>% 
  left_join(cbind.data.frame(setting = 1:24, sim_settings), by = "setting")

plt_fit_init2 = df_fit_init %>% 
  filter(SNR == 2) %>% 
  ggplot() +
  geom_line(aes(x = x, y = fit, group = replicate), alpha = 0.15) +
  geom_line(data = true_f %>% filter(SNR == 2),
            mapping = aes(x = x, y = y), col = "red") +
  facet_wrap(~setting, scales = "free", labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  labs(title = "Initial - SNR = 2")

plt_fit_init1 = df_fit_init %>% 
  filter(SNR == 1) %>% 
  ggplot() +
  geom_line(aes(x = x, y = fit, group = replicate), alpha = 0.15) +
  geom_line(data = true_f %>% filter(SNR == 1),
            mapping = aes(x = x, y = y), col = "red") +
  facet_wrap(~setting, scales = "free", labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  labs(title = "Initial - SNR = 1")

ggsave(filename = "fit.pdf",
       plot = marrangeGrob(list(plt_fit_init2, plt_fit_init1), nrow = 1, ncol = 1, top = NULL),
       path = PATH_IMG, width = 12, height = 9)

# DENSITY RESIDUALS ----
y_init = array(NA, c(nrow(sim_settings), n, 100),
               dimnames = list("setting" = 1:24, "unit" = 1:n, "replicate" = 1:100))

for (s in 1:nrow(sim_settings)){
  
  cat(s, "\t")
  
  I_i = input_init[[s]]
  
  y_i = sapply(I_i, function(I) I$y, simplify = "array")
  
  y_init[s, , ] = y_i
}

df_y_init = left_join(y_init %>% melt(value.name = "y"),
                      fit_init %>% melt(value.name = "fit"),
                      by = c("setting", "unit", "replicate")) %>% 
  left_join(cbind.data.frame(setting = 1:24, sim_settings), by = "setting")

plt_res_init2 = df_y_init %>%
  filter(SNR == 2) %>% 
  ggplot() +
  geom_density(aes(x = y - fit, group = replicate), color = alpha("black", alpha = 0.2)) +
  facet_wrap(~setting, scales = "free", labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  labs(title = "Initial - SNR = 2")

plt_res_init1 = df_y_init %>%
  filter(SNR == 1) %>%
  ggplot() +
  geom_density(aes(x = y - fit, group = replicate), color = alpha("black", alpha = 0.2)) +
  facet_wrap(~setting, scales = "free", labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  labs(title = "Initial - SNR = 1")

ggsave(filename = "residuals.pdf",
       plot = marrangeGrob(list(plt_res_init2, plt_res_init1), nrow = 1, ncol = 1),
       path = PATH_IMG, width = 12, height = 9)


# DIAGNOSTICS ----
## Traceplots coefficients ----
pdf(paste0(PATH_IMG, "traceplot_examples.pdf"), width = 12, height = 8)
for (s in 1:nrow(sim_settings)){
  I = input_init[[s]][[1]]
  
  # Add outer margin at the top for the page title
  par(mfrow = c(4, 6), mar = c(3, 2, 2, 1), oma = c(0, 0, 3, 0))
  
  for (j in 1:22){
    plot(I$draws[ , j], type = "l", 
         main = colnames(I$draws)[j], xlab = "MCMC Iteration", ylab = "")
    abline(h = 0, col = "red")
  }
  plot(I$draws_tau[ , 1], type = "l", main = "tau_lin")
  plot(I$draws_tau[ , 2], type = "l", main = "tau_nonlin")
  # err = colMeans(I$X_scaled %*% t(I$draws) - matrix(I$y, nrow = 1000, ncol = 1000, byrow = F))
  # plot(err, type = "l", main = "y - fit", xlab = "MCMC Iteration", ylab = "")
  # err_noInt = colMeans(I$X_scaled[ , -1] %*% t(I$draws[ , -1]) - matrix(I$y, nrow = 1000, ncol = 1000, byrow = F))
  # plot(err_noInt, type = "l", main = "y - fit (no beta0)", xlab = "MCMC Iteration", ylab = "")
  
  # Page title — written after all plots so it spans the full page
  page_title <- paste(paste0(sim_settings[s, ], collapse = " - "))
  mtext(page_title, outer = TRUE, cex = 1.2, font = 2, line = 1)
}
dev.off()
