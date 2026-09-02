rm(list = ls())

suppressPackageStartupMessages(library(tidyverse))
library(reshape2)
library(gridExtra)
source("src/01-utils.R")

input1000homo = readRDS("sim_study_nonC-GBI/splines_univ/ss_p1_n1000_homo.RDS")
input1000hetero = readRDS("sim_study_nonC-GBI/splines_univ/ss_p1_n1000_hetero.RDS")
input500homo = readRDS("sim_study_nonC-GBI/splines_univ/ss_p1_n500_homo.RDS")
input500hetero = readRDS("sim_study_nonC-GBI/splines_univ/ss_p1_n500_hetero.RDS")
input250homo = readRDS("sim_study_nonC-GBI/splines_univ/ss_p1_n250_homo.RDS")
input250hetero = readRDS("sim_study_nonC-GBI/splines_univ/ss_p1_n250_hetero.RDS")

summ1000homo = lapply(input1000homo$res,
                      function(I)
                        abind::abind(lapply(I, function(x) x$summ), along = 3))
summ1000hetero = lapply(input1000hetero$res,
                        function(I)
                          abind::abind(lapply(I, function(x) x$summ), along = 3))
summ500homo = lapply(input500homo$res,
                     function(I)
                       abind::abind(lapply(I, function(x) x$summ), along = 3))
summ500hetero = lapply(input500hetero$res,
                       function(I)
                         abind::abind(lapply(I, function(x) x$summ), along = 3))
summ250homo = lapply(input250homo$res,
                     function(I)
                       abind::abind(lapply(I, function(x) x$summ), along = 3))
summ250hetero = lapply(input250hetero$res,
                       function(I)
                         abind::abind(lapply(I, function(x) x$summ), along = 3))


PATH_IMG = "sim_study_nonC-GBI/splines_univ/"


# SIMULATION SETTINGS ----
error_distr = c("Gaussian", "Gamma")
cov_distr = c("Unif_rad2", "Beta_2_4")
transform_x = c("parabola", "cubic", "trigonometric")
SNR = c("High", "Low")
sim_settings = as.matrix(expand.grid(error_distr, cov_distr, transform_x, SNR))
colnames(sim_settings) = c("err_distr", "cov_distr", "transf_x", "SNR")

sim_settings_names = sim_settings
sim_settings_names[sim_settings_names[ , "transf_x"] == "trigonometric", "transf_x"] = "sin+exp"
sim_settings_names[sim_settings_names[ , "err_distr"] == "Gaussian", "err_distr"] = "Gauss"

# COLORS ----
colO1 = "#90e0ef"
colO2 = "#00b4d8"
colO3 = "#0077b6"
colE1 = "#ffba08"
colE2 = "#e85d04"
colE3 = "#d00000"
color_runs = c("Homo250" = colO1, "Homo500" = colO2, "Homo1000" = colO3,
               "Hetero250" = colE1, "Hetero500" = colE2, "Hetero1000" = colE3)


# CHECK FAILED JOBS ----
sapply(input1000homo$res, function(I) sum(sapply(I, function(x) is.null(x$draws))))
sapply(input1000hetero$res, function(I) sum(sapply(I, function(x) is.null(x$draws))))
sapply(input500homo$res, function(I) sum(sapply(I, function(x) is.null(x$draws))))
sapply(input500hetero$res, function(I) sum(sapply(I, function(x) is.null(x$draws))))
sapply(input250homo$res, function(I) sum(sapply(I, function(x) is.null(x$draws))))
sapply(input250hetero$res, function(I) sum(sapply(I, function(x) is.null(x$draws))))

# FIT vs TRUE FUNCTION vs OBSERVED Y----
## Fit ----
fit1000homo = fit1000hetero = array(NA, c(nrow(sim_settings), 1e3, 100),
                                    dimnames = list("setting" = 1:24, "unit" = 1:1e3, "replicate" = 1:100))
fit500homo = fit500hetero = array(NA, c(nrow(sim_settings), 5e2, 100),
                                  dimnames = list("setting" = 1:24, "unit" = 1:5e2, "replicate" = 1:100))
fit250homo = fit250hetero = array(NA, c(nrow(sim_settings), 2.5e2, 100),
                                  dimnames = list("setting" = 1:24, "unit" = 1:2.5e2, "replicate" = 1:100))


for (s in 1:nrow(sim_settings)){
  
  cat(s, "\t")
  
  I1000homo = input1000homo$res[[s]]
  I1000hetero = input1000hetero$res[[s]]
  I500homo = input500homo$res[[s]]
  I500hetero = input500hetero$res[[s]]
  I250homo = input250homo$res[[s]]
  I250hetero = input250hetero$res[[s]]
  
  fit1000homo_tmp = sapply(I1000homo, function(I) I$X_scaled %*% rowMeans(t(I$draws)))
  fit1000hetero_tmp = sapply(I1000hetero, function(I) I$X_scaled %*% rowMeans(t(I$draws)))
  fit500homo_tmp = sapply(I500homo, function(I) I$X_scaled %*% rowMeans(t(I$draws)))
  fit500hetero_tmp = sapply(I500hetero, function(I) I$X_scaled %*% rowMeans(t(I$draws)))
  fit250homo_tmp = sapply(I250homo, function(I) I$X_scaled %*% rowMeans(t(I$draws)))
  fit250hetero_tmp = sapply(I250hetero, function(I) I$X_scaled %*% rowMeans(t(I$draws)))
  
  fit1000homo[s, , 1:ncol(fit1000homo_tmp)] = fit1000homo_tmp
  fit1000hetero[s, , 1:ncol(fit1000hetero_tmp)] = fit1000hetero_tmp
  fit500homo[s, , 1:ncol(fit500homo_tmp)] = fit500homo_tmp
  fit500hetero[s, , 1:ncol(fit500hetero_tmp)] = fit500hetero_tmp
  fit250homo[s, , 1:ncol(fit250homo_tmp)] = fit250homo_tmp
  fit250hetero[s, , 1:ncol(fit250hetero_tmp)] = fit250hetero_tmp
}

## X ----
X1000homo = X1000hetero = array(NA, c(nrow(sim_settings), 1e3, 100),
                                dimnames = list("setting" = 1:24, "unit" = 1:1e3, "replicate" = 1:100))
X500homo = X500hetero = array(NA, c(nrow(sim_settings), 5e2, 100),
                              dimnames = list("setting" = 1:24, "unit" = 1:5e2, "replicate" = 1:100))
X250homo = X250hetero =array(NA, c(nrow(sim_settings), 2.5e2, 100),
                             dimnames = list("setting" = 1:24, "unit" = 1:2.5e2, "replicate" = 1:100))


for (s in 1:nrow(sim_settings)){
  
  cat(s, "\t")
  
  I1000homo = input1000homo$res[[s]]
  I1000hetero = input1000hetero$res[[s]]
  I500homo = input500homo$res[[s]]
  I500hetero = input500hetero$res[[s]]
  I250homo = input250homo$res[[s]]
  I250hetero = input250hetero$res[[s]]
  
  X1000homo_tmp = sapply(I1000homo, function(I) I$x1, simplify = "array")
  X1000hetero_tmp = sapply(I1000hetero, function(I) I$x1, simplify = "array")
  X500homo_tmp = sapply(I500homo, function(I) I$x1, simplify = "array")
  X500hetero_tmp = sapply(I500hetero, function(I) I$x1, simplify = "array")
  X250homo_tmp = sapply(I250homo, function(I) I$x1, simplify = "array")
  X250hetero_tmp = sapply(I250hetero, function(I) I$x1, simplify = "array")
  
  X1000homo[s, , 1:ncol(X1000homo_tmp)] = X1000homo_tmp
  X1000hetero[s, , 1:ncol(X1000hetero_tmp)] = X1000hetero_tmp
  X500homo[s, , 1:ncol(X500homo_tmp)] = X500homo_tmp
  X500hetero[s, , 1:ncol(X500hetero_tmp)] = X500hetero_tmp
  X250homo[s, , 1:ncol(X250homo_tmp)] = X250homo_tmp
  X250hetero[s, , 1:ncol(X250hetero_tmp)] = X250hetero_tmp
}


## True ----
true_f = data.frame(setting = integer(0L), x = numeric(0L), y = numeric(0L))
xmin = min(c(X1000homo, X1000hetero, X500homo, X500hetero, X250homo, X250hetero), na.rm = T)
xmax = max(c(X1000homo, X1000hetero, X500homo, X500hetero, X250homo, X250hetero), na.rm = T)
new_x = seq(xmin, xmax, length = 501)

for (s in 1:nrow(sim_settings)){
  
  transf = switch(sim_settings[s, "transf_x"],
                  "parabola" = function(x) x^2,
                  "cubic" = function(x) x^3,
                  "trigonometric" = function(x) sin(2.5*x) + 3*exp(-(5^2)*((x)^2)))
  
  y = transf(new_x)
  
  true_f = bind_rows(true_f, data.frame(setting = s, x = new_x, y = y))
}
true_f = left_join(true_f, cbind.data.frame(setting = 1:24, sim_settings), by = "setting")


## Observed y ----
y1000homo = y1000hetero = array(NA, c(nrow(sim_settings), 1e3, 100),
                                dimnames = list("setting" = 1:24, "unit" = 1:1e3, "replicate" = 1:100))
y500homo = y500hetero = array(NA, c(nrow(sim_settings), 5e2, 100),
                              dimnames = list("setting" = 1:24, "unit" = 1:5e2, "replicate" = 1:100))
y250homo = y250hetero = array(NA, c(nrow(sim_settings), 2.5e2, 100),
                              dimnames = list("setting" = 1:24, "unit" = 1:2.5e2, "replicate" = 1:100))


for (s in 1:nrow(sim_settings)){
  
  cat(s, "\t")
  
  I1000homo = input1000homo$res[[s]]
  I1000hetero = input1000hetero$res[[s]]
  I500homo = input500homo$res[[s]]
  I500hetero = input500hetero$res[[s]]
  I250homo = input250homo$res[[s]]
  I250hetero = input250hetero$res[[s]]
  
  y1000homo_tmp = sapply(I1000homo, function(I) I$y, simplify = "array")
  y1000hetero_tmp = sapply(I1000hetero, function(I) I$y, simplify = "array")
  y500homo_tmp = sapply(I500homo, function(I) I$y, simplify = "array")
  y500hetero_tmp = sapply(I500hetero, function(I) I$y, simplify = "array")
  y250homo_tmp = sapply(I250homo, function(I) I$y, simplify = "array")
  y250hetero_tmp = sapply(I250hetero, function(I) I$y, simplify = "array")
  
  y1000homo[s, , 1:ncol(X1000homo_tmp)] = y1000homo_tmp
  y1000hetero[s, , 1:ncol(X1000hetero_tmp)] = y1000hetero_tmp
  y500homo[s, , 1:ncol(X500homo_tmp)] = y500homo_tmp
  y500hetero[s, , 1:ncol(X500hetero_tmp)] = y500hetero_tmp
  y250homo[s, , 1:ncol(X250homo_tmp)] = y250homo_tmp
  y250hetero[s, , 1:ncol(X250hetero_tmp)] = y250hetero_tmp
}

## df (x+fit+y) ----
df_fit1000homo = left_join(fit1000homo %>% melt(value.name = "fit"),
                           X1000homo %>% melt(value.name = "x"),
                           by = c("setting", "unit", "replicate")) %>% 
  left_join(y1000homo %>% melt(value.name = "y"),
            by = c("setting", "unit", "replicate")) %>% 
  left_join(cbind.data.frame(setting = 1:24, sim_settings), by = "setting")
df_fit1000hetero = left_join(fit1000hetero %>% melt(value.name = "fit"),
                             X1000hetero %>% melt(value.name = "x"),
                             by = c("setting", "unit", "replicate")) %>% 
  left_join(y1000hetero %>% melt(value.name = "y"),
            by = c("setting", "unit", "replicate")) %>%
  left_join(cbind.data.frame(setting = 1:24, sim_settings), by = "setting")

df_fit500homo = left_join(fit500homo %>% melt(value.name = "fit"),
                          X500homo %>% melt(value.name = "x"),
                          y500homo %>% melt(value.name = "y"),
                          by = c("setting", "unit", "replicate")) %>%
  left_join(y500homo %>% melt(value.name = "y"),
            by = c("setting", "unit", "replicate")) %>% 
  left_join(cbind.data.frame(setting = 1:24, sim_settings), by = "setting")
df_fit500hetero = left_join(fit500hetero %>% melt(value.name = "fit"),
                            X500hetero %>% melt(value.name = "x"),
                            y500hetero %>% melt(value.name = "y"),
                            by = c("setting", "unit", "replicate")) %>%
  left_join(y500hetero %>% melt(value.name = "y"),
            by = c("setting", "unit", "replicate")) %>% 
  left_join(cbind.data.frame(setting = 1:24, sim_settings), by = "setting")

df_fit250homo = left_join(fit250homo %>% melt(value.name = "fit"),
                          X250homo %>% melt(value.name = "x"),
                          by = c("setting", "unit", "replicate")) %>%
  left_join(y250homo %>% melt(value.name = "y"),
            by = c("setting", "unit", "replicate")) %>%
  left_join(cbind.data.frame(setting = 1:24, sim_settings), by = "setting")
df_fit250hetero = left_join(fit250hetero %>% melt(value.name = "fit"),
                            X250hetero %>% melt(value.name = "x"),
                            by = c("setting", "unit", "replicate")) %>%
  left_join(y250hetero %>% melt(value.name = "y"),
            by = c("setting", "unit", "replicate")) %>%
  left_join(cbind.data.frame(setting = 1:24, sim_settings), by = "setting")

## y limits plots ----
y_limits <- list(
  scale_y_continuous(limits = c(-0.25, 2.5)), scale_y_continuous(limits = c(-0.25, 2.5)),
  scale_y_continuous(limits = c(-0.5, 3)),  scale_y_continuous(limits = c(-0.5, 3)),
  scale_y_continuous(limits = c(-3, 3)),  scale_y_continuous(limits = c(-3, 3)),
  scale_y_continuous(limits = c(-3, 3)),  scale_y_continuous(limits = c(-3, 3)),
  scale_y_continuous(limits = c(-1.5, 3.5)),  scale_y_continuous(limits = c(-1.5, 3.5)),
  scale_y_continuous(limits = c(-1.5, 3.5)),  scale_y_continuous(limits = c(-1.5, 3.5))
)

## Plot 1000 ----
plt_fit1000_SNRH = 
  bind_rows(df_fit1000homo %>% filter(SNR == "High") %>% mutate(err_var = "homoscedastic"),
            df_fit1000hetero %>% filter(SNR == "High") %>% mutate(err_var = "heteroscedastic")) %>% 
  ggplot() +
  geom_line(aes(x = x, y = fit, group = interaction(err_var, replicate), colour = err_var), alpha = 0.25) +
  stat_summary_bin(aes(x = x, y = y), fun = mean, bins = 35, geom = "point", size = 0.6, alpha = 1) +
  geom_line(data = true_f %>% filter(SNR == "High"),
            mapping = aes(x = x, y = y, col = "true")) +
  facet_wrap(~setting, scales = "free_y", drop = T,
             labeller = labeller(setting = setNames(apply(sim_settings_names, 1, paste, collapse = " - "), 1:24))) +
  ggh4x::facetted_pos_scales(y = y_limits) +
  scale_color_manual("", values = c("true" = "black", "homoscedastic" = colO3, "heteroscedastic" = colE3),
                     breaks = c("true", "homoscedastic", "heteroscedastic")) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5))) +
  labs(title = "n=1000, homoscedastic error [SNR = 5]") +
  theme(legend.position = "bottom")

plt_fit1000_SNRL =
  bind_rows(df_fit1000homo %>% filter(SNR == "Low") %>% mutate(err_var = "homoscedastic"),
            df_fit1000hetero %>% filter(SNR == "Low") %>% mutate(err_var = "heteroscedastic")) %>%
  ggplot() +
  geom_line(aes(x = x, y = fit, group = interaction(err_var, replicate), colour = err_var), alpha = 0.25) +
  stat_summary_bin(aes(x = x, y = y), fun = mean, bins = 35, geom = "point", size = 0.6, alpha = 1) +
  geom_line(data = true_f %>% filter(SNR == "Low"),
            mapping = aes(x = x, y = y, col = "true")) +
  facet_wrap(~setting, scales = "free_y", drop = T,
             labeller = labeller(setting = setNames(apply(sim_settings_names, 1, paste, collapse = " - "), 1:24))) +
  ggh4x::facetted_pos_scales(y = y_limits) +
  scale_color_manual("", values = c("true" = "black", "homoscedastic" = colO3, "heteroscedastic" = colE3),
                     breaks = c("true", "homoscedastic", "heteroscedastic")) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5))) +
  labs(title = "n=1000, homoscedastic error [SNR = 2]") +
  theme(legend.position = "bottom")


## Plot 500 ----
plt_fit500_SNRH = 
  bind_rows(df_fit500homo %>% filter(SNR == "High") %>% mutate(err_var = "homoscedastic"),
            df_fit500hetero %>% filter(SNR == "High") %>% mutate(err_var = "heteroscedastic")) %>% 
  ggplot() +
  geom_line(aes(x = x, y = fit, group = interaction(err_var, replicate), colour = err_var), alpha = 0.25) +
  stat_summary_bin(aes(x = x, y = y), fun = mean, bins = 35, geom = "point", size = 0.6, alpha = 1) +
  geom_line(data = true_f %>% filter(SNR == "High"),
            mapping = aes(x = x, y = y, col = "true")) +
  facet_wrap(~setting, scales = "free_y", drop = T,
             labeller = labeller(setting = setNames(apply(sim_settings_names, 1, paste, collapse = " - "), 1:24))) +
  ggh4x::facetted_pos_scales(y = y_limits) +
  scale_color_manual("", values = c("true" = "black", "homoscedastic" = colO2, "heteroscedastic" = colE2),
                     breaks = c("true", "homoscedastic", "heteroscedastic")) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5))) +
  labs(title = "n=500, homoscedastic error [SNR = 5]") +
  theme(legend.position = "bottom")

plt_fit500_SNRL =
  bind_rows(df_fit500homo %>% filter(SNR == "Low") %>% mutate(err_var = "homoscedastic"),
            df_fit500hetero %>% filter(SNR == "Low") %>% mutate(err_var = "heteroscedastic")) %>%
  ggplot() +
  geom_line(aes(x = x, y = fit, group = interaction(err_var, replicate), colour = err_var), alpha = 0.25) +
  stat_summary_bin(aes(x = x, y = y), fun = mean, bins = 35, geom = "point", size = 0.6, alpha = 1) +
  geom_line(data = true_f %>% filter(SNR == "Low"),
            mapping = aes(x = x, y = y, col = "true")) +
  facet_wrap(~setting, scales = "free_y", drop = T,
             labeller = labeller(setting = setNames(apply(sim_settings_names, 1, paste, collapse = " - "), 1:24))) +
  ggh4x::facetted_pos_scales(y = y_limits) +
  scale_color_manual("", values = c("true" = "black", "homoscedastic" = colO2, "heteroscedastic" = colE2),
                     breaks = c("true", "homoscedastic", "heteroscedastic")) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5))) +
  labs(title = "n=500, homoscedastic error [SNR = 2]") +
  theme(legend.position = "bottom")


## Plot 250 ----
plt_fit250_SNRH = 
  bind_rows(df_fit250homo %>% filter(SNR == "High") %>% mutate(err_var = "homoscedastic"),
            df_fit250hetero %>% filter(SNR == "High") %>% mutate(err_var = "heteroscedastic")) %>% 
  ggplot() +
  geom_line(aes(x = x, y = fit, group = interaction(err_var, replicate), colour = err_var), alpha = 0.25) +
  stat_summary_bin(aes(x = x, y = y), fun = mean, bins = 35, geom = "point", size = 0.6, alpha = 1) +
  geom_line(data = true_f %>% filter(SNR == "High"),
            mapping = aes(x = x, y = y, col = "true")) +
  facet_wrap(~setting, scales = "free_y", drop = T,
             labeller = labeller(setting = setNames(apply(sim_settings_names, 1, paste, collapse = " - "), 1:24))) +
  ggh4x::facetted_pos_scales(y = y_limits) +
  scale_color_manual("", values = c("true" = "black", "homoscedastic" = colO1, "heteroscedastic" = colE1),
                     breaks = c("true", "homoscedastic", "heteroscedastic")) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5))) +
  labs(title = "n=250, homoscedastic error [SNR = 5]") +
  theme(legend.position = "bottom")

plt_fit250_SNRL =
  bind_rows(df_fit250homo %>% filter(SNR == "Low") %>% mutate(err_var = "homoscedastic"),
            df_fit250hetero %>% filter(SNR == "Low") %>% mutate(err_var = "heteroscedastic")) %>% 
  ggplot() +
  geom_line(aes(x = x, y = fit, group = interaction(err_var, replicate), colour = err_var), alpha = 0.25) +
  stat_summary_bin(aes(x = x, y = y), fun = mean, bins = 35, geom = "point", size = 0.6, alpha = 1) +
  geom_line(data = true_f %>% filter(SNR == "Low"),
            mapping = aes(x = x, y = y, col = "true")) +
  facet_wrap(~setting, scales = "free_y", drop = T,
             labeller = labeller(setting = setNames(apply(sim_settings_names, 1, paste, collapse = " - "), 1:24))) +
  ggh4x::facetted_pos_scales(y = y_limits) +
  scale_color_manual("", values = c("true" = "black", "homoscedastic" = colO1, "heteroscedastic" = colE1),
                     breaks = c("true", "homoscedastic", "heteroscedastic")) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5))) +
  labs(title = "n=250, homoscedastic error [SNR = 2]") +
  theme(legend.position = "bottom")


ggsave(filename = "fit.pdf",
       plot = marrangeGrob(list(
         plt_fit1000_SNRH, plt_fit1000_SNRL,
         plt_fit500_SNRH, plt_fit500_SNRL,
         plt_fit250_SNRH, plt_fit250_SNRL 
       ), 
       nrow = 1, ncol = 1, top = NULL),
       path = PATH_IMG, width = 12, height = 9)


# FIT VS LINEAR COMPONENT ----
## Fit ----
fit_lin1000homo = fit_lin1000hetero = array(NA, c(nrow(sim_settings), 1e3, 100),
                                            dimnames = list("setting" = 1:24, "unit" = 1:1e3, "replicate" = 1:100))
fit_lin500homo = fit_lin500hetero = array(NA, c(nrow(sim_settings), 5e2, 100),
                                          dimnames = list("setting" = 1:24, "unit" = 1:5e2, "replicate" = 1:100))
fit_lin250homo = fit_lin250hetero = array(NA, c(nrow(sim_settings), 2.5e2, 100),
                                          dimnames = list("setting" = 1:24, "unit" = 1:2.5e2, "replicate" = 1:100))


for (s in 1:nrow(sim_settings)){
  
  cat(s, "\t")
  
  I1000homo = input1000homo$res[[s]]
  I1000hetero = input1000hetero$res[[s]]
  I500homo = input500homo$res[[s]]
  I500hetero = input500hetero$res[[s]]
  I250homo = input250homo$res[[s]]
  I250hetero = input250hetero$res[[s]]
  
  fit_lin1000homo_tmp = sapply(I1000homo, function(I) I$X_scaled[ , 1:2] %*% rowMeans(t(I$draws[ , 1:2])))
  fit_lin1000hetero_tmp = sapply(I1000hetero, function(I) I$X_scaled[ , 1:2] %*% rowMeans(t(I$draws[ , 1:2])))
  fit_lin500homo_tmp = sapply(I500homo, function(I) I$X_scaled[ , 1:2] %*% rowMeans(t(I$draws[ , 1:2])))
  fit_lin500hetero_tmp = sapply(I500hetero, function(I) I$X_scaled[ , 1:2] %*% rowMeans(t(I$draws[ , 1:2])))
  fit_lin250homo_tmp = sapply(I250homo, function(I) I$X_scaled[ , 1:2] %*% rowMeans(t(I$draws[ , 1:2])))
  fit_lin250hetero_tmp = sapply(I250hetero, function(I) I$X_scaled[ , 1:2] %*% rowMeans(t(I$draws[ , 1:2])))
  
  fit_lin1000homo[s, , 1:ncol(fit_lin1000homo_tmp)] = fit_lin1000homo_tmp
  fit_lin1000hetero[s, , 1:ncol(fit_lin1000hetero_tmp)] = fit_lin1000hetero_tmp
  fit_lin500homo[s, , 1:ncol(fit_lin500homo_tmp)] = fit_lin500homo_tmp
  fit_lin500hetero[s, , 1:ncol(fit_lin500hetero_tmp)] = fit_lin500hetero_tmp
  fit_lin250homo[s, , 1:ncol(fit_lin250homo_tmp)] = fit_lin250homo_tmp
  fit_lin250hetero[s, , 1:ncol(fit_lin250hetero_tmp)] = fit_lin250hetero_tmp
}

## df (x+fit) ----
df_lin1000homo = left_join(fit_lin1000homo %>% melt(value.name = "fit"),
                           X1000homo %>% melt(value.name = "x"),
                           by = c("setting", "unit", "replicate")) %>% 
  left_join(cbind.data.frame(setting = 1:24, sim_settings), by = "setting")
df_lin1000hetero = left_join(fit_lin1000hetero %>% melt(value.name = "fit"),
                             X1000hetero %>% melt(value.name = "x"),
                             by = c("setting", "unit", "replicate")) %>% 
  left_join(cbind.data.frame(setting = 1:24, sim_settings), by = "setting")

df_lin500homo = left_join(fit_lin500homo %>% melt(value.name = "fit"),
                          X500homo %>% melt(value.name = "x"),
                          by = c("setting", "unit", "replicate")) %>%
  left_join(cbind.data.frame(setting = 1:24, sim_settings), by = "setting")
df_lin500hetero = left_join(fit_lin500hetero %>% melt(value.name = "fit"),
                            X500hetero %>% melt(value.name = "x"),
                            by = c("setting", "unit", "replicate")) %>%
  left_join(cbind.data.frame(setting = 1:24, sim_settings), by = "setting")

df_lin250homo = left_join(fit_lin250homo %>% melt(value.name = "fit"),
                          X250homo %>% melt(value.name = "x"),
                          by = c("setting", "unit", "replicate")) %>%
  left_join(cbind.data.frame(setting = 1:24, sim_settings), by = "setting")
df_lin250hetero = left_join(fit_lin250hetero %>% melt(value.name = "fit"),
                            X250hetero %>% melt(value.name = "x"),
                            by = c("setting", "unit", "replicate")) %>%
  left_join(cbind.data.frame(setting = 1:24, sim_settings), by = "setting")


## True ----
true_lin = data.frame(setting = integer(0L), x = numeric(0L), y = numeric(0L))

for (s in 1:nrow(sim_settings)){
  
  s_s = sim_settings[s, ]
  
  if (s_s["transf_x"] == "parabola"){
    if(s_s["cov_distr"] == "Unif_rad2"){
      transf_lin = function(x) 2/3
    } else if(s_s["cov_distr"] == "Beta_2_4"){
      transf_lin = function(x) 1/7 - sqrt(2)/2*x
    }
  } else if (s_s["transf_x"] == "cubic"){
    if(s_s["cov_distr"] == "Unif_rad2"){
      transf_lin = function(x) 6/5*x
    } else if(s_s["cov_distr"] == "Beta_2_4"){
      transf_lin = function(x) sqrt(2)/21 + x
    }
  } else if (s_s["transf_x"] == "trigonometric"){
    if(s_s["cov_distr"] == "Unif_rad2"){
      transf_lin = function(x) 0.37599 + 0.488903*x
    } else if(s_s["cov_distr"] == "Beta_2_4"){
      transf_lin = function(x)  0.760395 + 1.52871*x
    }
  }
  
  
  y = transf_lin(new_x)
  
  true_lin = bind_rows(true_lin, data.frame(setting = s, x = new_x, y = y))
}
true_lin = left_join(true_lin, cbind.data.frame(setting = 1:24, sim_settings), by = "setting")

## y limits plots ----
y_limits <- list(
  scale_y_continuous(limits = c(0.25, 1)),  scale_y_continuous(limits = c(0.25, 1)),
  scale_y_continuous(limits = c(-1, 1.25)),  scale_y_continuous(limits = c(-1, 1.25)),
  scale_y_continuous(limits = c(-2, 2)),  scale_y_continuous(limits = c(-2, 2)),
  scale_y_continuous(limits = c(-2, 2)),  scale_y_continuous(limits = c(-2, 2)),
  scale_y_continuous(limits = c(-0.75, 1.5)), scale_y_continuous(limits = c(-0.75, 1.5)),
  scale_y_continuous(limits = c(-1.5, 3.25)),  scale_y_continuous(limits = c(-1.5, 3.25))
)

## Plot 1000 ----
plt_lin1000_SNRH = 
  bind_rows(df_lin1000homo %>% filter(SNR == "High") %>% mutate(err_var = "homoscedastic"),
            df_lin1000hetero %>% filter(SNR == "High") %>% mutate(err_var = "heteroscedastic")) %>% 
  ggplot() +
  geom_line(aes(x = x, y = fit, group = interaction(err_var, replicate), colour = err_var), alpha = 0.25) +
  geom_line(data = true_lin %>% filter(SNR == "High"),
            mapping = aes(x = x, y = y, col = "true")) +
  facet_wrap(~setting, scales = "free_y", drop = T,
             labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  ggh4x::facetted_pos_scales(y = y_limits) +
  scale_color_manual("", values = c("true" = "black", "homoscedastic" = colO3, "heteroscedastic" = colE3),
                     breaks = c("true", "homoscedastic", "heteroscedastic")) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5))) +
  labs(title = "n=1000, homoscedastic error [SNR = 5]") +
  theme(legend.position = "bottom")

plt_lin1000_SNRL =
  bind_rows(df_lin1000homo %>% filter(SNR == "Low") %>% mutate(err_var = "homoscedastic"),
            df_lin1000hetero %>% filter(SNR == "Low") %>% mutate(err_var = "heteroscedastic")) %>%
  ggplot() +
  geom_line(aes(x = x, y = fit, group = interaction(err_var, replicate), colour = err_var), alpha = 0.25) +
  geom_line(data = true_lin %>% filter(SNR == "Low"),
            mapping = aes(x = x, y = y, col = "true")) +
  facet_wrap(~setting, scales = "free_y", drop = T,
             labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  ggh4x::facetted_pos_scales(y = y_limits) +
  scale_color_manual("", values = c("true" = "black", "homoscedastic" = colO3, "heteroscedastic" = colE3),
                     breaks = c("true", "homoscedastic", "heteroscedastic")) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5))) +
  labs(title = "n=1000, homoscedastic error [SNR = 2]") +
  theme(legend.position = "bottom")


## Plot 500 ----
plt_lin500_SNRH = 
  bind_rows(df_lin500homo %>% filter(SNR == "High") %>% mutate(err_var = "homoscedastic"),
            df_lin500hetero %>% filter(SNR == "High") %>% mutate(err_var = "heteroscedastic")) %>% 
  ggplot() +
  geom_line(aes(x = x, y = fit, group = interaction(err_var, replicate), colour = err_var), alpha = 0.25) +
  geom_line(data = true_lin %>% filter(SNR == "High"),
            mapping = aes(x = x, y = y, col = "true")) +
  facet_wrap(~setting, scales = "free_y", drop = T,
             labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  ggh4x::facetted_pos_scales(y = y_limits) +
  scale_color_manual("", values = c("true" = "black", "homoscedastic" = colO2, "heteroscedastic" = colE2),
                     breaks = c("true", "homoscedastic", "heteroscedastic")) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5))) +
  labs(title = "n=500, homoscedastic error [SNR = 5]") +
  theme(legend.position = "bottom")

plt_lin500_SNRL =
  bind_rows(df_lin500homo %>% filter(SNR == "Low") %>% mutate(err_var = "homoscedastic"),
            df_lin500hetero %>% filter(SNR == "Low") %>% mutate(err_var = "heteroscedastic")) %>%
  ggplot() +
  geom_line(aes(x = x, y = fit, group = interaction(err_var, replicate), colour = err_var), alpha = 0.25) +
  geom_line(data = true_lin %>% filter(SNR == "Low"),
            mapping = aes(x = x, y = y, col = "true")) +
  facet_wrap(~setting, scales = "free_y", drop = T,
             labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  ggh4x::facetted_pos_scales(y = y_limits) +
  scale_color_manual("", values = c("true" = "black", "homoscedastic" = colO2, "heteroscedastic" = colE2),
                     breaks = c("true", "homoscedastic", "heteroscedastic")) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5))) +
  labs(title = "n=500, homoscedastic error [SNR = 2]") +
  theme(legend.position = "bottom")


## Plot 250 ----
plt_lin250_SNRH = 
  bind_rows(df_lin250homo %>% filter(SNR == "High") %>% mutate(err_var = "homoscedastic"),
            df_lin250hetero %>% filter(SNR == "High") %>% mutate(err_var = "heteroscedastic")) %>% 
  ggplot() +
  geom_line(aes(x = x, y = fit, group = interaction(err_var, replicate), colour = err_var), alpha = 0.25) +
  geom_line(data = true_lin %>% filter(SNR == "High"),
            mapping = aes(x = x, y = y, col = "true")) +
  facet_wrap(~setting, scales = "free_y", drop = T,
             labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  ggh4x::facetted_pos_scales(y = y_limits) +
  scale_color_manual("", values = c("true" = "black", "homoscedastic" = colO1, "heteroscedastic" = colE1),
                     breaks = c("true", "homoscedastic", "heteroscedastic")) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5))) +
  labs(title = "n=250, homoscedastic error [SNR = 5]") +
  theme(legend.position = "bottom")

plt_lin250_SNRL =
  bind_rows(df_lin250homo %>% filter(SNR == "Low") %>% mutate(err_var = "homoscedastic"),
            df_lin250hetero %>% filter(SNR == "Low") %>% mutate(err_var = "heteroscedastic")) %>% 
  ggplot() +
  geom_line(aes(x = x, y = fit, group = interaction(err_var, replicate), colour = err_var), alpha = 0.25) +
  geom_line(data = true_lin %>% filter(SNR == "Low"),
            mapping = aes(x = x, y = y, col = "true")) +
  facet_wrap(~setting, scales = "free_y", drop = T,
             labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  ggh4x::facetted_pos_scales(y = y_limits) +
  scale_color_manual("", values = c("true" = "black", "homoscedastic" = colO1, "heteroscedastic" = colE1),
                     breaks = c("true", "homoscedastic", "heteroscedastic")) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5))) +
  labs(title = "n=250, homoscedastic error [SNR = 2]") +
  theme(legend.position = "bottom")


ggsave(filename = "fit_lin.pdf",
       plot = marrangeGrob(list(
         plt_lin1000_SNRH, plt_lin1000_SNRL,
         plt_lin500_SNRH, plt_lin500_SNRL,
         plt_lin250_SNRH, plt_lin250_SNRL 
       ), 
       nrow = 1, ncol = 1, top = NULL),
       path = PATH_IMG, width = 12, height = 9)


# FIT VS NON-LINEAR COMPONENT ----
## Fit ----
fit_nonlin1000homo = fit1000homo - fit_lin1000homo
fit_nonlin1000hetero = fit1000hetero - fit_lin1000hetero
fit_nonlin500homo = fit500homo - fit_lin500homo
fit_nonlin500hetero = fit500hetero - fit_lin500hetero
fit_nonlin250homo = fit250homo - fit_lin250homo
fit_nonlin250hetero = fit250hetero - fit_lin250hetero

## df(x + fit) ----
df_nonlin1000homo = left_join(fit_nonlin1000homo %>% melt(value.name = "fit"),
                              X1000homo %>% melt(value.name = "x"),
                              by = c("setting", "unit", "replicate")) %>% 
  left_join(cbind.data.frame(setting = 1:24, sim_settings), by = "setting")
df_nonlin1000hetero = left_join(fit_nonlin1000hetero %>% melt(value.name = "fit"),
                                X1000hetero %>% melt(value.name = "x"),
                                by = c("setting", "unit", "replicate")) %>%
  left_join(cbind.data.frame(setting = 1:24, sim_settings), by = "setting")

df_nonlin500homo = left_join(fit_nonlin500homo %>% melt(value.name = "fit"),
                             X500homo %>% melt(value.name = "x"),
                             by = c("setting", "unit", "replicate")) %>%
  left_join(cbind.data.frame(setting = 1:24, sim_settings), by = "setting")
df_nonlin500hetero = left_join(fit_nonlin500hetero %>% melt(value.name = "fit"),
                               X500hetero %>% melt(value.name = "x"),
                               by = c("setting", "unit", "replicate")) %>%
  left_join(cbind.data.frame(setting = 1:24, sim_settings), by = "setting")

df_nonlin250homo = left_join(fit_nonlin250homo %>% melt(value.name = "fit"),
                             X250homo %>% melt(value.name = "x"),
                             by = c("setting", "unit", "replicate")) %>%
  left_join(cbind.data.frame(setting = 1:24, sim_settings), by = "setting")
df_nonlin250hetero = left_join(fit_nonlin250hetero %>% melt(value.name = "fit"),
                               X250hetero %>% melt(value.name = "x"),
                               by = c("setting", "unit", "replicate")) %>%
  left_join(cbind.data.frame(setting = 1:24, sim_settings), by = "setting")

## True ----
true_nonlin = true_f
true_nonlin$y = true_f$y - true_lin$y

## y limits plots ----
y_limits <- list(
  scale_y_continuous(limits = c(-0.75, 2)),  scale_y_continuous(limits = c(-0.75, 2)),
  scale_y_continuous(limits = c(-0.5, 3)),  scale_y_continuous(limits = c(-0.5, 3)),
  scale_y_continuous(limits = c(-2, 2)),  scale_y_continuous(limits = c(-2, 2)),
  scale_y_continuous(limits = c(-2, 2)),  scale_y_continuous(limits = c(-2, 2)),
  scale_y_continuous(limits = c(-2, 3)),  scale_y_continuous(limits = c(-2, 3)),
  scale_y_continuous(limits = c(-2, 3)),  scale_y_continuous(limits = c(-2, 3))
)

### Plot 1000 ----
plt_nonlin1000_SNRH = 
  bind_rows(df_nonlin1000homo %>% filter(SNR == "High") %>% mutate(err_var = "homoscedastic"),
            df_nonlin1000hetero %>% filter(SNR == "High") %>% mutate(err_var = "heteroscedastic")) %>% 
  ggplot() +
  geom_line(aes(x = x, y = fit, group = interaction(err_var, replicate), colour = err_var), alpha = 0.25) +
  geom_line(data = true_nonlin %>% filter(SNR == "High"),
            mapping = aes(x = x, y = y, col = "true")) +
  facet_wrap(~setting, scales = "free_y", drop = T,
             labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  ggh4x::facetted_pos_scales(y = y_limits) +
  scale_color_manual("", values = c("true" = "black", "homoscedastic" = colO3, "heteroscedastic" = colE3),
                     breaks = c("true", "homoscedastic", "heteroscedastic")) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5))) +
  labs(title = "n=1000, homoscedastic error [SNR = 5]") +
  theme(legend.position = "bottom")

plt_nonlin1000_SNRL =
  bind_rows(df_nonlin1000homo %>% filter(SNR == "Low") %>% mutate(err_var = "homoscedastic"),
            df_nonlin1000hetero %>% filter(SNR == "Low") %>% mutate(err_var = "heteroscedastic")) %>%
  ggplot() +
  geom_line(aes(x = x, y = fit, group = interaction(err_var, replicate), colour = err_var), alpha = 0.25) +
  geom_line(data = true_nonlin %>% filter(SNR == "Low"),
            mapping = aes(x = x, y = y, col = "true")) +
  facet_wrap(~setting, scales = "free_y", drop = T,
             labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  ggh4x::facetted_pos_scales(y = y_limits) +
  scale_color_manual("", values = c("true" = "black", "homoscedastic" = colO3, "heteroscedastic" = colE3),
                     breaks = c("true", "homoscedastic", "heteroscedastic")) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5))) +
  labs(title = "n=1000, homoscedastic error [SNR = 2]") +
  theme(legend.position = "bottom")


### Plot 500 ----
plt_nonlin500_SNRH = 
  bind_rows(df_nonlin500homo %>% filter(SNR == "High") %>% mutate(err_var = "homoscedastic"),
            df_nonlin500hetero %>% filter(SNR == "High") %>% mutate(err_var = "heteroscedastic")) %>% 
  ggplot() +
  geom_line(aes(x = x, y = fit, group = interaction(err_var, replicate), colour = err_var), alpha = 0.25) +
  geom_line(data = true_nonlin %>% filter(SNR == "High"),
            mapping = aes(x = x, y = y, col = "true")) +
  facet_wrap(~setting, scales = "free_y", drop = T,
             labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  ggh4x::facetted_pos_scales(y = y_limits) +
  scale_color_manual("", values = c("true" = "black", "homoscedastic" = colO2, "heteroscedastic" = colE2),
                     breaks = c("true", "homoscedastic", "heteroscedastic")) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5))) +
  labs(title = "n=500, homoscedastic error [SNR = 5]") +
  theme(legend.position = "bottom")

plt_nonlin500_SNRL =
  bind_rows(df_nonlin500homo %>% filter(SNR == "Low") %>% mutate(err_var = "homoscedastic"),
            df_nonlin500hetero %>% filter(SNR == "Low") %>% mutate(err_var = "heteroscedastic")) %>%
  ggplot() +
  geom_line(aes(x = x, y = fit, group = interaction(err_var, replicate), colour = err_var), alpha = 0.25) +
  geom_line(data = true_nonlin %>% filter(SNR == "Low"),
            mapping = aes(x = x, y = y, col = "true")) +
  facet_wrap(~setting, scales = "free_y", drop = T,
             labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  ggh4x::facetted_pos_scales(y = y_limits) +
  scale_color_manual("", values = c("true" = "black", "homoscedastic" = colO2, "heteroscedastic" = colE2),
                     breaks = c("true", "homoscedastic", "heteroscedastic")) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5))) +
  labs(title = "n=500, homoscedastic error [SNR = 2]") +
  theme(legend.position = "bottom")


### Plot 250 ----
plt_nonlin250_SNRH = 
  bind_rows(df_nonlin250homo %>% filter(SNR == "High") %>% mutate(err_var = "homoscedastic"),
            df_nonlin250hetero %>% filter(SNR == "High") %>% mutate(err_var = "heteroscedastic")) %>% 
  ggplot() +
  geom_line(aes(x = x, y = fit, group = interaction(err_var, replicate), colour = err_var), alpha = 0.25) +
  geom_line(data = true_nonlin %>% filter(SNR == "High"),
            mapping = aes(x = x, y = y, col = "true")) +
  facet_wrap(~setting, scales = "free_y", drop = T,
             labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  ggh4x::facetted_pos_scales(y = y_limits) +
  scale_color_manual("", values = c("true" = "black", "homoscedastic" = colO1, "heteroscedastic" = colE1),
                     breaks = c("true", "homoscedastic", "heteroscedastic")) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5))) +
  labs(title = "n=250, homoscedastic error [SNR = 5]") +
  theme(legend.position = "bottom")

plt_nonlin250_SNRL =
  bind_rows(df_nonlin250homo %>% filter(SNR == "Low") %>% mutate(err_var = "homoscedastic"),
            df_nonlin250hetero %>% filter(SNR == "Low") %>% mutate(err_var = "heteroscedastic")) %>% 
  ggplot() +
  geom_line(aes(x = x, y = fit, group = interaction(err_var, replicate), colour = err_var), alpha = 0.25) +
  geom_line(data = true_nonlin %>% filter(SNR == "Low"),
            mapping = aes(x = x, y = y, col = "true")) +
  facet_wrap(~setting, scales = "free_y", drop = T,
             labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  ggh4x::facetted_pos_scales(y = y_limits) +
  scale_color_manual("", values = c("true" = "black", "homoscedastic" = colO1, "heteroscedastic" = colE1),
                     breaks = c("true", "homoscedastic", "heteroscedastic")) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5))) +
  labs(title = "n=250, homoscedastic error [SNR = 2]") +
  theme(legend.position = "bottom")


ggsave(filename = "fit_nonlin.pdf",
       plot = marrangeGrob(list(
         plt_nonlin1000_SNRH, plt_nonlin1000_SNRL,
         plt_nonlin500_SNRH, plt_nonlin500_SNRL,
         plt_nonlin250_SNRH, plt_nonlin250_SNRL 
       ), 
       nrow = 1, ncol = 1, top = NULL),
       path = PATH_IMG, width = 12, height = 9)


# RMSE AND COVERAGE ----
cov_RMSE_fun <- function(II, new_x, true_f, true_lin, true_nonlin) {
  x <- II$x1
  DR = II$DR
  
  newB = mgcv::Predict.matrix(DR$smooth_object, data.frame(x = new_x))
  
  new_Z <- newB %*% DR$Trafo
  x_scaled <- (new_x - mean(x)) / sd(x)
  new_X <- cbind(1, x_scaled, new_Z)
  
  beta <- II$draws
  
  fit <- tcrossprod(beta, new_X)
  fit_lin <- tcrossprod(beta[, 1:2, drop = FALSE], new_X[, 1:2, drop = FALSE])
  fit_nonlin = fit - fit_lin
  
  q_fit <- matrixStats::colQuantiles(fit, probs = c(0.025, 0.975), drop = FALSE)
  q_lin <- matrixStats::colQuantiles(fit_lin, probs = c(0.025, 0.975), drop = FALSE)
  q_nonlin <- matrixStats::colQuantiles(fit_nonlin, probs = c(0.025, 0.975), drop = FALSE)
  
  fit025 <- matrixStats::colQuantiles(fit, probs = 0.025, drop = TRUE)
  fit975 <- matrixStats::colQuantiles(fit, probs = 0.975, drop = TRUE)
  
  fit_lin025 <- matrixStats::colQuantiles(fit_lin, probs = 0.025, drop = TRUE)
  fit_lin975 <- matrixStats::colQuantiles(fit_lin, probs = 0.975, drop = TRUE)
  
  fit_nonlin025 <- matrixStats::colQuantiles(fit_nonlin, probs = 0.025, drop = TRUE)
  fit_nonlin975 <- matrixStats::colQuantiles(fit_nonlin, probs = 0.975, drop = TRUE)
  
  fit_mean <- colMeans(fit)
  fit_lin_mean <- colMeans(fit_lin)
  fit_nonlin_mean <- fit_mean - fit_lin_mean
  
  cbind(
    cov = as.numeric(q_fit[, 1] <= true_f & q_fit[, 2] >= true_f),
    cov_lin = as.numeric(q_lin[, 1] <= true_lin & q_lin[, 2] >= true_lin),
    cov_nonlin = as.numeric(q_nonlin[, 1] <= true_nonlin & q_nonlin[, 2] >= true_nonlin),
    SE = (fit_mean - true_f)^2,
    SE_lin = (fit_lin_mean - true_lin)^2,
    SE_nonlin = (fit_nonlin_mean - true_nonlin)^2
  )
}

summarise_replicates <- function(Ilist, new_x, true_f, true_lin, true_nonlin) {
  
  tmp <- future.apply::future_lapply(
    Ilist,
    cov_RMSE_fun,
    new_x = new_x,
    true_f = true_f,
    true_lin = true_lin,
    true_nonlin = true_nonlin,
    future.seed = FALSE
  )
  
  nrep <- length(tmp)
  
  # Sum over replicates without creating a 3D array
  sums <- Reduce(`+`, tmp)
  
  out <- list(
    cov = sums[, "cov"] / nrep,
    cov_lin = sums[, "cov_lin"] / nrep,
    cov_nonlin = sums[, "cov_nonlin"] / nrep,
    RMSE = sqrt(sums[, "SE"] / nrep),
    RMSE_lin = sqrt(sums[, "SE_lin"] / nrep),
    RMSE_nonlin = sqrt(sums[, "SE_nonlin"] / nrep)
  )
  
  out
}


npoints = 500
new_x = seq(-sqrt(2), sqrt(2), length = npoints)
cov1000homo = cov1000hetero = cov500homo = cov500hetero = cov250homo = cov250hetero =
  array(NA, dim = c(nrow(sim_settings), npoints),
        dimnames = list("setting" = 1:24, "point" = 1:npoints))
cov_lin1000homo = cov_lin1000hetero = cov_lin500homo = cov_lin500hetero = cov_lin250homo = cov_lin250hetero =
  array(NA, dim = c(nrow(sim_settings), npoints),
        dimnames = list("setting" = 1:24, "point" = 1:npoints))
cov_nonlin1000homo = cov_nonlin1000hetero = cov_nonlin500homo = cov_nonlin500hetero = cov_nonlin250homo = cov_nonlin250hetero =
  array(NA, dim = c(nrow(sim_settings), npoints),
        dimnames = list("setting" = 1:24, "point" = 1:npoints))


RMSE1000homo = RMSE1000hetero = RMSE500homo = RMSE500hetero = RMSE250homo = RMSE250hetero = 
  array(NA, dim = c(nrow(sim_settings), npoints),
        dimnames = list("setting" = 1:24, "point" = 1:npoints))
RMSE_lin1000homo = RMSE_lin1000hetero = RMSE_lin500homo = RMSE_lin500hetero = RMSE_lin250homo = RMSE_lin250hetero =
  array(NA, dim = c(nrow(sim_settings), npoints),
        dimnames = list("setting" = 1:24, "point" = 1:npoints))
RMSE_nonlin1000homo = RMSE_nonlin1000hetero = RMSE_nonlin500homo = RMSE_nonlin500hetero = RMSE_nonlin250homo = RMSE_nonlin250hetero =
  array(NA, dim = c(nrow(sim_settings), npoints), 
        dimnames = list("setting" = 1:24, "point" = 1:npoints))


library(future)
library(future.apply)
## Choose number of workers
plan(multisession, workers = 8)

for (s in 1:nrow(sim_settings)){
  cat(s, "\t")
  
  s_s = sim_settings[s, ]
  
  I1000homo = input1000homo$res[[s]]
  I1000hetero = input1000hetero$res[[s]]
  I500homo = input500homo$res[[s]]
  I500hetero = input500hetero$res[[s]]
  I250homo = input250homo$res[[s]]
  I250hetero = input250hetero$res[[s]]
  
  transf = switch(sim_settings[s, "transf_x"],
                  "parabola" = function(x) x^2,
                  "cubic" = function(x) x^3,
                  "trigonometric" = function(x) sin(2.5*x) + 3*exp(-(5^2)*((x)^2))
  )
  
  # Linear transf
  if (s_s["transf_x"] == "parabola"){
    if(s_s["cov_distr"] == "Unif_rad2"){
      transf_lin = function(x) rep(2/3, length(x))
    } else if(s_s["cov_distr"] == "Beta_2_4"){
      transf_lin = function(x) 1/7 - sqrt(2)/2*x
    }
  } else if (s_s["transf_x"] == "cubic"){
    if(s_s["cov_distr"] == "Unif_rad2"){
      transf_lin = function(x) 6/5*x
    } else if(s_s["cov_distr"] == "Beta_2_4"){
      transf_lin = function(x) sqrt(2)/21 + x
    }
  } else if (s_s["transf_x"] == "trigonometric"){
    if(s_s["cov_distr"] == "Unif_rad2"){
      transf_lin = function(x) 0.37599 + 0.488903*x
    } else if(s_s["cov_distr"] == "Beta_2_4"){
      transf_lin = function(x)  0.760395 + 1.52871*x
    }
  }
  
  true_f = transf(new_x)
  true_lin = transf_lin(new_x)
  true_nonlin = true_f - true_lin
  
  res1000homo <- summarise_replicates(I1000homo, new_x, true_f, true_lin, true_nonlin)
  res1000hetero <- summarise_replicates(I1000hetero, new_x, true_f, true_lin, true_nonlin)
  res500homo <- summarise_replicates(I500homo, new_x, true_f, true_lin, true_nonlin)
  res500hetero <- summarise_replicates(I500hetero, new_x, true_f, true_lin, true_nonlin)
  res250homo <- summarise_replicates(I250homo, new_x, true_f, true_lin, true_nonlin)
  res250hetero <- summarise_replicates(I250hetero, new_x, true_f, true_lin, true_nonlin)
  
  cov1000homo[s, ] <- res1000homo$cov
  cov_lin1000homo[s, ] <- res1000homo$cov_lin
  cov_nonlin1000homo[s, ] <- res1000homo$cov_nonlin
  RMSE1000homo[s, ] <- res1000homo$RMSE
  RMSE_lin1000homo[s, ] <- res1000homo$RMSE_lin
  RMSE_nonlin1000homo[s, ] <- res1000homo$RMSE_nonlin
  
  cov500homo[s, ] <- res500homo$cov
  cov_lin500homo[s, ] <- res500homo$cov_lin
  cov_nonlin500homo[s, ] <- res500homo$cov_nonlin
  RMSE500homo[s, ] <- res500homo$RMSE
  RMSE_lin500homo[s, ] <- res500homo$RMSE_lin
  RMSE_nonlin500homo[s, ] <- res500homo$RMSE_nonlin
  
  cov250homo[s, ] <- res250homo$cov
  cov_lin250homo[s, ] <- res250homo$cov_lin
  cov_nonlin250homo[s, ] <- res250homo$cov_nonlin
  RMSE250homo[s, ] <- res250homo$RMSE
  RMSE_lin250homo[s, ] <- res250homo$RMSE_lin
  RMSE_nonlin250homo[s, ] <- res250homo$RMSE_nonlin
  
  cov1000hetero[s, ] <- res1000hetero$cov
  cov_lin1000hetero[s, ] <- res1000hetero$cov_lin
  cov_nonlin1000hetero[s, ] <- res1000hetero$cov_nonlin
  RMSE1000hetero[s, ] <- res1000hetero$RMSE
  RMSE_lin1000hetero[s, ] <- res1000hetero$RMSE_lin
  RMSE_nonlin1000hetero[s, ] <- res1000hetero$RMSE_nonlin
  
  cov500hetero[s, ] <- res500hetero$cov
  cov_lin500hetero[s, ] <- res500hetero$cov_lin
  cov_nonlin500hetero[s, ] <- res500hetero$cov_nonlin
  RMSE500hetero[s, ] <- res500hetero$RMSE
  RMSE_lin500hetero[s, ] <- res500hetero$RMSE_lin
  RMSE_nonlin500hetero[s, ] <- res500hetero$RMSE_nonlin
  
  cov250hetero[s, ] <- res250hetero$cov
  cov_lin250hetero[s, ] <- res250hetero$cov_lin
  cov_nonlin250hetero[s, ] <- res250hetero$cov_nonlin
  RMSE250hetero[s, ] <- res250hetero$RMSE
  RMSE_lin250hetero[s, ] <- res250hetero$RMSE_lin
  RMSE_nonlin250hetero[s, ] <- res250hetero$RMSE_nonlin
}
## Optional: return to sequential mode afterward
plan(sequential)

## Coverage ----
### complete function ----
df_cov = bind_rows(
  cov1000homo %>% melt(value.name = "cov") %>% mutate(Run = "Homo1000"),
  cov1000hetero %>% melt(value.name = "cov") %>% mutate(Run = "Hetero1000"),
  cov500homo %>% melt(value.name = "cov") %>% mutate(Run = "Homo500"),
  cov500hetero %>% melt(value.name = "cov") %>% mutate(Run = "Hetero500"),
  cov250homo %>% melt(value.name = "cov") %>% mutate(Run = "Homo250"),
  cov250hetero %>% melt(value.name = "cov") %>% mutate(Run = "Hetero250")) %>%
  mutate(Run = factor(Run, levels = c("Homo250", "Homo500", "Homo1000", "Hetero250", "Hetero500", "Hetero1000"))) %>% 
  left_join(cbind.data.frame(setting = 1:24, sim_settings), by = "setting") %>% 
  left_join(cbind.data.frame(point = 1:npoints, "x" = new_x), by = "point")


plt_cov_SNRH = 
  df_cov %>%
  filter(SNR == "High") %>%
  ggplot() +
  geom_hline(aes(col = "target", yintercept = 0.95), lty = 2) +
  geom_line(aes(x = x, y = cov, colour = Run)) +
  facet_wrap(~setting, scales = "fixed", drop = T,
             labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  scale_color_manual("", values = color_runs) +
  scale_y_continuous(limits = c(0, 1)) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5))) +
  labs(title = "Coverage (full function) (nominal 95%) - [SNR 5]") +
  theme(legend.position = "bottom", legend.byrow = T)

plt_cov_SNRL =
  df_cov %>%
  filter(SNR == "Low") %>%
  ggplot() +
  geom_hline(aes(col = "target", yintercept = 0.95), lty = 2) +
  geom_line(aes(x = x, y = cov, colour = Run)) +
  facet_wrap(~setting, scales = "fixed", drop = T,
             labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  scale_color_manual("", values = color_runs) +
  scale_y_continuous(limits = c(0, 1)) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5))) +
  labs(title = "Coverage (full function) (nominal 95%) - [SNR 2]") +
  theme(legend.position = "bottom", legend.byrow = T)

ggsave(filename = "coverage.pdf",
       plot = marrangeGrob(list(
         plt_cov_SNRH, plt_cov_SNRL
       ), 
       nrow = 1, ncol = 1, top = NULL),
       path = PATH_IMG, width = 12, height = 9)

### linear term ----
df_cov_lin = bind_rows(
  cov_lin1000homo %>% melt(value.name = "cov_lin") %>% mutate(Run = "Homo1000"),
  cov_lin1000hetero %>% melt(value.name = "cov_lin") %>% mutate(Run = "Hetero1000"),
  cov_lin500homo %>% melt(value.name = "cov_lin") %>% mutate(Run = "Homo500"),
  cov_lin500hetero %>% melt(value.name = "cov_lin") %>% mutate(Run = "Hetero500"),
  cov_lin250homo %>% melt(value.name = "cov_lin") %>% mutate(Run = "Homo250"),
  cov_lin250hetero %>% melt(value.name = "cov_lin") %>% mutate(Run = "Hetero250")) %>%
  mutate(Run = factor(Run, levels = c("Homo250", "Homo500", "Homo1000", "Hetero250", "Hetero500", "Hetero1000"))) %>% 
  left_join(cbind.data.frame(setting = 1:24, sim_settings), by = "setting") %>% 
  left_join(cbind.data.frame(point = 1:npoints, "x" = new_x), by = "point")


plt_cov_lin_SNRH = 
  df_cov_lin %>%
  filter(SNR == "High") %>%
  ggplot() +
  geom_hline(aes(col = "target", yintercept = 0.95), lty = 2) +
  geom_line(aes(x = x, y = cov_lin, colour = Run)) +
  facet_wrap(~setting, scales = "fixed", drop = T,
             labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  scale_color_manual("", values = color_runs) +
  scale_y_continuous(limits = c(0, 1)) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5))) +
  labs(title = "Coverage (linear term) (nominal 95%) - [SNR 5]") +
  theme(legend.position = "bottom", legend.byrow = T)

plt_cov_lin_SNRL =
  df_cov_lin %>%
  filter(SNR == "Low") %>%
  ggplot() +
  geom_hline(aes(col = "target", yintercept = 0.95), lty = 2) +
  geom_line(aes(x = x, y = cov_lin, colour = Run)) +
  facet_wrap(~setting, scales = "fixed", drop = T,
             labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  scale_color_manual("", values = color_runs) +
  scale_y_continuous(limits = c(0, 1)) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5))) +
  labs(title = "Coverage (linear term) (nominal 95%) - [SNR 2]") +
  theme(legend.position = "bottom", legend.byrow = T)

ggsave(filename = "coverage_lin.pdf",
       plot = marrangeGrob(list(
         plt_cov_lin_SNRH, plt_cov_lin_SNRL
       ), 
       nrow = 1, ncol = 1, top = NULL),
       path = PATH_IMG, width = 12, height = 9)


### nonlinear term ----
df_cov_nonlin = bind_rows(
  cov_nonlin1000homo %>% melt(value.name = "cov_nonlin") %>% mutate(Run = "Homo1000"),
  cov_nonlin1000hetero %>% melt(value.name = "cov_nonlin") %>% mutate(Run = "Hetero1000"),
  cov_nonlin500homo %>% melt(value.name = "cov_nonlin") %>% mutate(Run = "Homo500"),
  cov_nonlin500hetero %>% melt(value.name = "cov_nonlin") %>% mutate(Run = "Hetero500"),
  cov_nonlin250homo %>% melt(value.name = "cov_nonlin") %>% mutate(Run = "Homo250"),
  cov_nonlin250hetero %>% melt(value.name = "cov_nonlin") %>% mutate(Run = "Hetero250")) %>%
  mutate(Run = factor(Run, levels = c("Homo250", "Homo500", "Homo1000", "Hetero250", "Hetero500", "Hetero1000"))) %>% 
  left_join(cbind.data.frame(setting = 1:24, sim_settings), by = "setting") %>% 
  left_join(cbind.data.frame(point = 1:npoints, "x" = new_x), by = "point")


plt_cov_nonlin_SNRH = 
  df_cov_nonlin %>%
  filter(SNR == "High") %>%
  ggplot() +
  geom_hline(aes(col = "target", yintercept = 0.95), lty = 2) +
  geom_line(aes(x = x, y = cov_nonlin, colour = Run)) +
  facet_wrap(~setting, scales = "fixed", drop = T,
             labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  scale_color_manual("", values = color_runs) +
  scale_y_continuous(limits = c(0, 1)) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5))) +
  labs(title = "Coverage (non-linear term) (nominal 95%) - [SNR 5]") +
  theme(legend.position = "bottom", legend.byrow = T)

plt_cov_nonlin_SNRL =
  df_cov_nonlin %>%
  filter(SNR == "Low") %>%
  ggplot() +
  geom_hline(aes(col = "target", yintercept = 0.95), lty = 2) +
  geom_line(aes(x = x, y = cov_nonlin, colour = Run)) +
  facet_wrap(~setting, scales = "fixed", drop = T,
             labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  scale_color_manual("", values = color_runs) +
  scale_y_continuous(limits = c(0, 1)) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5))) +
  labs(title = "Coverage (non-linear term) (nominal 95%) - [SNR 2]") +
  theme(legend.position = "bottom", legend.byrow = T)

ggsave(filename = "coverage_nonlin.pdf",
       plot = marrangeGrob(list(
         plt_cov_nonlin_SNRH, plt_cov_nonlin_SNRL
       ), 
       nrow = 1, ncol = 1, top = NULL),
       path = PATH_IMG, width = 12, height = 9)

### intercept and linear coefficient (separately) ----
cov_lincoeffs1000homo = cov_lincoeffs1000hetero = 
  cov_lincoeffs500homo = cov_lincoeffs500hetero = 
  cov_lincoeffs250homo = cov_lincoeffs250hetero =  
  array(NA, dim = c(nrow(sim_settings), 2),
        dimnames = list("setting" = 1:24, "beta" = c("intercept", "linear")))

for (s in 1:nrow(sim_settings)){
  cat(s, "\t")
  
  s_s = sim_settings[s, ]
  
  I1000homo = input1000homo$res[[s]]
  I1000hetero = input1000hetero$res[[s]]
  I500homo = input500homo$res[[s]]
  I500hetero = input500hetero$res[[s]]
  I250homo = input250homo$res[[s]]
  I250hetero = input250hetero$res[[s]]
  
  # Linear transf
  if (s_s["transf_x"] == "parabola"){
    if(s_s["cov_distr"] == "Unif_rad2"){
      beta0 = 2/3
      beta1 = 0 
    } else if(s_s["cov_distr"] == "Beta_2_4"){
      beta0 = 10/21
      beta1 = -2*sqrt(14)/21
    }
  } else if (s_s["transf_x"] == "cubic"){
    if(s_s["cov_distr"] == "Unif_rad2"){
      beta0 = 0
      beta1 = 2*sqrt(6)/5
    } else if(s_s["cov_distr"] == "Beta_2_4"){
      beta0 = -2*sqrt(2)/7
      beta1 = 4*sqrt(7)/21
    }
  } else if (s_s["transf_x"] == "trigonometric"){
    if(s_s["cov_distr"] == "Unif_rad2"){
      beta0 = 0.37599
      beta1 = 0.39919
    } else if(s_s["cov_distr"] == "Beta_2_4"){
      beta0 = 0.03975
      beta1 = 0.770398
    }
  }
  
  
  tmp_lincoeffs1000hetero = sapply(I1000hetero, function(I) apply(I$draws[, 1:2], 2, quantile, probs = c(0.025, 0.975)))
  cov_lincoeffs1000hetero[s, "intercept"] = mean(tmp_lincoeffs1000hetero[1, ] < beta0 & tmp_lincoeffs1000hetero[2, ] >= beta0)
  cov_lincoeffs1000hetero[s, "linear"] = mean(tmp_lincoeffs1000hetero[3, ] < beta1 & tmp_lincoeffs1000hetero[4, ] >= beta1)
  
  tmp_lincoeffs1000homo = sapply(I1000homo, function(I) apply(I$draws[, 1:2], 2, quantile, probs = c(0.025, 0.975)))
  cov_lincoeffs1000homo[s, "intercept"] = mean(tmp_lincoeffs1000homo[1, ] < beta0 & tmp_lincoeffs1000homo[2, ] >= beta0)
  cov_lincoeffs1000homo[s, "linear"] = mean(tmp_lincoeffs1000homo[3, ] < beta1 & tmp_lincoeffs1000homo[4, ] >= beta1)
  
  tmp_lincoeffs500hetero = sapply(I500hetero, function(I) apply(I$draws[, 1:2], 2, quantile, probs = c(0.025, 0.975)))
  cov_lincoeffs500hetero[s, "intercept"] = mean(tmp_lincoeffs500hetero[1, ] < beta0 & tmp_lincoeffs500hetero[2, ] >= beta0)
  cov_lincoeffs500hetero[s, "linear"] = mean(tmp_lincoeffs500hetero[3, ] < beta1 & tmp_lincoeffs500hetero[4, ] >= beta1)
  
  tmp_lincoeffs500homo = sapply(I500homo, function(I) apply(I$draws[, 1:2], 2, quantile, probs = c(0.025, 0.975)))
  cov_lincoeffs500homo[s, "intercept"] = mean(tmp_lincoeffs500homo[1, ] < beta0 & tmp_lincoeffs500homo[2, ] >= beta0)
  cov_lincoeffs500homo[s, "linear"] = mean(tmp_lincoeffs500homo[3, ] < beta1 & tmp_lincoeffs500homo[4, ] >= beta1)
  
  tmp_lincoeffs250hetero = sapply(I250hetero, function(I) apply(I$draws[, 1:2], 2, quantile, probs = c(0.025, 0.975)))
  cov_lincoeffs250hetero[s, "intercept"] = mean(tmp_lincoeffs250hetero[1, ] < beta0 & tmp_lincoeffs250hetero[2, ] >= beta0)
  cov_lincoeffs250hetero[s, "linear"] = mean(tmp_lincoeffs250hetero[3, ] < beta1 & tmp_lincoeffs250hetero[4, ] >= beta1)
  
  tmp_lincoeffs250homo = sapply(I250homo, function(I) apply(I$draws[, 1:2], 2, quantile, probs = c(0.025, 0.975)))
  cov_lincoeffs250homo[s, "intercept"] = mean(tmp_lincoeffs250homo[1, ] < beta0 & tmp_lincoeffs250homo[2, ] >= beta0)
  cov_lincoeffs250homo[s, "linear"] = mean(tmp_lincoeffs250homo[3, ] < beta1 & tmp_lincoeffs250homo[4, ] >= beta1)
}

df_coverage_lincoeffs = bind_rows(
  cov_lincoeffs1000homo %>% melt(value.name = "cov") %>% mutate(Run = "Homo1000"),
  cov_lincoeffs1000hetero %>% melt(value.name = "cov") %>% mutate(Run = "Hetero1000"),
  cov_lincoeffs500homo %>% melt(value.name = "cov") %>% mutate(Run = "Homo500"),
  cov_lincoeffs500hetero %>% melt(value.name = "cov") %>% mutate(Run = "Hetero500"),
  cov_lincoeffs250homo %>% melt(value.name = "cov") %>% mutate(Run = "Homo250"),
  cov_lincoeffs250hetero %>% melt(value.name = "cov") %>% mutate(Run = "Hetero250")) %>%
  mutate(Run = factor(Run, levels = c("Homo250", "Homo500", "Homo1000", "Hetero250", "Hetero500", "Hetero1000"))) %>%
  left_join(cbind.data.frame(setting = 1:24, sim_settings), by = "setting")


plt_covlincoeff_SNRH = df_coverage_lincoeffs %>%
  filter(SNR == "High") %>% 
  ggplot() +
  geom_hline(aes(col = "target", yintercept = 0.95), lty = 2) +
  geom_point(aes(x = beta, y = cov, colour = Run)) +
  scale_color_manual("", values = color_runs) +
  scale_y_continuous(limits = c(0, 1)) +
  facet_wrap(~setting, scales = "fixed", drop = T, 
             labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5))) +
  labs(title = "Coverage of linear coefficient (nominal 95%)") +
  theme(legend.position = "bottom", legend.byrow = T)

plt_covlincoeff_SNRL = df_coverage_lincoeffs %>%
  filter(SNR == "Low") %>% 
  ggplot() +
  geom_hline(aes(col = "target", yintercept = 0.95), lty = 2) +
  geom_point(aes(x = beta, y = cov, colour = Run)) +
  scale_color_manual("", values = color_runs) +
  scale_y_continuous(limits = c(0, 1)) +
  facet_wrap(~setting, scales = "fixed", drop = T, 
             labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5))) +
  labs(title = "Coverage of linear coefficient (nominal 95%)") +
  theme(legend.position = "bottom", legend.byrow = T)

ggsave(filename = "coverage_lincoeff.pdf",
       plot = marrangeGrob(list(
         plt_covlincoeff_SNRH, plt_covlincoeff_SNRL
       ), 
       nrow = 1, ncol = 1, top = NULL),
       path = PATH_IMG, width = 12, height = 9)



## RMSE ----
### complete function ----
df_RMSE = bind_rows(
  RMSE1000homo %>% melt(value.name = "RMSE") %>% mutate(Run = "Homo1000"),
  RMSE1000hetero %>% melt(value.name = "RMSE") %>% mutate(Run = "Hetero1000"),
  RMSE500homo %>% melt(value.name = "RMSE") %>% mutate(Run = "Homo500"),
  RMSE500hetero %>% melt(value.name = "RMSE") %>% mutate(Run = "Hetero500"),
  RMSE250homo %>% melt(value.name = "RMSE") %>% mutate(Run = "Homo250"),
  RMSE250hetero %>% melt(value.name = "RMSE") %>% mutate(Run = "Hetero250")) %>%
  mutate(Run = factor(Run, levels = c("Homo250", "Homo500", "Homo1000", "Hetero250", "Hetero500", "Hetero1000"))) %>% 
  left_join(cbind.data.frame(setting = 1:24, sim_settings), by = "setting") %>% 
  left_join(cbind.data.frame(point = 1:npoints, "x" = new_x), by = "point")

plt_RMSE_SNRH = 
  df_RMSE %>%
  filter(SNR == "High") %>%
  ggplot() +
  geom_line(aes(x = x, y = RMSE, colour = Run)) +
  facet_wrap(~setting, scales = "free_y", drop = T,
             labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  scale_color_manual("", values = color_runs) +
  # scale_y_continuous(limits = c(0, 3.5)) + 
  labs(title = "RMSE (full function) [SNR = 5]") +
  theme(legend.position = "bottom", legend.byrow = T)

plt_RMSE_SNRL = 
  df_RMSE %>%
  filter(SNR == "Low") %>%
  ggplot() +
  geom_line(aes(x = x, y = RMSE, colour = Run)) +
  facet_wrap(~setting, scales = "free_y", drop = T,
             labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  scale_color_manual("", values = color_runs) +
  # scale_y_continuous(limits = c(0, 3.5)) + 
  labs(title = "RMSE (full function) [SNR = 2]") +
  theme(legend.position = "bottom", legend.byrow = T)

ggsave(filename = "RMSE_full.pdf",
       plot = marrangeGrob(list(
         plt_RMSE_SNRH, plt_RMSE_SNRL
       ), 
       nrow = 1, ncol = 1, top = NULL),
       path = PATH_IMG, width = 12, height = 9)

### linear term ----
df_RMSE_lin = bind_rows(
  RMSE_lin1000homo %>% melt(value.name = "RMSE") %>% mutate(Run = "Homo1000"),
  RMSE_lin1000hetero %>% melt(value.name = "RMSE") %>% mutate(Run = "Hetero1000"),
  RMSE_lin500homo %>% melt(value.name = "RMSE") %>% mutate(Run = "Homo500"),
  RMSE_lin500hetero %>% melt(value.name = "RMSE") %>% mutate(Run = "Hetero500"),
  RMSE_lin250homo %>% melt(value.name = "RMSE") %>% mutate(Run = "Homo250"),
  RMSE_lin250hetero %>% melt(value.name = "RMSE") %>% mutate(Run = "Hetero250")) %>%
  mutate(Run = factor(Run, levels = c("Homo250", "Homo500", "Homo1000", "Hetero250", "Hetero500", "Hetero1000"))) %>% 
  left_join(cbind.data.frame(setting = 1:24, sim_settings), by = "setting") %>% 
  left_join(cbind.data.frame(point = 1:npoints, "x" = new_x), by = "point")

plt_RMSE_lin_SNRH = 
  df_RMSE_lin %>%
  filter(SNR == "High") %>%
  ggplot() +
  geom_line(aes(x = x, y = RMSE, colour = Run)) +
  facet_wrap(~setting, scales = "free_y", drop = T,
             labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  scale_color_manual("", values = color_runs) +
  # scale_y_continuous(limits = c(0, 0.6)) +
  labs(title = "RMSE (linear term) [SNR = 5]") +
  theme(legend.position = "bottom", legend.byrow = T)

plt_RMSE_lin_SNRL = 
  df_RMSE_lin %>%
  filter(SNR == "Low") %>%
  ggplot() +
  geom_line(aes(x = x, y = RMSE, colour = Run)) +
  facet_wrap(~setting, scales = "free_y", drop = T,
             labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  scale_color_manual("", values = color_runs) +
  # scale_y_continuous(limits = c(0, 0.6)) +
  labs(title = "RMSE (linear term) [SNR = 2]") +
  theme(legend.position = "bottom", legend.byrow = T)


ggsave(filename = "RMSE_lin.pdf",
       plot = marrangeGrob(list(
         plt_RMSE_lin_SNRH, plt_RMSE_lin_SNRL
       ),
       nrow = 1, ncol = 1, top = NULL),
       path = PATH_IMG, width = 12, height = 9)

### non-linear term ----
df_RMSE_nonlin = bind_rows(
  RMSE_nonlin1000homo %>% melt(value.name = "RMSE") %>% mutate(Run = "Homo1000"),
  RMSE_nonlin1000hetero %>% melt(value.name = "RMSE") %>% mutate(Run = "Hetero1000"),
  RMSE_nonlin500homo %>% melt(value.name = "RMSE") %>% mutate(Run = "Homo500"),
  RMSE_nonlin500hetero %>% melt(value.name = "RMSE") %>% mutate(Run = "Hetero500"),
  RMSE_nonlin250homo %>% melt(value.name = "RMSE") %>% mutate(Run = "Homo250"),
  RMSE_nonlin250hetero %>% melt(value.name = "RMSE") %>% mutate(Run = "Hetero250")) %>%
  mutate(Run = factor(Run, levels = c("Homo250", "Homo500", "Homo1000", "Hetero250", "Hetero500", "Hetero1000"))) %>%
  left_join(cbind.data.frame(setting = 1:24, sim_settings), by = "setting") %>%
  left_join(cbind.data.frame(point = 1:npoints, "x" = new_x), by = "point")

plt_RMSE_nonlin_SNRH =
  df_RMSE_nonlin %>%
  filter(SNR == "High") %>%
  ggplot() +
  geom_line(aes(x = x, y = RMSE, colour = Run)) +
  facet_wrap(~setting, scales = "free_y", drop = T,
             labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  scale_color_manual("", values = color_runs) +
  # scale_y_continuous(limits = c(0, 3.5)) +
  labs(title = "RMSE (non-linear term) [SNR = 5]") +
  theme(legend.position = "bottom", legend.byrow = T)

plt_RMSE_nonlin_SNRL =
  df_RMSE_nonlin %>%
  filter(SNR == "Low") %>%
  ggplot() +
  geom_line(aes(x = x, y = RMSE, colour = Run)) +
  facet_wrap(~setting, scales = "free_y", drop = T,
             labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  scale_color_manual("", values = color_runs) +
  # scale_y_continuous(limits = c(0, 3.5)) +
  labs(title = "RMSE (non-linear term) [SNR = 2]") +
  theme(legend.position = "bottom", legend.byrow = T)

ggsave(filename = "RMSE_nonlin.pdf",
       plot = marrangeGrob(list(
         plt_RMSE_nonlin_SNRH, plt_RMSE_nonlin_SNRL
       ),
       nrow = 1, ncol = 1, top = NULL),
       path = PATH_IMG, width = 12, height = 9)

### bias interc. and lin. coeff. ----
bias_lincoeffs1000homo = bias_lincoeffs1000hetero =
  bias_lincoeffs500homo = bias_lincoeffs500hetero = 
  bias_lincoeffs250homo = bias_lincoeffs250hetero =  
  array(NA, dim = c(nrow(sim_settings), 2),
        dimnames = list("setting" = 1:24, "beta" = c("intercept", "linear")))  


for (s in 1:nrow(sim_settings)){
  
  s_s = sim_settings[s, ]
  
  I1000homo = input1000homo$res[[s]]
  I1000hetero = input1000hetero$res[[s]]
  I500homo = input500homo$res[[s]]
  I500hetero = input500hetero$res[[s]]
  I250homo = input250homo$res[[s]]
  I250hetero = input250hetero$res[[s]]
  
  # Linear transf
  if (s_s["transf_x"] == "parabola"){
    if(s_s["cov_distr"] == "Unif_rad2"){
      beta0 = 2/3
      beta1 = 0 
    } else if(s_s["cov_distr"] == "Beta_2_4"){
      beta0 = 10/21
      beta1 = -2*sqrt(14)/21
    }
  } else if (s_s["transf_x"] == "cubic"){
    if(s_s["cov_distr"] == "Unif_rad2"){
      beta0 = 0
      beta1 = 2*sqrt(6)/5
    } else if(s_s["cov_distr"] == "Beta_2_4"){
      beta0 = -2*sqrt(2)/7
      beta1 = 4*sqrt(7)/21
    }
  } else if (s_s["transf_x"] == "trigonometric"){
    if(s_s["cov_distr"] == "Unif_rad2"){
      beta0 = 0.37599
      beta1 = 0.39919
    } else if(s_s["cov_distr"] == "Beta_2_4"){
      beta0 = 0.03975
      beta1 = 0.770398
    }
  }
  
  bias_lincoeffs1000hetero[s, ] = rowMeans(sapply(I1000hetero, function(I) colMeans(I$draws[, 1:2])) - c(beta0, beta1))
  bias_lincoeffs1000homo[s, ] = rowMeans(sapply(I1000homo, function(I) colMeans(I$draws[, 1:2])) - c(beta0, beta1))
  
  bias_lincoeffs500hetero[s, ] = rowMeans(sapply(I500hetero, function(I) colMeans(I$draws[, 1:2])) - c(beta0, beta1))
  bias_lincoeffs500homo[s, ] = rowMeans(sapply(I500homo, function(I) colMeans(I$draws[, 1:2])) - c(beta0, beta1))
  
  bias_lincoeffs250hetero[s, ] = rowMeans(sapply(I250hetero, function(I) colMeans(I$draws[, 1:2])) - c(beta0, beta1))
  bias_lincoeffs250homo[s, ] = rowMeans(sapply(I250homo, function(I) colMeans(I$draws[, 1:2])) - c(beta0, beta1))
}


df_bias_lincoeffs = bind_rows(
  bias_lincoeffs1000homo %>% melt(value.name = "bias") %>% mutate(Run = "Homo1000"),
  bias_lincoeffs1000hetero %>% melt(value.name = "bias") %>% mutate(Run = "Hetero1000"),
  bias_lincoeffs500homo %>% melt(value.name = "bias") %>% mutate(Run = "Homo500"),
  bias_lincoeffs500hetero %>% melt(value.name = "bias") %>% mutate(Run = "Hetero500"),
  bias_lincoeffs250homo %>% melt(value.name = "bias") %>% mutate(Run = "Homo250"),
  bias_lincoeffs250hetero %>% melt(value.name = "bias") %>% mutate(Run = "Hetero250")) %>%
  mutate(Run = factor(Run, levels = c("Homo250", "Homo500", "Homo1000", "Hetero250", "Hetero500", "Hetero1000"))) %>%
  left_join(cbind.data.frame(setting = 1:24, sim_settings), by = "setting")


plt_biaslincoeff_SNRH = df_bias_lincoeffs %>%
  filter(SNR == "High") %>% 
  ggplot() +
  geom_hline(aes(col = "target", yintercept = 0), lty = 2) +
  geom_point(aes(x = beta, y = bias, colour = Run)) +
  scale_color_manual("", values = color_runs) +
  facet_wrap(~setting, scales = "fixed", drop = T, 
             labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5))) +
  labs(title = "Bias of linear coefficient [SNR = 5]") +
  theme(legend.position = "bottom", legend.byrow = T)

plt_biaslincoeff_SNRL = df_bias_lincoeffs %>%
  filter(SNR == "Low") %>%
  ggplot() +
  geom_hline(aes(col = "target", yintercept = 0), lty = 2) +
  geom_point(aes(x = beta, y = bias, colour = Run)) +
  scale_color_manual("", values = color_runs) +
  facet_wrap(~setting, scales = "fixed", drop = T,
             labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5))) +
  labs(title = "Bias of linear coefficient [SNR = 2]") +
  theme(legend.position = "bottom", legend.byrow = T)


ggsave(filename = "bias_lincoeff.pdf",
       plot = marrangeGrob(list(
         plt_biaslincoeff_SNRH, plt_biaslincoeff_SNRL
       ),
       nrow = 1, ncol = 1, top = NULL),
       path = PATH_IMG, width = 12, height = 9)


# ESS ----
df_ESS1000homo = lapply(input1000homo$res, function(I1) sapply(I1, function(I2) I2$diagn["ess_thin10", ])) %>% simplify2array() %>% 
  melt() %>% rename(Beta = Var1, Replicate = Var2, ESS = value, Setting = Var3)
df_ESS1000hetero = lapply(input1000hetero$res, function(I1) sapply(I1, function(I2) I2$diagn["ess_thin10", ])) %>% simplify2array() %>% 
  melt() %>% rename(Beta = Var1, Replicate = Var2, ESS = value, Setting = Var3)
df_ESS500homo = lapply(input500homo$res, function(I1) sapply(I1, function(I2) I2$diagn["ess_thin10", ])) %>% simplify2array() %>% 
  melt() %>% rename(Beta = Var1, Replicate = Var2, ESS = value, Setting = Var3)
df_ESS500hetero = lapply(input500hetero$res, function(I1) sapply(I1, function(I2) I2$diagn["ess_thin10", ])) %>% simplify2array() %>% 
  melt() %>% rename(Beta = Var1, Replicate = Var2, ESS = value, Setting = Var3)
df_ESS250homo = lapply(input250homo$res, function(I1) sapply(I1, function(I2) I2$diagn["ess_thin10", ])) %>% simplify2array() %>% 
  melt() %>% rename(Beta = Var1, Replicate = Var2, ESS = value, Setting = Var3)
df_ESS250hetero = lapply(input250hetero$res, function(I1) sapply(I1, function(I2) I2$diagn["ess_thin10", ])) %>% simplify2array() %>% 
  melt() %>% rename(Beta = Var1, Replicate = Var2, ESS = value, Setting = Var3)

df_ESS = bind_rows(
  df_ESS1000homo %>% mutate(Run = "Homo1000"),
  df_ESS1000hetero %>% mutate(Run = "Hetero1000"),
  df_ESS500homo %>% mutate(Run = "Homo500"),
  df_ESS500hetero %>% mutate(Run = "Hetero500"),
  df_ESS250homo %>% mutate(Run = "Homo250"),
  df_ESS250hetero %>% mutate(Run = "Hetero250")) %>% 
  mutate(ESS = if_else(is.na(ESS), 0, ESS),
         Run = factor(Run, levels = c("Homo250", "Homo500", "Homo1000", "Hetero250", "Hetero500", "Hetero1000")))


plt_ESS1000 = bind_rows(
  df_ESS1000homo %>% mutate(Run = "Homo1000"),
  df_ESS1000hetero %>% mutate(Run = "Hetero1000")) %>% 
  ggplot() + 
  geom_boxplot(aes(y = ESS, x = Beta, col = Run)) +
  geom_hline(aes(yintercept = 1000), col = 2, lty = 2) +
  facet_wrap(~Setting, labeller = labeller(Setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24)),
             nrow = 4, dir = "v") +
  scale_color_manual("", values = color_runs) +
  scale_y_continuous(limits = c(0, 1500)) +
  labs(title = "ESS - n=1000") +
  theme(legend.position = "bottom",
        axis.text.x = element_blank(),
        text = element_text(size = 14))

plt_ESS500 = bind_rows(
  df_ESS500homo %>% mutate(Run = "Homo500"),
  df_ESS500hetero %>% mutate(Run = "Hetero500")) %>% 
  ggplot() + 
  geom_boxplot(aes(y = ESS, x = Beta, col = Run)) +
  geom_hline(aes(yintercept = 1000), col = 2, lty = 2) +
  facet_wrap(~Setting, labeller = labeller(Setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24)),
             nrow = 4, dir = "v") +
  scale_color_manual("", values = color_runs) +
  scale_y_continuous(limits = c(0, 1500)) +
  labs(title = "ESS - n=500") +
  theme(legend.position = "bottom",
        axis.text.x = element_blank(),
        text = element_text(size = 14))

plt_ESS250 = bind_rows(
  df_ESS250homo %>% mutate(Run = "Homo250"),
  df_ESS250hetero %>% mutate(Run = "Hetero250")) %>% 
  ggplot() + 
  geom_boxplot(aes(y = ESS, x = Beta, col = Run)) +
  geom_hline(aes(yintercept = 1000), col = 2, lty = 2) +
  facet_wrap(~Setting, labeller = labeller(Setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24)),
             nrow = 4, dir = "v") +
  scale_color_manual("", values = color_runs) +
  scale_y_continuous(limits = c(0, 1500)) +
  labs(title = "ESS - n=250") +
  theme(legend.position = "bottom",
        axis.text.x = element_blank(),
        text = element_text(size = 14))


ggsave(filename = "ESS_n.pdf",
       plot = marrangeGrob(list(plt_ESS1000, plt_ESS500, plt_ESS250), nrow = 1, ncol = 1, top = NULL),
       path = PATH_IMG, width = 18, height = 10)


# Rhat ----
df_RHAT1000homo = lapply(input1000homo$res, function(I1) sapply(I1, function(I2) I2$diagn["rhat_thin10", ])) %>% simplify2array() %>% 
  melt() %>% rename(Beta = Var1, Replicate = Var2, RHAT = value, setting = Var3)
df_RHAT1000hetero = lapply(input1000hetero$res, function(I1) sapply(I1, function(I2) I2$diagn["rhat_thin10", ])) %>% simplify2array() %>% 
  melt() %>% rename(Beta = Var1, Replicate = Var2, RHAT = value, setting = Var3)
df_RHAT500homo = lapply(input500homo$res, function(I1) sapply(I1, function(I2) I2$diagn["rhat_thin10", ])) %>% simplify2array() %>% 
  melt() %>% rename(Beta = Var1, Replicate = Var2, RHAT = value, setting = Var3)
df_RHAT500hetero = lapply(input500hetero$res, function(I1) sapply(I1, function(I2) I2$diagn["rhat_thin10", ])) %>% simplify2array() %>% 
  melt() %>% rename(Beta = Var1, Replicate = Var2, RHAT = value, setting = Var3)
df_RHAT250homo = lapply(input250homo$res, function(I1) sapply(I1, function(I2) I2$diagn["rhat_thin10", ])) %>% simplify2array() %>% 
  melt() %>% rename(Beta = Var1, Replicate = Var2, RHAT = value, setting = Var3)
df_RHAT250hetero = lapply(input250hetero$res, function(I1) sapply(I1, function(I2) I2$diagn["rhat_thin10", ])) %>% simplify2array() %>% 
  melt() %>% rename(Beta = Var1, Replicate = Var2, RHAT = value, setting = Var3)

df_RHAT = bind_rows(
  df_RHAT1000homo %>% mutate(Run = "Homo1000"),
  df_RHAT1000hetero %>% mutate(Run = "Hetero1000"),
  df_RHAT500homo %>% mutate(Run = "Homo500"),
  df_RHAT500hetero %>% mutate(Run = "Hetero500"),
  df_RHAT250homo %>% mutate(Run = "Homo250"),
  df_RHAT250hetero %>% mutate(Run = "Hetero250")) %>%
  mutate(Run = factor(Run, levels = c("Homo1000", "Hetero1000", "Homo500",  "Hetero500", "Homo250", "Hetero250"))) %>% 
  left_join(cbind.data.frame(setting = 1:24, sim_settings), by = "setting")

df_RHAT %>% 
  group_by(Run, setting) %>%
  summarise(below1.01 = mean(RHAT < 1.01, na.rm = T)) %>% 
  pivot_wider(names_from = Run, values_from = below1.01)

df_RHAT %>% 
  group_by(Beta, setting) %>%
  summarise(below1.01 = mean(RHAT < 1.01, na.rm = T)) %>% 
  pivot_wider(names_from = Beta, values_from = below1.01) %>% print(n = Inf, width = Inf)


plt_RHAT_SNRH = df_RHAT %>% 
  filter(SNR == "High") %>% 
  group_by(setting, Beta, Run) %>%
  summarise(below1.01 = mean(RHAT < 1.01, na.rm = T)) %>% 
  ggplot() +
  geom_point(aes(x = Beta, y = below1.01, col = Run)) +
  facet_wrap(~setting, labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24)),
             ncol = 4) +
  scale_color_manual("", values = color_runs) +
  labs(title = "Proportion of replicates with Rhat < 1.01 [SNR = 5]") +
  theme(legend.position = "bottom",
        axis.text.x = element_blank(),
        text = element_text(size = 14))

plt_RHAT_SNRL = df_RHAT %>% 
  filter(SNR == "Low") %>% 
  group_by(setting, Beta, Run) %>%
  summarise(below1.01 = mean(RHAT < 1.01, na.rm = T)) %>% 
  ggplot() +
  geom_point(aes(x = Beta, y = below1.01, col = Run)) +
  facet_wrap(~setting, labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24)),
             ncol = 4) +
  scale_color_manual("", values = color_runs) +
  labs(title = "Proportion of replicates with Rhat < 1.01 [SNR = 2]") +
  theme(legend.position = "bottom",
        axis.text.x = element_blank(),
        text = element_text(size = 14))

ggsave(filename = "Rhat.pdf",
       plot = marrangeGrob(list(plt_RHAT_SNRH, plt_RHAT_SNRL), nrow = 1, ncol = 1, top = NULL),
       path = PATH_IMG, width = 12, height = 9)



# W (learning rate) ----
df_w1000homo = sapply(input1000homo$res, function(I1) sapply(I1, function(I2) I2$w)) %>%
  melt() %>% rename(Replicate = Var1, w = value, Setting = Var2)
df_w1000hetero = sapply(input1000hetero$res, function(I1) sapply(I1, function(I2) I2$w)) %>%
  melt() %>% rename(Replicate = Var1, w = value, Setting = Var2)
df_w500homo = sapply(input500homo$res, function(I1) sapply(I1, function(I2) I2$w)) %>%
  melt() %>% rename(Replicate = Var1, w = value, Setting = Var2)
df_w500hetero = sapply(input500hetero$res, function(I1) sapply(I1, function(I2) I2$w)) %>%
  melt() %>% rename(Replicate = Var1, w = value, Setting = Var2)
df_w250homo = sapply(input250homo$res, function(I1) sapply(I1, function(I2) I2$w)) %>%
  melt() %>% rename(Replicate = Var1, w = value, Setting = Var2)
df_w250hetero = sapply(input250hetero$res, function(I1) sapply(I1, function(I2) I2$w)) %>%
  melt() %>% rename(Replicate = Var1, w = value, Setting = Var2)

df_w = bind_rows(
  df_w1000homo %>% mutate(Run = "Homo1000"),
  df_w1000hetero %>% mutate(Run = "Hetero1000"),
  df_w500homo %>% mutate(Run = "Homo500"),
  df_w500hetero %>% mutate(Run = "Hetero500"),
  df_w250homo %>% mutate(Run = "Homo250"),
  df_w250hetero %>% mutate(Run = "Hetero250")) %>% 
  mutate(Run = factor(Run, levels = c("Homo250", "Homo500", "Homo1000", "Hetero250", "Hetero500", "Hetero1000")))


plt_w = df_w %>% 
  ggplot() + 
  geom_boxplot(aes(y = w, x = Run, col = Run)) +
  facet_wrap(~Setting, scales = "fixed",
             labeller = labeller(Setting = setNames(apply(sim_settings_names, 1, paste, collapse = " - "), 1:24)),
             nrow = 4, dir = "v", ) +
  scale_color_manual("", values = color_runs) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.25))) +
  labs(title = "w (learning rate)") +
  theme(legend.position = "bottom", legend.byrow = T,
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = 14))

ggsave(filename = "w.pdf",
       plot = plt_w,
       path = PATH_IMG, width = 18, height = 10)


# k (flatness of loss) ----
df_k1000homo = lapply(input1000homo$res, function(I1) sapply(I1, function(I2) I2$k)) %>% simplify2array() %>%
  melt() %>% rename(Side = Var1, Replicate = Var2, k = value, Setting = Var3)
df_k1000hetero = lapply(input1000hetero$res, function(I1) sapply(I1, function(I2) I2$k)) %>% simplify2array() %>%
  melt() %>% rename(Side = Var1, Replicate = Var2, k = value, Setting = Var3)
df_k500homo = lapply(input500homo$res, function(I1) sapply(I1, function(I2) I2$k)) %>% simplify2array() %>%
  melt() %>% rename(Side = Var1, Replicate = Var2, k = value, Setting = Var3)
df_k500hetero = lapply(input500hetero$res, function(I1) sapply(I1, function(I2) I2$k)) %>% simplify2array() %>%
  melt() %>% rename(Side = Var1, Replicate = Var2, k = value, Setting = Var3)
df_k250homo = lapply(input250homo$res, function(I1) sapply(I1, function(I2) I2$k)) %>% simplify2array() %>%
  melt() %>% rename(Side = Var1, Replicate = Var2, k = value, Setting = Var3)
df_k250hetero = lapply(input250hetero$res, function(I1) sapply(I1, function(I2) I2$k)) %>% simplify2array() %>%
  melt() %>% rename(Side = Var1, Replicate = Var2, k = value, Setting = Var3)

df_k = bind_rows(
  df_k1000homo %>% mutate(Run = "Homo1000"),
  df_k1000hetero %>% mutate(Run = "Hetero1000"),
  df_k500homo %>% mutate(Run = "Homo500"),
  df_k500hetero %>% mutate(Run = "Hetero500"),
  df_k250homo %>% mutate(Run = "Homo250"),
  df_k250hetero %>% mutate(Run = "Hetero250")) %>% 
  mutate(Run = factor(Run, levels = c("Homo250", "Homo500", "Homo1000", "Hetero250", "Hetero500", "Hetero1000")),
         Side = factor(Side, levels = c(1, 2), labels = c("Left", "Right")))


plt_k = df_k %>% 
  ggplot() + 
  geom_boxplot(aes(y = k, x = interaction(Side, Run), col = Run)) +
  facet_wrap(~Setting, scales = "fixed",
             labeller = labeller(Setting = setNames(apply(sim_settings_names, 1, paste, collapse = " - "), 1:24)),
             nrow = 4, dir = "v", ) +
  scale_color_manual("", values = color_runs) +
  scale_x_discrete("Side", labels = function(x) str_split_i(x, "\\.", 1)) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.25))) +
  labs(title = "k (flatness of loss)") +
  theme(legend.position = "bottom", legend.byrow = T,
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 14))

ggsave(filename = "k.pdf",
       plot = plt_k,
       path = PATH_IMG, width = 18, height = 10)


# Frequentist RMSE ----
RMSEfreq1000homo = RMSEfreq1000hetero = RMSEfreq500homo = RMSEfreq500hetero = RMSEfreq250homo = RMSEfreq250hetero = 
  array(NA, dim = c(nrow(sim_settings), npoints),
        dimnames = list("setting" = 1:24, "point" = 1:npoints))
RMSEfreq_lin1000homo = RMSEfreq_lin1000hetero = RMSEfreq_lin500homo = RMSEfreq_lin500hetero = RMSEfreq_lin250homo = RMSEfreq_lin250hetero =
  array(NA, dim = c(nrow(sim_settings), npoints),
        dimnames = list("setting" = 1:24, "point" = 1:npoints))
RMSEfreq_nonlin1000homo = RMSEfreq_nonlin1000hetero = RMSEfreq_nonlin500homo = RMSEfreq_nonlin500hetero = RMSEfreq_nonlin250homo = RMSEfreq_nonlin250hetero =
  array(NA, dim = c(nrow(sim_settings), npoints), 
        dimnames = list("setting" = 1:24, "point" = 1:npoints))


cov_RMSE_fun_freq = function(I, new_x, true_f, true_lin, true_nonlin){
  opt = which.min(apply(I$summ[10:14, ], 1, function(b) loss_asymm_pop(b, I$y, I$X_scaled, I$k, 1e-1, g = 1)))
  b_opt = I$summ[9+opt, ]
  DR = I$DR
  newB = mgcv::Predict.matrix(DR$smooth_object, data.frame(x = new_x))
  new_Z <- newB %*% DR$Trafo
  x_scaled <- (new_x - mean(x)) / sd(x)
  new_X <- cbind(1, x_scaled, new_Z)
  fit = drop(new_X %*% b_opt)
  fit_lin = drop(new_X[, 1:2] %*% b_opt[1:2])
  fit_nonlin = drop(fit - fit_lin)
  SE = (fit - true_f)^2
  SE_lin = (fit_lin - true_lin)^2
  SE_nonlin = (fit_nonlin - (true_f - true_nonlin))^2
  
  out = cbind(SE = SE, SE_lin = SE_lin, SE_nonlin = SE_nonlin)
  
  return(out)
}


library(future)
library(future.apply)
## Choose number of workers
plan(multisession, workers = 8)

for (s in 1:nrow(sim_settings)){
  cat(s, "\t")
  
  s_s = sim_settings[s, ]
  
  I1000homo = input1000homo$res[[s]]
  I1000hetero = input1000hetero$res[[s]]
  I500homo = input500homo$res[[s]]
  I500hetero = input500hetero$res[[s]]
  I250homo = input250homo$res[[s]]
  I250hetero = input250hetero$res[[s]]
  
  transf = switch(sim_settings[s, "transf_x"],
                  "parabola" = function(x) x^2,
                  "cubic" = function(x) x^3,
                  "trigonometric" = function(x) sin(2.5*x) + 3*exp(-(5^2)*((x)^2))
  )
  
  # Linear transf
  if (s_s["transf_x"] == "parabola"){
    if(s_s["cov_distr"] == "Unif_rad2"){
      transf_lin = function(x) rep(2/3, length(x))
    } else if(s_s["cov_distr"] == "Beta_2_4"){
      transf_lin = function(x) 1/7 - sqrt(2)/2*x
    }
  } else if (s_s["transf_x"] == "cubic"){
    if(s_s["cov_distr"] == "Unif_rad2"){
      transf_lin = function(x) 6/5*x
    } else if(s_s["cov_distr"] == "Beta_2_4"){
      transf_lin = function(x) sqrt(2)/21 + x
    }
  } else if (s_s["transf_x"] == "trigonometric"){
    if(s_s["cov_distr"] == "Unif_rad2"){
      transf_lin = function(x) 0.37599 + 0.488903*x
    } else if(s_s["cov_distr"] == "Beta_2_4"){
      transf_lin = function(x)  0.760395 + 1.52871*x
    }
  }
  
  true_f = transf(new_x)
  true_lin = transf_lin(new_x)
  true_nonlin = true_f - true_lin
  
  tmp1000homo = sqrt(rowMeans(future_lapply(I1000homo, function(I) cov_RMSE_fun_freq(I, new_x, true_f, true_lin, true_nonlin)) %>% simplify2array(), dims = 2))
  tmp1000hetero = sqrt(rowMeans(future_lapply(I1000hetero, function(I) cov_RMSE_fun_freq(I, new_x, true_f, true_lin, true_nonlin)) %>% simplify2array(), dims = 2))
  tmp500homo = sqrt(rowMeans(future_lapply(I500homo, function(I) cov_RMSE_fun_freq(I, new_x, true_f, true_lin, true_nonlin)) %>% simplify2array(), dims = 2))
  tmp500hetero = sqrt(rowMeans(future_lapply(I500hetero, function(I) cov_RMSE_fun_freq(I, new_x, true_f, true_lin, true_nonlin)) %>% simplify2array(), dims = 2))
  tmp250homo = sqrt(rowMeans(future_lapply(I250homo, function(I) cov_RMSE_fun_freq(I, new_x, true_f, true_lin, true_nonlin)) %>% simplify2array(), dims = 2))
  tmp250hetero = sqrt(rowMeans(future_lapply(I250hetero, function(I) cov_RMSE_fun_freq(I, new_x, true_f, true_lin, true_nonlin)) %>% simplify2array(), dims = 2))
  
  RMSEfreq1000homo[s, ] = tmp1000homo[, "SE"]
  RMSEfreq1000hetero[s, ] = tmp1000hetero[, "SE"]
  RMSEfreq500homo[s, ] = tmp500homo[, "SE"]
  RMSEfreq500hetero[s, ] = tmp500hetero[, "SE"]
  RMSEfreq250homo[s, ] = tmp250homo[, "SE"]
  RMSEfreq250hetero[s, ] = tmp250hetero[, "SE"]
  
  RMSEfreq_lin1000homo[s, ] = tmp1000homo[, "SE_lin"]
  RMSEfreq_lin1000hetero[s, ] = tmp1000hetero[, "SE_lin"]
  RMSEfreq_lin500homo[s, ] = tmp500homo[, "SE_lin"]
  RMSEfreq_lin500hetero[s, ] = tmp500hetero[, "SE_lin"]
  RMSEfreq_lin250homo[s, ] = tmp250homo[, "SE_lin"]
  RMSEfreq_lin250hetero[s, ] = tmp250hetero[, "SE_lin"]
  
  RMSEfreq_nonlin1000homo[s, ] = tmp1000homo[, "SE_nonlin"]
  RMSEfreq_nonlin1000hetero[s, ] = tmp1000hetero[, "SE_nonlin"]
  RMSEfreq_nonlin500homo[s, ] = tmp500homo[, "SE_nonlin"]
  RMSEfreq_nonlin500hetero[s, ] = tmp500hetero[, "SE_nonlin"]
  RMSEfreq_nonlin250homo[s, ] = tmp250homo[, "SE_nonlin"]
  RMSEfreq_nonlin250hetero[s, ] = tmp250hetero[, "SE_nonlin"]
}
plan(sequential)

df_RMSEfreq = bind_rows(
  RMSEfreq1000homo %>% melt(value.name = "RMSE_freq") %>% mutate(Run = "Homo1000"),
  RMSEfreq1000hetero %>% melt(value.name = "RMSE_freq") %>% mutate(Run = "Hetero1000"),
  RMSEfreq500homo %>% melt(value.name = "RMSE_freq") %>% mutate(Run = "Homo500"),
  RMSEfreq500hetero %>% melt(value.name = "RMSE_freq") %>% mutate(Run = "Hetero500"),
  RMSEfreq250homo %>% melt(value.name = "RMSE_freq") %>% mutate(Run = "Homo250"),
  RMSEfreq250hetero %>% melt(value.name = "RMSE_freq") %>% mutate(Run = "Hetero250")) %>%
  mutate(Run = factor(Run, levels = c("Homo250", "Homo500", "Homo1000", "Hetero250", "Hetero500", "Hetero1000"))) %>%
  left_join(cbind.data.frame(setting = 1:24, sim_settings), by = "setting")


plt_RMSEfreq_SNRH = 
  df_RMSEfreq %>%
  filter(SNR == "High") %>%
  ggplot() +
  geom_line(aes(x = x, y = RMSE_freq, colour = Run)) +
  facet_wrap(~setting, scales = "free_y", drop = T,
             labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  scale_color_manual("", values = color_runs) +
  # scale_y_continuous(limits = c(0, 3.5)) + 
  labs(title = "RMSE freq. (full function) [SNR = 5]") +
  theme(legend.position = "bottom", legend.byrow = T)

plt_RMSEfreq_SNRL = 
  df_RMSEfreq %>%
  filter(SNR == "Low") %>%
  ggplot() +
  geom_line(aes(x = x, y = RMSE_freq, colour = Run)) +
  facet_wrap(~setting, scales = "free_y", drop = T,
             labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  scale_color_manual("", values = color_runs) +
  # scale_y_continuous(limits = c(0, 3.5)) + 
  labs(title = "RMSE freq. (full function) [SNR = 2]") +
  theme(legend.position = "bottom", legend.byrow = T)

ggsave(filename = "RMSE_full.pdf",
       plot = marrangeGrob(list(
         plt_RMSE_SNRH, plt_RMSE_SNRL
       ), 
       nrow = 1, ncol = 1, top = NULL),
       path = PATH_IMG, width = 12, height = 9)

full_join(df_RMSEfreq, df_RMSE, by = c("setting", "point", "Run", "SNR", "cov_distr", "transf_x")) %>% 
  group_by(setting) %>% summarize(Pearson = cor(RMSE, RMSE_freq) %>% round(2),
                                  Spearman = cor(RMSE, RMSE_freq, method = "spearman") %>% round(2)) %>% 
  bind_cols(sim_settings) %>% 
  write_csv(file.path(PATH_IMG, "cor_RMSE_freq.csv"), col_names = T)

full_join(df_RMSE, df_cov, by = c("setting", "point", "Run", "SNR", "cov_distr", "transf_x")) %>% 
  group_by(setting) %>% summarize(Pearson = cor(RMSE, cov) %>% round(2),
                                  Spearman = cor(RMSE, cov, method = "spearman") %>% round(2)) %>% 
  bind_cols(sim_settings) %>% 
  write_csv(file.path(PATH_IMG, "cor_RMSE_cov.csv"), col_names = T)
