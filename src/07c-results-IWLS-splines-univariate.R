rm(list = ls())

suppressPackageStartupMessages(library(tidyverse))
library(reshape2)
library(gridExtra)

input1 = readRDS("sim_study_nonC-GBI/splines_univ/sim_study_nonCalibrated_asymm_TL_init_PropFlatter.RDS")
input2 = readRDS("sim_study_nonC-GBI/splines_univ/sim_study_nonCalibrated_asymm_TL_init_ksqrt_c1e-1_PropFC.RDS")
input3 = readRDS("sim_study_nonC-GBI/splines_univ/sim_study_nonCalibrated_asymm_TL_init_ksqrt_c1e-2_PropFlatter.RDS")
inputD1 = readRDS("sim_study_nonC-GBI/splines_univ/sim_study_nonCalibrated_asymm_TL_init_diagH_Propc1e-1.RDS")$res
inputD2 = readRDS("sim_study_nonC-GBI/splines_univ/sim_study_nonCalibrated_asymm_TL_init_diagH_sqrtK_c1e-1_PropFC.RDS")$res

summ1 = lapply(input1,
               function(I)
                 abind::abind(lapply(I, function(x) x$summ), along = 3))
summ2 = lapply(input2,
               function(I)
                 abind::abind(lapply(I, function(x) x$summ), along = 3))

summ3 = lapply(input3,
               function(I)
                 abind::abind(lapply(I, function(x) x$summ), along = 3))

summD1 = lapply(inputD1,
                function(I)
                  abind::abind(lapply(I, function(x) x$summ), along = 3))

summD2 = lapply(inputD2,
                function(I)
                  abind::abind(lapply(I, function(x) x$summ), along = 3))

PATH_IMG = "sim_study_nonC-GBI/splines_univ/Together/"

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
col1 = "#8ecae6"
col2 = "#219ebc"
col3 = "#023047"
colD1 = "#ffd60a"
colD2 = "#fb8500"

# CHECK FAILED JOBS ----
input1_wfail = input1
input2_wfail = input2
input3_wfail = input3
inputD1_wfail = inputD1
inputD2_wfail = inputD2
sapply(input1_wfail, function(I) sum(sapply(I, function(x) is.null(x$draws))))
sapply(input2_wfail, function(I) sum(sapply(I, function(x) is.null(x$draws))))
sapply(input3_wfail, function(I) sum(sapply(I, function(x) is.null(x$draws))))
sapply(inputD1_wfail, function(I) sum(sapply(I, function(x) is.null(x$draws))))
sapply(inputD2_wfail, function(I) sum(sapply(I, function(x) is.null(x$draws))))


# Remove failed jobs
for (s in 1:nrow(sim_settings)){
  failed1 = which(sapply(input1[[s]], function(I) is.null(I$draws)))
  failed2 = which(sapply(input2[[s]], function(I) is.null(I$draws)))
  failed3 = which(sapply(input3[[s]], function(I) is.null(I$draws)))
  
  input1[[s]][failed1] = NULL
  input2[[s]][failed2] = NULL
  input3[[s]][failed3] = NULL
}


# FIT vs TRUE FUNCTION ----
## Fit ----
fit1 = fit2 = fit3 = 
  fitD1 = fitD2 = 
  array(NA, c(nrow(sim_settings), n, 100),
        dimnames = list("setting" = 1:24, "unit" = 1:n, "replicate" = 1:100))

for (s in 1:nrow(sim_settings)){
  
  cat(s, "\t")
  
  I1 = input1[[s]]
  I2 = input2[[s]]
  I3 = input3[[s]]
  ID1 = inputD1[[s]]
  ID2 = inputD2[[s]]
  
  fit1_tmp = sapply(I1, function(I) rowMeans(I$X_scaled %*% t(I$draws)))
  fit2_tmp = sapply(I2, function(I) rowMeans(I$X_scaled %*% t(I$draws)))
  fit3_tmp = sapply(I3, function(I) rowMeans(I$X_scaled %*% t(I$draws)))
  fitD1_tmp = sapply(ID1, function(I) rowMeans(I$X_scaled %*% t(I$draws)))
  fitD2_tmp = sapply(ID2, function(I) rowMeans(I$X_scaled %*% t(I$draws)))
  
  fit1[s, , 1:ncol(fit1_tmp)] = fit1_tmp
  fit2[s, , 1:ncol(fit2_tmp)] = fit2_tmp
  fit3[s, , 1:ncol(fit3_tmp)] = fit3_tmp
  fitD1[s, , 1:ncol(fitD1_tmp)] = fitD1_tmp
  fitD2[s, , 1:ncol(fitD2_tmp)] = fitD2_tmp
}

## X ----
X1 = X2 = X3 = 
  XD1 = XD2 = 
  array(NA, c(nrow(sim_settings), n, 100),
        dimnames = list("setting" = 1:24, "unit" = 1:n, "replicate" = 1:100))

for (s in 1:nrow(sim_settings)){
  
  cat(s, "\t")
  
  I1 = input1[[s]]
  I2 = input2[[s]]
  I3 = input3[[s]]
  ID1 = inputD1[[s]]
  ID2 = inputD2[[s]]
  
  X1_tmp = sapply(I1, function(I) I$x1, simplify = "array")
  X2_tmp = sapply(I2, function(I) I$x1, simplify = "array")
  X3_tmp = sapply(I3, function(I) I$x1, simplify = "array")
  XD1_tmp = sapply(ID1, function(I) I$x1, simplify = "array")
  XD2_tmp = sapply(ID2, function(I) I$x1, simplify = "array")
  
  X1[s, , 1:ncol(X1_tmp)] = X1_tmp
  X2[s, , 1:ncol(X2_tmp)] = X2_tmp
  X3[s, , 1:ncol(X3_tmp)] = X3_tmp
  XD1[s, , 1:ncol(XD1_tmp)] = XD1_tmp
  XD2[s, , 1:ncol(XD2_tmp)] = XD2_tmp
}


## True functions ----
true_f = data.frame(setting = integer(0L), x = numeric(0L), y = numeric(0L))
for (s in 1:nrow(sim_settings)){
  xmin = min(c(X1[s, , ], X2[s, , ], X3[s, , ], XD1[s, , ], XD2[s, , ]), na.rm = T)
  xmax = max(c(X1[s, , ], X2[s, , ], X3[s, , ], XD1[s, , ], XD2[s, , ]), na.rm = T)
  
  transf = switch(sim_settings[s, "transf_x"],
                  "parabola" = function(x) x^2,
                  "cubic" = function(x) x^3,
                  "trigonometric" = function(x) sin(2.5*x) + 3*exp(-(5^2)*((x)^2)))
  
  x = seq(xmin, xmax, length = 501)
  y = transf(x)
  
  true_f = bind_rows(true_f, data.frame(setting = s, x = x, y = y))
}
true_f = left_join(true_f, cbind.data.frame(setting = 1:24, sim_settings), by = "setting")


df_fit1 = left_join(fit1 %>% melt(value.name = "fit"),
                    X1 %>% melt(value.name = "x"),
                    by = c("setting", "unit", "replicate")) %>% 
  left_join(cbind.data.frame(setting = 1:24, sim_settings), by = "setting")

df_fit2 = left_join(fit2 %>% melt(value.name = "fit"),
                    X2 %>% melt(value.name = "x"),
                    by = c("setting", "unit", "replicate")) %>% 
  left_join(cbind.data.frame(setting = 1:24, sim_settings), by = "setting")

df_fit3 = left_join(fit3 %>% melt(value.name = "fit"),
                    X3 %>% melt(value.name = "x"),
                    by = c("setting", "unit", "replicate")) %>%
  left_join(cbind.data.frame(setting = 1:24, sim_settings), by = "setting")

df_fitD1 = left_join(fitD1 %>% melt(value.name = "fit"),
                     XD1 %>% melt(value.name = "x"),
                     by = c("setting", "unit", "replicate")) %>%
  left_join(cbind.data.frame(setting = 1:24, sim_settings), by = "setting")

df_fitD2 = left_join(fitD2 %>% melt(value.name = "fit"),
                     XD2 %>% melt(value.name = "x"),
                     by = c("setting", "unit", "replicate")) %>%
  left_join(cbind.data.frame(setting = 1:24, sim_settings), by = "setting")

y_limits <- list(
  scale_y_continuous(limits = c(-0.5, 3)),  scale_y_continuous(limits = c(-0.5, 3)),
  scale_y_continuous(limits = c(-0.5, 3)),  scale_y_continuous(limits = c(-0.5, 3)),
  scale_y_continuous(limits = c(-3, 3)),  scale_y_continuous(limits = c(-3, 3)),
  scale_y_continuous(limits = c(-3, 3)),  scale_y_continuous(limits = c(-3, 3)),
  scale_y_continuous(limits = c(-1.5, 3.5)),  scale_y_continuous(limits = c(-1.5, 3.5)),
  scale_y_continuous(limits = c(-1.5, 3.5)),  scale_y_continuous(limits = c(-1.5, 3.5))
)



plt_fit1_SNR2 = df_fit1%>% 
  filter(SNR == 2) %>% 
  ggplot() +
  geom_line(aes(x = x, y = fit, group = replicate), alpha = 0.15) +
  geom_line(data = true_f %>% filter(SNR == 2),
            mapping = aes(x = x, y = y), col = "red") +
  facet_wrap(~setting, scales = "free", labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  ggh4x::facetted_pos_scales(y = y_limits) +
  labs(title = "k = init, k_deriv = sqrt(init), c = 1e-3, c_deriv = 1e-1 [SNR = 2]")
plt_fit1_SNR1 = df_fit1 %>%
  filter(SNR == 1) %>% 
  ggplot() +
  geom_line(aes(x = x, y = fit, group = replicate), alpha = 0.15) +
  geom_line(data = true_f %>% filter(SNR == 1),
            mapping = aes(x = x, y = y), col = "red") +
  facet_wrap(~setting, scales = "free", labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  ggh4x::facetted_pos_scales(y = y_limits) +
  labs(title = "k = init, k_deriv = sqrt(init), c = 1e-3, c_deriv = 1e-1 [SNR = 1]")


plt_fit2_SNR2 = df_fit2 %>% 
  filter(SNR == 2) %>% 
  ggplot() +
  geom_line(aes(x = x, y = fit, group = replicate), alpha = 0.15) +
  geom_line(data = true_f %>% filter(SNR == 2),
            mapping = aes(x = x, y = y), col = "red") +
  facet_wrap(~setting, scales = "free", labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  ggh4x::facetted_pos_scales(y = y_limits) +
  labs(title = "k = sqrt(init), k_deriv = sqrt(init), c = 1e-1, c_deriv = 1e-1 [SNR = 2]")
plt_fit2_SNR1 = df_fit2 %>%
  filter(SNR == 1) %>% 
  ggplot() +
  geom_line(aes(x = x, y = fit, group = replicate), alpha = 0.15) +
  geom_line(data = true_f %>% filter(SNR == 1),
            mapping = aes(x = x, y = y), col = "red") +
  facet_wrap(~setting, scales = "free", labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  ggh4x::facetted_pos_scales(y = y_limits) +
  labs(title = "k = sqrt(init), k_deriv = sqrt(init), c = 1e-1, c_deriv = 1e-1 [SNR = 2]")

plt_fit3_SNR2 = df_fit3 %>% 
  filter(SNR == 2) %>% 
  ggplot() +
  geom_line(aes(x = x, y = fit, group = replicate), alpha = 0.15) +
  geom_line(data = true_f %>% filter(SNR == 2),
            mapping = aes(x = x, y = y), col = "red") +
  facet_wrap(~setting, scales = "free", labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  ggh4x::facetted_pos_scales(y = y_limits) +
  labs(title = "k = sqrt(init), k_deriv = (init)^(1/4), c = 1e-2, c_deriv = 1e-2 [SNR = 2]")
plt_fit3_SNR1 = df_fit3 %>%
  filter(SNR == 1) %>% 
  ggplot() +
  geom_line(aes(x = x, y = fit, group = replicate), alpha = 0.15) +
  geom_line(data = true_f %>% filter(SNR == 1),
            mapping = aes(x = x, y = y), col = "red") +
  facet_wrap(~setting, scales = "free", labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  ggh4x::facetted_pos_scales(y = y_limits) +
  labs(title = "k = sqrt(init), k_deriv = (init)^(1/4), c = 1e-2, c_deriv = 1e-2 [SNR = 1]")


plt_fitD1_SNR2 = df_fitD1 %>% 
  filter(SNR == 2) %>% 
  ggplot() +
  geom_line(aes(x = x, y = fit, group = replicate), alpha = 0.15) +
  geom_line(data = true_f %>% filter(SNR == 2),
            mapping = aes(x = x, y = y), col = "red") +
  facet_wrap(~setting, scales = "free", labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  ggh4x::facetted_pos_scales(y = y_limits) +
  labs(title = "diag H | k = init, k_deriv = init, c = 1e-3, c_deriv = 1e-1, [SNR = 2]")
plt_fitD1_SNR1 = df_fitD1 %>%
  filter(SNR == 1) %>% 
  ggplot() +
  geom_line(aes(x = x, y = fit, group = replicate), alpha = 0.15) +
  geom_line(data = true_f %>% filter(SNR == 1),
            mapping = aes(x = x, y = y), col = "red") +
  facet_wrap(~setting, scales = "free", labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  ggh4x::facetted_pos_scales(y = y_limits) +
  labs(title = "diag H | k = init, k_deriv = init, c = 1e-3, c_deriv = 1e-1, [SNR = 1]")

plt_fitD2_SNR2 = df_fitD2 %>% 
  filter(SNR == 2) %>% 
  ggplot() +
  geom_line(aes(x = x, y = fit, group = replicate), alpha = 0.15) +
  geom_line(data = true_f %>% filter(SNR == 2),
            mapping = aes(x = x, y = y), col = "red") +
  facet_wrap(~setting, scales = "free", labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  ggh4x::facetted_pos_scales(y = y_limits) +
  labs(title = "diag H | k = sqrt(init), k_deriv = sqrt(init), c = 1e-1, c_deriv = 1e-1 [SNR = 2]")
plt_fitD2_SNR1 = df_fitD2 %>%
  filter(SNR == 1) %>% 
  ggplot() +
  geom_line(aes(x = x, y = fit, group = replicate), alpha = 0.15) +
  geom_line(data = true_f %>% filter(SNR == 1),
            mapping = aes(x = x, y = y), col = "red") +
  facet_wrap(~setting, scales = "free", labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  ggh4x::facetted_pos_scales(y = y_limits) +
  labs(title = "diag H | k = sqrt(init), k_deriv = sqrt(init), c = 1e-1, c_deriv = 1e-1 [SNR = 1]")

ggsave(filename = "fit.pdf",
       plot = marrangeGrob(list(
         plt_fit1_SNR2, plt_fit1_SNR1, 
         plt_fit2_SNR2, plt_fit2_SNR1, 
         plt_fit3_SNR2, plt_fit3_SNR1,
         plt_fitD1_SNR2, plt_fitD1_SNR1,
         plt_fitD2_SNR2, plt_fitD2_SNR1
       ), 
       nrow = 1, ncol = 1, top = NULL),
       path = PATH_IMG, width = 12, height = 9)

# DENSITY RESIDUALS ----
y1 = y2 = y3 = 
  yD1 = yD2 = 
  array(NA, c(nrow(sim_settings), n, 100),
        dimnames = list("setting" = 1:24, "unit" = 1:n, "replicate" = 1:100))

for (s in 1:nrow(sim_settings)){
  
  cat(s, "\t")
  
  I1 = input1[[s]]
  I2 = input2[[s]]
  I3 = input3[[s]]
  ID1 = inputD1[[s]]
  ID2 = inputD2[[s]]
  
  y1_tmp = sapply(I1, function(I) I$y, simplify = "array")
  y2_tmp = sapply(I2, function(I) I$y, simplify = "array")
  y3_tmp = sapply(I3, function(I) I$y, simplify = "array")
  yD1_tmp = sapply(ID1, function(I) I$y, simplify = "array")
  yD2_tmp = sapply(ID2, function(I) I$y, simplify = "array")
  
  y1[s, , 1:ncol(y1_tmp)] = y1_tmp
  y2[s, , 1:ncol(y2_tmp)] = y2_tmp
  y3[s, , 1:ncol(y3_tmp)] = y3_tmp
  yD1[s, , 1:ncol(yD1_tmp)] = yD1_tmp
  yD2[s, , 1:ncol(yD2_tmp)] = yD2_tmp
}


df_y1 = left_join(y1 %>% melt(value.name = "y"),
                      fit1 %>% melt(value.name = "fit"),
                      by = c("setting", "unit", "replicate")) %>% 
  left_join(cbind.data.frame(setting = 1:24, sim_settings), by = "setting")
df_y2 = left_join(y2 %>% melt(value.name = "y"),
                      fit2 %>% melt(value.name = "fit"),
                      by = c("setting", "unit", "replicate")) %>%
  left_join(cbind.data.frame(setting = 1:24, sim_settings), by = "setting")
df_y3 = left_join(y3 %>% melt(value.name = "y"),
                      fit3 %>% melt(value.name = "fit"),
                      by = c("setting", "unit", "replicate")) %>%
  left_join(cbind.data.frame(setting = 1:24, sim_settings), by = "setting")
df_yD1 = left_join(yD1 %>% melt(value.name = "y"),
                       fitD1 %>% melt(value.name = "fit"),
                       by = c("setting", "unit", "replicate")) %>%
  left_join(cbind.data.frame(setting = 1:24, sim_settings), by = "setting")
df_yD2 = left_join(yD2 %>% melt(value.name = "y"),
                       fitD2 %>% melt(value.name = "fit"),
                       by = c("setting", "unit", "replicate")) %>%
  left_join(cbind.data.frame(setting = 1:24, sim_settings), by = "setting")



plt_res1_SNR2 = df_y1 %>%
  filter(SNR == 2) %>% 
  ggplot() +
  geom_density(aes(x = y - fit, group = replicate), color = alpha("black", alpha = 0.2)) +
  facet_wrap(~setting, scales = "free", labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  labs(title = "k = init, k_deriv = sqrt(init), c = 1e-3, c_deriv = 1e-1 [SNR = 2]")
plt_res1_SNR1 = df_y1 %>%
  filter(SNR == 1) %>% 
  ggplot() +
  geom_density(aes(x = y - fit, group = replicate), color = alpha("black", alpha = 0.2)) +
  facet_wrap(~setting, scales = "free", labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  labs(title = "k = init, k_deriv = sqrt(init), c = 1e-3, c_deriv = 1e-1 [SNR = 1]")

plt_res2_SNR2 = df_y2 %>%
  filter(SNR == 2) %>% 
  ggplot() +
  geom_density(aes(x = y - fit, group = replicate), color = alpha("black", alpha = 0.2)) +
  facet_wrap(~setting, scales = "free", labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  labs(title = "k = sqrt(init), k_deriv = sqrt(init), c = 1e-1, c_deriv = 1e-1 [SNR = 2]")
plt_res2_SNR1 = df_y2 %>%
  filter(SNR == 1) %>% 
  ggplot() +
  geom_density(aes(x = y - fit, group = replicate), color = alpha("black", alpha = 0.2)) +
  facet_wrap(~setting, scales = "free", labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  labs(title = "k = sqrt(init), k_deriv = sqrt(init), c = 1e-1, c_deriv = 1e-1 [SNR = 1]")

plt_res3_SNR2 = df_y3 %>%
  filter(SNR == 2) %>% 
  ggplot() +
  geom_density(aes(x = y - fit, group = replicate), color = alpha("black", alpha = 0.2)) +
  facet_wrap(~setting, scales = "free", labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  labs(title = "k = sqrt(init), k_deriv = (init)^(1/4), c = 1e-2, c_deriv = 1e-2 [SNR = 2]")
plt_res3_SNR1 = df_y3 %>%
  filter(SNR == 1) %>% 
  ggplot() +
  geom_density(aes(x = y - fit, group = replicate), color = alpha("black", alpha = 0.2)) +
  facet_wrap(~setting, scales = "free", labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  labs(title = "k = sqrt(init), k_deriv = (init)^(1/4), c = 1e-2, c_deriv = 1e-2 [SNR = 1]")

plt_resD1_SNR2 = df_yD1 %>%
  filter(SNR == 2) %>% 
  ggplot() +
  geom_density(aes(x = y - fit, group = replicate), color = alpha("black", alpha = 0.2)) +
  facet_wrap(~setting, scales = "free", labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  labs(title = "diag H | k = init, k_deriv = init, c = 1e-3, c_deriv = 1e-1, [SNR = 2]")
plt_resD1_SNR1 = df_yD1 %>%
  filter(SNR == 1) %>% 
  ggplot() +
  geom_density(aes(x = y - fit, group = replicate), color = alpha("black", alpha = 0.2)) +
  facet_wrap(~setting, scales = "free", labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  labs(title = "diag H | k = init, k_deriv = init, c = 1e-3, c_deriv = 1e-1, [SNR = 1]")

plt_resD2_SNR2 = df_yD2 %>%
  filter(SNR == 2) %>% 
  ggplot() +
  geom_density(aes(x = y - fit, group = replicate), color = alpha("black", alpha = 0.2)) +
  facet_wrap(~setting, scales = "free", labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  labs(title = "diag H | k = sqrt(init), k_deriv = sqrt(init), c = 1e-1, c_deriv = 1e-1 [SNR = 2]")

ggsave(filename = "residuals.pdf",
       plot = marrangeGrob(list(
         plt_res1_SNR2, plt_res1_SNR1, 
         plt_res2_SNR2, plt_res2_SNR1, 
         plt_res3_SNR2, plt_res3_SNR1,
         plt_resD1_SNR2, plt_resD1_SNR1,
         plt_resD2_SNR2
       ), nrow = 1, ncol = 1),
       path = PATH_IMG, width = 12, height = 9)


# DIAGNOSTICS ----
## Traceplots coefficients ----
pdf(paste0(PATH_IMG, "traceplot_examples1.pdf"), width = 12, height = 8)
for (s in 1:nrow(sim_settings)){
  I = input1[[s]][[1]]
  
  # Add outer margin at the top for the page title
  par(mfrow = c(4, 6), mar = c(3, 2, 2, 1), oma = c(0, 0, 3, 0))
  
  for (j in 1:22){
    plot(I$draws[ , j], type = "l", 
         main = colnames(I$draws)[j], xlab = "MCMC Iteration", ylab = "")
    abline(h = 0, col = "red")
  }
  plot(I$draws_tau[ , 1], type = "l", main = "tau_lin")
  plot(I$draws_tau[ , 2], type = "l", main = "tau_nonlin")
  
  # Page title — written after all plots so it spans the full page
  page_title <- paste(paste0(sim_settings[s, ], collapse = " - "))
  mtext(page_title, outer = TRUE, cex = 1.2, font = 2, line = 1)
}
dev.off()

pdf(paste0(PATH_IMG, "traceplot_examples2.pdf"), width = 12, height = 8)
for (s in 1:nrow(sim_settings)){
  I = input2[[s]][[1]]
  
  # Add outer margin at the top for the page title
  par(mfrow = c(4, 6), mar = c(3, 2, 2, 1), oma = c(0, 0, 3, 0))
  
  for (j in 1:22){
    plot(I$draws[ , j], type = "l", 
         main = colnames(I$draws)[j], xlab = "MCMC Iteration", ylab = "")
    abline(h = 0, col = "red")
  }
  plot(I$draws_tau[ , 1], type = "l", main = "tau_lin")
  plot(I$draws_tau[ , 2], type = "l", main = "tau_nonlin")
  
  # Page title — written after all plots so it spans the full page
  page_title <- paste(paste0(sim_settings[s, ], collapse = " - "))
  mtext(page_title, outer = TRUE, cex = 1.2, font = 2, line = 1)
}
dev.off()

pdf(paste0(PATH_IMG, "traceplot_examples3.pdf"), width = 12, height = 8)
for (s in 1:nrow(sim_settings)){
  I = input3[[s]][[1]]
  
  # Add outer margin at the top for the page title
  par(mfrow = c(4, 6), mar = c(3, 2, 2, 1), oma = c(0, 0, 3, 0))
  
  for (j in 1:22){
    plot(I$draws[ , j], type = "l", 
         main = colnames(I$draws)[j], xlab = "MCMC Iteration", ylab = "")
    abline(h = 0, col = "red")
  }
  plot(I$draws_tau[ , 1], type = "l", main = "tau_lin")
  plot(I$draws_tau[ , 2], type = "l", main = "tau_nonlin")
  
  # Page title — written after all plots so it spans the full page
  page_title <- paste(paste0(sim_settings[s, ], collapse = " - "))
  mtext(page_title, outer = TRUE, cex = 1.2, font = 2, line = 1)
}
dev.off()


pdf(paste0(PATH_IMG, "traceplot_examplesD1.pdf"), width = 12, height = 8)
for (s in 1:nrow(sim_settings)){
  I = inputD1[[s]][[1]]
  
  # Add outer margin at the top for the page title
  par(mfrow = c(4, 6), mar = c(3, 2, 2, 1), oma = c(0, 0, 3, 0))
  
  for (j in 1:22){
    plot(I$draws[ , j], type = "l", 
         main = colnames(I$draws)[j], xlab = "MCMC Iteration", ylab = "")
    abline(h = 0, col = "red")
  }
  plot(I$draws_tau[ , 1], type = "l", main = "tau_lin")
  plot(I$draws_tau[ , 2], type = "l", main = "tau_nonlin")
  
  # Page title — written after all plots so it spans the full page
  page_title <- paste(paste0(sim_settings[s, ], collapse = " - "))
  mtext(page_title, outer = TRUE, cex = 1.2, font = 2, line = 1)
}
dev.off()

pdf(paste0(PATH_IMG, "traceplot_examplesD2.pdf"), width = 12, height = 8)
for (s in 1:nrow(sim_settings)){
  I = inputD2[[s]][[1]]
  
  # Add outer margin at the top for the page title
  par(mfrow = c(4, 6), mar = c(3, 2, 2, 1), oma = c(0, 0, 3, 0))
  
  for (j in 1:22){
    plot(I$draws[ , j], type = "l", 
         main = colnames(I$draws)[j], xlab = "MCMC Iteration", ylab = "")
    abline(h = 0, col = "red")
  }
  plot(I$draws_tau[ , 1], type = "l", main = "tau_lin")
  plot(I$draws_tau[ , 2], type = "l", main = "tau_nonlin")
  
  # Page title — written after all plots so it spans the full page
  page_title <- paste(paste0(sim_settings[s, ], collapse = " - "))
  mtext(page_title, outer = TRUE, cex = 1.2, font = 2, line = 1)
}
dev.off()





## ESS ----
df_ESS1 = sapply(input1, function(I1) sapply(I1, function(I2) I2$diagn["ess_thin10", ])) %>% 
  melt() %>% rename(Beta = Var1, Replicate = Var2, ESS = value, Setting = L1)
df_ESS2 = sapply(input2, function(I1) sapply(I1, function(I2) I2$diagn["ess_thin10", ])) %>%
  melt() %>% rename(Beta = Var1, Replicate = Var2, ESS = value, Setting = L1)
df_ESS3 = sapply(input3, function(I1) sapply(I1, function(I2) I2$diagn["ess_thin10", ])) %>%
  melt() %>% rename(Beta = Var1, Replicate = Var2, ESS = value, Setting = L1)
df_ESSD1 = lapply(inputD1, function(I1) sapply(I1, function(I2) I2$diagn["ess_thin10", ])) %>% simplify2array() %>% 
  melt() %>% rename(Beta = Var1, Replicate = Var2, ESS = value, Setting = Var3)
df_ESSD2 = lapply(inputD2, function(I1) sapply(I1, function(I2) I2$diagn["ess_thin10", ])) %>% simplify2array() %>% 
  melt() %>% rename(Beta = Var1, Replicate = Var2, ESS = value, Setting = Var3)
df_ESS = bind_rows(
  df_ESS1 %>% mutate(Approach = "1"),
  df_ESS2 %>% mutate(Approach = "2"),
  df_ESS3 %>% mutate(Approach = "3"),
  df_ESSD1 %>% mutate(Approach = "D1"),
  df_ESSD2 %>% mutate(Approach = "D2")) %>%
  mutate(ESS = if_else(is.na(ESS), 0, ESS),
         Approach = factor(Approach, levels = c("1", "2", "3", "D1", "D2")))

plot_ESS_init_SNR2 = df_ESS %>% filter(Setting <= 12) %>% 
  mutate(Setting = factor(Setting)) %>% 
  ggplot() +
  geom_boxplot(aes(y = ESS, col = Approach)) +
  geom_hline(aes(yintercept = 1000), col = 2, lty = 2) + 
  facet_wrap(~Setting, labeller = labeller(Setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24)),
             nrow = 4, dir = "v") + 
  scale_color_manual("",
                     values = c("1" = col1, "2" = col2, "3" = col3, "D1" = colD1, "D2" = colD2),
                     labels = c("k = init, k_deriv = sqrt(init), c = 1e-3, c_deriv = 1e-1",
                                "k = sqrt(init), k_deriv = sqrt(init), c = 1e-1, c_deriv = 1e-1",
                                "k = sqrt(init), k_deriv = (init)^(1/4), c = 1e-2, c_deriv = 1e-2", 
                                "diag H | k = init, k_deriv = init, c = 1e-3, c_deriv = 1e-1",
                                "diag H | k = sqrt(init), k_deriv = sqrt(init), c = 1e-1, c_deriv = 1e-1")) +
  scale_y_continuous(limits = c(0, 1500)) +
  guides(color = guide_legend(nrow = 3, byrow = FALSE)) +
  labs(title = "ESS - SNR=2") + 
  theme(legend.position = "bottom",
        axis.text.x = element_blank(),
        text = element_text(size = 16))

plot_ESS_init_SNR1 =  df_ESS %>% filter(Setting > 12) %>% 
  mutate(Setting = factor(Setting)) %>% 
  ggplot() +
  geom_boxplot(aes(y = ESS, col = Approach)) +
  geom_hline(aes(yintercept = 1000), col = 2, lty = 2) + 
  facet_wrap(~Setting, labeller = labeller(Setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24)),
             nrow = 4, dir = "v") + 
  scale_color_manual("", 
                     values = c("1" = col1, "2" = col2, "3" = col3, "D1" = colD1, "D2" = colD2),
                     labels = c("k = init, k_deriv = sqrt(init), c = 1e-3, c_deriv = 1e-1",
                                "k = sqrt(init), k_deriv = sqrt(init), c = 1e-1, c_deriv = 1e-1",
                                "k = sqrt(init), k_deriv = (init)^(1/4), c = 1e-2, c_deriv = 1e-2", 
                                "diag H | k = init, k_deriv = init, c = 1e-3, c_deriv = 1e-1",
                                "diag H | k = sqrt(init), k_deriv = sqrt(init), c = 1e-1, c_deriv = 1e-1")) +
  scale_y_continuous(limits = c(0, 1500)) +
  guides(color = guide_legend(nrow = 3, byrow = FALSE)) +
  labs(title = "ESS - SNR=1") + 
  theme(legend.position = "bottom",
        axis.text.x = element_blank(),
        text = element_text(size = 16))

ggsave(filename = "ESS.pdf",
       plot = marrangeGrob(list(plot_ESS_init_SNR2, plot_ESS_init_SNR1), nrow = 1, ncol = 1, top = NULL),
       path = PATH_IMG, width = 12, height = 9)

