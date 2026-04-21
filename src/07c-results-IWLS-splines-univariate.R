rm(list = ls())

suppressPackageStartupMessages(library(tidyverse))
library(reshape2)
library(gridExtra)

input_init = readRDS("sim_study_nonC-GBI/splines_univ/sim_study_nonCalibrated_asymm_TL_init.RDS")
input_flat = readRDS("sim_study_nonC-GBI/splines_univ/sim_study_nonCalibrated_asymm_TL_init_smallClargeK.RDS")

summ_init = lapply(input_init,
                   function(I)
                     abind::abind(lapply(I, function(x) x$summ), along = 3))
summ_flat = lapply(input_flat,
                   function(I)
                     abind::abind(lapply(I, function(x) x$summ), along = 3))


PATH_IMG = "sim_study_nonC-GBI/splines_univ/"

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
col_flat = "#fb8500"
col_init = "#588157"
col_mid = "#8338ec"

# CHECK FAILED JOBS ----
input_flat_wfail = input_flat
input_init_wfail = input_init

sapply(input_init_wfail, function(I) sum(sapply(I, function(x) is.null(x$draws))))
sapply(input_flat_wfail, function(I) sum(sapply(I, function(x) is.null(x$draws))))

# Remove failed jobs
for (s in 1:nrow(sim_settings)){
  failed_flat = which(sapply(input_flat[[s]], function(I) is.null(I$draws)))
  failed_init = which(sapply(input_init[[s]], function(I) is.null(I$draws)))
  
  input_flat[[s]][failed_flat] = NULL
  input_init[[s]][failed_init] = NULL
}


# FIT vs TRUE FUNCTION ----
## Fit ----
fit_flat = 
  fit_init = array(NA, c(nrow(sim_settings), n, 100),
                   dimnames = list("setting" = 1:24, "unit" = 1:n, "replicate" = 1:100))

for (s in 1:nrow(sim_settings)){
  
  cat(s, "\t")
  
  I_f = input_flat[[s]]
  I_i = input_init[[s]]
  
  fit_f = sapply(I_f, function(I) rowMeans(I$X_scaled %*% t(I$draws)))
  fit_i = sapply(I_i, function(I) rowMeans(I$X_scaled %*% t(I$draws)))
  
  fit_flat[s, , 1:ncol(fit_f)] = fit_f
  fit_init[s, , 1:ncol(fit_i)] = fit_i
}

## X ----
X_flat = 
  X_init = array(NA, c(nrow(sim_settings), n, 100),
                 dimnames = list("setting" = 1:24, "unit" = 1:n, "replicate" = 1:100))

for (s in 1:nrow(sim_settings)){
  
  cat(s, "\t")
  
  I_f = input_flat[[s]]
  I_i = input_init[[s]]
  
  X_f = sapply(I_f, function(I) I$x1, simplify = "array")
  
  X_i = sapply(I_i, function(I) I$x1, simplify = "array")
  
  X_flat[s, , 1:ncol(X_f)] = X_f
  X_init[s, , 1:ncol(X_i)] = X_i
}


## True functions ----
true_f = data.frame(setting = integer(0L), x = numeric(0L), y = numeric(0L))
for (s in 1:nrow(sim_settings)){
  xmin = min(c(X_flat[s, , ], X_init[s, , ]), na.rm = T)
  xmax = max(c(X_flat[s, , ], X_init[s, , ]), na.rm = T)
  
  transf = switch(sim_settings[s, "transf_x"],
                  "parabola" = function(x) x^2,
                  "cubic" = function(x) x^3,
                  "trigonometric" = function(x) sin(2.5*x) + 2*exp(-(5^2)*((x)^2)))
  
  x = seq(xmin, xmax, length = 100)
  y = transf(x)
  
  true_f = bind_rows(true_f, data.frame(setting = s, x = x, y = y))
}
true_f = left_join(true_f, cbind.data.frame(setting = 1:24, sim_settings), by = "setting")


df_fit_flat = left_join(fit_flat %>% melt(value.name = "fit"),
                        X_flat %>% melt(value.name = "x"),
                        by = c("setting", "unit", "replicate")) %>% 
  left_join(cbind.data.frame(setting = 1:24, sim_settings), by = "setting")

df_fit_init = left_join(fit_init %>% melt(value.name = "fit"),
                        X_init %>% melt(value.name = "x"),
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

plt_fit_flat2 = df_fit_flat %>% 
  filter(SNR == 2) %>% 
  ggplot() +
  geom_line(aes(x = x, y = fit, group = replicate), alpha = 0.15) +
  geom_line(data = true_f %>% filter(SNR == 2),
            mapping = aes(x = x, y = y), col = "red") +
  facet_wrap(~setting, scales = "free", labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  ggh4x::facetted_pos_scales(y = y_limits) +
  labs(title = "Flat - SNR = 2")

plt_fit_init2 = df_fit_init %>% 
  filter(SNR == 2) %>% 
  ggplot() +
  geom_line(aes(x = x, y = fit, group = replicate), alpha = 0.15) +
  geom_line(data = true_f %>% filter(SNR == 2),
            mapping = aes(x = x, y = y), col = "red") +
  facet_wrap(~setting, scales = "free", labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  ggh4x::facetted_pos_scales(y = y_limits) +
  labs(title = "Initial - SNR = 2")

plt_fit_flat1 = df_fit_flat %>% 
  filter(SNR == 1) %>% 
  ggplot() +
  geom_line(aes(x = x, y = fit, group = replicate), alpha = 0.15) +
  geom_line(data = true_f %>% filter(SNR == 1),
            mapping = aes(x = x, y = y), col = "red") +
  facet_wrap(~setting, scales = "free", labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  ggh4x::facetted_pos_scales(y = y_limits) +
  labs(title = "Flat - SNR = 1")

plt_fit_init1 = df_fit_init %>% 
  filter(SNR == 1) %>% 
  ggplot() +
  geom_line(aes(x = x, y = fit, group = replicate), alpha = 0.15) +
  geom_line(data = true_f %>% filter(SNR == 1),
            mapping = aes(x = x, y = y), col = "red") +
  facet_wrap(~setting, scales = "free", labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  ggh4x::facetted_pos_scales(y = y_limits) +
  labs(title = "Initial - SNR = 1")

ggsave(filename = "fit.pdf",
       plot = marrangeGrob(list(plt_fit_flat2, plt_fit_init2, plt_fit_flat1, plt_fit_init1), nrow = 1, ncol = 1, top = NULL),
       path = PATH_IMG, width = 12, height = 9)

# DENSITY RESIDUALS ----
y_init = 
  y_flat = array(NA, c(nrow(sim_settings), n, 100),
                 dimnames = list("setting" = 1:24, "unit" = 1:n, "replicate" = 1:100))

for (s in 1:nrow(sim_settings)){
  
  cat(s, "\t")
  
  I_f = input_flat[[s]]
  I_i = input_init[[s]]
  
  y_f = sapply(I_f, function(I) I$y, simplify = "array")
  
  y_i = sapply(I_i, function(I) I$y, simplify = "array")
  
  y_flat[s, , 1:ncol(y_f)] = y_f
  y_init[s, , 1:ncol(y_i)] = y_i
}

df_y_init = left_join(y_init %>% melt(value.name = "y"),
                      fit_init %>% melt(value.name = "fit"),
                      by = c("setting", "unit", "replicate")) %>% 
  left_join(cbind.data.frame(setting = 1:24, sim_settings), by = "setting")

df_y_flat = left_join(y_flat %>% melt(value.name = "y"),
                      fit_flat %>% melt(value.name = "fit"),
                      by = c("setting", "unit", "replicate")) %>% 
  left_join(cbind.data.frame(setting = 1:24, sim_settings), by = "setting")


plt_res_flat2 = df_y_flat %>%
  filter(SNR == 2) %>% 
  ggplot() +
  geom_density(aes(x = y - fit, group = replicate), color = alpha("black", alpha = 0.2)) +
  facet_wrap(~setting, scales = "free", labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  labs(title = "Flat - SNR = 2")

plt_res_init2 = df_y_init %>%
  filter(SNR == 2) %>% 
  ggplot() +
  geom_density(aes(x = y - fit, group = replicate), color = alpha("black", alpha = 0.2)) +
  facet_wrap(~setting, scales = "free", labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  labs(title = "Initial - SNR = 2")

plt_res_flat1 = df_y_flat %>%
  filter(SNR == 1) %>% 
  ggplot() +
  geom_density(aes(x = y - fit, group = replicate), color = alpha("black", alpha = 0.2)) +
  facet_wrap(~setting, scales = "free", labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  labs(title = "Flat - SNR = 1")

plt_res_init1 = df_y_init %>%
  filter(SNR == 1) %>%
  ggplot() +
  geom_density(aes(x = y - fit, group = replicate), color = alpha("black", alpha = 0.2)) +
  facet_wrap(~setting, scales = "free", labeller = labeller(setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24))) +
  labs(title = "Initial - SNR = 1")

ggsave(filename = "residuals.pdf",
       plot = marrangeGrob(list(plt_res_flat2, plt_res_init2, plt_res_flat1, plt_res_init1), nrow = 1, ncol = 1),
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
  # plot(I$draws_tau[ , 1], type = "l", main = "tau_lin")
  # plot(I$draws_tau[ , 2], type = "l", main = "tau_nonlin")
  # err = colMeans(I$X_scaled %*% t(I$draws) - matrix(I$y, nrow = 1000, ncol = 1000, byrow = F))
  # plot(err, type = "l", main = "y - fit", xlab = "MCMC Iteration", ylab = "")
  # err_noInt = colMeans(I$X_scaled[ , -1] %*% t(I$draws[ , -1]) - matrix(I$y, nrow = 1000, ncol = 1000, byrow = F))
  # plot(err_noInt, type = "l", main = "y - fit (no beta0)", xlab = "MCMC Iteration", ylab = "")
  
  # Page title — written after all plots so it spans the full page
  page_title <- paste(paste0(sim_settings[s, ], collapse = " - "))
  mtext(page_title, outer = TRUE, cex = 1.2, font = 2, line = 1)
}
dev.off()

## Traceplots coefficients flattened ----
pdf(paste0(PATH_IMG, "traceplot_examples_flat.pdf"), width = 12, height = 8)
for (s in 1:nrow(sim_settings)){
  I = input_flat[[s]][[1]]
  
  # Add outer margin at the top for the page title
  par(mfrow = c(4, 6), mar = c(3, 2, 2, 1), oma = c(0, 0, 3, 0))
  for (j in 1:22){
    plot(I$draws[ , j], type = "l", 
         main = colnames(I$draws)[j], xlab = "MCMC Iteration", ylab = "")
    abline(h = 0, col = "red")
  }
  # plot(I$draws_tau[ , 1], type = "l", main = "tau_lin")
  # plot(I$draws_tau[ , 2], type = "l", main = "tau_nonlin")
  # err = colMeans(I$X_scaled %*% t(I$draws) - matrix(I$y, nrow = 1000, ncol = 1000, byrow = F))
  # plot(err, type = "l", main = "y - fit", xlab = "MCMC Iteration", ylab = "")
  # err_noInt = colMeans(I$X_scaled[ , -1] %*% t(I$draws[ , -1]) - matrix(I$y, nrow = 1000, ncol = 1000, byrow = F))
  # plot(err_noInt, type = "l", main = "y - fit (no beta0)", xlab = "MCMC Iteration", ylab = "")
  
  # Page title — written after all plots so it spans the full page
  page_title <- paste(paste0(sim_settings[s, ], collapse = " - "))
  mtext(page_title, outer = TRUE, cex = 1.2, font = 2, line = 1)
}
dev.off()


## ESS ----
input_mid = readRDS("sim_study_nonC-GBI/splines_univ/sim_study_nonCalibrated_asymm_TL_init_flattened.RDS")
# Init
df_ESS_init = sapply(input_init, function(I1) sapply(I1, function(I2) I2$diagn["ess_thin10", ])) %>% 
  melt() %>% rename(Beta = Var1, Replicate = Var2, ESS = value, Setting = L1)
df_ESS_flat = sapply(input_flat, function(I1) sapply(I1, function(I2) I2$diagn["ess_thin10", ])) %>% 
  melt() %>% rename(Beta = Var1, Replicate = Var2, ESS = value, Setting = L1)
df_ESS_mid = lapply(input_mid, function(I1) sapply(I1, function(I2) {
  if (!is.null(I2$diagn)){
    return(I2$diagn["ess_thin10", ])
  } else {
    return(rep(NA, 22))
  }})) %>% melt() %>% rename(Beta = Var1, Replicate = Var2, ESS = value, Setting = L1)
df_ESS = bind_rows(
  df_ESS_init %>% mutate(Approach = "Init"),
  df_ESS_flat %>% mutate(Approach = "Flat"),
  df_ESS_mid %>% mutate(Approach = "Mid")) %>% 
  mutate(ESS = if_else(is.na(ESS), 0, ESS),
         Approach = factor(Approach, levels = c("Init", "Mid", "Flat")))

plot_ESS_init_SNR2 = df_ESS %>% filter(Setting<=12) %>% 
  mutate(Setting = factor(Setting)) %>% 
  ggplot() +
  geom_boxplot(aes(x = Beta, y = ESS, col = Approach)) +
  facet_wrap(~Setting, labeller = labeller(Setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24)),
             nrow = 4, dir = "v") + 
  scale_color_manual(values = c("Init" = col_init, "Mid" = col_mid, "Flat" = col_flat), 
                     labels = c("k=init, c=0.001 [flatter proposal]", "k=sqrt(init), c=0.01 [flatter proposal]",
                                "k=sqrt(init), c=0.1 [proposal=full-cond]")) +
  scale_y_continuous(limits = c(0, 1000)) +
  labs(title = "ESS - SNR=2") + 
  theme(legend.position = "bottom",
        axis.text.x = element_blank(),
        text = element_text(size = 16))

plot_ESS_init_SNR1 =  df_ESS %>% filter(Setting<=12) %>% 
  mutate(Setting = factor(Setting)) %>% 
  ggplot() +
  geom_boxplot(aes(x = Beta, y = ESS, col = Approach)) +
  facet_wrap(~Setting, labeller = labeller(Setting = setNames(apply(sim_settings, 1, paste, collapse = " - "), 1:24)),
             nrow = 4, dir = "v") + 
  scale_color_manual(values = c("Init" = col_init, "Mid" = col_mid, "Flat" = col_flat), 
                     labels = c("k=init, c=0.001 [flatter proposal]", "k=sqrt(init), c=0.01 [flatter proposal]",
                                "k=sqrt(init), c=0.1 [proposal=full-cond]")) +
  scale_y_continuous(limits = c(0, 1000)) +
  labs(title = "ESS - SNR=1") + 
  theme(legend.position = "bottom",
        axis.text.x = element_blank(), 
        text = element_text(size = 16))

ggsave(filename = "ESS.pdf",
       plot = marrangeGrob(list(plot_ESS_init_SNR2, plot_ESS_init_SNR1), nrow = 1, ncol = 1, top = NULL),
       path = PATH_IMG, width = 24, height = 16)

