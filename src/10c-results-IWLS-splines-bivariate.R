rm(list = ls())

suppressPackageStartupMessages(library(tidyverse))
library(reshape2)
library(gridExtra)
source("src/01-utils.R")

input1000homo = readRDS("sim_study_nonC-GBI/splines_biv/ss_p2_n1000_homo.RDS")
input1000hetero = readRDS("sim_study_nonC-GBI/splines_biv/ss_p2_n1000_hetero.RDS")
input500homo = readRDS("sim_study_nonC-GBI/splines_biv/ss_p2_n500_homo.RDS")
input500hetero = readRDS("sim_study_nonC-GBI/splines_biv/ss_p2_n500_hetero.RDS")
input250homo = readRDS("sim_study_nonC-GBI/splines_biv/ss_p2_n250_homo.RDS")
input250hetero = readRDS("sim_study_nonC-GBI/splines_biv/ss_p2_n250_hetero.RDS")


PATH_IMG = "sim_study_nonC-GBI/splines_biv/"


# SIMULATION SETTINGS ----
error_distr = c("Gaussian", "Gamma")
cov_distr = "Unif_rad2" #c("Unif_rad2", "Beta_2_4") #"Gamma_1.5_1.5") #, "BVN")
transform_x1 = factor(c("parabola", "cubic", "trigonometric"), levels = c("parabola", "cubic", "trigonometric"), ordered = T)
transform_x2 = factor(c("parabola", "cubic", "trigonometric"), levels = c("parabola", "cubic", "trigonometric"), ordered = T)
SNR = c("High", "Low")
sim_settings = expand.grid(error_distr, cov_distr, transform_x1, transform_x2, SNR)
colnames(sim_settings) = c("err_distr", "cov_distr", "transf_x1", "transf_x2", "SNR")
sim_settings = sim_settings[sim_settings$transf_x1<sim_settings$transf_x2, ]  
sim_settings = as.matrix(sim_settings)

sim_settings_names = sim_settings
sim_settings_names[sim_settings_names[ , "transf_x1"] == "trigonometric", "transf_x1"] = "sin+exp"
sim_settings_names[sim_settings_names[ , "transf_x2"] == "trigonometric", "transf_x2"] = "sin+exp"
sim_settings_names[sim_settings_names[ , "transf_x1"] == "parabola", "transf_x1"] = "parab"
sim_settings_names[sim_settings_names[ , "transf_x2"] == "parabola", "transf_x2"] = "parab"
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


# LINEAR TRANSFORMATIONS ----
true_linear_components <- data.frame(
  transf = c("parabola", "parabola",
             "cubic", "cubic",
             "trigonometric", "trigonometric"),
  cov_distr = c("Unif_rad2", "Beta_2_4",
                "Unif_rad2", "Beta_2_4",
                "Unif_rad2", "Beta_2_4"),
  intercept = c(2/3, 10/21,
                0, -2*sqrt(2)/7,
                0.37599, 0.03975),
  slope = c(0, -2*sqrt(14)/21,
            2*sqrt(6)/5, 4*sqrt(7)/21,
            0.39919, 0.770398)
)


# RMSE AND COVERAGE ----
cov_RMSE_fun <- function(II, new_x, 
                         true_int1, true_int2, true_slope1, true_slope2,
                         true_nonlin1, true_nonlin2, true_f1, true_f2) {
  
  true_full1 = true_f1 - true_int1
  true_full2 = true_f2 - true_int2
  true_f = true_f1 + true_f2
  true_int = true_int1 + true_int2
  
  DR1 = II$DR1
  DR2 = II$DR2
  
  newB1 = mgcv::Predict.matrix(DR1$smooth_object, data.frame(x = new_x))
  newB2 = mgcv::Predict.matrix(DR2$smooth_object, data.frame(x = new_x))
  
  new_Z1 <- newB1 %*% DR1$Trafo
  new_Z2 <- newB2 %*% DR2$Trafo
  
  new_xtilde1 = (new_x - mean(II$x1))/sd(II$x1)
  new_xtilde2 = (new_x - mean(II$x2))/sd(II$x2)
  
  new_X = cbind(1, new_xtilde1, new_xtilde2, new_Z1, new_Z2)
  
  beta <- II$draws
  
  fit <- tcrossprod(beta, new_X)
  fit_int <- beta[ , 1]
  fit_slope1 <- beta[ , 2, drop = F]
  fit_slope2 <- beta[ , 3, drop = F]
  fit_nonlin1 <- tcrossprod(beta[, 4:23, drop = FALSE], new_Z1)
  fit_nonlin2 <- tcrossprod(beta[, 24:43, drop = FALSE], new_Z2)
  fit_full1 <- tcrossprod(fit_slope1, new_xtilde1) + fit_nonlin1
  fit_full2 <- tcrossprod(fit_slope2, new_xtilde2) + fit_nonlin2
  
  q_fit <- matrixStats::colQuantiles(fit, probs = c(0.025, 0.975), drop = FALSE)
  q_int <- quantile(fit_int, probs = c(0.025, 0.975))
  q_slope1 <- quantile(fit_slope1, probs = c(0.025, 0.975))
  q_slope2 <- quantile(fit_slope2, probs = c(0.025, 0.975))
  q_nonlin1 <- matrixStats::colQuantiles(fit_nonlin1, probs = c(0.025, 0.975), drop = FALSE)
  q_nonlin2 <- matrixStats::colQuantiles(fit_nonlin2, probs = c(0.025, 0.975), drop = FALSE)
  q_full1  <- matrixStats::colQuantiles(fit_full1, probs = c(0.025, 0.975), drop = FALSE)
  q_full2  <- matrixStats::colQuantiles(fit_full2, probs = c(0.025, 0.975), drop = FALSE)
  
  fit_mean <- colMeans(fit)
  fit_int_mean <- mean(fit_int)
  fit_nonlin1_mean <- colMeans(fit_nonlin1)
  fit_nonlin2_mean <- colMeans(fit_nonlin2)
  fit_slope1_mean <- mean(fit_slope1)
  fit_slope2_mean <- mean(fit_slope2)
  fit_full1_mean <- colMeans(fit_full1)
  fit_full2_mean <- colMeans(fit_full2)
  
  out_cov = rbind(
    cov_fit = as.numeric(q_fit[, 1] <= true_f & q_fit[, 2] >= true_f),
    cov_int = as.numeric(q_int[1] <= true_int & q_int[2] >= true_int),
    cov_slope1 = as.numeric(q_slope1[1] <= true_slope1 & q_slope1[2] >= true_slope1),
    cov_slope2 = as.numeric(q_slope2[1] <= true_slope2 & q_slope2[2] >= true_slope2),
    cov_nonlin1 = as.numeric(q_nonlin1[, 1] <= true_nonlin1 & q_nonlin1[, 2] >= true_nonlin1),
    cov_nonlin2 = as.numeric(q_nonlin2[, 1] <= true_nonlin2 & q_nonlin2[, 2] >= true_nonlin2),
    cov_full1 = as.numeric(q_full1[, 1] <= true_full1 & q_full1[, 2] >= true_full1),
    cov_full2 = as.numeric(q_full2[, 1] <= true_full2 & q_full2[, 2] >= true_full2))
  
  out_SE = rbind(SE_fit = (fit_mean - true_f)^2,
                 SE_int = rep((fit_int_mean - true_int)^2, npoints),
                 SE_slope1 = (fit_slope1_mean - true_slope1)^2,
                 SE_slope2 = (fit_slope2_mean - true_slope2)^2,
                 SE_nonlin1 = (fit_nonlin1_mean - true_nonlin1)^2,
                 SE_nonlin2 = (fit_nonlin2_mean - true_nonlin2)^2,
                 SE_full1 = (fit_full1_mean - true_full1)^2,
                 SE_full2 = (fit_full2_mean - true_full2)^2
  )
  
  out_fit = rbind(
    f1pf2 = fit_mean, int = rep(fit_int_mean, npoints), 
    slope1 = rep(fit_slope1_mean, npoints), slope2 = rep(fit_slope2_mean, npoints),
    nonlin1 = fit_nonlin1_mean, nonlin2 = fit_nonlin2_mean,
    full1 = fit_full1_mean, full2 = fit_full2_mean
  )
  
  return(list(cov = out_cov, SE = out_SE, fit = out_fit))
}

npoints = 101
new_x = seq(-sqrt(2), sqrt(2), length = npoints)


cov1000homo = cov500homo = cov250homo =
   cov1000hetero = cov500hetero = cov250hetero =
  array(NA, dim = c(8, nrow(sim_settings), npoints),
        dimnames = list("component" = c("f1+f2", "int", "slope1", "slope2", 
                                        "nonlin1", "nonlin2", "full1", "full2"), 
                        "setting" = 1:12, "point" = 1:npoints))

RMSE1000homo = RMSE500homo = RMSE250homo =
   RMSE1000hetero = RMSE500hetero = RMSE250hetero =
  array(NA, dim = c(8, nrow(sim_settings), npoints),
        dimnames = list("component" = c("f1+f2", "int", "slope1", "slope2", 
                                        "nonlin1", "nonlin2", "full1", "full2"), 
                        "setting" = 1:12, "point" = 1:npoints))

fit1000homo = fit500homo = fit250homo =
   fit1000hetero = fit500hetero = fit250hetero =
  array(NA, dim = c(8, nrow(sim_settings), npoints, 100),
        dimnames = list("component" = c("f1+f2", "int", "slope1", "slope2", 
                                        "nonlin1", "nonlin2", "full1", "full2"), 
                        "setting" = 1:12, "point" = 1:npoints, "replicate" = 1:100))

true_functions = 
  array(NA, dim = c(8, nrow(sim_settings), npoints),
        dimnames = list("component" = c("f1+f2", "int", "slope1", "slope2", 
                                        "nonlin1", "nonlin2", "full1", "full2"), 
                        "setting" = 1:12, "point" = 1:npoints))




for (s in 1:nrow(sim_settings)){
  cat(s, "\t")
  
  s_s = sim_settings[s, ]
  
  I1000homo = input1000homo$res[[s]]
  I1000hetero = input1000hetero$res[[s]]
  I500homo = input500homo$res[[s]]
  I500hetero = input500hetero$res[[s]]
  I250homo = input250homo$res[[s]]
  I250hetero = input250hetero$res[[s]]
  
  transf1 = switch(s_s["transf_x1"],
                   "parabola" = function(x) x^2,
                   "cubic" = function(x) x^3,
                   "trigonometric" = function(x) sin(2.5*x) + 3*exp(-(5^2)*((x)^2))
  )
  
  transf2 = switch(s_s["transf_x2"],
                   "parabola" = function(x) x^2,
                   "cubic" = function(x) x^3,
                   "trigonometric" = function(x) sin(2.5*x) + 3*exp(-(5^2)*((x)^2))
  )
  
  true_f1 = transf1(new_x)
  true_f2 = transf2(new_x)
  true_f = true_f1 + true_f2
  
  true_int1 = true_linear_components[true_linear_components$transf == s_s["transf_x1"] & true_linear_components$cov_distr == "Unif_rad2", 3]
  true_int2 = true_linear_components[true_linear_components$transf == s_s["transf_x2"] & true_linear_components$cov_distr == "Unif_rad2", 3]
  
  true_slope1 = true_linear_components[true_linear_components$transf == s_s["transf_x1"] & true_linear_components$cov_distr == "Unif_rad2", 4]
  true_slope2 = true_linear_components[true_linear_components$transf == s_s["transf_x2"] & true_linear_components$cov_distr == "Unif_rad2", 4]
  
  true_nonlin1 = true_f1 - true_int1 - true_slope1*new_x*sqrt(3/2) 
  true_nonlin2 = true_f2 - true_int2 - true_slope2*new_x*sqrt(3/2) 
  
  true_full1 = true_f1 - true_int1
  true_full2 = true_f2 - true_int2
  
  true_functions[ , s, ] <- rbind(true_f, true_int1+true_int2, true_slope1, true_slope2,
                                  true_nonlin1, true_nonlin2, true_full1, true_full2)
  
  
  res1000homo = lapply(I1000homo, function(II) 
    cov_RMSE_fun(II, new_x, 
                 true_int1 = true_int1, true_int2 = true_int2, 
                 true_slope1 = true_slope1, true_slope2 = true_slope2, 
                 true_nonlin1 = true_nonlin1, true_nonlin2 = true_nonlin2,
                 true_f1 = true_f1, true_f2 = true_f2)) 
  res1000hetero = lapply(I1000hetero, function(II)
    cov_RMSE_fun(II, new_x,
                 true_int1 = true_int1, true_int2 = true_int2,
                 true_slope1 = true_slope1, true_slope2 = true_slope2,
                 true_nonlin1 = true_nonlin1, true_nonlin2 = true_nonlin2,
                 true_f1 = true_f1, true_f2 = true_f2))
  res500homo = lapply(I500homo, function(II) 
    cov_RMSE_fun(II, new_x, 
                 true_int1 = true_int1, true_int2 = true_int2, 
                 true_slope1 = true_slope1, true_slope2 = true_slope2, 
                 true_nonlin1 = true_nonlin1, true_nonlin2 = true_nonlin2,
                 true_f1 = true_f1, true_f2 = true_f2))
  res500hetero = lapply(I500hetero, function(II)
    cov_RMSE_fun(II, new_x,
                 true_int1 = true_int1, true_int2 = true_int2,
                 true_slope1 = true_slope1, true_slope2 = true_slope2,
                 true_nonlin1 = true_nonlin1, true_nonlin2 = true_nonlin2,
                 true_f1 = true_f1, true_f2 = true_f2))
  res250homo = lapply(I250homo, function(II) 
    cov_RMSE_fun(II, new_x, 
                 true_int1 = true_int1, true_int2 = true_int2, 
                 true_slope1 = true_slope1, true_slope2 = true_slope2, 
                 true_nonlin1 = true_nonlin1, true_nonlin2 = true_nonlin2,
                 true_f1 = true_f1, true_f2 = true_f2))
  res250hetero = lapply(I250hetero, function(II)
    cov_RMSE_fun(II, new_x,
                 true_int1 = true_int1, true_int2 = true_int2,
                 true_slope1 = true_slope1, true_slope2 = true_slope2,
                 true_nonlin1 = true_nonlin1, true_nonlin2 = true_nonlin2,
                 true_f1 = true_f1, true_f2 = true_f2))
  
  cov1000homo[ , s, ] <- lapply(res1000homo, function(x) x$cov) %>% simplify2array() %>% apply(., 1:2, mean)
  RMSE1000homo[ , s, ] <- lapply(res1000homo, function(x) x$SE) %>% simplify2array() %>% apply(., 1:2, mean)
  fit1000homo[ , s, , ] <- lapply(res1000homo, function(x) x$fit) %>% simplify2array()
  cov500homo[ , s, ] <- lapply(res500homo, function(x) x$cov) %>% simplify2array() %>% apply(., 1:2, mean)
  RMSE500homo[ , s, ] <- lapply(res500homo, function(x) x$SE) %>% simplify2array() %>% apply(., 1:2, mean)
  fit500homo[ , s, , ] <- lapply(res500homo, function(x) x$fit) %>% simplify2array()
  cov250homo[ , s, ] <- lapply(res250homo, function(x) x$cov) %>% simplify2array() %>% apply(., 1:2, mean)
  RMSE250homo[ , s, ] <- lapply(res250homo, function(x) x$SE) %>% simplify2array() %>% apply(., 1:2, mean)
  fit250homo[ , s, , ] <- lapply(res250homo, function(x) x$fit) %>% simplify2array()
  
  cov1000hetero[ , s, ] <- lapply(res1000hetero, function(x) x$cov) %>% simplify2array() %>% apply(., 1:2, mean)
  RMSE1000hetero[ , s, ] <- lapply(res1000hetero, function(x) x$SE) %>% simplify2array() %>% apply(., 1:2, mean)
  fit1000hetero[ , s, , ] <- lapply(res1000hetero, function(x) x$fit) %>% simplify2array()
  cov500hetero[ , s, ] <- lapply(res500hetero, function(x) x$cov) %>% simplify2array() %>% apply(., 1:2, mean)
  RMSE500hetero[ , s, ] <- lapply(res500hetero, function(x) x$SE) %>% simplify2array() %>% apply(., 1:2, mean)
  fit500hetero[ , s, , ] <- lapply(res500hetero, function(x) x$fit) %>% simplify2array()
  cov250hetero[ , s, ] <- lapply(res250hetero, function(x) x$cov) %>% simplify2array() %>% apply(., 1:2, mean)
  RMSE250hetero[ , s, ] <- lapply(res250hetero, function(x) x$SE) %>% simplify2array() %>% apply(., 1:2, mean)
  fit250hetero[ , s, , ] <- lapply(res250hetero, function(x) x$fit) %>% simplify2array()
}


## Coverage ----
df_cov = bind_rows(
  cov1000homo %>% melt(value.name = "cov") %>% mutate(Run = "Homo1000"),
  cov1000hetero %>% melt(value.name = "cov") %>% mutate(Run = "Hetero1000"),
  cov500homo %>% melt(value.name = "cov") %>% mutate(Run = "Homo500"),
  cov500hetero %>% melt(value.name = "cov") %>% mutate(Run = "Hetero500"),
  cov250homo %>% melt(value.name = "cov") %>% mutate(Run = "Homo250"),
  cov250hetero %>% melt(value.name = "cov") %>% mutate(Run = "Hetero250")
) %>%
  mutate(Run = factor(Run, levels = c("Homo250", "Homo500", "Homo1000", "Hetero250", "Hetero500", "Hetero1000"))) %>% 
  left_join(cbind.data.frame(setting = 1:12, sim_settings), by = "setting") %>% 
  left_join(cbind.data.frame(point = 1:npoints, "x" = new_x), by = "point") %>% 
  mutate(transf_x1 = factor(transf_x1, levels = c("parabola", "cubic", "trigonometric"), ordered = T))


plt_cov_SNRH = 
  df_cov %>%
  filter(SNR == "High") %>%
  ggplot() +
  geom_hline(aes(col = "target", yintercept = 0.95), lty = 2) +
  geom_line(aes(x = x, y = cov, colour = Run)) +
  facet_grid(component ~ setting, scales = "fixed",
             labeller = labeller(setting = setNames(apply(sim_settings_names[ , c(1, 3:4)], 1, paste, collapse = " - "), 1:12))) +
  scale_color_manual("", values = color_runs) +
  scale_y_continuous(limits = c(0, 1)) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5))) +
  labs(title = "Coverage (nominal 95%) - [SNR 5]") +
  theme(legend.position = "bottom", legend.byrow = T)

plt_cov_SNRL =
  df_cov %>%
  filter(SNR == "Low") %>%
  ggplot() +
  geom_hline(aes(col = "target", yintercept = 0.95), lty = 2) +
  geom_line(aes(x = x, y = cov, colour = Run)) +
  facet_grid(component ~ setting, scales = "fixed",
             labeller = labeller(setting = setNames(apply(sim_settings_names[ , c(1, 3:4)], 1, paste, collapse = " - "), 1:12))) +
  scale_color_manual("", values = color_runs) +
  scale_y_continuous(limits = c(0, 1)) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5))) +
  labs(title = "Coverage (nominal 95%) - [SNR 2]") +
  theme(legend.position = "bottom", legend.byrow = T)

ggsave(filename = "coverage.pdf",
       plot = marrangeGrob(list(
         plt_cov_SNRH, plt_cov_SNRL
       ), 
       nrow = 1, ncol = 1, top = NULL),
       path = PATH_IMG, width = 12, height = 9)

## RMSE ----
df_RMSE = bind_rows(
  RMSE1000homo %>% melt(value.name = "RMSE") %>% mutate(Run = "Homo1000"),
  RMSE1000hetero %>% melt(value.name = "RMSE") %>% mutate(Run = "Hetero1000"),
  RMSE500homo %>% melt(value.name = "RMSE") %>% mutate(Run = "Homo500"),
  RMSE500hetero %>% melt(value.name = "RMSE") %>% mutate(Run = "Hetero500"),
  RMSE250homo %>% melt(value.name = "RMSE") %>% mutate(Run = "Homo250"),
  RMSE250hetero %>% melt(value.name = "RMSE") %>% mutate(Run = "Hetero250")
) %>%
  mutate(Run = factor(Run, levels = c("Homo250", "Homo500", "Homo1000", "Hetero250", "Hetero500", "Hetero1000"))) %>% 
  left_join(cbind.data.frame(setting = 1:12, sim_settings), by = "setting") %>% 
  left_join(cbind.data.frame(point = 1:npoints, "x" = new_x), by = "point") %>% 
  mutate(transf_x1 = factor(transf_x1, levels = c("parabola", "cubic", "trigonometric"), ordered = T))


plt_RMSE_SNRH = 
  df_RMSE %>%
  filter(SNR == "High") %>%
  ggplot() +
  geom_line(aes(x = x, y = RMSE, colour = Run)) +
  facet_grid(component ~ setting, scales = "free_y",
             labeller = labeller(setting = setNames(apply(sim_settings_names[ , c(1, 3:4)], 1, paste, collapse = " - "), 1:12))) +
  scale_color_manual("", values = color_runs) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5))) +
  labs(title = "RMSE - [SNR 5]") +
  theme(legend.position = "bottom", legend.byrow = T)

plt_RMSE_SNRL =
  df_RMSE %>%
  filter(SNR == "Low") %>%
  ggplot() +
  geom_line(aes(x = x, y = RMSE, colour = Run)) +
  facet_grid(component ~ setting, scales = "free_y",
             labeller = labeller(setting = setNames(apply(sim_settings_names[ , c(1, 3:4)], 1, paste, collapse = " - "), 1:12))) +
  scale_color_manual("", values = color_runs) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5))) +
  labs(title = "RMSE - [SNR 2]") +
  theme(legend.position = "bottom", legend.byrow = T)

ggsave(filename = "RMSE.pdf",
       plot = marrangeGrob(list(
         plt_RMSE_SNRH, plt_RMSE_SNRL
       ), 
       nrow = 1, ncol = 1, top = NULL),
       path = PATH_IMG, width = 12, height = 9)


## Fit plot ----
df_fit1000 = bind_rows(
  fit1000homo %>% melt(value.name = "fit") %>% mutate(Run = "Homo1000"),
  fit1000hetero %>% melt(value.name = "fit") %>% mutate(Run = "Hetero1000")
) %>%
  mutate(Run = factor(Run, levels = c("Homo250", "Homo500", "Homo1000", "Hetero250", "Hetero500", "Hetero1000"))) %>% 
  left_join(cbind.data.frame(setting = 1:12, sim_settings), by = "setting") %>% 
  left_join(cbind.data.frame(point = 1:npoints, "x" = new_x), by = "point") %>% 
  mutate(transf_x1 = factor(transf_x1, levels = c("parabola", "cubic", "trigonometric"), ordered = T))

df_fit500 = bind_rows(
  fit500homo %>% melt(value.name = "fit") %>% mutate(Run = "Homo500"),
  fit500hetero %>% melt(value.name = "fit") %>% mutate(Run = "Hetero500")
) %>%
  mutate(Run = factor(Run, levels = c("Homo250", "Homo500", "Homo1000", "Hetero250", "Hetero500", "Hetero1000"))) %>% 
  left_join(cbind.data.frame(setting = 1:12, sim_settings), by = "setting") %>% 
  left_join(cbind.data.frame(point = 1:npoints, "x" = new_x), by = "point") %>% 
  mutate(transf_x1 = factor(transf_x1, levels = c("parabola", "cubic", "trigonometric"), ordered = T))

df_fit250 = bind_rows(
  fit250homo %>% melt(value.name = "fit") %>% mutate(Run = "Homo250"),
  fit250hetero %>% melt(value.name = "fit") %>% mutate(Run = "Hetero250")
) %>%
  mutate(Run = factor(Run, levels = c("Homo250", "Homo500", "Homo1000", "Hetero250", "Hetero500", "Hetero1000"))) %>% 
  left_join(cbind.data.frame(setting = 1:12, sim_settings), by = "setting") %>% 
  left_join(cbind.data.frame(point = 1:npoints, "x" = new_x), by = "point") %>% 
  mutate(transf_x1 = factor(transf_x1, levels = c("parabola", "cubic", "trigonometric"), ordered = T))

df_true = left_join(true_functions %>% melt(value.name = "true"),
                    data.frame(point = 1:npoints, x = new_x),
                    by = "point") %>% 
  left_join(cbind.data.frame(setting = 1:12, sim_settings), by = "setting") %>% 
  mutate(transf_x1 = factor(transf_x1, levels = c("parabola", "cubic", "trigonometric"), ordered = T))

### 1000 ----
plt_fit1000_SNRH = 
  df_fit1000 %>%
  filter(SNR == "High") %>%
  ggplot() +
  geom_line(aes(x = x, y = fit, colour = Run, group = interaction(Run, replicate)), alpha = 0.25) +
  geom_line(data = df_true %>% filter(SNR == "High", transf_x1 < transf_x2),
            mapping = aes(x = x, y = true), col = "black") +
  facet_grid(component ~ setting, scales = "free_y",
             labeller = labeller(setting = setNames(apply(sim_settings_names[ , c(1, 3:4)], 1, paste, collapse = " - "), 1:12))) +
  scale_color_manual("", values = color_runs) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5))) +
  labs(title = "fit vs true - n = 1000 [SNR 5]") +
  theme(legend.position = "bottom", legend.byrow = T)

plt_fit1000_SNRL = 
  df_fit1000 %>%
  filter(SNR == "Low") %>%
  ggplot() +
  geom_line(aes(x = x, y = fit, colour = Run, group = interaction(Run, replicate)), alpha = 0.25) +
  geom_line(data = df_true %>% filter(SNR == "Low", transf_x1 < transf_x2),
            mapping = aes(x = x, y = true), col = "black") +
  facet_grid(component ~ setting, scales = "free_y",
             labeller = labeller(setting = setNames(apply(sim_settings_names[ , c(1, 3:4)], 1, paste, collapse = " - "), 1:12))) +
  scale_color_manual("", values = color_runs) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5))) +
  labs(title = "fit vs true - n = 1000 [SNR 2]") +
  theme(legend.position = "bottom", legend.byrow = T)

### 500 ----
plt_fit500_SNRH = 
  df_fit500 %>%
  filter(SNR == "High") %>%
  ggplot() +
  geom_line(aes(x = x, y = fit, colour = Run, group = interaction(Run, replicate)), alpha = 0.25) +
  geom_line(data = df_true %>% filter(SNR == "High", transf_x1 < transf_x2),
            mapping = aes(x = x, y = true), col = "black") +
  facet_grid(component ~ setting, scales = "free_y",
             labeller = labeller(setting = setNames(apply(sim_settings_names[ , c(1, 3:4)], 1, paste, collapse = " - "), 1:12))) +
  scale_color_manual("", values = color_runs) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5))) +
  labs(title = "fit vs true - n = 500 [SNR 5]") +
  theme(legend.position = "bottom", legend.byrow = T)

plt_fit500_SNRL = 
  df_fit1000 %>%
  filter(SNR == "Low", transf_x1 < transf_x2) %>%
  ggplot() +
  geom_line(aes(x = x, y = fit, colour = Run, group = interaction(Run, replicate)), alpha = 0.25) +
  geom_line(data = df_true %>% filter(SNR == "Low", transf_x1 < transf_x2),
            mapping = aes(x = x, y = true), col = "black") +
  facet_grid(component ~ setting, scales = "free_y",
             labeller = labeller(setting = setNames(apply(sim_settings_names[ , c(1, 3:4)], 1, paste, collapse = " - "), 1:12))) +
  scale_color_manual("", values = color_runs) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5))) +
  labs(title = "fit vs true - n = 1000 [SNR 2]") +
  theme(legend.position = "bottom", legend.byrow = T)

### 250 ----
plt_fit250_SNRH = 
  df_fit250 %>%
  filter(SNR == "High", transf_x1 < transf_x2) %>%
  ggplot() +
  geom_line(aes(x = x, y = fit, colour = Run, group = interaction(Run, replicate)), alpha = 0.25) +
  geom_line(data = df_true %>% filter(SNR == "High", transf_x1 < transf_x2),
            mapping = aes(x = x, y = true), col = "black") +
  facet_grid(component ~ setting, scales = "free_y",
             labeller = labeller(setting = setNames(apply(sim_settings_names[ , c(1, 3:4)], 1, paste, collapse = " - "), 1:12))) +
  scale_color_manual("", values = color_runs) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5))) +
  labs(title = "fit vs true - n = 250 [SNR 5]") +
  theme(legend.position = "bottom", legend.byrow = T)

plt_fit250_SNRL = 
  df_fit250 %>%
  filter(SNR == "Low", transf_x1 < transf_x2) %>%
  ggplot() +
  geom_line(aes(x = x, y = fit, colour = Run, group = interaction(Run, replicate)), alpha = 0.25) +
  geom_line(data = df_true %>% filter(SNR == "Low", transf_x1 < transf_x2),
            mapping = aes(x = x, y = true), col = "black") +
  facet_grid(component ~ setting, scales = "free_y",
             labeller = labeller(setting = setNames(apply(sim_settings_names[ , c(1, 3:4)], 1, paste, collapse = " - "), 1:12))) +
  scale_color_manual("", values = color_runs) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5))) +
  labs(title = "fit vs true - n = 250 [SNR 2]") +
  theme(legend.position = "bottom", legend.byrow = T)

### save ----
ggsave(filename = "fit.pdf",
       plot = marrangeGrob(list(
         plt_fit1000_SNRH, plt_fit1000_SNRL,
         plt_fit500_SNRH, plt_fit500_SNRL,
         plt_fit250_SNRH, plt_fit250_SNRL 
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
  df_ESS250hetero %>% mutate(Run = "Hetero250")
) %>% 
  mutate(ESS = if_else(is.na(ESS), 0, ESS),
         Run = factor(Run, levels = c("Homo1000", "Hetero1000", "Homo500", "Hetero500", "Homo250", "Hetero250"))) %>% 
  left_join(cbind.data.frame(Setting = 1:12, sim_settings), by = "Setting")


plt_ESS_SNRH = 
  df_ESS %>% filter(SNR == "High") %>%
  ggplot() + 
  geom_boxplot(aes(y = ESS, x = Beta, col = Run)) +
  geom_hline(aes(yintercept = 1000), col = 2, lty = 2) +
  facet_grid(Run~Setting, labeller = labeller(Setting = setNames(apply(sim_settings_names[ , c(1, 3:4)], 1, paste, collapse = " - "), 1:12))) +
  scale_color_manual("", values = color_runs) +
  scale_y_continuous(limits = c(0, 1500)) +
  labs(title = "ESS - [SNR = 5]") +
  theme(legend.position = "bottom",
        axis.text.x = element_blank(),
        text = element_text(size = 14))

plt_ESS_SNRL = 
  df_ESS %>% filter(SNR == "Low") %>%
  ggplot() + 
  geom_boxplot(aes(y = ESS, x = Beta, col = Run)) +
  geom_hline(aes(yintercept = 1000), col = 2, lty = 2) +
  facet_grid(Run~Setting, labeller = labeller(Setting = setNames(apply(sim_settings_names[ , c(1, 3:4)], 1, paste, collapse = " - "), 1:12))) +
  scale_color_manual("", values = color_runs) +
  scale_y_continuous(limits = c(0, 1500)) +
  labs(title = "ESS [SNR = 2]") +
  theme(legend.position = "bottom",
        axis.text.x = element_blank(),
        text = element_text(size = 14))


ggsave(filename = "ESS_n.pdf",
       plot = marrangeGrob(list(plt_ESS_SNRH, plt_ESS_SNRL), nrow = 1, ncol = 1, top = NULL),
       path = PATH_IMG, width = 15, height = 15)


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
  df_RHAT250hetero %>% mutate(Run = "Hetero250")
) %>%
  mutate(Run = factor(Run, levels = c("Homo1000", "Hetero1000", "Homo500",  "Hetero500", "Homo250", "Hetero250"))) %>% 
  left_join(cbind.data.frame(setting = 1:12, sim_settings), by = "setting")

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
  facet_wrap(~setting, labeller = labeller(setting = setNames(apply(sim_settings[ , -c(2, 5)], 1, paste, collapse = " - "), 1:12)),
             ncol = 3, dir = "h") +
  scale_color_manual("", values = color_runs) +
  scale_y_continuous(limits = c(0, 1)) +
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
  facet_wrap(~setting, labeller = labeller(setting = setNames(apply(sim_settings[ , -c(2, 5)], 1, paste, collapse = " - "), 1:12)),
             ncol = 3, dir = "h") +
  scale_color_manual("", values = color_runs) +
  scale_y_continuous(limits = c(0, 1)) +
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
  df_w250hetero %>% mutate(Run = "Hetero250")
) %>% 
  mutate(Run = factor(Run, levels = c("Homo250", "Homo500", "Homo1000", "Hetero250", "Hetero500", "Hetero1000")))


plt_w = df_w %>% 
  ggplot() + 
  geom_boxplot(aes(y = w, x = Run, col = Run)) +
  facet_wrap(~Setting, scales = "fixed",
             labeller = labeller(Setting = setNames(apply(sim_settings_names[ , -4], 1, paste, collapse = " - "), 1:12)),
             nrow = 4, dir = "h") +
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
  df_k250hetero %>% mutate(Run = "Hetero250")
) %>% 
  mutate(Run = factor(Run, levels = c("Homo250", "Homo500", "Homo1000", "Hetero250", "Hetero500", "Hetero1000")),
         Side = factor(Side, levels = c(1, 2), labels = c("Left", "Right")))


plt_k = df_k %>% 
  ggplot() + 
  geom_boxplot(aes(y = k, x = interaction(Side, Run), col = Run)) +
  facet_wrap(~Setting, scales = "fixed",
             labeller = labeller(Setting = setNames(apply(sim_settings_names[ , -2], 1, paste, collapse = " - "), 1:12)),
             nrow = 4, dir = "h") +
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
