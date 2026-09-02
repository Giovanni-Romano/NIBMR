rm(list = ls())

suppressPackageStartupMessages(library(tidyverse))

input_TL = readRDS("sim_study_nonC-GBI/02-w_unit_loss_SNR/sim_study_nonCalibrated_multiv_asymm_TL.RDS")
input_sd = readRDS("sim_study_nonC-GBI/02-w_unit_loss_SNR/sim_study_nonCalibrated_multiv_asymm_sd.RDS")
input_TL_init = readRDS("sim_study_nonC-GBI/02-w_unit_loss_SNR/sim_study_nonCalibrated_multiv_asymm_TL_init.RDS")
input_sd_init = readRDS("sim_study_nonC-GBI/02-w_unit_loss_SNR/sim_study_nonCalibrated_multiv_asymm_sd_init.RDS")

results_TL = lapply(input_TL,
                    function(I)
                      abind::abind(lapply(I, function(x) x$summ), along = 3))
results_sd = lapply(input_sd,
                    function(I)
                      abind::abind(lapply(I, function(x) x$summ), along = 3))
apply(input_TL_init,
      function(I)
        abind::abind(lapply(I, function(x) x$summ), along = 3))
results_TL_init = lapply(input_TL_init,
                         function(I)
                           abind::abind(lapply(I, function(x) x$summ), along = 3))
results_sd_init = lapply(input_sd_init,
                         function(I)
                           abind::abind(lapply(I, function(x) x$summ), along = 3))


PATH_IMG = "sim_study_nonC-GBI/02-w_unit_loss_SNR/"

# TRUE BETA AND SAMPLE SIZE ----
n = 1e3
true_beta = c(0, -1, +1.5)

# SIMULATION SETTINGS ----
idx_iter = seq(10, 10^4, by = 10)
error_distr = c("Gaussian", "Gamma_2_rad2", "Gaussian_rad2", "Gamma_2_1")
covar_distr = c("Unif_rad2", "Gamma_1.5_1.5", "BVN")
SNR = paste0("SNR", rep(c(2, 2, 1, 1), 3))
sim_settings = cbind(as.matrix(expand.grid(error_distr, covar_distr)), SNR)
colnames(sim_settings) = c("err_distr", "cov_distr", "SNR")
order_runs = c(1, 2, 5, 6, 9, 10, 3, 4, 7, 8, 11, 12)

# Colors
col_sd = "#d00000"
col_TL = "#023e8a"


# Plot point estimates ----
## Mean ----
# First regression coefficient
pdf(paste0(PATH_IMG, "multiv_beta1.pdf"), width = 9, height = 6)
par(mfrow = c(3, 4), mar = c(3, 2, 2, 1))
for (j in order_runs){
  min_x = -2
  max_x = 0
  plot(1, xlim = c(min_x, max_x), ylim = c(0, 12), type = "n", main = paste0("b1 | k final | ", paste(sim_settings[j, ], collapse = " - ")),
       xlab = "", ylab = "", cex.main = 0.75)
  lines(density(results_sd[[j]]["mean", 2, ]), col = col_sd)
  lines(density(results_sd[[j]]["SANN", 2, ]), col = col_sd, lty = 2)
  lines(density(results_TL[[j]]["mean", 2, ]), col = col_TL)
  lines(density(results_TL[[j]]["SANN", 2, ]), col = col_TL, lty = 2)
  abline(v = true_beta[2], col = "gold")
  legend("topright", col = c(col_sd, col_sd, "gold"), lty = c(1, 2, 1),
         legend = c("Bayes sd", "Freq sd", "True"), cex = 0.8)
  legend("topleft", col = c(col_TL), lty = c(1, 2),
         legend = c("Bayes TL", "Freq TL"), cex = 0.8)
  
  plot(1, xlim = c(min_x, max_x), ylim = c(0, 12), type = "n", main = paste0("b1 | k init | ", paste(sim_settings[j, ], collapse = " - ")),
       xlab = "", ylab = "", cex.main = 0.75)
  lines(density(results_sd_init[[j]]["mean", 2, ]), col = col_sd)
  lines(density(results_sd_init[[j]]["SANN", 2, ]), col = col_sd, lty = 2)
  lines(density(results_TL_init[[j]]["mean", 2, ]), col = col_TL)
  lines(density(results_TL_init[[j]]["SANN", 2, ]), col = col_TL, lty = 2)
  abline(v = true_beta[2], col = "gold")
  legend("topright", col = c(col_sd, col_sd, "gold"), lty = c(1, 2, 1),
         legend = c("Bayes sd", "Freq sd init", "True"), cex = 0.8)
  legend("topleft", col = c(col_TL), lty = c(1, 2),
         legend = c("Bayes TL", "Freq TL"), cex = 0.8)
}
dev.off()

# Second regression coefficient
pdf(paste0(PATH_IMG, "multiv_beta2.pdf"), width = 9, height = 6)
par(mfrow = c(3, 4), mar = c(3, 2, 2, 1))
for (j in order_runs){
  min_x = 0.5
  max_x = 2.5
  plot(1, xlim = c(min_x, max_x), ylim = c(0, 12), type = "n", main = paste0("b2 | k final | ", paste(sim_settings[j, ], collapse = " - ")),
       xlab = "", ylab = "", cex.main = 0.75)
  lines(density(results_sd[[j]]["mean", 3, ]), col = col_sd)
  lines(density(results_sd[[j]]["SANN", 3, ]), col = col_sd, lty = 2)
  lines(density(results_TL[[j]]["mean", 3, ]), col = col_TL)
  lines(density(results_TL[[j]]["SANN", 3, ]), col = col_TL, lty = 2)
  abline(v = true_beta[3], col = "gold")
  legend("topright", col = c(col_sd, col_sd, "gold"), lty = c(1, 2, 1),
         legend = c("Bayes sd", "Freq sd", "True"), cex = 0.8)
  legend("topleft", col = c(col_TL), lty = c(1, 2),
         legend = c("Bayes TL", "Freq TL"), cex = 0.8)
  
  plot(1, xlim = c(min_x, max_x), ylim = c(0, 12), type = "n", main = paste0("b2 | k init | ", paste(sim_settings[j, ], collapse = " - ")),
       xlab = "", ylab = "", cex.main = 0.75)
  lines(density(results_sd_init[[j]]["mean", 3, ]), col = col_sd)
  lines(density(results_sd_init[[j]]["SANN", 3, ]), col = col_sd, lty = 2)
  lines(density(results_TL_init[[j]]["mean", 3, ]), col = col_TL)
  lines(density(results_TL_init[[j]]["SANN", 3, ]), col = col_TL, lty = 2)
  abline(v = true_beta[3], col = "gold")
  legend("topright", col = c(col_sd, col_sd, "gold"), lty = c(1, 2, 1),
         legend = c("Bayes sd", "Freq sd init", "True"), cex = 0.8)
  legend("topleft", col = c(col_TL), lty = c(1, 2),
         legend = c("Bayes TL", "Freq TL"), cex = 0.8)
}
dev.off()

# Intercept
pdf(paste0(PATH_IMG, "multiv_beta0.pdf"), width = 9, height = 6)
par(mfrow = c(3, 4), mar = c(3, 2, 2, 1))
for (j in order_runs){
  min_x = -1
  max_x = +1
  plot(1, xlim = c(min_x, max_x), ylim = c(0, 12), type = "n", main = paste0("b0 | k final | ", paste(sim_settings[j, ], collapse = " - ")),
       xlab = "", ylab = "", cex.main = 0.75)
  lines(density(results_sd[[j]]["mean", 1, ]), col = col_sd)
  lines(density(results_sd[[j]]["SANN", 1, ]), col = col_sd, lty = 2)
  lines(density(results_TL[[j]]["mean", 1, ]), col = col_TL)
  lines(density(results_TL[[j]]["SANN", 1, ]), col = col_TL, lty = 2)
  abline(v = true_beta[1], col = "gold")
  legend("topright", col = c(col_sd, col_sd, "gold"), lty = c(1, 2, 1),
         legend = c("Bayes sd", "Freq sd", "True"), cex = 0.8)
  legend("topleft", col = c(col_TL), lty = c(1, 2),
         legend = c("Bayes TL", "Freq TL"), cex = 0.8)
  
  plot(1, xlim = c(min_x, max_x), ylim = c(0, 12), type = "n", main = paste0("b0 | k init | ", paste(sim_settings[j, ], collapse = " - ")),
       xlab = "", ylab = "", cex.main = 0.75)
  lines(density(results_sd_init[[j]]["mean", 1, ]), col = col_sd)
  lines(density(results_sd_init[[j]]["SANN", 1, ]), col = col_sd, lty = 2)
  lines(density(results_TL_init[[j]]["mean", 1, ]), col = col_TL)
  lines(density(results_TL_init[[j]]["SANN", 1, ]), col = col_TL, lty = 2)
  abline(v = true_beta[1], col = "gold")
  legend("topright", col = c(col_sd, col_sd, "gold"), lty = c(1, 2, 1),
         legend = c("Bayes sd", "Freq sd init", "True"), cex = 0.8)
  legend("topleft", col = c(col_TL), lty = c(1, 2),
         legend = c("Bayes TL", "Freq TL"), cex = 0.8)
}
dev.off()


# Credible intervals frequentist coverage ----
## First regression coefficient ----
cov_beta1_quant_sd = t(sapply(results_sd, function(R) {
  I = R[c("q025", "q975"), 2, ]
  rowMeans(apply(I, 2, function(i) unname(c(true_beta[2] >= i[1] & true_beta[2] <= i[2], i[2]-i[1]))))
}))
colnames(cov_beta1_quant_sd) = c("cov_q_sd", "len_q_sd")
cov_beta1_quant_TL = t(sapply(results_TL, function(R) {
  I = R[c("q025", "q975"), 2, ]
  rowMeans(apply(I, 2, function(i) unname(c(true_beta[2] >= i[1] & true_beta[2] <= i[2], i[2]-i[1]))))
}))
colnames(cov_beta1_quant_TL) = c("cov_q_TL", "len_q_TL")

cov_beta1_quant_sd_init = t(sapply(results_sd_init, function(R) {
  I = R[c("q025", "q975"), 2, ]
  rowMeans(apply(I, 2, function(i) unname(c(true_beta[2] >= i[1] & true_beta[2] <= i[2], i[2]-i[1]))))
}))
colnames(cov_beta1_quant_sd_init) = c("cov_q_sd", "len_q_sd")
cov_beta1_quant_TL_init = t(sapply(results_TL_init, function(R) {
  I = R[c("q025", "q975"), 2, ]
  rowMeans(apply(I, 2, function(i) unname(c(true_beta[2] >= i[1] & true_beta[2] <= i[2], i[2]-i[1]))))
}))
colnames(cov_beta1_quant_TL_init) = c("cov_q_TL", "len_q_TL")

cov_beta1_HDI_sd = t(sapply(results_sd, function(R) {
  I1 = R[c("HDI_L1", "HDI_U1"), 2, ]; I2 = R[c("HDI_L2", "HDI_U2"), 2, ]
  tmp1 = apply(I1, 2, function(i) true_beta[2] >= i[1] & true_beta[2] <= i[2])
  tmp2 = apply(I2, 2, function(i) true_beta[2] >= i[1] & true_beta[2] <= i[2])
  tmp2[is.na(tmp2)] = FALSE
  len1 = diff(I1); len2 = ifelse(is.na(diff(I2)), 0, diff(I2))
  c(mean(tmp1 | tmp2), mean(len1+len2))
}))
colnames(cov_beta1_HDI_sd) = c("cov_HDI_sd", "len_HDI_sd")
cov_beta1_HDI_TL = t(sapply(results_TL, function(R) {
  I1 = R[c("HDI_L1", "HDI_U1"), 2, ]; I2 = R[c("HDI_L2", "HDI_U2"), 2, ]
  tmp1 = apply(I1, 2, function(i) true_beta[2] >= i[1] & true_beta[2] <= i[2])
  tmp2 = apply(I2, 2, function(i) true_beta[2] >= i[1] & true_beta[2] <= i[2])
  tmp2[is.na(tmp2)] = FALSE
  len1 = diff(I1); len2 = ifelse(is.na(diff(I2)), 0, diff(I2))
  c(mean(tmp1 | tmp2), mean(len1+len2))
}))
colnames(cov_beta1_HDI_TL) = c("cov_HDI_TL", "len_HDI_TL")

cov_beta1_HDI_sd_init = t(sapply(results_sd_init, function(R) {
  I1 = R[c("HDI_L1", "HDI_U1"), 2, ]; I2 = R[c("HDI_L2", "HDI_U2"), 2, ]
  tmp1 = apply(I1, 2, function(i) true_beta[2] >= i[1] & true_beta[2] <= i[2])
  tmp2 = apply(I2, 2, function(i) true_beta[2] >= i[1] & true_beta[2] <= i[2])
  tmp2[is.na(tmp2)] = FALSE
  len1 = diff(I1); len2 = ifelse(is.na(diff(I2)), 0, diff(I2))
  c(mean(tmp1 | tmp2), mean(len1+len2))
}))
colnames(cov_beta1_HDI_sd_init) = c("cov_HDI_sd", "len_HDI_sd")
cov_beta1_HDI_TL_init = t(sapply(results_TL_init, function(R) {
  I1 = R[c("HDI_L1", "HDI_U1"), 2, ]; I2 = R[c("HDI_L2", "HDI_U2"), 2, ]
  tmp1 = apply(I1, 2, function(i) true_beta[2] >= i[1] & true_beta[2] <= i[2])
  tmp2 = apply(I2, 2, function(i) true_beta[2] >= i[1] & true_beta[2] <= i[2])
  tmp2[is.na(tmp2)] = FALSE
  len1 = diff(I1); len2 = ifelse(is.na(diff(I2)), 0, diff(I2))
  c(mean(tmp1 | tmp2), mean(len1+len2))
}))
colnames(cov_beta1_HDI_TL_init) = c("cov_HDI_TL", "len_HDI_TL")

df_beta1 = data.frame(sim_settings, 
                      cov_beta1_quant_sd, cov_beta1_quant_TL,
                      cov_beta1_HDI_sd, cov_beta1_HDI_TL)
df_beta1_init = data.frame(sim_settings, 
                           cov_beta1_quant_sd_init, cov_beta1_quant_TL_init,
                           cov_beta1_HDI_sd_init, cov_beta1_HDI_TL_init)

tab = knitr::kable(df_beta1[order_runs , c(1:3, 4, 6, 5, 7, 8, 10, 9, 11)], format = "pipe", digits = 3)
tab_init = knitr::kable(df_beta1_init[order_runs , c(1:3, 4, 6, 5, 7, 8, 10, 9, 11)], format = "pipe", digits = 3)
writeLines(c("k final", tab, "\n\n k init", tab_init), paste0(PATH_IMG, "multiv_beta1.txt"))

## Second regression coefficient ----
cov_beta2_quant_sd= t(sapply(results_sd, function(R) {
  I = R[c("q025", "q975"), 3, ]
  rowMeans(apply(I, 2, function(i) unname(c(true_beta[3] >= i[1] & true_beta[3] <= i[2], i[2]-i[1]))))
}))
colnames(cov_beta2_quant_sd) = c("cov_q_sd", "len_q_sd")
cov_beta2_quant_TL = t(sapply(results_TL, function(R) {
  I = R[c("q025", "q975"), 3, ]
  rowMeans(apply(I, 2, function(i) unname(c(true_beta[3] >= i[1] & true_beta[3] <= i[2], i[2]-i[1]))))
}))
colnames(cov_beta2_quant_TL) = c("cov_q_TL", "len_q_TL")
cov_beta2_quant_sd_init = t(sapply(results_sd_init, function(R) {
  I = R[c("q025", "q975"), 3, ]
  rowMeans(apply(I, 2, function(i) unname(c(true_beta[3] >= i[1] & true_beta[3] <= i[2], i[2]-i[1]))))
}))
colnames(cov_beta2_quant_sd_init) = c("cov_q_sd", "len_q_sd")
cov_beta2_quant_TL_init = t(sapply(results_TL_init, function(R) {
  I = R[c("q025", "q975"), 3, ]
  rowMeans(apply(I, 2, function(i) unname(c(true_beta[3] >= i[1] & true_beta[3] <= i[2], i[2]-i[1]))))
}))
colnames(cov_beta2_quant_TL_init) = c("cov_q_TL", "len_q_TL")

cov_beta2_HDI_sd = t(sapply(results_sd, function(R) {
  I1 = R[c("HDI_L1", "HDI_U1"), 3, ]; I2 = R[c("HDI_L2", "HDI_U2"), 3, ]
  tmp1 = apply(I1, 2, function(i) true_beta[3] >= i[1] & true_beta[3] <= i[2])
  tmp2 = apply(I2, 2, function(i) true_beta[3] >= i[1] & true_beta[3] <= i[2])
  tmp2[is.na(tmp2)] = FALSE
  len1 = diff(I1); len2 = ifelse(is.na(diff(I2)), 0, diff(I2))
  c(mean(tmp1 | tmp2), mean(len1+len2))
}))
colnames(cov_beta2_HDI_sd) = c("cov_HDI_sd", "len_HDI_sd")
cov_beta2_HDI_TL = t(sapply(results_TL, function(R) {
  I1 = R[c("HDI_L1", "HDI_U1"), 3, ]; I2 = R[c("HDI_L2", "HDI_U2"), 3, ]
  tmp1 = apply(I1, 2, function(i) true_beta[3] >= i[1] & true_beta[3] <= i[2])
  tmp2 = apply(I2, 2, function(i) true_beta[3] >= i[1] & true_beta[3] <= i[2])
  tmp2[is.na(tmp2)] = FALSE
  len1 = diff(I1); len2 = ifelse(is.na(diff(I2)), 0, diff(I2))
  c(mean(tmp1 | tmp2), mean(len1+len2))
}))
colnames(cov_beta2_HDI_TL) = c("cov_HDI_TL", "len_HDI_TL")
cov_beta2_HDI_sd_init = t(sapply(results_sd_init, function(R) {
  I1 = R[c("HDI_L1", "HDI_U1"), 3, ]; I2 = R[c("HDI_L2", "HDI_U2"), 3, ]
  tmp1 = apply(I1, 2, function(i) true_beta[3] >= i[1] & true_beta[3] <= i[2])
  tmp2 = apply(I2, 2, function(i) true_beta[3] >= i[1] & true_beta[3] <= i[2])
  tmp2[is.na(tmp2)] = FALSE
  len1 = diff(I1); len2 = ifelse(is.na(diff(I2)), 0, diff(I2))
  c(mean(tmp1 | tmp2), mean(len1+len2))
}))
colnames(cov_beta2_HDI_sd_init) = c("cov_HDI_sd", "len_HDI_sd")
cov_beta2_HDI_TL_init = t(sapply(results_TL_init, function(R) {
  I1 = R[c("HDI_L1", "HDI_U1"), 3, ]; I2 = R[c("HDI_L2", "HDI_U2"), 3, ]
  tmp1 = apply(I1, 2, function(i) true_beta[3] >= i[1] & true_beta[3] <= i[2])
  tmp2 = apply(I2, 2, function(i) true_beta[3] >= i[1] & true_beta[3] <= i[2])
  tmp2[is.na(tmp2)] = FALSE
  len1 = diff(I1); len2 = ifelse(is.na(diff(I2)), 0, diff(I2))
  c(mean(tmp1 | tmp2), mean(len1+len2))
}))
colnames(cov_beta2_HDI_TL_init) = c("cov_HDI_TL", "len_HDI_TL")

df_beta2 = data.frame(sim_settings, 
                      cov_beta2_quant_sd, cov_beta2_quant_TL,
                      cov_beta2_HDI_sd, cov_beta2_HDI_TL)
df_beta2_init = data.frame(sim_settings, 
                           cov_beta2_quant_sd_init, cov_beta2_quant_TL_init,
                           cov_beta2_HDI_sd_init, cov_beta2_HDI_TL_init)

tab = knitr::kable(df_beta2[order_runs ,  c(1:3, 4, 6, 5, 7, 8, 10, 9, 11)], format = "pipe", digits = 3)
tab_init = knitr::kable(df_beta2_init[order_runs ,  c(1:3, 4, 6, 5, 7, 8, 10, 9, 11)], format = "pipe", digits = 3)
writeLines(c("k final", tab, "\n\n k init", tab_init), paste0(PATH_IMG, "multiv_beta2.txt"))

## Intercept ----
cov_int_quant_sd = t(sapply(results_sd, function(R) {
  I = R[c("q025", "q975"), 1, ]
  rowMeans(apply(I, 2, function(i) unname(c(true_beta[1] >= i[1] & true_beta[1] <= i[2], i[2]-i[1]))))
}))
colnames(cov_int_quant_sd) = c("cov_q_sd", "len_q_sd")
cov_int_quant_TL = t(sapply(results_TL, function(R) {
  I = R[c("q025", "q975"), 1, ]
  rowMeans(apply(I, 2, function(i) unname(c(true_beta[1] >= i[1] & true_beta[1] <= i[2], i[2]-i[1]))))
}))
colnames(cov_int_quant_TL) = c("cov_q_TL", "len_q_TL")
cov_int_quant_sd_init = t(sapply(results_sd_init, function(R) {
  I = R[c("q025", "q975"), 1, ]
  rowMeans(apply(I, 2, function(i) unname(c(true_beta[1] >= i[1] & true_beta[1] <= i[2], i[2]-i[1]))))
}))
colnames(cov_int_quant_sd_init) = c("cov_q_sd", "len_q_sd")
cov_int_quant_TL_init = t(sapply(results_TL_init, function(R) {
  I = R[c("q025", "q975"), 1, ]
  rowMeans(apply(I, 2, function(i) unname(c(true_beta[1] >= i[1] & true_beta[1] <= i[2], i[2]-i[1]))))
}))
colnames(cov_int_quant_TL_init) = c("cov_q_TL", "len_q_TL")


cov_int_HDI_sd = t(sapply(results_sd, function(R) {
  I1 = R[c("HDI_L1", "HDI_U1"), 1, ]; I2 = R[c("HDI_L2", "HDI_U2"), 1, ]
  tmp1 = apply(I1, 2, function(i) true_beta[1] >= i[1] & true_beta[1] <= i[2])
  tmp2 = apply(I2, 2, function(i) true_beta[1] >= i[1] & true_beta[1] <= i[2])
  tmp2[is.na(tmp2)] = FALSE
  len1 = diff(I1); len2 = ifelse(is.na(diff(I2)), 0, diff(I2))
  c(mean(tmp1 | tmp2), mean(len1+len2))
}))
colnames(cov_int_HDI_sd) = c("cov_HDI_sd", "len_HDI_sd")
cov_int_HDI_TL = t(sapply(results_TL, function(R) {
  I1 = R[c("HDI_L1", "HDI_U1"), 1, ]; I2 = R[c("HDI_L2", "HDI_U2"), 1, ]
  tmp1 = apply(I1, 2, function(i) true_beta[1] >= i[1] & true_beta[1] <= i[2])
  tmp2 = apply(I2, 2, function(i) true_beta[1] >= i[1] & true_beta[1] <= i[2])
  tmp2[is.na(tmp2)] = FALSE
  len1 = diff(I1); len2 = ifelse(is.na(diff(I2)), 0, diff(I2))
  c(mean(tmp1 | tmp2), mean(len1+len2))
}))
colnames(cov_int_HDI_TL) = c("cov_HDI_TL", "len_HDI_TL")
cov_int_HDI_sd_init = t(sapply(results_sd_init, function(R) {
  I1 = R[c("HDI_L1", "HDI_U1"), 1, ]; I2 = R[c("HDI_L2", "HDI_U2"), 1, ]
  tmp1 = apply(I1, 2, function(i) true_beta[1] >= i[1] & true_beta[1] <= i[2])
  tmp2 = apply(I2, 2, function(i) true_beta[1] >= i[1] & true_beta[1] <= i[2])
  tmp2[is.na(tmp2)] = FALSE
  len1 = diff(I1); len2 = ifelse(is.na(diff(I2)), 0, diff(I2))
  c(mean(tmp1 | tmp2), mean(len1+len2))
}))
colnames(cov_int_HDI_sd_init) = c("cov_HDI_sd", "len_HDI_sd")
cov_int_HDI_TL_init = t(sapply(results_TL_init, function(R) {
  I1 = R[c("HDI_L1", "HDI_U1"), 1, ]; I2 = R[c("HDI_L2", "HDI_U2"), 1, ]
  tmp1 = apply(I1, 2, function(i) true_beta[1] >= i[1] & true_beta[1] <= i[2])
  tmp2 = apply(I2, 2, function(i) true_beta[1] >= i[1] & true_beta[1] <= i[2])
  tmp2[is.na(tmp2)] = FALSE
  len1 = diff(I1); len2 = ifelse(is.na(diff(I2)), 0, diff(I2))
  c(mean(tmp1 | tmp2), mean(len1+len2))
}))
colnames(cov_int_HDI_TL_init) = c("cov_HDI_TL", "len_HDI_TL")


df_beta0 = data.frame(sim_settings, 
                      cov_int_quant_sd, cov_int_quant_TL,
                      cov_int_HDI_sd, cov_int_HDI_TL)
df_beta0_init = data.frame(sim_settings, 
                           cov_int_quant_sd_init, cov_int_quant_TL_init,
                           cov_int_HDI_sd_init, cov_int_HDI_TL_init)

tab = knitr::kable(df_beta0[order_runs , c(1:3, 4, 6, 5, 7, 8, 10, 9, 11)], format = "pipe", digits = 3)
tab_init = knitr::kable(df_beta0_init[order_runs , c(1:3, 4, 6, 5, 7, 8, 10, 9, 11)], format = "pipe", digits = 3)
writeLines(c("k final", tab, "\n\n k init", tab_init), paste0(PATH_IMG, "multiv_beta0.txt"))

# Bias ----
bias_sd = lapply(results_sd, function(R) R["mean", , ] - true_beta) %>% abind::abind(along = 3) %>% aperm(c(3, 1, 2))
bias_TL = lapply(results_TL, function(R) R["mean", , ] - true_beta) %>% abind::abind(along = 3) %>% aperm(c(3, 1, 2))
bias_sd_init = lapply(results_sd_init, function(R) R["mean", , ] - true_beta) %>% abind::abind(along = 3) %>% aperm(c(3, 1, 2))
bias_TL_init = lapply(results_TL_init, function(R) R["mean", , ] - true_beta) %>% abind::abind(along = 3) %>% aperm(c(3, 1, 2))
dimnames(bias_sd) = dimnames(bias_TL) = 
  dimnames(bias_sd_init) = dimnames(bias_TL_init) = 
  list("sim_setting" = 1:12,
       "param" = c("beta0", "beta1", "beta2"),
       "replicate" = 1:100)

df_bias = bind_rows(bias_sd %>% reshape2::melt() %>% mutate(type = "sd", k = "final"),
                    bias_TL %>% reshape2::melt() %>% mutate(type = "TL", k = "final"),
                    bias_sd_init %>% reshape2::melt() %>% mutate(type = "sd", k = "init"),
                    bias_TL_init %>% reshape2::melt() %>% mutate(type = "TL", k = "init")) %>% 
  left_join(cbind.data.frame(sim_settings, sim_setting = 1:12), by = "sim_setting") %>% 
  mutate(err_distr = factor(err_distr, levels = error_distr),
         cov_distr = factor(cov_distr, levels = covar_distr),
         k = factor(k, levels = c("final", "init")))

plot_bias_beta1 = df_bias %>% filter(param == "beta1") %>% 
  ggplot(aes(x = k, y = value, color = type)) + 
  geom_hline(aes(yintercept = 0)) +
  geom_boxplot() +
  facet_grid(err_distr ~ cov_distr, scales = "free_y") +
  scale_color_manual("Type Asymm.", values = c("sd" = col_sd, "TL" = col_TL)) +
  labs(title = "Bias beta1 across replicates (rows = err_distr, cols = cov_distr)") +
  theme(axis.title.y = element_blank())

plot_bias_beta2 = df_bias %>% filter(param == "beta2") %>%
  ggplot(aes(x = k, y = value, color = type)) + 
  geom_hline(aes(yintercept = 0)) +
  geom_boxplot() +
  facet_grid(err_distr ~ cov_distr, scales = "free_y") +
  scale_color_manual("Type Asymm.", values = c("sd" = col_sd, "TL" = col_TL)) +
  labs(title = "Bias beta2 across replicates (rows = err_distr, cols = cov_distr)") +
  theme(axis.title.y = element_blank())

plot_bias_beta0 = df_bias %>% filter(param == "beta0") %>%
  ggplot(aes(x = k, y = value, color = type)) + 
  geom_hline(aes(yintercept = 0)) +
  geom_boxplot() +
  facet_grid(err_distr ~ cov_distr, scales = "free_y") +
  scale_color_manual("Type Asymm.", values = c("sd" = col_sd, "TL" = col_TL)) +
  labs(title = "Bias beta0 across replicates (rows = err_distr, cols = cov_distr)") +
  theme(axis.title.y = element_blank())

ggsave(plot = gridExtra::marrangeGrob(list(plot_bias_beta0, plot_bias_beta1, plot_bias_beta2), 
                                      nrow = 1, ncol = 1, top = NULL),
       filename = "bias.pdf", path = PATH_IMG,
       width = 9, height = 6)

# Analysis w TL ----
w_TL = t(sapply(input_TL, function(I) sapply(I, function(x) x$w)))
w_sd = t(sapply(input_sd, function(I) sapply(I, function(x) x$w)))
w_TL_init = t(sapply(input_TL_init, function(I) sapply(I, function(x) x$w)))
w_sd_init = t(sapply(input_sd_init, function(I) sapply(I, function(x) x$w)))
colnames(w_TL) = colnames(w_sd) = 
  colnames(w_TL_init) = colnames(w_sd_init) =
  paste0("replicate", 1:100)

plot_w = bind_rows(cbind.data.frame(sim_settings, w_TL) %>% mutate(type = "sd", k = "final") %>% 
                     pivot_longer(replicate1:replicate100, names_to = "replicate", values_to = "w"),
                   cbind.data.frame(sim_settings, w_sd) %>% mutate(type = "TL", k = "final") %>%
                     pivot_longer(replicate1:replicate100, names_to = "replicate", values_to = "w"),
                   cbind.data.frame(sim_settings, w_TL_init) %>% mutate(type = "sd", k = "init") %>%
                     pivot_longer(replicate1:replicate100, names_to = "replicate", values_to = "w"),
                   cbind.data.frame(sim_settings, w_sd_init) %>% mutate(type = "TL", k = "init") %>%
                     pivot_longer(replicate1:replicate100, names_to = "replicate", values_to = "w")) %>% 
  mutate(err_distr = factor(err_distr, levels = error_distr),
         cov_distr = factor(cov_distr, levels = covar_distr),
         k = factor(k, levels = c("final", "init"))) %>% 
  ggplot(aes(x = k, y = w, color = type)) + 
  geom_boxplot() +
  facet_grid(err_distr ~ cov_distr, scales = "free_y") +
  scale_color_manual("Type Asymm.", values = c("sd" = col_sd, "TL" = col_TL)) +
  labs(title = "w values across replicates (rows = err_distr, cols = cov_distr)") +
  theme(axis.title = element_blank())

ggsave(plot = plot_w, 
       filename = "w.pdf", path = PATH_IMG,
       width = 9, height = 9)

# Average density TL ----
draws_TL = lapply(input_TL, 
                  function(I) 
                    lapply(I, function(x) x$draws[seq(10, 10^4, by = 10), ]) %>% abind::abind(along = 3)) %>% 
  abind::abind(along = 4) %>% aperm(c(4, 3, 2, 1))
draws_sd = lapply(input_sd, 
                  function(I) 
                    lapply(I, function(x) x$draws[seq(10, 10^4, by = 10), ]) %>% abind::abind(along = 3)) %>% 
  abind::abind(along = 4) %>% aperm(c(4, 3, 2, 1))
draws_TL_init = lapply(input_TL_init, 
                       function(I) 
                         lapply(I, function(x) x$draws[seq(10, 10^4, by = 10), ]) %>% abind::abind(along = 3)) %>% 
  abind::abind(along = 4) %>% aperm(c(4, 3, 2, 1))
draws_sd_init = lapply(input_sd_init, 
                       function(I) 
                         lapply(I, function(x) x$draws[seq(10, 10^4, by = 10), ]) %>% abind::abind(along = 3)) %>% 
  abind::abind(along = 4) %>% aperm(c(4, 3, 2, 1))


min_b0 = quantile(c(draws_TL[ , , 1, ], draws_sd[ , , 1, ], draws_TL_init[ , , 1, ], draws_sd_init[ , , 1, ]), 0.01)
max_b0 = quantile(c(draws_TL[ , , 1, ], draws_sd[ , , 1, ], draws_TL_init[ , , 1, ], draws_sd_init[ , , 1, ]), 0.99)
min_b1 = quantile(c(draws_TL[ , , 2, ], draws_sd[ , , 2, ], draws_TL_init[ , , 2, ], draws_sd_init[ , , 2, ]), 0.01)
max_b1 = quantile(c(draws_TL[ , , 2, ], draws_sd[ , , 2, ], draws_TL_init[ , , 2, ], draws_sd_init[ , , 2, ]), 0.99)
min_b2 = quantile(c(draws_TL[ , , 3, ], draws_sd[ , , 3, ], draws_TL_init[ , , 3, ], draws_sd_init[ , , 3, ]), 0.01)
max_b2 = quantile(c(draws_TL[ , , 3, ], draws_sd[ , , 3, ], draws_TL_init[ , , 3, ], draws_sd_init[ , , 3, ]), 0.99)

dens_b0_TL = apply(draws_TL[ , , 1, ], 1,
                   function(D1) apply(D1, 1,
                                      function(D2){
                                        tmp = density(D2, from = min_b0, to = max_b0, n = 100)
                                        cbind(tmp$x, tmp$y)
                                      }, simplify = F) %>% abind::abind(along = 3),
                   simplify = F) %>% abind::abind(along = 4) %>% 
  aperm(c(4, 3, 2, 1))
dens_b0_sd = apply(draws_sd[ , , 1, ], 1,
                   function(D1) apply(D1, 1,
                                      function(D2){
                                        tmp = density(D2, from = min_b0, to = max_b0, n = 100)
                                        cbind(tmp$x, tmp$y)
                                      }, simplify = F) %>% abind::abind(along = 3),
                   simplify = F) %>% abind::abind(along = 4) %>% 
  aperm(c(4, 3, 2, 1))
dens_b0_TL_init = apply(draws_TL_init[ , , 1, ], 1,
                        function(D1) apply(D1, 1,
                                           function(D2){
                                             tmp = density(D2, from = min_b0, to = max_b0, n = 100)
                                             cbind(tmp$x, tmp$y)
                                           }, simplify = F) %>% abind::abind(along = 3),
                        simplify = F) %>% abind::abind(along = 4) %>% 
  aperm(c(4, 3, 2, 1))
dens_b0_sd_init = apply(draws_sd_init[ , , 1, ], 1,
                        function(D1) apply(D1, 1,
                                           function(D2){
                                             tmp = density(D2, from = min_b0, to = max_b0, n = 100)
                                             cbind(tmp$x, tmp$y)
                                           }, simplify = F) %>% abind::abind(along = 3),
                        simplify = F) %>% abind::abind(along = 4) %>%
  aperm(c(4, 3, 2, 1))

dens_b1_TL = apply(draws_TL[ , , 2, ], 1,
                   function(D1) apply(D1, 1,
                                      function(D2){
                                        tmp = density(D2, from = min_b1, to = max_b1, n = 100)
                                        cbind(tmp$x, tmp$y)
                                      }, simplify = F) %>% abind::abind(along = 3),
                   simplify = F) %>% abind::abind(along = 4) %>% 
  aperm(c(4, 3, 2, 1))
dens_b1_sd = apply(draws_sd[ , , 2, ], 1,
                   function(D1) apply(D1, 1,
                                      function(D2){
                                        tmp = density(D2, from = min_b1, to = max_b1, n = 100)
                                        cbind(tmp$x, tmp$y)
                                      }, simplify = F) %>% abind::abind(along = 3),
                   simplify = F) %>% abind::abind(along = 4) %>% 
  aperm(c(4, 3, 2, 1))
dens_b1_TL_init = apply(draws_TL_init[ , , 2, ], 1,
                        function(D1) apply(D1, 1,
                                           function(D2){
                                             tmp = density(D2, from = min_b1, to = max_b1, n = 100)
                                             cbind(tmp$x, tmp$y)
                                           }, simplify = F) %>% abind::abind(along = 3),
                        simplify = F) %>% abind::abind(along = 4) %>% 
  aperm(c(4, 3, 2, 1))
dens_b1_sd_init = apply(draws_sd_init[ , , 2, ], 1,
                        function(D1) apply(D1, 1,
                                           function(D2){
                                             tmp = density(D2, from = min_b1, to = max_b1, n = 100)
                                             cbind(tmp$x, tmp$y)
                                           }, simplify = F) %>% abind::abind(along = 3),
                        simplify = F) %>% abind::abind(along = 4) %>% 
  aperm(c(4, 3, 2, 1))

dens_b2_TL = apply(draws_TL[ , , 3, ], 1,
                   function(D1) apply(D1, 1,
                                      function(D2){
                                        tmp = density(D2, from = min_b2, to = max_b2, n = 100)
                                        cbind(tmp$x, tmp$y)
                                      }, simplify = F) %>% abind::abind(along = 3),
                   simplify = F) %>% abind::abind(along = 4) %>% 
  aperm(c(4, 3, 2, 1))
dens_b2_sd = apply(draws_sd[ , , 3, ], 1,
                   function(D1) apply(D1, 1,
                                      function(D2){
                                        tmp = density(D2, from = min_b2, to = max_b2, n = 100)
                                        cbind(tmp$x, tmp$y)
                                      }, simplify = F) %>% abind::abind(along = 3),
                   simplify = F) %>% abind::abind(along = 4) %>% 
  aperm(c(4, 3, 2, 1))
dens_b2_TL_init = apply(draws_TL_init[ , , 3, ], 1,
                        function(D1) apply(D1, 1,
                                           function(D2){
                                             tmp = density(D2, from = min_b2, to = max_b2, n = 100)
                                             cbind(tmp$x, tmp$y)
                                           }, simplify = F) %>% abind::abind(along = 3),
                        simplify = F) %>% abind::abind(along = 4) %>%
  aperm(c(4, 3, 2, 1))
dens_b2_sd_init = apply(draws_sd_init[ , , 3, ], 1,
                        function(D1) apply(D1, 1,
                                           function(D2){
                                             tmp = density(D2, from = min_b2, to = max_b2, n = 100)
                                             cbind(tmp$x, tmp$y)
                                           }, simplify = F) %>% abind::abind(along = 3),
                        simplify = F) %>% abind::abind(along = 4) %>%
  aperm(c(4, 3, 2, 1))

dimnames(dens_b0_TL) = dimnames(dens_b1_TL) = dimnames(dens_b2_TL) = 
  dimnames(dens_b0_sd) = dimnames(dens_b1_sd) = dimnames(dens_b2_sd) =
  dimnames(dens_b0_TL_init) = dimnames(dens_b1_TL_init) = dimnames(dens_b2_TL_init) =
  dimnames(dens_b0_sd_init) = dimnames(dens_b1_sd_init) = dimnames(dens_b2_sd_init) =
  list("setting" = 1:12, "replicate" = 1:100, "axis" = c("x", "y"), "knot" = 1:100) 




# Plot density b0 ----
df_dens_b0 = 
  bind_rows(
    left_join(reshape2::melt(dens_b0_TL), cbind.data.frame(setting = 1:12, sim_settings), by = "setting") %>% 
      mutate(type = "TL", k = "final"),
    left_join(reshape2::melt(dens_b0_TL_init), cbind.data.frame(setting = 1:12, sim_settings), by = "setting") %>%
      mutate(type = "TL", k = "init"),
    left_join(reshape2::melt(dens_b0_sd), cbind.data.frame(setting = 1:12, sim_settings), by = "setting") %>% 
      mutate(type = "sd", k = "final"),
    left_join(reshape2::melt(dens_b0_sd_init), cbind.data.frame(setting = 1:12, sim_settings), by = "setting") %>%
      mutate(type = "sd", k = "init")
  ) %>% 
  pivot_wider(names_from = axis, values_from = value) %>%
  group_by(setting, knot, type, .drop = F) %>% 
  mutate(mean_y = mean(y)) %>% 
  mutate(err_distr = factor(err_distr, levels = error_distr),
         cov_distr = factor(cov_distr, levels = covar_distr))

plt_dens_b0_SNR2 = df_dens_b0 %>% filter(SNR == "SNR2") %>% 
  ggplot(aes(x = x)) +
  geom_vline(xintercept = true_beta[1], col = "gold") +
  geom_line(aes(y = y, group = interaction(replicate, type), col = type), alpha = 0.1) +
  scale_colour_manual("Type Asymm.", values = colorspace::lighten(c("sd" = col_sd, "TL" = col_TL), amount = 0.5)) +
  ggnewscale::new_scale_colour() +
  scale_color_manual("Type Asymm.", values = c("sd" = col_sd, "TL" = col_TL)) +
  geom_line(aes(y = mean_y, col = type)) + 
  ggh4x::facet_nested(cov_distr ~ err_distr + k) +
  labs(title = "SNR = 2 | Density beta0 [light = per replica; dark = average] (rows = err_distr, cols = cov_distr)", 
       y = "density", x = "beta2") +
  scale_x_continuous(limits = true_beta[1] + c(-.5, +.5)) +
  theme(axis.title.y = element_blank())

plt_dens_b0_SNR1 = df_dens_b0 %>% filter(SNR == "SNR1") %>% 
  ggplot(aes(x = x)) +
  geom_vline(xintercept = true_beta[1], col = "gold") +
  geom_line(aes(y = y, group = interaction(replicate, type), col = type), alpha = 0.1) +
  scale_colour_manual("Type Asymm.", values = colorspace::lighten(c("sd" = col_sd, "TL" = col_TL), amount = 0.5)) +
  ggnewscale::new_scale_colour() +
  scale_color_manual("Type Asymm.", values = c("sd" = col_sd, "TL" = col_TL)) +
  geom_line(aes(y = mean_y, col = type)) + 
  ggh4x::facet_nested(cov_distr ~ err_distr + k) +
  labs(title = "SNR = 1 | Density beta0 [light = per replica; dark = average] (rows = err_distr, cols = cov_distr)", 
       y = "density", x = "beta2") +
  scale_x_continuous(limits = true_beta[1] + c(-.5, +.5)) +
  theme(axis.title.y = element_blank())

ggsave(plot = gridExtra::marrangeGrob(list(plt_dens_b0_SNR2, plt_dens_b0_SNR1), nrow = 1, ncol = 1, top = NULL),
       filename = "density_b0.pdf", path = PATH_IMG,
       width = 9, height = 5)

# Plot density b1 ----
df_dens_b1 = 
  bind_rows(
    left_join(reshape2::melt(dens_b1_TL), cbind.data.frame(setting = 1:12, sim_settings), by = "setting") %>% 
      mutate(type = "TL", k = "final"),
    left_join(reshape2::melt(dens_b1_TL_init), cbind.data.frame(setting = 1:12, sim_settings), by = "setting") %>%
      mutate(type = "TL", k = "init"),
    left_join(reshape2::melt(dens_b1_sd), cbind.data.frame(setting = 1:12, sim_settings), by = "setting") %>% 
      mutate(type = "sd", k = "final"),
    left_join(reshape2::melt(dens_b1_sd_init), cbind.data.frame(setting = 1:12, sim_settings), by = "setting") %>%
      mutate(type = "sd", k = "init")
  ) %>% 
  pivot_wider(names_from = axis, values_from = value) %>%
  group_by(setting, knot, type, .drop = F) %>% 
  mutate(mean_y = mean(y)) %>% 
  mutate(err_distr = factor(err_distr, levels = error_distr),
         cov_distr = factor(cov_distr, levels = covar_distr))

plt_dens_b1_SNR2 = df_dens_b1 %>% filter(SNR == "SNR2") %>% 
  ggplot(aes(x = x)) +
  geom_vline(xintercept = true_beta[2], col = "gold") +
  geom_line(aes(y = y, group = interaction(replicate, type), col = type), alpha = 0.1) +
  scale_colour_manual("Type Asymm.", values = colorspace::lighten(c("sd" = col_sd, "TL" = col_TL), amount = 0.5)) +
  ggnewscale::new_scale_colour() +
  scale_color_manual("Type Asymm.", values = c("sd" = col_sd, "TL" = col_TL)) +
  geom_line(aes(y = mean_y, col = type)) + 
  ggh4x::facet_nested(cov_distr ~ err_distr + k) +
  labs(title = "SNR = 2 | Density beta1 [light = per replica; dark = average] (rows = err_distr, cols = cov_distr)", 
       y = "density", x = "beta2") +
  scale_x_continuous(limits = true_beta[2] + c(-.5, +.5)) +
  theme(axis.title.y = element_blank())

plt_dens_b1_SNR1 = df_dens_b1 %>% filter(SNR == "SNR1") %>% 
  ggplot(aes(x = x)) +
  geom_vline(xintercept = true_beta[2], col = "gold") +
  geom_line(aes(y = y, group = interaction(replicate, type), col = type), alpha = 0.1) +
  scale_colour_manual("Type Asymm.", values = colorspace::lighten(c("sd" = col_sd, "TL" = col_TL), amount = 0.5)) +
  ggnewscale::new_scale_colour() +
  scale_color_manual("Type Asymm.", values = c("sd" = col_sd, "TL" = col_TL)) +
  geom_line(aes(y = mean_y, col = type)) + 
  ggh4x::facet_nested(cov_distr ~ err_distr + k) +
  labs(title = "SNR = 1 | Density beta1 [light = per replica; dark = average] (rows = err_distr, cols = cov_distr)", 
       y = "density", x = "beta2") +
  scale_x_continuous(limits = true_beta[2] + c(-.5, +.5)) +
  theme(axis.title.y = element_blank())

ggsave(plot = gridExtra::marrangeGrob(list(plt_dens_b1_SNR2, plt_dens_b1_SNR1), nrow = 1, ncol = 1, top = NULL),
       filename = "density_b1.pdf", path = PATH_IMG,
       width = 9, height = 5)

# Plot density beta2 ----
df_dens_b2 = 
  bind_rows(
    left_join(reshape2::melt(dens_b2_TL), cbind.data.frame(setting = 1:12, sim_settings), by = "setting") %>% 
      mutate(type = "TL", k = "final"),
    left_join(reshape2::melt(dens_b2_TL_init), cbind.data.frame(setting = 1:12, sim_settings), by = "setting") %>%
      mutate(type = "TL", k = "init"),
    left_join(reshape2::melt(dens_b2_sd), cbind.data.frame(setting = 1:12, sim_settings), by = "setting") %>% 
      mutate(type = "sd", k = "final"),
    left_join(reshape2::melt(dens_b2_sd_init), cbind.data.frame(setting = 1:12, sim_settings), by = "setting") %>%
      mutate(type = "sd", k = "init")
  ) %>% 
  pivot_wider(names_from = axis, values_from = value) %>%
  group_by(setting, knot, type, .drop = F) %>% 
  mutate(mean_y = mean(y)) %>% 
  mutate(err_distr = factor(err_distr, levels = error_distr),
         cov_distr = factor(cov_distr, levels = covar_distr))


plt_dens_b2_SNR2 = df_dens_b2 %>% filter(SNR == "SNR2") %>%
  ggplot(aes(x = x)) +
  geom_vline(xintercept = true_beta[3], col = "gold") +
  geom_line(aes(y = y, group = interaction(replicate, type), col = type), alpha = 0.1) +
  scale_colour_manual("Type Asymm.", values = colorspace::lighten(c("sd" = col_sd, "TL" = col_TL), amount = 0.5)) +
  ggnewscale::new_scale_colour() +
  scale_color_manual("Type Asymm.", values = c("sd" = col_sd, "TL" = col_TL)) +
  geom_line(aes(y = mean_y, col = type)) + 
  ggh4x::facet_nested(cov_distr ~ err_distr + k) +
  labs(title = "SNR = 2 | Density beta2 [light = per replica; dark = average] (rows = err_distr, cols = cov_distr)", 
       y = "density", x = "beta2") +
  scale_x_continuous(limits = true_beta[3] + c(-.5, +.5)) +
  theme(axis.title.y = element_blank())

plt_dens_b2_SNR1 = df_dens_b2 %>% filter(SNR == "SNR1") %>%
  ggplot(aes(x = x)) +
  geom_vline(xintercept = true_beta[3], col = "gold") +
  geom_line(aes(y = y, group = interaction(replicate, type), col = type), alpha = 0.1) +
  scale_colour_manual("Type Asymm.", values = colorspace::lighten(c("sd" = col_sd, "TL" = col_TL), amount = 0.5)) +
  ggnewscale::new_scale_colour() +
  scale_color_manual("Type Asymm.", values = c("sd" = col_sd, "TL" = col_TL)) +
  geom_line(aes(y = mean_y, col = type)) + 
  ggh4x::facet_nested(cov_distr ~ err_distr + k) +
  labs(title = "SNR = 1 | Density beta2 [light = per replica; dark = average] (rows = err_distr, cols = cov_distr)", 
       y = "density", x = "beta2") +
  scale_x_continuous(limits = true_beta[3] + c(-.5, +.5)) +
  theme(axis.title.y = element_blank())

ggsave(plot = gridExtra::marrangeGrob(list(plt_dens_b2_SNR2, plt_dens_b2_SNR1), nrow = 1, ncol = 1, top = NULL),
       filename = "density_b2.pdf", path = PATH_IMG,
       width = 9, height = 5)


# # Analysis k_multi from TL ----
# k_TL = lapply(input_TL, function(I) sapply(I, function(x) x$k)) %>% abind::abind(along = 3)
# dimnames(k_TL) = list("side" = c("neg", "pos"), 
#                       "replicate" = 1:100, 
#                       "sim_setting" = 1:12)
# 
# df_k_TL = k_TL %>% reshape2::melt() %>% 
#   left_join(cbind.data.frame(sim_setting = 1:12, sim_settings),
#             by = "sim_setting") %>% 
#   mutate(err_distr = factor(err_distr, levels = error_distr),
#          cov_distr = factor(cov_distr, levels = covar_distr),
#          side = factor(side, levels = c("neg", "pos"), labels = c("Left (<0)", "Right (>0)")))
# 
# plt_k_TL = df_k_TL %>%
#   ggplot() + 
#   geom_boxplot(aes(x = side, y = value)) +
#   facet_grid(err_distr ~ cov_distr, scales = "free_y") +
#   labs(title = "k (k1 and k2) values with Tower Law (rows = err_distr, cols = cov_distr)") +
#   theme(axis.title = element_blank())
# 
# ggsave(plot = plt_k_TL, 
#        filename = "k_TL.pdf", path = PATH_IMG,
#        width = 6, height = 6)
# 
# 
# # Analysis k_multi from sd ----
# k_sd = lapply(input_sd, function(I) sapply(I, function(x) x$k)) %>% abind::abind(along = 3)
# dimnames(k_sd) = list("side" = c("neg", "pos"), 
#                       "replicate" = 1:100, 
#                       "sim_setting" = 1:12)
# 
# df_k_sd = k_sd %>% reshape2::melt() %>% 
#   left_join(cbind.data.frame(sim_setting = 1:12, sim_settings),
#             by = "sim_setting") %>% 
#   mutate(err_distr = factor(err_distr, levels = error_distr),
#          cov_distr = factor(cov_distr, levels = covar_distr),
#          side = factor(side, levels = c("neg", "pos"), labels = c("Left (<0)", "Right (>0)")))
# 
# plt_k_sd = df_k_sd %>%
#   ggplot() + 
#   geom_boxplot(aes(x = side, y = value)) +
#   facet_grid(err_distr ~ cov_distr, scales = "free_y") +
#   labs(title = "k (k1 and k2) values with SD (rows = err_distr, cols = cov_distr)") +
#   theme(axis.title = element_blank())
# 
# ggsave(plot = plt_k_sd,
#        filename = "k_TL.pdf", path = PATH_IMG,
#        width = 6, height = 6)
# 
# 
# ## Plot together ----
# df_k_complete = bind_rows(df_k_TL %>% mutate(type = "TL"), df_k_sd %>% mutate(type = "sd")) %>% 
#   mutate(type = factor(type, levels = c("sd", "TL"))) %>% 
#   mutate(err_distr = factor(err_distr, levels = error_distr),
#          cov_distr = factor(cov_distr, levels = covar_distr))
# 
# plt_k_together = df_k_complete %>% 
#   ggplot() + 
#   geom_boxplot(aes(x = side, y = value, col = type)) +
#   facet_grid(err_distr ~ cov_distr, scales = "free_y") +
#   scale_color_manual("Type Asymm.", values = c("sd" = col_sd, "TL" = col_TL)) +
#   labs(title = "k (k1 and k2) values across replicates (rows = err_distr, cols = cov_distr)") +
#   theme(axis.title = element_blank())
# 
# ggsave(plot = plt_k_together,
#        filename = "k_together.pdf", path = PATH_IMG,
#        width = 9, height = 6)


# Predictor ----
err_pred_sd = lapply(input_sd, function(I1) sapply(I1, function(I2) {
  b = I2$draws[idx_iter, ]; X = I2$X
  colMeans(b %*% t(X)) - true_beta %*% t(X)
})) %>% abind::abind(along = 3) %>% aperm(c(3, 2, 1))
dimnames(err_pred_sd) = list("sim_setting" = 1:12, "replicate" = 1:100, "unit" = 1:1000)
err_pred_TL = lapply(input_TL, function(I1) sapply(I1, function(I2) {
  b = I2$draws[idx_iter, ]; X_scaled = I2$X_scaled; X = I2$X
  colMeans(b %*% t(X)) - true_beta %*% t(X)
})) %>% abind::abind(along = 3) %>% aperm(c(3, 2, 1))
dimnames(err_pred_TL) = list("sim_setting" = 1:12, "replicate" = 1:100, "unit" = 1:1000)
err_pred_sd_init = lapply(input_sd_init, function(I1) sapply(I1, function(I2) {
  b = I2$draws[idx_iter, ]; X = I2$X
  colMeans(b %*% t(X)) - true_beta %*% t(X)
})) %>% abind::abind(along = 3) %>% aperm(c(3, 2, 1))
dimnames(err_pred_sd_init) = list("sim_setting" = 1:12, "replicate" = 1:100, "unit" = 1:1000)
err_pred_TL_init = lapply(input_TL_init, function(I1) sapply(I1, function(I2) {
  b = I2$draws[idx_iter, ]; X_scaled = I2$X_scaled; X = I2$X
  colMeans(b %*% t(X)) - true_beta %*% t(X)
})) %>% abind::abind(along = 3) %>% aperm(c(3, 2, 1))
dimnames(err_pred_TL_init) = list("sim_setting" = 1:12, "replicate" = 1:100, "unit" = 1:1000)

dens_err_pred_sd = apply(err_pred_sd, 1, function(E1) apply(E1, 1, function(E) {
  tmp = density(E, from = -1.5, to = 1.5, n = 25)
  cbind(tmp$x, tmp$y)
}, simplify = F) %>% abind::abind(along = 3), simplify = F) %>% 
  abind::abind(along = 4)  %>% aperm(c(4, 3, 2, 1))
dens_err_pred_TL = apply(err_pred_TL, 1, function(E1) apply(E1, 1, function(E) {
  tmp = density(E, from = -1.5, to = 1.5, n = 25)
  cbind(tmp$x, tmp$y)
}, simplify = F) %>% abind::abind(along = 3), simplify = F) %>% 
  abind::abind(along = 4)  %>% aperm(c(4, 3, 2, 1))
dens_err_pred_sd_init = apply(err_pred_sd_init, 1, function(E1) apply(E1, 1, function(E) {
  tmp = density(E, from = -1.5, to = 1.5, n = 25)
  cbind(tmp$x, tmp$y)
}, simplify = F) %>% abind::abind(along = 3), simplify = F) %>% 
  abind::abind(along = 4)  %>% aperm(c(4, 3, 2, 1))
dens_err_pred_TL_init = apply(err_pred_TL_init, 1, function(E1) apply(E1, 1, function(E) {
  tmp = density(E, from = -1.5, to = 1.5, n = 25)
  cbind(tmp$x, tmp$y)
}, simplify = F) %>% abind::abind(along = 3), simplify = F) %>% 
  abind::abind(along = 4)  %>% aperm(c(4, 3, 2, 1))

dimnames(dens_err_pred_sd) = dimnames(dens_err_pred_TL) =
  dimnames(dens_err_pred_sd_init) = dimnames(dens_err_pred_TL_init) =
  list("setting" = 1:12, "replicate" = 1:100, "axis" = c("x", "y"), "knot" = 1:25) 

df_dens_err_pred = 
  left_join(
    bind_rows(
      reshape2::melt(dens_err_pred_sd) %>% mutate(type = "sd", k = "final"),
      reshape2::melt(dens_err_pred_TL) %>% mutate(type = "TL", k = "final"),
      reshape2::melt(dens_err_pred_sd_init) %>% mutate(type = "sd", k = "init"),
      reshape2::melt(dens_err_pred_TL_init) %>% mutate(type = "TL", k = "init")),
    cbind.data.frame(setting = 1:12, sim_settings), by = "setting") %>% 
  pivot_wider(names_from = axis, values_from = value) %>%
  group_by(setting, err_distr, cov_distr, knot, type, x, .drop = F) %>% 
  mutate(mean_y = mean(y)) %>% 
  mutate(err_distr = factor(err_distr, levels = error_distr),
         cov_distr = factor(cov_distr, levels = covar_distr))


plt_dens_err_pred_SNR2 = df_dens_err_pred %>% filter(SNR == "SNR2") %>% 
  ggplot(aes(x = x)) +
  geom_vline(aes(xintercept = 0), col = "gold") +
  geom_line(aes(y = y, group = interaction(replicate, type), col = type), alpha = 0.1) +
  scale_colour_manual("Type Asymm.", values = colorspace::lighten(c("sd" = col_sd, "TL" = col_TL), amount = 0.5)) +
  ggnewscale::new_scale_color() +
  geom_line(aes(y = mean_y, col = type)) + 
  scale_colour_manual("Type Asymm.", values = c("sd" = col_sd, "TL" = col_TL)) +
  scale_y_continuous(limits = c(0, 5)) +
  ggh4x::facet_nested(cov_distr ~ err_distr + k) +
  labs(title = "SNR = 2 | Error on linear predictor [light = per replicate, dark = average] (rows = err_distr, cols = cov_distr)", 
       y = "density", x = "error linear predictor") +
  theme(axis.title.x = element_blank())

plt_dens_err_pred_SNR1 = df_dens_err_pred %>% filter(SNR == "SNR1") %>%
  ggplot(aes(x = x)) +
  geom_vline(aes(xintercept = 0), col = "gold") +
  geom_line(aes(y = y, group = interaction(replicate, type), col = type), alpha = 0.1) +
  scale_colour_manual("Type Asymm.", values = colorspace::lighten(c("sd" = col_sd, "TL" = col_TL), amount = 0.5)) +
  ggnewscale::new_scale_color() +
  geom_line(aes(y = mean_y, col = type)) + 
  scale_colour_manual("Type Asymm.", values = c("sd" = col_sd, "TL" = col_TL)) +
  scale_y_continuous(limits = c(0, 5)) +
  ggh4x::facet_nested(cov_distr ~ err_distr + k) +
  labs(title = "SNR = 1 | Error on linear predictor [light = per replicate, dark = average] (rows = err_distr, cols = cov_distr)", 
       y = "density", x = "error linear predictor") +
  theme(axis.title.x = element_blank())

ggsave(plot = gridExtra::marrangeGrob(list(plt_dens_err_pred_SNR2, plt_dens_err_pred_SNR1), nrow = 1, ncol = 1, top = NULL),
       filename = "err_lin_pred.pdf", path = PATH_IMG,
       width = 9, height = 6)


# Residuals ----
residuals_sd = lapply(input_sd, function(I1) sapply(I1, function(I2) {
  b = I2$draws[idx_iter, ]; X = I2$X; y = I2$y
  y - colMeans(b %*% t(X))
})) %>% abind::abind(along = 3) %>% aperm(c(3, 2, 1))
dimnames(err_pred_sd) = list("sim_setting" = 1:12, "replicate" = 1:100, "unit" = 1:1000)
residuals_TL = lapply(input_TL, function(I1) sapply(I1, function(I2) {
  b = I2$draws[idx_iter, ]; X = I2$X; y = I2$y
  y - colMeans(b %*% t(X))
})) %>% abind::abind(along = 3) %>% aperm(c(3, 2, 1))
dimnames(err_pred_TL) = list("sim_setting" = 1:12, "replicate" = 1:100, "unit" = 1:1000)
residuals_sd_init = lapply(input_sd_init, function(I1) sapply(I1, function(I2) {
  b = I2$draws[idx_iter, ]; X = I2$X; y = I2$y
  y - colMeans(b %*% t(X))
})) %>% abind::abind(along = 3) %>% aperm(c(3, 2, 1))
dimnames(err_pred_sd_init) = list("sim_setting" = 1:12, "replicate" = 1:100, "unit" = 1:1000)
residuals_TL_init = lapply(input_TL_init, function(I1) sapply(I1, function(I2) {
  b = I2$draws[idx_iter, ]; X = I2$X; y = I2$y
  y - colMeans(b %*% t(X))
})) %>% abind::abind(along = 3) %>% aperm(c(3, 2, 1))
dimnames(err_pred_TL_init) = list("sim_setting" = 1:12, "replicate" = 1:100, "unit" = 1:1000)

dens_residuals_sd = apply(residuals_sd, 1, function(R1) apply(R1, 1, function(R) {
  tmp = density(R, from = -3, to = 3, n = 51)
  cbind(tmp$x, tmp$y)
}, simplify = F) %>% abind::abind(along = 3), simplify = F) %>% 
  abind::abind(along = 4)  %>% aperm(c(4, 3, 2, 1))
dens_residuals_TL = apply(residuals_TL, 1, function(R1) apply(R1, 1, function(R) {
  tmp = density(R, from = -3, to = 3, n = 51)
  cbind(tmp$x, tmp$y)
}, simplify = F) %>% abind::abind(along = 3), simplify = F) %>% 
  abind::abind(along = 4)  %>% aperm(c(4, 3, 2, 1))
dens_residuals_sd_init = apply(residuals_sd_init, 1, function(R1) apply(R1, 1, function(R) {
  tmp = density(R, from = -3, to = 3, n = 51)
  cbind(tmp$x, tmp$y)
}, simplify = F) %>% abind::abind(along = 3), simplify = F) %>% 
  abind::abind(along = 4)  %>% aperm(c(4, 3, 2, 1))
dens_residuals_TL_init = apply(residuals_TL_init, 1, function(R1) apply(R1, 1, function(R) {
  tmp = density(R, from = -3, to = 3, n = 51)
  cbind(tmp$x, tmp$y)
}, simplify = F) %>% abind::abind(along = 3), simplify = F) %>% 
  abind::abind(along = 4)  %>% aperm(c(4, 3, 2, 1))

dimnames(dens_residuals_sd) = dimnames(dens_residuals_TL) =
  dimnames(dens_residuals_sd_init) = dimnames(dens_residuals_TL_init) =
  list("setting" = 1:12, "replicate" = 1:100, "axis" = c("x", "y"), "knot" = 1:51) 

df_residuals = 
  left_join(
    bind_rows(
      reshape2::melt(dens_residuals_sd) %>% mutate(type = "sd", k = "final"),
      reshape2::melt(dens_residuals_TL) %>% mutate(type = "TL", k = "final"),
      reshape2::melt(dens_residuals_sd_init) %>% mutate(type = "sd", k = "init"),
      reshape2::melt(dens_residuals_TL_init) %>% mutate(type = "TL", k = "init")),
    cbind.data.frame(setting = 1:12, sim_settings), by = "setting"
  ) %>% 
  pivot_wider(names_from = axis, values_from = value) %>%
  group_by(setting, knot, type, .drop = F) %>% 
  mutate(mean_y = mean(y)) %>% 
  mutate(err_distr = factor(err_distr, levels = error_distr),
         cov_distr = factor(cov_distr, levels = covar_distr))

df_true_density_final = data.frame(x = seq(-3, 3, length = 101),
                             y = c(dnorm(seq(-3, 3, length = 101), sd = 1),
                                   dgamma(seq(-3, 3, length = 101) + 1/sqrt(2), 2, rate = sqrt(2)),
                                   dnorm(seq(-3, 3, length = 101), sd = sqrt(2)),
                                   dgamma(seq(-3, 3, length = 101) + 1, 2, rate = 1)),
                             err_distr = rep(error_distr, each = 101),
                             SNR = rep(c("SNR2", "SNR2", "SNR1", "SNR1"), each = 101)) %>% 
  mutate(err_distr = factor(err_distr, levels = error_distr),
         k = factor("final", levels = c("final", "init")))
df_true_density_init = data.frame(x = seq(-3, 3, length = 101),
                             y = c(dnorm(seq(-3, 3, length = 101), sd = 1),
                                   dgamma(seq(-3, 3, length = 101) + 1/sqrt(2), 2, rate = sqrt(2)),
                                   dnorm(seq(-3, 3, length = 101), sd = sqrt(2)),
                                   dgamma(seq(-3, 3, length = 101) + 1, 2, rate = 1)),
                             err_distr = rep(error_distr, each = 101),
                             SNR = rep(c("SNR2", "SNR2", "SNR1", "SNR1"), each = 101)) %>% 
  mutate(err_distr = factor(err_distr, levels = error_distr),
         k = factor("init", levels = c("final", "init")))
df_true_density = bind_rows(df_true_density_final, df_true_density_init)

plt_dens_residuals_SNR2 = df_residuals %>% filter(SNR == "SNR2") %>%
  ggplot(aes(x = x)) +
  geom_line(aes(y = y, group = interaction(replicate, type), col = type), alpha = 0.1) +
  scale_colour_manual("Type Asymm.", values = colorspace::lighten(c("sd" = col_sd, "TL" = col_TL), amount = 0.5)) +
  ggnewscale::new_scale_color() +
  geom_line(aes(y = mean_y, col = type)) + 
  scale_color_manual("Type Asymm.", values = c("sd" = col_sd, "TL" = col_TL)) +
  geom_line(data = df_true_density %>% filter(SNR == "SNR2"),
            aes(x = x, y = y),
            col = "gold", lty = 2) +
  ggh4x::facet_nested(cov_distr ~ err_distr + k, drop = T) +
  labs(title = "SNR = 2 | Residuals [light = per replicate, dark = average] (rows = err_distr, cols = cov_distr)", 
       y = "density", x = "residual") + 
  theme(axis.title.x = element_blank())

plt_dens_residuals_SNR1 = df_residuals %>% filter(SNR == "SNR1") %>%
  ggplot(aes(x = x)) +
  geom_line(aes(y = y, group = interaction(replicate, type), col = type), alpha = 0.1) +
  scale_colour_manual("Type Asymm.", values = colorspace::lighten(c("sd" = col_sd, "TL" = col_TL), amount = 0.5)) +
  ggnewscale::new_scale_color() +
  geom_line(aes(y = mean_y, col = type)) + 
  scale_color_manual("Type Asymm.", values = c("sd" = col_sd, "TL" = col_TL)) +
  geom_line(data = df_true_density %>% filter(SNR == "SNR1"),
            aes(x = x, y = y),
            col = "gold", lty = 2) +
  ggh4x::facet_nested(cov_distr ~ err_distr + k) +
  labs(title = "SNR = 1 | Residuals [light = per replicate, dark = average] (rows = err_distr, cols = cov_distr)", 
       y = "density", x = "residual") + 
  theme(axis.title.x = element_blank())

ggsave(plot = gridExtra::marrangeGrob(list(plt_dens_residuals_SNR2, plt_dens_residuals_SNR1), nrow = 1, ncol = 1, top = NULL),
       filename = "residuals_dens.pdf", path = PATH_IMG,
       width = 9, height = 6)


# Bimodality ----
LB_bimod_TL = cbind(sim_settings, 
                    sapply(results_TL, function(R) apply(R["HDI_L2", , ], 1, function(x) 1-mean(is.na(x)))) %>% t())
LB_bimod_sd = cbind(sim_settings,
                    sapply(results_sd, function(R) apply(R["HDI_L2", , ], 1, function(x) 1-mean(is.na(x)))) %>% t())
colnames(LB_bimod_TL)[4:6] = colnames(LB_bimod_sd)[4:6] = c("beta0", "beta1", "beta2")
LB_bimod_TL
LB_bimod_sd

# Diagnostics ----
diagn_TL = lapply(input_TL, function(I1) 
  lapply(I1, function(I2) I2$diagn) %>% abind::abind(along = 3)) %>% 
  abind::abind(along = 4)
diagn_sd = lapply(input_sd, function(I1)
  lapply(I1, function(I2) I2$diagn) %>% abind::abind(along = 3)) %>% 
  abind::abind(along = 4)

dimnames(diagn_sd) = 
  dimnames(diagn_TL) =
  list("index" = dimnames(diagn_sd)[[1]], 
       "param" = dimnames(diagn_sd)[[2]],
       "replicate" = 1:100,
       "sim_setting" = 1:12)

df_diagn_sd = diagn_sd %>% reshape2::melt() %>% mutate(type = "sd") %>% 
  left_join(cbind.data.frame(sim_setting = 1:12, sim_settings),
            by = "sim_setting") %>% 
  separate(index, into = c("index", "thin")) %>% 
  pivot_wider(names_from = index, values_from = value) %>% 
  mutate(err_distr = factor(err_distr, levels = error_distr),
         cov_distr = factor(cov_distr, levels = covar_distr))

df_diagn_TL = diagn_TL %>% reshape2::melt() %>% mutate(type = "TL") %>% 
  left_join(cbind.data.frame(sim_setting = 1:12, sim_settings),
            by = "sim_setting") %>% 
  separate(index, into = c("index", "thin")) %>% 
  pivot_wider(names_from = index, values_from = value) %>% 
  mutate(err_distr = factor(err_distr, levels = error_distr),
         cov_distr = factor(cov_distr, levels = covar_distr))

df_diagn = bind_rows(df_diagn_sd, df_diagn_TL)

df_failed = df_diagn %>% 
  group_by(type, sim_setting, err_distr, cov_distr, thin) %>% 
  summarize(failed_101 = mean(rhat > 1.01),
            failed_11 = mean(rhat > 1.1)) 

df_diagn_sd %>% 
  ggplot() +
  geom_boxplot(aes(x = err_distr, y = rhat)) +
  geom_hline(aes(yintercept = 1.01), col = "red") +
  facet_grid(k ~ unif_cov, scales = "free") +
  scale_y_continuous(transform = "log") +
  labs(title = "Rhat in asymmetric sd")

df_diagn_TL %>% 
  ggplot() +
  geom_boxplot(aes(x = err_distr, y = rhat)) +
  geom_hline(aes(yintercept = 1.01), col = "red") +
  facet_grid(k ~ unif_cov, scales = "free") +
  scale_y_continuous(transform = "log")

df_failed %>% 
  ggplot() +
  geom_point(aes(x = unif_cov, y = failed_101, col = type, shape = thin)) +
  facet_grid(k ~ err_distr, scales = "free")

df_failed %>% 
  ggplot() +
  geom_point(aes(x = unif_cov, y = failed_11, col = type, shape = thin)) +
  facet_grid(k ~ err_distr, scales = "free")
