rm(list =)
results_symm = readRDS("sim_study_nonC-GBI/sim_study_nonCalibrated.RDS")
results_asymm = readRDS("sim_study_nonC-GBI/sim_study_nonCalibrated_asymm.RDS")

# TRUE BETA AND SAMPLE SIZE ----
n = 1e3
true_beta = c(0, -0.5)

# SIMULATION SETTINGS ----
error_distr = c("Gaussian", "Gaussian_0.5", "Gamma_2_2")
covariate = c("0.5", "1", "asymm")
k = c("0.5", "init", "final")
sim_settings = as.matrix(expand.grid(error_distr, covariate, k))
colnames(sim_settings) = c("err_distr", "unif_cov", "k")

# Plot point estimates ----
# First regression coefficient
pdf("sim_study_nonC-GBI/univ_beta1.pdf", width = 9, height = 6)
par(mfrow = c(3, 3), mar = c(3, 2, 2, 1))
for (j in seq(1, 27, by = 3)){
  min_x = -1 #min(c(results[[j]]["mean", 2, ], results[[j]]["BFGS", 2, ])) 
  max_x = 0 #max(c(results[[j]]["mean", 2, ], results[[j]]["BFGS", 2, ])) 
  plot(1, xlim = c(min_x, max_x), ylim = c(0, 12), type = "n", main = paste(sim_settings[j, ], collapse = " - "),
       xlab = "", ylab = "")
  lines(density(results_symm[[j]]["mean", 2, ]), col = "blue")
  lines(density(results_symm[[j]]["BFGS", 2, ]), col = "blue", lty = 2)
  lines(density(results_asymm[[j]]["mean", 2, ]), col = "red")
  lines(density(results_asymm[[j]]["BFGS", 2, ]), col = "red", lty = 2)
  abline(v = true_beta[2], col = "gold")
  legend("topright", col = c(4, 4, "gold"), lty = c(1, 2, 1),
         legend = c("Bayes Symm", "Freq Symm", "True"))
  legend("topleft", col = c(2, 2), lty = c(1, 2),
         legend = c("Bayes Asymm", "Freq Asymm"))
}

for (j in seq(2, 27, by = 3)){
  min_x = -1 #min(c(results[[j]]["mean", 2, ], results[[j]]["BFGS", 2, ])) 
  max_x = 0 #max(c(results[[j]]["mean", 2, ], results[[j]]["BFGS", 2, ])) 
  plot(1, xlim = c(min_x, max_x), ylim = c(0, 12), type = "n", main = paste(sim_settings[j, ], collapse = " - "),
       xlab = "", ylab = "")
  lines(density(results_symm[[j]]["mean", 2, ]), col = "blue")
  lines(density(results_symm[[j]]["BFGS", 2, ]), col = "blue", lty = 2)
  lines(density(results_asymm[[j]]["mean", 2, ]), col = "red")
  lines(density(results_asymm[[j]]["BFGS", 2, ]), col = "red", lty = 2)
  abline(v = true_beta[2], col = "gold")
  legend("topright", col = c(4, 4, "gold"), lty = c(1, 2, 1),
         legend = c("Bayes Symm", "Freq Symm", "True"))
  legend("topleft", col = c(2, 2), lty = c(1, 2),
         legend = c("Bayes Asymm", "Freq Asymm"))
}

for (j in seq(3, 27, by = 3)){
  min_x = -1 #min(c(results[[j]]["mean", 2, ], results[[j]]["BFGS", 2, ])) 
  max_x = 0 #max(c(results[[j]]["mean", 2, ], results[[j]]["BFGS", 2, ])) 
  plot(1, xlim = c(min_x, max_x), ylim = c(0, 12), type = "n", main = paste(sim_settings[j, ], collapse = " - "),
       xlab = "", ylab = "")
  lines(density(results_symm[[j]]["mean", 2, ]), col = "blue")
  lines(density(results_symm[[j]]["BFGS", 2, ]), col = "blue", lty = 2)
  lines(density(results_asymm[[j]]["mean", 2, ]), col = "red")
  lines(density(results_asymm[[j]]["BFGS", 2, ]), col = "red", lty = 2)
  abline(v = true_beta[2], col = "gold")
  legend("topright", col = c(4, 4, "gold"), lty = c(1, 2, 1),
         legend = c("Bayes Symm", "Freq Symm", "True"))
  legend("topleft", col = c(2, 2), lty = c(1, 2),
         legend = c("Bayes Asymm", "Freq Asymm"))
}
dev.off()


# Intercept
pdf("sim_study_nonC-GBI/univ_beta0.pdf", width = 9, height = 6)
par(mfrow = c(3, 3), mar = c(3, 2, 2, 1))
for (j in seq(1, 27, by = 3)){
  min_x = -0.25 #min(c(results[[j]]["mean", 2, ], results[[j]]["BFGS", 2, ])) 
  max_x = +0.25 #max(c(results[[j]]["mean", 2, ], results[[j]]["BFGS", 2, ])) 
  plot(1, xlim = c(min_x, max_x), ylim = c(0, 12), type = "n", main = paste(sim_settings[j, ], collapse = " - "),
       xlab = "", ylab = "")
  lines(density(results_symm[[j]]["mean", 1, ]), col = "blue")
  lines(density(results_symm[[j]]["BFGS", 1, ]), col = "blue", lty = 2)
  lines(density(results_asymm[[j]]["mean", 1, ]), col = "red")
  lines(density(results_asymm[[j]]["BFGS", 1, ]), col = "red", lty = 2)
  abline(v = true_beta[1], col = "gold")
  legend("topright", col = c(4, 4, "gold"), lty = c(1, 2, 1),
         legend = c("Bayes Symm", "Freq Symm", "True"))
  legend("topleft", col = c(2, 2), lty = c(1, 2),
         legend = c("Bayes Asymm", "Freq Asymm"))
}

for (j in seq(2, 27, by = 3)){
  min_x = -0.25 #min(c(results[[j]]["mean", 2, ], results[[j]]["BFGS", 2, ])) 
  max_x = +0.25 #max(c(results[[j]]["mean", 2, ], results[[j]]["BFGS", 2, ])) 
  plot(1, xlim = c(min_x, max_x), ylim = c(0, 12), type = "n", main = paste(sim_settings[j, ], collapse = " - "),
       xlab = "", ylab = "")
  lines(density(results_symm[[j]]["mean", 1, ]), col = "blue")
  lines(density(results_symm[[j]]["BFGS", 1, ]), col = "blue", lty = 2)
  lines(density(results_asymm[[j]]["mean", 1, ]), col = "red")
  lines(density(results_asymm[[j]]["BFGS", 1, ]), col = "red", lty = 2)
  abline(v = true_beta[1], col = "gold")
  legend("topright", col = c(4, 4, "gold"), lty = c(1, 2, 1),
         legend = c("Bayes Symm", "Freq Symm", "True"))
  legend("topleft", col = c(2, 2), lty = c(1, 2),
         legend = c("Bayes Asymm", "Freq Asymm"))
}

for (j in seq(3, 27, by = 3)){
  min_x = -0.25 #min(c(results[[j]]["mean", 2, ], results[[j]]["BFGS", 2, ])) 
  max_x = +0.25 #max(c(results[[j]]["mean", 2, ], results[[j]]["BFGS", 2, ])) 
  plot(1, xlim = c(min_x, max_x), ylim = c(0, 12), type = "n", main = paste(sim_settings[j, ], collapse = " - "),
       xlab = "", ylab = "")
  lines(density(results_symm[[j]]["mean", 1, ]), col = "blue")
  lines(density(results_symm[[j]]["BFGS", 1, ]), col = "blue", lty = 2)
  lines(density(results_asymm[[j]]["mean", 1, ]), col = "red")
  lines(density(results_asymm[[j]]["BFGS", 1, ]), col = "red", lty = 2)
  abline(v = true_beta[1], col = "gold")
  legend("topright", col = c(4, 4, "gold"), lty = c(1, 2, 1),
         legend = c("Bayes Symm", "Freq Symm", "True"))
  legend("topleft", col = c(2, 2), lty = c(1, 2),
         legend = c("Bayes Asymm", "Freq Asymm"))
}
dev.off()

# Credible intervals frequentist coverage ----
# First regression coefficient
cov_beta1_quant_symm = t(sapply(results_symm, function(R) {
  I = R[c("q025", "q0975"), 2, ]
  rowMeans(apply(I, 2, function(i) unname(c(true_beta[2] >= i[1] & true_beta[2] <= i[2], i[2]-i[1]))))
}))
colnames(cov_beta1_quant_symm) = c("cov_quant_symm", "length_quant_symm")

cov_beta1_quant_asymm = t(sapply(results_asymm, function(R) {
  I = R[c("q025", "q0975"), 2, ]
  rowMeans(apply(I, 2, function(i) unname(c(true_beta[2] >= i[1] & true_beta[2] <= i[2], i[2]-i[1]))))
}))
colnames(cov_beta1_quant_asymm) = c("cov_quant_asymm", "length_quant_asymm")

cov_beta1_HDI_symm = t(sapply(results_symm, function(R) {
  I1 = R[c("HDI_L1", "HDI_U1"), 2, ]
  I2 = R[c("HDI_L2", "HDI_U2"), 2, ]
  tmp1 = apply(I1, 2, function(i) true_beta[2] >= i[1] & true_beta[2] <= i[2])
  tmp2 = apply(I2, 2, function(i) true_beta[2] >= i[1] & true_beta[2] <= i[2])
  tmp2[is.na(tmp2)] = FALSE
  len1 = diff(I1)
  len2 = ifelse(is.na(diff(I2)), 0, diff(I2))
  c(mean(tmp1 | tmp2), mean(len1+len2))
}))
colnames(cov_beta1_HDI_symm) = c("cov_HDI_symm", "length_HDI_symm")

cov_beta1_HDI_asymm = t(sapply(results_asymm, function(R) {
  I1 = R[c("HDI_L1", "HDI_U1"), 2, ]
  I2 = R[c("HDI_L2", "HDI_U2"), 2, ]
  tmp1 = apply(I1, 2, function(i) true_beta[2] >= i[1] & true_beta[2] <= i[2])
  tmp2 = apply(I2, 2, function(i) true_beta[2] >= i[1] & true_beta[2] <= i[2])
  tmp2[is.na(tmp2)] = FALSE
  len1 = diff(I1)
  len2 = ifelse(is.na(diff(I2)), 0, diff(I2))
  c(mean(tmp1 | tmp2), mean(len1+len2))
}))
colnames(cov_beta1_HDI_asymm) = c("cov_HDI_asymm", "length_HDI_asymm")

df_beta1 = data.frame(sim_settings, 
                      cov_beta1_quant_symm, cov_beta1_quant_asymm, 
                      cov_beta1_HDI_symm, cov_beta1_HDI_asymm)

tab = knitr::kable(df_beta1[ , c(1:3, 4, 6, 5, 7, 8, 10, 9, 11)], format = "pipe", digits = 3)
writeLines(tab, "sim_study_nonC-GBI/univ_beta1.txt")


# Intercept
cov_int_quant_symm = t(sapply(results_symm, function(R) {
  I = R[c("q025", "q0975"), 1, ]
  rowMeans(apply(I, 2, function(i) unname(c(true_beta[1] >= i[1] & true_beta[1] <= i[2], i[2]-i[1]))))
}))
colnames(cov_int_quant_symm) = c("cov_quant_symm", "length_quant_symm")

cov_int_quant_asymm = t(sapply(results_asymm, function(R) {
  I = R[c("q025", "q0975"), 1, ]
  rowMeans(apply(I, 2, function(i) unname(c(true_beta[1] >= i[1] & true_beta[1] <= i[2], i[2]-i[1]))))
}))
colnames(cov_int_quant_asymm) = c("cov_quant_asymm", "length_quant_asymm")

cov_int_HDI_symm = t(sapply(results_symm, function(R) {
  I1 = R[c("HDI_L1", "HDI_U1"), 1, ]
  I2 = R[c("HDI_L2", "HDI_U2"), 1, ]
  tmp1 = apply(I1, 2, function(i) true_beta[1] >= i[1] & true_beta[1] <= i[2])
  tmp2 = apply(I2, 2, function(i) true_beta[1] >= i[1] & true_beta[1] <= i[2])
  tmp2[is.na(tmp2)] = FALSE
  len1 = diff(I1)
  len2 = ifelse(is.na(diff(I2)), 0, diff(I2))
  c(mean(tmp1 | tmp2), mean(len1+len2))
}))
colnames(cov_int_HDI_symm) = c("cov_HDI_symm", "length_HDI_symm")

cov_int_HDI_asymm = t(sapply(results_asymm, function(R) {
  I1 = R[c("HDI_L1", "HDI_U1"), 1, ]
  I2 = R[c("HDI_L2", "HDI_U2"), 1, ]
  tmp1 = apply(I1, 2, function(i) true_beta[1] >= i[1] & true_beta[1] <= i[2])
  tmp2 = apply(I2, 2, function(i) true_beta[1] >= i[1] & true_beta[1] <= i[2])
  tmp2[is.na(tmp2)] = FALSE
  len1 = diff(I1)
  len2 = ifelse(is.na(diff(I2)), 0, diff(I2))
  c(mean(tmp1 | tmp2), mean(len1+len2))
}))
colnames(cov_int_HDI_asymm) = c("cov_HDI_asymm", "length_HDI_asymm")

df_beta0 = data.frame(sim_settings, 
                      cov_int_quant_symm, cov_int_quant_asymm, 
                      cov_int_HDI_symm, cov_int_HDI_asymm)

tab = knitr::kable(df_beta0[ , c(1:3, 4, 6, 5, 7, 8, 10, 9, 11)], format = "pipe", digits = 3)
writeLines(tab, "sim_study_nonC-GBI/univ_beta0.txt")


# Regression coefficient
mean_mean_symm = sapply(results_symm, function(R) rowMeans(R["mean", , ]))
mean_mean_asymm = sapply(results_asymm, function(R) rowMeans(R["mean", , ]))

bias_symm = t(mean_mean_symm - true_beta); colnames(bias_symm) = c("beta0_symm", "beta1_symm")
bias_asymm = t(mean_mean_asymm - true_beta); colnames(bias_asymm) = c("beta0_asymm", "beta1_asymm")

df_bias = data.frame(sim_settings, bias_symm, bias_asymm)[c(1:3, 4, 6, 5, 7)]

tab = knitr::kable(df_bias, format = "pipe", digits = 3)
writeLines(tab, "sim_study_nonC-GBI/univ_bias.txt")
