rm(list = ls())

suppressPackageStartupMessages(library(tidyverse))

input_asymm_TL = readRDS("sim_study_nonC-GBI/02-w_unit_loss/sim_study_nonCalibrated_multiv_asymm_TL.RDS")

results_symm = readRDS("sim_study_nonC-GBI/02-w_unit_loss/sim_study_nonCalibrated_multiv.RDS")
results_asymm = readRDS("sim_study_nonC-GBI/02-w_unit_loss/sim_study_nonCalibrated_multiv_asymm.RDS")
results_asymm_TL = lapply(input_asymm_TL,
                          function(I)
                            abind::abind(lapply(I, function(x) x$summ), along = 3))


# TRUE BETA AND SAMPLE SIZE ----
n = 1e3
true_beta = c(0, -0.5, +0.7)

# SIMULATION SETTINGS ----
error_distr = c("Gaussian", "Gaussian_0.5", "Gamma_2_2")
covariate = c("0.5", "1", "asymm")
k = c("0.5", "init", "final")
sim_settings = as.matrix(expand.grid(error_distr, covariate, k))
colnames(sim_settings) = c("err_distr", "unif_cov", "k")

# Plot point estimates ----
## Mean ----
# First regression coefficient
pdf("sim_study_nonC-GBI/multiv_beta1.pdf", width = 9, height = 6)
par(mfrow = c(3, 3), mar = c(3, 2, 2, 1))
for (j in seq(1, 27, by = 3)){
  min_x = -1
  max_x = 0
  plot(1, xlim = c(min_x, max_x), ylim = c(0, 12), type = "n", main = paste(sim_settings[j, ], collapse = " - "),
       xlab = "", ylab = "")
  lines(density(results_symm[[j]]["mean", 2, ]), col = "blue")
  lines(density(results_symm[[j]]["BFGS", 2, ]), col = "blue", lty = 2)
  lines(density(results_asymm[[j]]["mean", 2, ]), col = "red")
  lines(density(results_asymm[[j]]["BFGS", 2, ]), col = "red", lty = 2)
  lines(density(results_asymm_TL[[j]]["mean", 2, ]), col = "green3")
  lines(density(results_asymm_TL[[j]]["BFGS", 2, ]), col = "green3", lty = 2)
  abline(v = true_beta[2], col = "gold")
  legend("topright", col = c(4, 4, "gold"), lty = c(1, 2, 1),
         legend = c("Bayes Symm", "Freq Symm", "True"), cex = 0.8)
  legend("topleft", col = c(2, 2, "green3", "green3"), lty = c(1, 2),
         legend = c("Bayes Asymm", "Freq Asymm", "Bayes TL", "Freq TL"), cex = 0.8)
}

for (j in seq(2, 27, by = 3)){
  min_x = -1
  max_x = 0
  plot(1, xlim = c(min_x, max_x), ylim = c(0, 12), type = "n", main = paste(sim_settings[j, ], collapse = " - "),
       xlab = "", ylab = "")
  lines(density(results_symm[[j]]["mean", 2, ]), col = "blue")
  lines(density(results_symm[[j]]["BFGS", 2, ]), col = "blue", lty = 2)
  lines(density(results_asymm[[j]]["mean", 2, ]), col = "red")
  lines(density(results_asymm[[j]]["BFGS", 2, ]), col = "red", lty = 2)
  lines(density(results_asymm_TL[[j]]["mean", 2, ]), col = "green3")
  lines(density(results_asymm_TL[[j]]["BFGS", 2, ]), col = "green3", lty = 2)
  abline(v = true_beta[2], col = "gold")
  legend("topright", col = c(4, 4, "gold"), lty = c(1, 2, 1),
         legend = c("Bayes Symm", "Freq Symm", "True"), cex = 0.8)
  legend("topleft", col = c(2, 2, "green3", "green3"), lty = c(1, 2),
         legend = c("Bayes Asymm", "Freq Asymm", "Bayes TL", "Freq TL"), cex = 0.8)
}

for (j in seq(3, 27, by = 3)){
  min_x = -1
  max_x = 0
  plot(1, xlim = c(min_x, max_x), ylim = c(0, 12), type = "n", main = paste(sim_settings[j, ], collapse = " - "),
       xlab = "", ylab = "")
  lines(density(results_symm[[j]]["mean", 2, ]), col = "blue")
  lines(density(results_symm[[j]]["BFGS", 2, ]), col = "blue", lty = 2)
  lines(density(results_asymm[[j]]["mean", 2, ]), col = "red")
  lines(density(results_asymm[[j]]["BFGS", 2, ]), col = "red", lty = 2)
  lines(density(results_asymm_TL[[j]]["mean", 2, ]), col = "green3")
  lines(density(results_asymm_TL[[j]]["BFGS", 2, ]), col = "green3", lty = 2)
  abline(v = true_beta[2], col = "gold")
  legend("topright", col = c(4, 4, "gold"), lty = c(1, 2, 1),
         legend = c("Bayes Symm", "Freq Symm", "True"), cex = 0.8)
  legend("topleft", col = c(2, 2, "green3", "green3"), lty = c(1, 2),
         legend = c("Bayes Asymm", "Freq Asymm", "Bayes TL", "Freq TL"), cex = 0.8)
}
dev.off()

# Second regression coefficient
pdf("sim_study_nonC-GBI/multiv_beta2.pdf", width = 9, height = 6)
par(mfrow = c(3, 3), mar = c(3, 2, 2, 1))
for (j in seq(1, 27, by = 3)){
  min_x = 0.25
  max_x = 1.25
  plot(1, xlim = c(min_x, max_x), ylim = c(0, 12), type = "n", main = paste(sim_settings[j, ], collapse = " - "),
       xlab = "", ylab = "")
  lines(density(results_symm[[j]]["mean", 3, ]), col = "blue")
  lines(density(results_symm[[j]]["BFGS", 3, ]), col = "blue", lty = 2)
  lines(density(results_asymm[[j]]["mean", 3, ]), col = "red")
  lines(density(results_asymm[[j]]["BFGS", 3, ]), col = "red", lty = 2)
  lines(density(results_asymm_TL[[j]]["mean", 3, ]), col = "green3")
  lines(density(results_asymm_TL[[j]]["BFGS", 3, ]), col = "green3", lty = 2)
  abline(v = true_beta[3], col = "gold")
  legend("topright", col = c(4, 4, "gold"), lty = c(1, 2, 1),
         legend = c("Bayes Symm", "Freq Symm", "True"), cex = 0.8)
  legend("topleft", col = c(2, 2, "green3", "green3"), lty = c(1, 2),
         legend = c("Bayes Asymm", "Freq Asymm", "Bayes TL", "Freq TL"), cex = 0.8)
}

for (j in seq(2, 27, by = 3)){
  min_x = 0.25
  max_x = 1.25
  plot(1, xlim = c(min_x, max_x), ylim = c(0, 12), type = "n", main = paste(sim_settings[j, ], collapse = " - "),
       xlab = "", ylab = "")
  lines(density(results_symm[[j]]["mean", 3, ]), col = "blue")
  lines(density(results_symm[[j]]["BFGS", 3, ]), col = "blue", lty = 2)
  lines(density(results_asymm[[j]]["mean", 3, ]), col = "red")
  lines(density(results_asymm[[j]]["BFGS", 3, ]), col = "red", lty = 2)
  lines(density(results_asymm_TL[[j]]["mean", 3, ]), col = "green3")
  lines(density(results_asymm_TL[[j]]["BFGS", 3, ]), col = "green3", lty = 2)
  abline(v = true_beta[3], col = "gold")
  legend("topright", col = c(4, 4, "gold"), lty = c(1, 2, 1),
         legend = c("Bayes Symm", "Freq Symm", "True"), cex = 0.8)
  legend("topleft", col = c(2, 2, "green3", "green3"), lty = c(1, 2),
         legend = c("Bayes Asymm", "Freq Asymm", "Bayes TL", "Freq TL"), cex = 0.8)
}

for (j in seq(3, 27, by = 3)){
  min_x = 0.25
  max_x = 1.25
  plot(1, xlim = c(min_x, max_x), ylim = c(0, 12), type = "n", main = paste(sim_settings[j, ], collapse = " - "),
       xlab = "", ylab = "")
  lines(density(results_symm[[j]]["mean", 3, ]), col = "blue")
  lines(density(results_symm[[j]]["BFGS", 3, ]), col = "blue", lty = 2)
  lines(density(results_asymm[[j]]["mean", 3, ]), col = "red")
  lines(density(results_asymm[[j]]["BFGS", 3, ]), col = "red", lty = 2)
  lines(density(results_asymm_TL[[j]]["mean", 3, ]), col = "green3")
  lines(density(results_asymm_TL[[j]]["BFGS", 3, ]), col = "green3", lty = 2)
  abline(v = true_beta[3], col = "gold")
  legend("topright", col = c(4, 4, "gold"), lty = c(1, 2, 1),
         legend = c("Bayes Symm", "Freq Symm", "True"), cex = 0.8)
  legend("topleft", col = c(2, 2, "green3", "green3"), lty = c(1, 2),
         legend = c("Bayes Asymm", "Freq Asymm", "Bayes TL", "Freq TL"), cex = 0.8)
}
dev.off()

# Intercept
pdf("sim_study_nonC-GBI/multiv_beta0.pdf", width = 9, height = 6)
par(mfrow = c(3, 3), mar = c(3, 2, 2, 1))
for (j in seq(1, 27, by = 3)){
  min_x = -0.25
  max_x = +0.7
  plot(1, xlim = c(min_x, max_x), ylim = c(0, 12), type = "n", main = paste(sim_settings[j, ], collapse = " - "),
       xlab = "", ylab = "")
  lines(density(results_symm[[j]]["mean", 1, ]), col = "blue")
  lines(density(results_symm[[j]]["BFGS", 1, ]), col = "blue", lty = 2)
  lines(density(results_asymm[[j]]["mean", 1, ]), col = "red")
  lines(density(results_asymm[[j]]["BFGS", 1, ]), col = "red", lty = 2)
  lines(density(results_asymm_TL[[j]]["mean", 1, ]), col = "green3")
  lines(density(results_asymm_TL[[j]]["BFGS", 1, ]), col = "green3", lty = 2)
  abline(v = true_beta[1], col = "gold")
  legend("topright", col = c(4, 4, "gold"), lty = c(1, 2, 1),
         legend = c("Bayes Symm", "Freq Symm", "True"), cex = 0.8)
  legend("topleft", col = c(2, 2, "green3", "green3"), lty = c(1, 2),
         legend = c("Bayes Asymm", "Freq Asymm", "Bayes TL", "Freq TL"), cex = 0.8)
}

for (j in seq(2, 27, by = 3)){
  min_x = -0.25
  max_x = +0.7
  plot(1, xlim = c(min_x, max_x), ylim = c(0, 12), type = "n", main = paste(sim_settings[j, ], collapse = " - "),
       xlab = "", ylab = "")
  lines(density(results_symm[[j]]["mean", 1, ]), col = "blue")
  lines(density(results_symm[[j]]["BFGS", 1, ]), col = "blue", lty = 2)
  lines(density(results_asymm[[j]]["mean", 1, ]), col = "red")
  lines(density(results_asymm[[j]]["BFGS", 1, ]), col = "red", lty = 2)
  lines(density(results_asymm_TL[[j]]["mean", 1, ]), col = "green3")
  lines(density(results_asymm_TL[[j]]["BFGS", 1, ]), col = "green3", lty = 2)
  abline(v = true_beta[1], col = "gold")
  legend("topright", col = c(4, 4, "gold"), lty = c(1, 2, 1),
         legend = c("Bayes Symm", "Freq Symm", "True"), cex = 0.8)
  legend("topleft", col = c(2, 2, "green3", "green3"), lty = c(1, 2),
         legend = c("Bayes Asymm", "Freq Asymm", "Bayes TL", "Freq TL"), cex = 0.8)
}

for (j in seq(3, 27, by = 3)){
  min_x = -0.25 #min(c(results[[j]]["mean", 2, ], results[[j]]["BFGS", 2, ])) 
  max_x = +0.7 #max(c(results[[j]]["mean", 2, ], results[[j]]["BFGS", 2, ])) 
  plot(1, xlim = c(min_x, max_x), ylim = c(0, 12), type = "n", main = paste(sim_settings[j, ], collapse = " - "),
       xlab = "", ylab = "")
  lines(density(results_symm[[j]]["mean", 1, ]), col = "blue")
  lines(density(results_symm[[j]]["BFGS", 1, ]), col = "blue", lty = 2)
  lines(density(results_asymm[[j]]["mean", 1, ]), col = "red")
  lines(density(results_asymm[[j]]["BFGS", 1, ]), col = "red", lty = 2)
  lines(density(results_asymm_TL[[j]]["mean", 1, ]), col = "green3")
  lines(density(results_asymm_TL[[j]]["BFGS", 1, ]), col = "green3", lty = 2)
  abline(v = true_beta[1], col = "gold")
  legend("topright", col = c(4, 4, "gold"), lty = c(1, 2, 1),
         legend = c("Bayes Symm", "Freq Symm", "True"), cex = 0.8)
  legend("topleft", col = c(2, 2, "green3", "green3"), lty = c(1, 2),
         legend = c("Bayes Asymm", "Freq Asymm", "Bayes TL", "Freq TL"), cex = 0.8)
}
dev.off()


## Median ----
# First regression coefficient
pdf("sim_study_nonC-GBI/multiv_beta1_median.pdf", width = 9, height = 6)
par(mfrow = c(3, 3), mar = c(3, 2, 2, 1))
for (j in seq(1, 27, by = 3)){
  min_x = -1
  max_x = 0
  plot(1, xlim = c(min_x, max_x), ylim = c(0, 12), type = "n", main = paste(sim_settings[j, ], collapse = " - "),
       xlab = "", ylab = "")
  lines(density(results_symm[[j]]["median", 2, ]), col = "blue")
  lines(density(results_symm[[j]]["BFGS", 2, ]), col = "blue", lty = 2)
  lines(density(results_asymm[[j]]["median", 2, ]), col = "red")
  lines(density(results_asymm[[j]]["BFGS", 2, ]), col = "red", lty = 2)
  lines(density(results_asymm_TL[[j]]["median", 2, ]), col = "green3")
  lines(density(results_asymm_TL[[j]]["BFGS", 2, ]), col = "green3", lty = 2)
  abline(v = true_beta[2], col = "gold")
  legend("topright", col = c(4, 4, "gold"), lty = c(1, 2, 1),
         legend = c("Bayes Symm", "Freq Symm", "True"), cex = 0.8)
  legend("topleft", col = c(2, 2, "green3", "green3"), lty = c(1, 2),
         legend = c("Bayes Asymm", "Freq Asymm", "Bayes TL", "Freq TL"), cex = 0.8)
}

for (j in seq(2, 27, by = 3)){
  min_x = -1
  max_x = 0
  plot(1, xlim = c(min_x, max_x), ylim = c(0, 12), type = "n", main = paste(sim_settings[j, ], collapse = " - "),
       xlab = "", ylab = "")
  lines(density(results_symm[[j]]["median", 2, ]), col = "blue")
  lines(density(results_symm[[j]]["BFGS", 2, ]), col = "blue", lty = 2)
  lines(density(results_asymm[[j]]["median", 2, ]), col = "red")
  lines(density(results_asymm[[j]]["BFGS", 2, ]), col = "red", lty = 2)
  lines(density(results_asymm_TL[[j]]["median", 2, ]), col = "green3")
  lines(density(results_asymm_TL[[j]]["BFGS", 2, ]), col = "green3", lty = 2)
  abline(v = true_beta[2], col = "gold")
  legend("topright", col = c(4, 4, "gold"), lty = c(1, 2, 1),
         legend = c("Bayes Symm", "Freq Symm", "True"), cex = 0.8)
  legend("topleft", col = c(2, 2, "green3", "green3"), lty = c(1, 2),
         legend = c("Bayes Asymm", "Freq Asymm", "Bayes TL", "Freq TL"), cex = 0.8)
}

for (j in seq(3, 27, by = 3)){
  min_x = -1
  max_x = 0
  plot(1, xlim = c(min_x, max_x), ylim = c(0, 12), type = "n", main = paste(sim_settings[j, ], collapse = " - "),
       xlab = "", ylab = "")
  lines(density(results_symm[[j]]["median", 2, ]), col = "blue")
  lines(density(results_symm[[j]]["BFGS", 2, ]), col = "blue", lty = 2)
  lines(density(results_asymm[[j]]["median", 2, ]), col = "red")
  lines(density(results_asymm[[j]]["BFGS", 2, ]), col = "red", lty = 2)
  lines(density(results_asymm_TL[[j]]["median", 2, ]), col = "green3")
  lines(density(results_asymm_TL[[j]]["BFGS", 2, ]), col = "green3", lty = 2)
  abline(v = true_beta[2], col = "gold")
  legend("topright", col = c(4, 4, "gold"), lty = c(1, 2, 1),
         legend = c("Bayes Symm", "Freq Symm", "True"), cex = 0.8)
  legend("topleft", col = c(2, 2, "green3", "green3"), lty = c(1, 2),
         legend = c("Bayes Asymm", "Freq Asymm", "Bayes TL", "Freq TL"), cex = 0.8)
}
dev.off()

# Second regression coefficient
pdf("sim_study_nonC-GBI/multiv_beta2_median.pdf", width = 9, height = 6)
par(mfrow = c(3, 3), mar = c(3, 2, 2, 1))
for (j in seq(1, 27, by = 3)){
  min_x = 0.25
  max_x = 1.25
  plot(1, xlim = c(min_x, max_x), ylim = c(0, 12), type = "n", main = paste(sim_settings[j, ], collapse = " - "),
       xlab = "", ylab = "")
  lines(density(results_symm[[j]]["median", 3, ]), col = "blue")
  lines(density(results_symm[[j]]["BFGS", 3, ]), col = "blue", lty = 2)
  lines(density(results_asymm[[j]]["median", 3, ]), col = "red")
  lines(density(results_asymm[[j]]["BFGS", 3, ]), col = "red", lty = 2)
  lines(density(results_asymm_TL[[j]]["median", 3, ]), col = "green3")
  lines(density(results_asymm_TL[[j]]["BFGS", 3, ]), col = "green3", lty = 2)
  abline(v = true_beta[3], col = "gold")
  legend("topright", col = c(4, 4, "gold"), lty = c(1, 2, 1),
         legend = c("Bayes Symm", "Freq Symm", "True"), cex = 0.8)
  legend("topleft", col = c(2, 2, "green3", "green3"), lty = c(1, 2),
         legend = c("Bayes Asymm", "Freq Asymm", "Bayes TL", "Freq TL"), cex = 0.8)
}

for (j in seq(2, 27, by = 3)){
  min_x = 0.25
  max_x = 1.25
  plot(1, xlim = c(min_x, max_x), ylim = c(0, 12), type = "n", main = paste(sim_settings[j, ], collapse = " - "),
       xlab = "", ylab = "")
  lines(density(results_symm[[j]]["median", 3, ]), col = "blue")
  lines(density(results_symm[[j]]["BFGS", 3, ]), col = "blue", lty = 2)
  lines(density(results_asymm[[j]]["median", 3, ]), col = "red")
  lines(density(results_asymm[[j]]["BFGS", 3, ]), col = "red", lty = 2)
  lines(density(results_asymm_TL[[j]]["median", 3, ]), col = "green3")
  lines(density(results_asymm_TL[[j]]["BFGS", 3, ]), col = "green3", lty = 2)
  abline(v = true_beta[3], col = "gold")
  legend("topright", col = c(4, 4, "gold"), lty = c(1, 2, 1),
         legend = c("Bayes Symm", "Freq Symm", "True"), cex = 0.8)
  legend("topleft", col = c(2, 2, "green3", "green3"), lty = c(1, 2),
         legend = c("Bayes Asymm", "Freq Asymm", "Bayes TL", "Freq TL"), cex = 0.8)
}

for (j in seq(3, 27, by = 3)){
  min_x = 0.25
  max_x = 1.25
  plot(1, xlim = c(min_x, max_x), ylim = c(0, 12), type = "n", main = paste(sim_settings[j, ], collapse = " - "),
       xlab = "", ylab = "")
  lines(density(results_symm[[j]]["median", 3, ]), col = "blue")
  lines(density(results_symm[[j]]["BFGS", 3, ]), col = "blue", lty = 2)
  lines(density(results_asymm[[j]]["median", 3, ]), col = "red")
  lines(density(results_asymm[[j]]["BFGS", 3, ]), col = "red", lty = 2)
  lines(density(results_asymm_TL[[j]]["median", 3, ]), col = "green3")
  lines(density(results_asymm_TL[[j]]["BFGS", 3, ]), col = "green3", lty = 2)
  abline(v = true_beta[3], col = "gold")
  legend("topright", col = c(4, 4, "gold"), lty = c(1, 2, 1),
         legend = c("Bayes Symm", "Freq Symm", "True"), cex = 0.8)
  legend("topleft", col = c(2, 2, "green3", "green3"), lty = c(1, 2),
         legend = c("Bayes Asymm", "Freq Asymm", "Bayes TL", "Freq TL"), cex = 0.8)
}
dev.off()

# Intercept
pdf("sim_study_nonC-GBI/multiv_beta0_median.pdf", width = 9, height = 6)
par(mfrow = c(3, 3), mar = c(3, 2, 2, 1))
for (j in seq(1, 27, by = 3)){
  min_x = -0.25
  max_x = +0.7
  plot(1, xlim = c(min_x, max_x), ylim = c(0, 12), type = "n", main = paste(sim_settings[j, ], collapse = " - "),
       xlab = "", ylab = "")
  lines(density(results_symm[[j]]["median", 1, ]), col = "blue")
  lines(density(results_symm[[j]]["BFGS", 1, ]), col = "blue", lty = 2)
  lines(density(results_asymm[[j]]["median", 1, ]), col = "red")
  lines(density(results_asymm[[j]]["BFGS", 1, ]), col = "red", lty = 2)
  lines(density(results_asymm_TL[[j]]["median", 1, ]), col = "green3")
  lines(density(results_asymm_TL[[j]]["BFGS", 1, ]), col = "green3", lty = 2)
  abline(v = true_beta[1], col = "gold")
  legend("topright", col = c(4, 4, "gold"), lty = c(1, 2, 1),
         legend = c("Bayes Symm", "Freq Symm", "True"), cex = 0.8)
  legend("topleft", col = c(2, 2, "green3", "green3"), lty = c(1, 2),
         legend = c("Bayes Asymm", "Freq Asymm", "Bayes TL", "Freq TL"), cex = 0.8)
}

for (j in seq(2, 27, by = 3)){
  min_x = -0.25
  max_x = +0.7
  plot(1, xlim = c(min_x, max_x), ylim = c(0, 12), type = "n", main = paste(sim_settings[j, ], collapse = " - "),
       xlab = "", ylab = "")
  lines(density(results_symm[[j]]["median", 1, ]), col = "blue")
  lines(density(results_symm[[j]]["BFGS", 1, ]), col = "blue", lty = 2)
  lines(density(results_asymm[[j]]["median", 1, ]), col = "red")
  lines(density(results_asymm[[j]]["BFGS", 1, ]), col = "red", lty = 2)
  lines(density(results_asymm_TL[[j]]["median", 1, ]), col = "green3")
  lines(density(results_asymm_TL[[j]]["BFGS", 1, ]), col = "green3", lty = 2)
  abline(v = true_beta[1], col = "gold")
  legend("topright", col = c(4, 4, "gold"), lty = c(1, 2, 1),
         legend = c("Bayes Symm", "Freq Symm", "True"), cex = 0.8)
  legend("topleft", col = c(2, 2, "green3", "green3"), lty = c(1, 2),
         legend = c("Bayes Asymm", "Freq Asymm", "Bayes TL", "Freq TL"), cex = 0.8)
}

for (j in seq(3, 27, by = 3)){
  min_x = -0.25
  max_x = +0.7
  plot(1, xlim = c(min_x, max_x), ylim = c(0, 12), type = "n", main = paste(sim_settings[j, ], collapse = " - "),
       xlab = "", ylab = "")
  lines(density(results_symm[[j]]["median", 1, ]), col = "blue")
  lines(density(results_symm[[j]]["BFGS", 1, ]), col = "blue", lty = 2)
  lines(density(results_asymm[[j]]["median", 1, ]), col = "red")
  lines(density(results_asymm[[j]]["BFGS", 1, ]), col = "red", lty = 2)
  lines(density(results_asymm_TL[[j]]["median", 1, ]), col = "green3")
  lines(density(results_asymm_TL[[j]]["BFGS", 1, ]), col = "green3", lty = 2)
  abline(v = true_beta[1], col = "gold")
  legend("topright", col = c(4, 4, "gold"), lty = c(1, 2, 1),
         legend = c("Bayes Symm", "Freq Symm", "True"), cex = 0.8)
  legend("topleft", col = c(2, 2, "green3", "green3"), lty = c(1, 2),
         legend = c("Bayes Asymm", "Freq Asymm", "Bayes TL", "Freq TL"), cex = 0.8)
}
dev.off()

# Credible intervals frequentist coverage ----
# First regression coefficient
cov_beta1_quant_symm = t(sapply(results_symm, function(R) {
  I = R[c("q025", "q975"), 2, ]
  rowMeans(apply(I, 2, function(i) unname(c(true_beta[2] >= i[1] & true_beta[2] <= i[2], i[2]-i[1]))))
}))
colnames(cov_beta1_quant_symm) = c("cov_q_symm", "len_q_symm")

cov_beta1_quant_asymm = t(sapply(results_asymm, function(R) {
  I = R[c("q025", "q975"), 2, ]
  rowMeans(apply(I, 2, function(i) unname(c(true_beta[2] >= i[1] & true_beta[2] <= i[2], i[2]-i[1]))))
}))
colnames(cov_beta1_quant_asymm) = c("cov_q_asymm", "len_q_asymm")

cov_beta1_quant_asymm_TL = t(sapply(results_asymm_TL, function(R) {
  I = R[c("q025", "q975"), 2, ]
  rowMeans(apply(I, 2, function(i) unname(c(true_beta[2] >= i[1] & true_beta[2] <= i[2], i[2]-i[1]))))
}))
colnames(cov_beta1_quant_asymm_TL) = c("cov_q_asymm_TL", "len_q_asymm_TL")


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
colnames(cov_beta1_HDI_symm) = c("cov_HDI_symm", "len_HDI_symm")

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
colnames(cov_beta1_HDI_asymm) = c("cov_HDI_asymm", "len_HDI_asymm")

cov_beta1_HDI_asymm_TL = t(sapply(results_asymm_TL, function(R) {
  I1 = R[c("HDI_L1", "HDI_U1"), 2, ]
  I2 = R[c("HDI_L2", "HDI_U2"), 2, ]
  tmp1 = apply(I1, 2, function(i) true_beta[2] >= i[1] & true_beta[2] <= i[2])
  tmp2 = apply(I2, 2, function(i) true_beta[2] >= i[1] & true_beta[2] <= i[2])
  tmp2[is.na(tmp2)] = FALSE
  len1 = diff(I1)
  len2 = ifelse(is.na(diff(I2)), 0, diff(I2))
  c(mean(tmp1 | tmp2), mean(len1+len2))
}))
colnames(cov_beta1_HDI_asymm_TL) = c("cov_HDI_asymm_TL", "len_HDI_asymm_TL")

df_beta1 = data.frame(sim_settings, 
                      cov_beta1_quant_symm, cov_beta1_quant_asymm, cov_beta1_quant_asymm_TL,
                      cov_beta1_HDI_symm, cov_beta1_HDI_asymm, cov_beta1_HDI_asymm_TL)

# tab = knitr::kable(df_beta1[ , c(1:3, 4, 6, 5, 7, 8, 10, 9, 11)], format = "pipe", digits = 3)
tab = knitr::kable(df_beta1[ , c(1:3, 4, 6, 8, 5, 7, 9, 10, 12, 14, 11, 13, 15)], format = "pipe", digits = 3)
writeLines(tab, "sim_study_nonC-GBI/multiv_beta1.txt")

# Second regression coefficient
cov_beta2_quant_symm = t(sapply(results_symm, function(R) {
  I = R[c("q025", "q975"), 3, ]
  rowMeans(apply(I, 2, function(i) unname(c(true_beta[3] >= i[1] & true_beta[3] <= i[2], i[2]-i[1]))))
}))
colnames(cov_beta2_quant_symm) = c("cov_q_symm", "len_q_symm")

cov_beta2_quant_asymm = t(sapply(results_asymm, function(R) {
  I = R[c("q025", "q975"), 3, ]
  rowMeans(apply(I, 2, function(i) unname(c(true_beta[3] >= i[1] & true_beta[3] <= i[2], i[2]-i[1]))))
}))
colnames(cov_beta2_quant_asymm) = c("cov_q_asymm", "len_q_asymm")

cov_beta2_quant_asymm_TL = t(sapply(results_asymm_TL, function(R) {
  I = R[c("q025", "q975"), 3, ]
  rowMeans(apply(I, 2, function(i) unname(c(true_beta[3] >= i[1] & true_beta[3] <= i[2], i[2]-i[1]))))
}))
colnames(cov_beta2_quant_asymm_TL) = c("cov_q_asymm_TL", "len_q_asymm_TL")

cov_beta2_HDI_symm = t(sapply(results_symm, function(R) {
  I1 = R[c("HDI_L1", "HDI_U1"), 3, ]
  I2 = R[c("HDI_L2", "HDI_U2"), 3, ]
  tmp1 = apply(I1, 2, function(i) true_beta[3] >= i[1] & true_beta[3] <= i[2])
  tmp2 = apply(I2, 2, function(i) true_beta[3] >= i[1] & true_beta[3] <= i[2])
  tmp2[is.na(tmp2)] = FALSE
  len1 = diff(I1)
  len2 = ifelse(is.na(diff(I2)), 0, diff(I2))
  c(mean(tmp1 | tmp2), mean(len1+len2))
}))
colnames(cov_beta2_HDI_symm) = c("cov_HDI_symm", "len_HDI_symm")

cov_beta2_HDI_asymm = t(sapply(results_asymm, function(R) {
  I1 = R[c("HDI_L1", "HDI_U1"), 3, ]
  I2 = R[c("HDI_L2", "HDI_U2"), 3, ]
  tmp1 = apply(I1, 2, function(i) true_beta[3] >= i[1] & true_beta[3] <= i[2])
  tmp2 = apply(I2, 2, function(i) true_beta[3] >= i[1] & true_beta[3] <= i[2])
  tmp2[is.na(tmp2)] = FALSE
  len1 = diff(I1)
  len2 = ifelse(is.na(diff(I2)), 0, diff(I2))
  c(mean(tmp1 | tmp2), mean(len1+len2))
}))
colnames(cov_beta2_HDI_asymm) = c("cov_HDI_asymm", "len_HDI_asymm")

cov_beta2_HDI_asymm_TL = t(sapply(results_asymm_TL, function(R) {
  I1 = R[c("HDI_L1", "HDI_U1"), 3, ]
  I2 = R[c("HDI_L2", "HDI_U2"), 3, ]
  tmp1 = apply(I1, 2, function(i) true_beta[3] >= i[1] & true_beta[3] <= i[2])
  tmp2 = apply(I2, 2, function(i) true_beta[3] >= i[1] & true_beta[3] <= i[2])
  tmp2[is.na(tmp2)] = FALSE
  len1 = diff(I1)
  len2 = ifelse(is.na(diff(I2)), 0, diff(I2))
  c(mean(tmp1 | tmp2), mean(len1+len2))
}))
colnames(cov_beta2_HDI_asymm_TL) = c("cov_HDI_asymm_TL", "len_HDI_asymm_TL")

df_beta2 = data.frame(sim_settings, 
                      cov_beta2_quant_symm, cov_beta2_quant_asymm, cov_beta2_quant_asymm_TL,
                      cov_beta2_HDI_symm, cov_beta2_HDI_asymm, cov_beta2_HDI_asymm_TL)

# tab = knitr::kable(df_beta2[ , c(1:3, 4, 6, 5, 7, 8, 10, 9, 11)], format = "pipe", digits = 3)
tab = knitr::kable(df_beta2[ , c(1:3, 4, 6, 8, 5, 7, 9, 10, 12, 14, 11, 13, 15)], format = "pipe", digits = 3)
writeLines(tab, "sim_study_nonC-GBI/multiv_beta2.txt")

# Intercept
cov_int_quant_symm = t(sapply(results_symm, function(R) {
  I = R[c("q025", "q975"), 1, ]
  rowMeans(apply(I, 2, function(i) unname(c(true_beta[1] >= i[1] & true_beta[1] <= i[2], i[2]-i[1]))))
}))
colnames(cov_int_quant_symm) = c("cov_q_symm", "len_q_symm")

cov_int_quant_asymm = t(sapply(results_asymm, function(R) {
  I = R[c("q025", "q975"), 1, ]
  rowMeans(apply(I, 2, function(i) unname(c(true_beta[1] >= i[1] & true_beta[1] <= i[2], i[2]-i[1]))))
}))
colnames(cov_int_quant_asymm) = c("cov_q_asymm", "len_q_asymm")

cov_int_quant_asymm_TL = t(sapply(results_asymm_TL, function(R) {
  I = R[c("q025", "q975"), 1, ]
  rowMeans(apply(I, 2, function(i) unname(c(true_beta[1] >= i[1] & true_beta[1] <= i[2], i[2]-i[1]))))
}))
colnames(cov_int_quant_asymm_TL) = c("cov_q_asymm_TL", "len_q_asymm_TL")


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
colnames(cov_int_HDI_symm) = c("cov_HDI_symm", "len_HDI_symm")

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
colnames(cov_int_HDI_asymm) = c("cov_HDI_asymm", "len_HDI_asymm")

cov_int_HDI_asymm_TL = t(sapply(results_asymm_TL, function(R) {
  I1 = R[c("HDI_L1", "HDI_U1"), 1, ]
  I2 = R[c("HDI_L2", "HDI_U2"), 1, ]
  tmp1 = apply(I1, 2, function(i) true_beta[1] >= i[1] & true_beta[1] <= i[2])
  tmp2 = apply(I2, 2, function(i) true_beta[1] >= i[1] & true_beta[1] <= i[2])
  tmp2[is.na(tmp2)] = FALSE
  len1 = diff(I1)
  len2 = ifelse(is.na(diff(I2)), 0, diff(I2))
  c(mean(tmp1 | tmp2), mean(len1+len2))
}))
colnames(cov_int_HDI_asymm_TL) = c("cov_HDI_asymm_TL", "len_HDI_asymm_TL")


df_beta0 = data.frame(sim_settings, 
                      cov_int_quant_symm, cov_int_quant_asymm, cov_int_quant_asymm_TL,
                      cov_int_HDI_symm, cov_int_HDI_asymm, cov_int_HDI_asymm_TL)

# tab = knitr::kable(df_beta0[ , c(1:3, 4, 6, 5, 7, 8, 10, 9, 11)], format = "pipe", digits = 3)
tab = knitr::kable(df_beta0[ , c(1:3, 4, 6, 8, 5, 7, 9, 10, 12, 14, 11, 13, 15)], format = "pipe", digits = 3)
writeLines(tab, "sim_study_nonC-GBI/multiv_beta0.txt")

# Regression coefficient
mean_mean_symm = sapply(results_symm, function(R) rowMeans(R["mean", , ]))
mean_mean_asymm = sapply(results_asymm, function(R) rowMeans(R["mean", , ]))
mean_mean_asymm_TL = sapply(results_asymm_TL, function(R) rowMeans(R["mean", , ]))

bias_symm = t(mean_mean_symm - true_beta); colnames(bias_symm) = c("beta0_symm", "beta1_symm", "beta2_symm")
bias_asymm = t(mean_mean_asymm - true_beta); colnames(bias_asymm) = c("beta0_asymm", "beta1_asymm", "beta2_asymm")
bias_asymm_TL = t(mean_mean_asymm_TL - true_beta); colnames(bias_asymm_TL) = c("beta0_asymm_TL", "beta1_asymm_TL", "beta2_asymm_TL")

# df_bias = data.frame(sim_settings, bias_symm, bias_asymm)[c(1:3, 4, 7, 5, 8, 6, 9)]
df_bias = data.frame(sim_settings, bias_symm, bias_asymm, bias_asymm_TL)[c(1:3, 4, 7, 10, 5, 8, 11, 6, 9, 12)]

tab = knitr::kable(df_bias, format = "pipe", digits = 3)
writeLines(tab, "sim_study_nonC-GBI/multiv_bias.txt")


# Analysis w TL ----
w_TL = t(sapply(input_asymm_TL, function(I) sapply(I, function(x) x$w)))
colnames(w_TL) = paste0("replicate", 1:100)
plot_w = cbind.data.frame(sim_settings, w_TL) %>% 
  pivot_longer(replicate1:replicate100, names_to = "replicate", values_to = "w") %>% 
  mutate(k = factor(k, levels = c("0.5", "init", "final"))) %>% 
  ggplot(aes(x = k, y = w)) + 
  geom_boxplot() +
  facet_grid(err_distr ~ unif_cov)
ggsave(plot = plot_w, 
       filename = "boxplot_w_replicates.pdf", path = "sim_study_nonC-GBI/02-w_unit_loss/",
       width = 5, height = 5)

# Average density TL ----
draws_TL = lapply(input_asymm_TL, 
                  function(I) 
                    lapply(I, function(x) x$draws[seq(10, 10^4, by = 10), ]) %>% abind::abind(along = 3)) %>% 
  abind::abind(along = 4) %>% aperm(c(4, 3, 2, 1))

draws_TL[ , 1] = drawsTL[ , 1] - draws_TL[ , 2]*mu_x1/sd_x1 - draws_TL[ , 3]*mu_x2/sd_x2
draws[ , 2] = draws_TL[ , 2]/sd_x1
draws[ , 3] = draws_TL[ , 3]/sd_x2

min_b0 = min(draws_TL[ , , 1, ]); max_b0 = max(draws_TL[ , , 1, ])
min_b1 = min(draws_TL[ , , 2, ]); max_b1 = max(draws_TL[ , , 2, ])
min_b2 = min(draws_TL[ , , 3, ]); max_b2 = max(draws_TL[ , , 3, ])
bw_b0 = (max_b0*1.05 - min_b0*1.05)/99
bw_b1 = (max_b1*1.05 - min_b1*1.05)/99
bw_b2 = (max_b2*1.05 - min_b2*1.05)/99

dens_b0 = apply(draws_TL[ , , 1, ], 1,
                function(D1) apply(D1, 1,
                                   function(D2){
                                     tmp = density(D2, from = min_b0, to = max_b0, n = 100)
                                     cbind(tmp$x, tmp$y)
                                   }, simplify = F) %>% abind::abind(along = 3),
                simplify = F) %>% abind::abind(along = 4) %>% 
  aperm(c(4, 3, 2, 1))

dens_b1 = apply(draws_TL[ , , 2, ], 1,
                function(D1) apply(D1, 1,
                                   function(D2){
                                     tmp = density(D2, from = min_b0, to = max_b0, n = 100)
                                     cbind(tmp$x, tmp$y)
                                   }, simplify = F) %>% abind::abind(along = 3),
                simplify = F) %>% abind::abind(along = 4) %>% 
  aperm(c(4, 3, 2, 1))
dens_b2 = apply(draws_TL[ , , 3, ], 1,
                function(D1) apply(D1, 1,
                                   function(D2){
                                     tmp = density(D2, from = min_b0, to = max_b0, n = 100)
                                     cbind(tmp$x, tmp$y)
                                   }, simplify = F) %>% abind::abind(along = 3),
                simplify = F) %>% abind::abind(along = 4) %>% 
  aperm(c(4, 3, 2, 1))

dimnames(dens_b0) = dimnames(dens_b1) = dimnames(dens_b2) = 
  list("setting" = 1:27, "replicate" = 1:100, "axis" = c("x", "y"), "knot" = 1:100) 

df_dens_b0 = left_join(reshape2::melt(dens_b0),
          cbind.data.frame(setting = 1:27, sim_settings),
          by = "setting") %>% 
  pivot_wider(names_from = axis, values_from = value) %>% 
  group_by(setting, knot, .drop = F) %>% 
  mutate(mean_y = mean(y),
         k = factor(k, levels = c("0.5", "init", "final")))

# Plot b0
df_dens_b0 %>% 
  filter(err_distr == "Gaussian") %>% 
  ggplot(aes(x = x, y = mean_y)) +
  geom_vline(xintercept = true_beta[1], col = "gold") +
  geom_line(aes(y = y, group = replicate), col = "grey", alpha = 0.2) +
  geom_line() + 
  facet_grid(k ~ unif_cov) +
  labs(title = "Gaussian error distribution")

df_dens_b0 %>% 
  filter(err_distr == "Gaussian_0.5") %>% 
  ggplot(aes(x = x, y = mean_y)) +
  geom_vline(xintercept = true_beta[1], col = "gold") +
  geom_line(aes(y = y, group = replicate), col = "grey", alpha = 0.2) +
  geom_line() + 
  facet_grid(k ~ unif_cov) +
  labs(title = "Gaussian error distribution with sd = 0.5")

df_dens_b0 %>% 
  filter(err_distr == "Gamma_2_2") %>% 
  ggplot(aes(x = x, y = mean_y)) +
  geom_vline(xintercept = true_beta[1], col = "gold") +
  geom_line(aes(y = y, group = replicate), col = "grey", alpha = 0.2) +
  geom_line() + 
  facet_grid(k ~ unif_cov) +
  labs(title = "Gamma error distribution with shape = 2 and scale = 2")

