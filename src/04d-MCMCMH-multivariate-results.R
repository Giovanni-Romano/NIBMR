rm(list = ls())

suppressPackageStartupMessages(library(tidyverse))

input_asymm_TL = readRDS("sim_study_nonC-GBI/02-w_unit_loss/sim_study_nonCalibrated_multiv_asymm_TL.RDS")
input_asymm_sd = readRDS("sim_study_nonC-GBI/02-w_unit_loss/sim_study_nonCalibrated_multiv_asymm_sd.RDS")

results_symm = readRDS("sim_study_nonC-GBI/02-w_unit_loss/sim_study_nonCalibrated_multiv.RDS")
# results_asymm = readRDS("sim_study_nonC-GBI/02-w_unit_loss/sim_study_nonCalibrated_multiv_asymm.RDS")
results_asymm = lapply(input_asymm_sd,
                       function(I)
                         abind::abind(lapply(I, function(x) x$summ), along = 3))
results_asymm_TL = lapply(input_asymm_TL,
                          function(I)
                            abind::abind(lapply(I, function(x) x$summ), along = 3))


# TRUE BETA AND SAMPLE SIZE ----
n = 1e3
true_beta = c(0, -0.5, +0.7)

# SIMULATION SETTINGS ----
idx_iter = seq(10, 10^4, by = 10)
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
  max_x = +0.5
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
  max_x = +0.5
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
  min_x = -0.25
  max_x = +0.5
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
  max_x = +0.5
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
  max_x = +0.5
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
  max_x = +0.5
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
w_sd = t(sapply(input_asymm_sd, function(I) sapply(I, function(x) x$w)))
colnames(w_TL) = colnames(w_sd) = paste0("replicate", 1:100)
plot_w = cbind.data.frame(sim_settings, w_TL) %>% 
  pivot_longer(replicate1:replicate100, names_to = "replicate", values_to = "TL") %>% 
  left_join(cbind.data.frame(sim_settings, w_sd) %>% 
              pivot_longer(replicate1:replicate100, names_to = "replicate", values_to = "sd"),
            by = c("err_distr", "unif_cov", "k", "replicate")) %>% 
  mutate(k = factor(k, levels = c("0.5", "init", "final")),
         err_distr = factor(err_distr, levels = c("Gaussian", "Gaussian_0.5", "Gamma_2_2"),
                            labels = c("Gauss_1", "Gauss_0.5", "Gamma_2_2"))) %>% 
  pivot_longer(c("TL", "sd"), names_to = "type", values_to = "w") %>% 
  ggplot(aes(x = err_distr, y = w, color = type)) + 
  geom_boxplot() +
  facet_grid(k ~ unif_cov, scales = "free_y") +
  scale_color_manual(values = c(3, 4)) +
  labs(title = "w values across replicates")
ggsave(plot = plot_w, 
       filename = "w_together.pdf", path = "sim_study_nonC-GBI/",
       width = 9, height = 6)

# Average density TL ----
draws_TL = lapply(input_asymm_TL, 
                  function(I) 
                    lapply(I, function(x) x$draws[seq(10, 10^4, by = 10), ]) %>% abind::abind(along = 3)) %>% 
  abind::abind(along = 4) %>% aperm(c(4, 3, 2, 1))
draws_sd = lapply(input_asymm_sd, 
                  function(I) 
                    lapply(I, function(x) x$draws[seq(10, 10^4, by = 10), ]) %>% abind::abind(along = 3)) %>% 
  abind::abind(along = 4) %>% aperm(c(4, 3, 2, 1))

min_b0 = quantile(c(draws_TL[ , , 1, ], draws_sd[ , , 1, ]), 0.01)
max_b0 = quantile(c(draws_TL[ , , 1, ], draws_sd[ , , 1, ]), 0.99)
min_b1 = quantile(c(draws_TL[ , , 2, ], draws_sd[ , , 2, ]), 0.01)
max_b1 = quantile(c(draws_TL[ , , 2, ], draws_sd[ , , 2, ]), 0.99)
min_b2 = quantile(c(draws_TL[ , , 3, ], draws_sd[ , , 3, ]), 0.01)
max_b2 = quantile(c(draws_TL[ , , 3, ], draws_sd[ , , 3, ]), 0.99)
bw_b0 = (max_b0 - min_b0)/99
bw_b1 = (max_b1 - min_b1)/99
bw_b2 = (max_b2 - min_b2)/99

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

dimnames(dens_b0_TL) = dimnames(dens_b1_TL) = dimnames(dens_b2_TL) = 
  dimnames(dens_b0_sd) = dimnames(dens_b1_sd) = dimnames(dens_b2_sd) =
  list("setting" = 1:27, "replicate" = 1:100, "axis" = c("x", "y"), "knot" = 1:100) 




# Plot b0 ----
df_dens_b0 = left_join(reshape2::melt(dens_b0_TL),
                       cbind.data.frame(setting = 1:27, sim_settings),
                       by = "setting") %>% 
  pivot_wider(names_from = axis, values_from = value) %>%
  mutate(type = "TL") %>% 
  bind_rows(left_join(reshape2::melt(dens_b0_sd),
                      cbind.data.frame(setting = 1:27, sim_settings),
                      by = "setting") %>% 
              pivot_wider(names_from = axis, values_from = value) %>% 
              mutate(type = "sd")) %>% 
  group_by(setting, knot, type, .drop = F) %>% 
  mutate(mean_y = mean(y),
         k = factor(k, levels = c("0.5", "init", "final")))

plt_dens_b0_1 = df_dens_b0 %>%
  filter(err_distr == "Gaussian") %>% 
  ggplot(aes(x = x)) +
  geom_vline(xintercept = true_beta[1], col = "gold") +
  geom_line(aes(y = y, group = interaction(replicate, type), col = type), alpha = 0.1) +
  ggnewscale::new_scale_colour() +
  scale_colour_manual(values = c("sd" = "red3", "TL" = "blue3")) +
  geom_line(aes(y = mean_y, col = type)) + 
  facet_grid(k ~ unif_cov) +
  labs(title = "Density beta0 [light = per replica; dark = average] - Gaussian(0, sd = 1) error distribution", 
       y = "density", x = "beta2") +
  scale_x_continuous(limits = true_beta[1] + c(-.5, +.5))

plt_dens_b0_2 = df_dens_b0 %>%
  filter(err_distr == "Gaussian_0.5") %>% 
  ggplot(aes(x = x)) +
  geom_vline(xintercept = true_beta[1], col = "gold") +
  geom_line(aes(y = y, group = interaction(replicate, type), col = type), alpha = 0.1) +
  ggnewscale::new_scale_colour() +
  scale_colour_manual(values = c("sd" = "red3", "TL" = "blue3")) +
  geom_line(aes(y = mean_y, col = type)) + 
  facet_grid(k ~ unif_cov) +
  labs(title = "Density beta0 [light = per replica; dark = average] - Gaussian(0, sd = 0.5) error distribution", 
       y = "density", x = "beta2") +
  scale_x_continuous(limits = true_beta[1] + c(-.5, +.5))

plt_dens_b0_3 = df_dens_b0 %>%
  filter(err_distr == "Gamma_2_2") %>% 
  ggplot(aes(x = x)) +
  geom_vline(xintercept = true_beta[1], col = "gold") +
  geom_line(aes(y = y, group = interaction(replicate, type), col = type), alpha = 0.1) +
  ggnewscale::new_scale_colour() +
  scale_colour_manual(values = c("sd" = "red3", "TL" = "blue3")) +
  geom_line(aes(y = mean_y, col = type)) + 
  facet_grid(k ~ unif_cov) +
  labs(title = "Density beta0 [light = per replica; dark = average] - Gamma(2, sd = 2) error distribution", 
       y = "density", x = "beta2") +
  scale_x_continuous(limits = true_beta[1] + c(-.5, +.5))

ggsave(plot = gridExtra::marrangeGrob(grobs = list(plt_dens_b0_1, plt_dens_b0_2, plt_dens_b0_3), 
                                      nrow = 1, ncol = 1, top = NULL), 
       filename = "density_b0.pdf", path = "sim_study_nonC-GBI/",
       width = 9, height = 5)


# Plot density b1 ----
df_dens_b1 = left_join(reshape2::melt(dens_b1_TL),
                       cbind.data.frame(setting = 1:27, sim_settings),
                       by = "setting") %>% 
  pivot_wider(names_from = axis, values_from = value) %>%
  mutate(type = "TL") %>% 
  bind_rows(left_join(reshape2::melt(dens_b1_sd),
                      cbind.data.frame(setting = 1:27, sim_settings),
                      by = "setting") %>% 
              pivot_wider(names_from = axis, values_from = value) %>% 
              mutate(type = "sd")) %>% 
  group_by(setting, knot, type, .drop = F) %>% 
  mutate(mean_y = mean(y),
         k = factor(k, levels = c("0.5", "init", "final")))

plt_dens_b1_1 = df_dens_b1 %>%
  filter(err_distr == "Gaussian") %>% 
  ggplot(aes(x = x)) +
  geom_vline(xintercept = true_beta[2], col = "gold") +
  geom_line(aes(y = y, group = interaction(replicate, type), col = type), alpha = 0.1) +
  ggnewscale::new_scale_colour() +
  scale_colour_manual(values = c("sd" = "red3", "TL" = "blue3")) +
  geom_line(aes(y = mean_y, col = type)) + 
  facet_grid(k ~ unif_cov) +
  labs(title = "Density beta1 [light = per replica; dark = average] - Gaussian(0, sd = 1) error distribution", 
       y = "density", x = "beta2") +
  scale_x_continuous(limits = true_beta[2] + c(-.75, +.75))

plt_dens_b1_2 = df_dens_b1 %>%
  filter(err_distr == "Gaussian_0.5") %>% 
  ggplot(aes(x = x)) +
  geom_vline(xintercept = true_beta[2], col = "gold") +
  geom_line(aes(y = y, group = interaction(replicate, type), col = type), alpha = 0.1) +
  ggnewscale::new_scale_colour() +
  scale_colour_manual(values = c("sd" = "red3", "TL" = "blue3")) +
  geom_line(aes(y = mean_y, col = type)) + 
  facet_grid(k ~ unif_cov) +
  labs(title = "Density beta1 [light = per replica; dark = average] - Gaussian(0, sd = 0.5) error distribution", 
       y = "density", x = "beta2") +
  scale_x_continuous(limits = true_beta[2] + c(-.75, +.75))

plt_dens_b1_3 = df_dens_b1 %>%
  filter(err_distr == "Gamma_2_2") %>% 
  ggplot(aes(x = x)) +
  geom_vline(xintercept = true_beta[2], col = "gold") +
  geom_line(aes(y = y, group = interaction(replicate, type), col = type), alpha = 0.1) +
  ggnewscale::new_scale_colour() +
  scale_colour_manual(values = c("sd" = "red3", "TL" = "blue3")) +
  geom_line(aes(y = mean_y, col = type)) + 
  facet_grid(k ~ unif_cov) +
  labs(title = "Density beta1 [light = per replica; dark = average] - Gamma(2, sd = 2) error distribution", 
       y = "density", x = "beta2") +
  scale_x_continuous(limits = true_beta[2] + c(-.75, +.75))

ggsave(plot = gridExtra::marrangeGrob(grobs = list(plt_dens_b1_1, plt_dens_b1_2, plt_dens_b1_3), 
                                      nrow = 1, ncol = 1, top = NULL), 
       filename = "density_b1.pdf", path = "sim_study_nonC-GBI/",
       width = 9, height = 5)


# Plot density b2 ----
df_dens_b2 = left_join(reshape2::melt(dens_b2_TL),
                       cbind.data.frame(setting = 1:27, sim_settings),
                       by = "setting") %>% 
  pivot_wider(names_from = axis, values_from = value) %>%
  mutate(type = "TL") %>% 
  bind_rows(left_join(reshape2::melt(dens_b2_sd),
                      cbind.data.frame(setting = 1:27, sim_settings),
                      by = "setting") %>% 
              pivot_wider(names_from = axis, values_from = value) %>% 
              mutate(type = "sd")) %>% 
  group_by(setting, knot, type, .drop = F) %>% 
  mutate(mean_y = mean(y),
         k = factor(k, levels = c("0.5", "init", "final")))

plt_dens_b2_1 = df_dens_b2 %>%
  filter(err_distr == "Gaussian") %>% 
  ggplot(aes(x = x)) +
  geom_vline(xintercept = true_beta[3], col = "gold") +
  geom_line(aes(y = y, group = interaction(replicate, type), col = type), alpha = 0.1) +
  ggnewscale::new_scale_colour() +
  scale_colour_manual(values = c("sd" = "red3", "TL" = "blue3")) +
  geom_line(aes(y = mean_y, col = type)) + 
  facet_grid(k ~ unif_cov) +
  labs(title = "Density beta2 [light = per replica; dark = average] - Gaussian(0, sd = 1) error distribution", 
       y = "density", x = "beta2") +
  scale_x_continuous(limits = true_beta[3] + c(-.75, +.75))

plt_dens_b2_2 = df_dens_b2 %>%
  filter(err_distr == "Gaussian_0.5") %>% 
  ggplot(aes(x = x)) +
  geom_vline(xintercept = true_beta[3], col = "gold") +
  geom_line(aes(y = y, group = interaction(replicate, type), col = type), alpha = 0.1) +
  ggnewscale::new_scale_colour() +
  scale_colour_manual(values = c("sd" = "red3", "TL" = "blue3")) +
  geom_line(aes(y = mean_y, col = type)) + 
  facet_grid(k ~ unif_cov) +
  labs(title = "Density beta2 [light = per replica; dark = average] - Gaussian(0, sd = 0.5) error distribution", 
       y = "density", x = "beta2") +
  scale_x_continuous(limits = true_beta[3] + c(-.75, +.75))

plt_dens_b2_3 = df_dens_b2 %>%
  filter(err_distr == "Gamma_2_2") %>% 
  ggplot(aes(x = x)) +
  geom_vline(xintercept = true_beta[3], col = "gold") +
  geom_line(aes(y = y, group = interaction(replicate, type), col = type), alpha = 0.1) +
  ggnewscale::new_scale_colour() +
  scale_colour_manual(values = c("sd" = "red3", "TL" = "blue3")) +
  geom_line(aes(y = mean_y, col = type)) + 
  facet_grid(k ~ unif_cov) +
  labs(title = "Density beta2 [light = per replica; dark = average] - Gamma(2, sd = 2) error distribution", 
       y = "density", x = "beta2") +
  scale_x_continuous(limits = true_beta[3] + c(-.75, +.75))

ggsave(plot = gridExtra::marrangeGrob(grobs = list(plt_dens_b2_1, plt_dens_b2_2, plt_dens_b2_3), 
                                      nrow = 1, ncol = 1, top = NULL), 
       filename = "density_b2.pdf", path = "sim_study_nonC-GBI/",
       width = 9, height = 5)

# Analysis k_multi from TL ----
k_TL = lapply(input_asymm_TL, function(I) sapply(I, function(x) x$k)) %>% abind::abind(along = 3)
dimnames(k_TL) = list("side" = c("neg", "pos"), 
                      "replicate" = 1:100, 
                      "sim_setting" = 1:27)

df_k_TL = k_TL %>% reshape2::melt() %>% 
  left_join(cbind.data.frame(sim_setting = 1:27, sim_settings),
            by = "sim_setting") %>% 
  mutate(k = factor(k, levels = c("0.5", "init", "final")))

plt_k_TL_Gauss05 = df_k_TL %>%  filter(err_distr == "Gaussian_0.5") %>% 
  ggplot() + 
  geom_boxplot(aes(x = side, y = value)) +
  facet_grid(k ~ unif_cov, scales = "free_y") +
  labs(title = "k values with Tower Law - Gaussian(0, sd = 0.5) error distribution")

plt_k_TL_Gauss1 = df_k_TL %>%  filter(err_distr == "Gaussian") %>% 
  ggplot() + 
  geom_boxplot(aes(x = side, y = value)) +
  facet_grid(k ~ unif_cov, scales = "free_y") +
  labs(title = "k values with Tower Law - Gaussian(0, sd = 1) error distribution")

plt_k_TL_Gamma = df_k_TL %>%  filter(err_distr == "Gamma_2_2") %>% 
  ggplot() + 
  geom_boxplot(aes(x = side, y = value)) +
  facet_grid(k ~ unif_cov, scales = "free_y") +
  labs(title = "k values with Tower Law - Gamma(2, 2) error distribution")

ggsave(plot = gridExtra::marrangeGrob(grobs = list(plt_k_TL_Gauss1, plt_k_TL_Gauss05, plt_k_TL_Gamma), 
                                      nrow = 1, ncol = 1, top = NULL), 
       filename = "k_TL.pdf", path = "sim_study_nonC-GBI/",
       width = 6, height = 6)


# Analysis k_multi from sd ----
k_sd = lapply(input_asymm_sd, function(I) sapply(I, function(x) x$k)) %>% abind::abind(along = 3)
dimnames(k_sd) = list("side" = c("neg", "pos"), 
                      "replicate" = 1:100, 
                      "sim_setting" = 1:27)

df_k_sd = k_sd %>% reshape2::melt() %>% 
  left_join(cbind.data.frame(sim_setting = 1:27, sim_settings),
            by = "sim_setting") %>% 
  mutate(k = factor(k, levels = c("0.5", "init", "final")))

plt_k_sd_Gauss05 = df_k_sd %>%  filter(err_distr == "Gaussian_0.5") %>% 
  ggplot() + 
  geom_boxplot(aes(x = side, y = value)) +
  facet_grid(k ~ unif_cov, scales = "free_y") +
  labs(title = "k values with Tower Law - Gaussian(0, sd = 0.5) error distribution")

plt_k_sd_Gauss1 = df_k_sd %>%  filter(err_distr == "Gaussian") %>% 
  ggplot() + 
  geom_boxplot(aes(x = side, y = value)) +
  facet_grid(k ~ unif_cov, scales = "free_y") +
  labs(title = "k values with Tower Law - Gaussian(0, sd = 1) error distribution")

plt_k_sd_Gamma = df_k_sd %>%  filter(err_distr == "Gamma_2_2") %>% 
  ggplot() + 
  geom_boxplot(aes(x = side, y = value)) +
  facet_grid(k ~ unif_cov, scales = "free_y") +
  labs(title = "k values with Tower Law - Gamma(2, 2) error distribution")

ggsave(plot = gridExtra::marrangeGrob(grobs = list(plt_k_TL_Gauss1, plt_k_TL_Gauss05, plt_k_TL_Gamma), 
                                      nrow = 1, ncol = 1, top = NULL), 
       filename = "k_TL.pdf", path = "sim_study_nonC-GBI/",
       width = 6, height = 6)



# Small simulation to get k_symm ----
k_symm = matrix(NA, nrow = 100, ncol = 27)
dimnames(k_symm) = list("replicate" = 1:100, "sim_setting" = 1:27)
nrep = 100
for (s in 1:nrow(sim_settings)){
  
  s_s = sim_settings[s, ]
  ed_s = unname(s_s[1])
  cd_s = unname(s_s[2])
  k_s = unname(s_s[3])
  
  res_tmp = sapply(1:nrep, function(x) {
    
    # DEFINE ERROR SAMPLER
    rerr = switch(ed_s,
                  Gaussian = function(ndraws) rnorm(ndraws, 0, 1),
                  Gaussian_0.5 = function(ndraws) rnorm(ndraws, 0, 0.5),
                  Gamma_2_2 = function(ndraws) rgamma(ndraws, 2, 2),
                  Gamma_1.5_1.5 = function(ndraws) rgamma(ndraws, 1.5, 1.5),
                  Gamma_1.5_3 = function(ndraws) rgamma(ndraws, 1.5, 3))
    
    mode_shift = switch(ed_s,
                        Gaussian = 0,
                        Gaussian_0.5 = 0,
                        Gamma_2_2 = (2-1)/2,
                        Gamma_1.5_1.5 = (1.5-1)/1.5,
                        Gamma_1.5_3 = (1.5-1)/3)
    
    rcov = switch(cd_s,
                  "0.5" = function(ndraws) runif(ndraws, -0.5, 0.5),
                  "1" = function(ndraws) runif(ndraws, -1, 1),
                  asymm = function(ndraws) runif(ndraws, 0, 2))
    
    # DRAW COVARIATES
    x1 = rcov(n)
    x2 = rcov(n)
    X = cbind(1, x1, x2)
    mu_x1 = mean(x1); sd_x1 = sd(x1)
    mu_x2 = mean(x2); sd_x2 = sd(x2)
    X_scaled = cbind(1, (x1-mu_x1)/sd_x1, (x2-mu_x2)/sd_x2)
    eta = X%*%true_beta
    
    # DRAW RESPONSE
    y = eta + rerr(n) - mode_shift
    
    
    # OPTIMAL K
    const = (n^(5/12)) * (log(n)^(7/12)) / (n^(7/12))
    kinit = (n/const)^(1/7)
    res_qr = quantreg::rq(y ~ -1 + X_scaled, tau = 0.5)$residuals
    kfinal = kinit * (n^(1/7))
    k_symm = switch(k_s,
                    "0.5" = 0.5,
                    "init" = kinit,
                    "final" = kfinal/sd(res_qr))
    
    
    k_symm
  })
  
  k_symm[ , s] = res_tmp
}


## Plot symm ----
df_k_symm = cbind.data.frame(t(k_symm), sim_settings) %>% 
  pivot_longer(1:100, names_to = "replicate", values_to = "value") %>% 
  mutate(k = factor(k, levels = c("0.5", "init", "final")),
         err_distr = factor(err_distr, levels = c("Gaussian", "Gaussian_0.5", "Gamma_2_2")),
         replicate = as.integer(replicate),
         side = "symm")

plt_k_symm = df_k_symm %>%
  ggplot() + 
  geom_boxplot(aes(x = unif_cov, y = value)) +
  facet_grid(k ~ err_distr, scales = "free_y")

ggsave(plot = plt_k_symm, 
       filename = "k_symm.pdf", path = "sim_study_nonC-GBI/",
       width = 6, height = 6)


## Plot together ----
df_k_complete = bind_rows(df_k_TL %>% mutate(type = "TL"), df_k_asymm %>% mutate(type = "asymm"), df_k_symm %>% mutate(type = "symm")) %>% 
  mutate(type = factor(type, levels = c("symm", "asymm", "TL"), labels = c("symm", "asymm sd", "asymm TL")))

plt_k_Gauss05 = df_k_complete %>% 
  filter(err_distr == "Gaussian_0.5") %>%
  ggplot() + 
  geom_boxplot(aes(x = side, y = value, col = type)) +
  facet_grid(k ~ unif_cov, scales = "free_y") +
  labs(title = "k values across replicates - Gaussian(0, sd = 0.5) error distribution")

plt_k_Gauss1 = df_k_complete %>% 
  filter(err_distr == "Gaussian") %>%
  ggplot() + 
  geom_boxplot(aes(x = side, y = value, col = type)) +
  facet_grid(k ~ unif_cov, scales = "free_y") +
  labs(title = "k values across replicates - Gaussian(0, sd = 1) error distribution")

plt_k_Gamma = df_k_complete %>% 
  filter(err_distr == "Gamma_2_2") %>%
  ggplot() + 
  geom_boxplot(aes(x = side, y = value, col = type)) +
  facet_grid(k ~ unif_cov, scales = "free_y") +
  labs(title = "k values across replicates - Gamma(2, 2) error distribution")

ggsave(plot = gridExtra::marrangeGrob(grobs = list(plt_k_Gauss1, plt_k_Gauss05, plt_k_Gamma), 
                                      nrow = 1, ncol = 1, top = NULL), 
       filename = "k_together.pdf", path = "sim_study_nonC-GBI/",
       width = 9, height = 6)


# Predictor ----
err_pred_sd = lapply(input_asymm_sd, function(I1) sapply(I1, function(I2) {
  b = I2$draws[idx_iter, ]; X_scaled = I2$X_scaled; X = I2$X
  colMeans(b %*% t(X_scaled)) - true_beta %*% t(X)
})) %>% abind::abind(along = 3) %>% aperm(c(3, 2, 1))
dimnames(err_pred_sd) = list("sim_setting" = 1:27, "replicate" = 1:100, "unit" = 1:1000)

err_pred_TL = lapply(input_asymm_TL, function(I1) sapply(I1, function(I2) {
  b = I2$draws[idx_iter, ]; X_scaled = I2$X_scaled; X = I2$X
  colMeans(b %*% t(X_scaled)) - true_beta %*% t(X)
})) %>% abind::abind(along = 3) %>% aperm(c(3, 2, 1))
dimnames(err_pred_TL) = list("sim_setting" = 1:27, "replicate" = 1:100, "unit" = 1:1000)


dens_err_pred_sd = apply(err_pred_sd, 1, function(E1) apply(E1, 1, function(E) {
  tmp = density(E, from = -1.5, to = 1.5, n = 100)
  cbind(tmp$x, tmp$y)
}, simplify = F) %>% abind::abind(along = 3), simplify = F) %>% 
  abind::abind(along = 4)  %>% aperm(c(4, 3, 2, 1))

dens_err_pred_TL = apply(err_pred_TL, 1, function(E1) apply(E1, 1, function(E) {
  tmp = density(E, from = -1.5, to = 1.5, n = 100)
  cbind(tmp$x, tmp$y)
}, simplify = F) %>% abind::abind(along = 3), simplify = F) %>% 
  abind::abind(along = 4)  %>% aperm(c(4, 3, 2, 1))

dimnames(dens_err_pred_sd) = dimnames(dens_err_pred_TL) =
  list("setting" = 1:27, "replicate" = 1:100, "axis" = c("x", "y"), "knot" = 1:100) 

df_dens_err_pred = left_join(reshape2::melt(dens_err_pred_sd),
                       cbind.data.frame(setting = 1:27, sim_settings),
                       by = "setting") %>% 
  pivot_wider(names_from = axis, values_from = value) %>%
  mutate(type = "sd") %>% 
  bind_rows(left_join(reshape2::melt(dens_err_pred_TL),
                      cbind.data.frame(setting = 1:27, sim_settings),
                      by = "setting") %>% 
              pivot_wider(names_from = axis, values_from = value) %>% 
              mutate(type = "TL")) %>% 
  group_by(setting, knot, type, .drop = F) %>% 
  mutate(mean_y = mean(y),
         k = factor(k, levels = c("0.5", "init", "final")))


plt_dens_err_pred_Gauss05 = df_dens_err_pred %>%
  filter(err_distr == "Gaussian_0.5") %>% 
  ggplot(aes(x = x)) +
  geom_vline(aes(xintercept = 0), col = "gold") + 
  geom_line(aes(y = y, group = interaction(replicate, type), col = type), alpha = 0.1) +
  ggnewscale::new_scale_color() +
  geom_line(aes(y = mean_y, col = type)) + 
  scale_colour_manual(values = c("sd" = "red3", "TL" = "blue3")) +
  scale_y_continuous(limits = c(0, 5)) +
  ggh4x::facet_nested(k ~ unif_cov) +
  labs(title = "Error on linear predictor [light = per replicate, dark = average] - Gaussian(0, sd = 0.5) error distribution", 
       y = "density", x = "error linear predictor")

plt_dens_err_pred_Gauss1 = df_dens_err_pred %>%
  filter(err_distr == "Gaussian") %>% 
  ggplot(aes(x = x)) +
  geom_vline(aes(xintercept = 0), col = "gold") + 
  geom_line(aes(y = y, group = interaction(replicate, type), col = type), alpha = 0.1) +
  ggnewscale::new_scale_color() +
  geom_line(aes(y = mean_y, col = type)) + 
  scale_colour_manual(values = c("sd" = "red3", "TL" = "blue3")) +
  scale_y_continuous(limits = c(0, 5)) +
  ggh4x::facet_nested(k ~ unif_cov) +
  labs(title = "Error on linear predictor [light = per replicate, dark = average] - Gaussian(0, sd = 1) error distribution", 
       y = "density", x = "error linear predictor")

plt_dens_err_pred_Gamma = df_dens_err_pred %>%
  filter(err_distr == "Gamma_2_2") %>% 
  ggplot(aes(x = x)) +
  geom_vline(aes(xintercept = 0), col = "gold") + 
  geom_line(aes(y = y, group = interaction(replicate, type), col = type), alpha = 0.1) +
  ggnewscale::new_scale_color() +
  geom_line(aes(y = mean_y, col = type)) + 
  scale_colour_manual(values = c("sd" = "red3", "TL" = "blue3")) +
  scale_y_continuous(limits = c(0, 5)) +
  ggh4x::facet_nested(k ~ unif_cov) +
  labs(title = "Error on linear predictor [light = per replicate, dark = average] - Gamma(2, 2) error distribution", 
       y = "density", x = "error linear predictor")

ggsave(plot = gridExtra::marrangeGrob(grobs = list(plt_dens_err_pred_Gauss1, plt_dens_err_pred_Gauss05, plt_dens_err_pred_Gamma),
                                      nrow = 1, ncol = 1, top = NULL), 
       filename = "err_lin_pred.pdf", path = "sim_study_nonC-GBI/",
       width = 9, height = 6)


# Residuals ----
residuals_sd = lapply(input_asymm_sd, function(I1) sapply(I1, function(I2) {
  b = I2$draws[idx_iter, ]; X_scaled = I2$X_scaled; y = I2$y
  y - colMeans(b %*% t(X_scaled))
})) %>% abind::abind(along = 3) %>% aperm(c(3, 2, 1))
dimnames(err_pred_sd) = list("sim_setting" = 1:27, "replicate" = 1:100, "unit" = 1:1000)

residuals_TL = lapply(input_asymm_TL, function(I1) sapply(I1, function(I2) {
  b = I2$draws[idx_iter, ]; X_scaled = I2$X_scaled; y = I2$y
  y - colMeans(b %*% t(X_scaled))
})) %>% abind::abind(along = 3) %>% aperm(c(3, 2, 1))
dimnames(err_pred_TL) = list("sim_setting" = 1:27, "replicate" = 1:100, "unit" = 1:1000)


dens_residuals_sd = apply(residuals_sd, 1, function(R1) apply(R1, 1, function(R) {
  tmp = density(R, from = -3, to = 3, n = 100)
  cbind(tmp$x, tmp$y)
}, simplify = F) %>% abind::abind(along = 3), simplify = F) %>% 
  abind::abind(along = 4)  %>% aperm(c(4, 3, 2, 1))

dens_residuals_TL = apply(residuals_TL, 1, function(R1) apply(R1, 1, function(R) {
  tmp = density(R, from = -3, to = 3, n = 100)
  cbind(tmp$x, tmp$y)
}, simplify = F) %>% abind::abind(along = 3), simplify = F) %>% 
  abind::abind(along = 4)  %>% aperm(c(4, 3, 2, 1))

dimnames(dens_residuals_sd) = dimnames(dens_residuals_TL) =
  list("setting" = 1:27, "replicate" = 1:100, "axis" = c("x", "y"), "knot" = 1:100) 

df_residuals = left_join(reshape2::melt(dens_residuals_sd),
                             cbind.data.frame(setting = 1:27, sim_settings),
                             by = "setting") %>% 
  pivot_wider(names_from = axis, values_from = value) %>%
  mutate(type = "sd") %>% 
  bind_rows(left_join(reshape2::melt(dens_residuals_TL),
                      cbind.data.frame(setting = 1:27, sim_settings),
                      by = "setting") %>% 
              pivot_wider(names_from = axis, values_from = value) %>% 
              mutate(type = "TL")) %>% 
  group_by(setting, knot, type, .drop = F) %>% 
  mutate(mean_y = mean(y),
         k = factor(k, levels = c("0.5", "init", "final")))


plt_dens_residuals_Gauss05 = df_residuals %>%
  filter(err_distr == "Gaussian_0.5") %>% 
  ggplot(aes(x = x)) +
  geom_line(aes(y = y, group = interaction(replicate, type), col = type), alpha = 0.1) +
  geom_line(data = data.frame(x = seq(-3, 3, length = 101),
                              y = dnorm(seq(-3, 3, length = 101), sd = 0.5)),
            aes(x = x, y = y),
            col = "gold", lty = 1) +
  ggnewscale::new_scale_color() +
  geom_line(aes(y = mean_y, col = type)) + 
  scale_colour_manual(values = c("sd" = "red3", "TL" = "blue3")) +
  ggh4x::facet_nested(k ~ unif_cov) +
  labs(title = "Residuals [light = per replicate, dark = average] - Gaussian(0, sd = 0.5) error distribution", 
       y = "density", x = "error linear predictor")

plt_dens_residuals_Gauss1 = df_residuals %>%
  filter(err_distr == "Gaussian") %>% 
  ggplot(aes(x = x)) +
  geom_line(aes(y = y, group = interaction(replicate, type), col = type), alpha = 0.1) +
  geom_line(data = data.frame(x = seq(-3, 3, length = 101),
                              y = dnorm(seq(-3, 3, length = 101), sd = 1)),
            aes(x = x, y = y),
            col = "gold2", lty = 1) +
  ggnewscale::new_scale_color() +
  geom_line(aes(y = mean_y, col = type)) + 
  scale_colour_manual(values = c("sd" = "red3", "TL" = "blue3")) +
  ggh4x::facet_nested(k ~ unif_cov) +
  labs(title = "Residuals [light = per replicate, dark = average] - Gaussian(0, sd = 1) error distribution", 
       y = "density", x = "error linear predictor")

plt_dens_residuals_Gamma = df_residuals %>%
  filter(err_distr == "Gamma_2_2") %>% 
  ggplot(aes(x = x)) +
  geom_line(aes(y = y, group = interaction(replicate, type), col = type), alpha = 0.1) +
  geom_line(data = data.frame(x = seq(-3, 3, length = 101),
                              y = dgamma(seq(-3, 3, length = 101) + 0.5, shape = 2, rate = 2)),
            aes(x = x, y = y),
            col = "gold", lty = 1) +
  ggnewscale::new_scale_color() +
  geom_line(aes(y = mean_y, col = type)) + 
  scale_colour_manual(values = c("sd" = "red3", "TL" = "blue3")) +
  ggh4x::facet_nested(k ~ unif_cov) +
  labs(title = "Residuals [light = per replicate, dark = average] - Gamma(2, 2) error distribution", 
       y = "density", x = "error linear predictor")

ggsave(plot = gridExtra::marrangeGrob(grobs = list(plt_dens_residuals_Gauss1, plt_dens_residuals_Gauss05, plt_dens_residuals_Gamma),
                                      nrow = 1, ncol = 1, top = NULL), 
       filename = "residuals_dens.pdf", path = "sim_study_nonC-GBI/",
       width = 9, height = 6)

# Diagnostics ----
diagn_TL = lapply(input_asymm_TL, function(I1) 
  lapply(I1, function(I2) I2$diagn) %>% abind::abind(along = 3)) %>% 
  abind::abind(along = 4)
diagn_sd = lapply(input_asymm_sd, function(I1)
  lapply(I1, function(I2) I2$diagn) %>% abind::abind(along = 3)) %>% 
  abind::abind(along = 4)

dimnames(diagn_sd) = 
  dimnames(diagn_TL) =
  list("index" = dimnames(diagn_sd)[[1]], 
                       "param" = dimnames(diagn_sd)[[2]],
                       "replicate" = 1:100,
                       "sim_setting" = 1:27)

df_diagn_sd = diagn_sd %>% reshape2::melt() %>% mutate(type = "sd") %>% 
  left_join(cbind.data.frame(sim_setting = 1:27, sim_settings),
            by = "sim_setting") %>% 
  separate(index, into = c("index", "thin")) %>% 
  pivot_wider(names_from = index, values_from = value) %>% 
  mutate(k = factor(k, levels = c("0.5", "init", "final")),
         err_distr = factor(err_distr, levels = c("Gaussian", "Gaussian_0.5", "Gamma_2_2")))

df_diagn_TL = diagn_TL %>% reshape2::melt() %>% mutate(type = "TL") %>% 
  left_join(cbind.data.frame(sim_setting = 1:27, sim_settings),
            by = "sim_setting") %>% 
  separate(index, into = c("index", "thin")) %>% 
  pivot_wider(names_from = index, values_from = value) %>% 
  mutate(k = factor(k, levels = c("0.5", "init", "final")),
         err_distr = factor(err_distr, levels = c("Gaussian", "Gaussian_0.5", "Gamma_2_2")))

df_diagn = bind_rows(df_diagn_sd, df_diagn_TL)

df_failed = df_diagn %>% 
  group_by(type, sim_setting, err_distr, unif_cov, k, thin) %>% 
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
