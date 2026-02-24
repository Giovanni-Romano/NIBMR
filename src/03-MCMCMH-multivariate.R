rm(list = ls())
source("src/01-utils.R")

# TRUE BETA AND SAMPLE SIZE ----
n = 1e3
true_beta = c(0, -0.5, +0.7)

# SIMULATION SETTINGS ----
error_distr = c("Gaussian", "Gaussian_0.5", "Gamma_2_2")
covariate = c("1", "asymm")
k = c("0.5", "init", "final")
sim_settings = as.matrix(expand.grid(error_distr, covariate, k))
colnames(sim_settings) = c("err_distr", "unif_cov", "k")

results = vector("list", nrow(sim_settings))
niter = 5e3
nrep = 100

# PARALLEL BACKEND ----
library(doParallel); library(foreach)
ncores = 5 #max(1, parallel::detectCores() - 1)
cl = makeCluster(ncores)
registerDoParallel(cl)

cat("Starting time: ", format(Sys.time()), "\n", sep = "")

for (s in 1:nrow(sim_settings)){
  
  s_s = sim_settings[s, ]
  ed_s = unname(s_s[1])
  cd_s = unname(s_s[2])
  k_s = unname(s_s[3])
  
  res_tmp = foreach(
    rep      = 1:nrep,
    .packages = c("quantreg", "HDInterval")) %dopar% {
      
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
      sd_quant = sd(quantreg::rq(y ~ -1 + X_scaled, tau = 0.5)$residuals)
      kfinal = kinit * (n^(1/7)) / sd_quant
      k = switch(k_s,
                 "0.5" = 0.5,
                 "init" = kinit,
                 "final" = kfinal) 
      
      out_optim = sapply(c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"), 
                         function(m){
                           est = optim(par = c(0.1, -0.25, +0.5), 
                                       fn = function(b) sum(loss(y - X_scaled%*%b, k, 1e-3)),
                                       method = m)$par
                           b0 = est[1]; b1 = est[2]; b2 = est[3]
                           c(b0-b1*mu_x1/sd_x1-b2*mu_x2/sd_x2, b1/sd_x1, b2/sd_x2)
                         }
                         
      )
      
      out_MCMC = MCMC_MH(niter = niter, X = X_scaled, y = y, k = k, c = 1e-3, 
                         p = 2,  beta0 = c(0.1, -0.25), prior_sd = 1/3)
      
      draws_tmp = out_MCMC$beta[-(1:1000), ]
      draws = draws_tmp
      draws[ , 1] = draws_tmp[ , 1] - draws_tmp[ , 2]*mu_x/sd_x
      draws[ , 2] = draws_tmp[ , 2]/sd_x
      
      rbind(apply(draws, 2,
                  function(D){
                    mean_b = mean(D)
                    med_b = median(D)
                    sd_b = sd(D)
                    q025 = unname(quantile(D, probs = c(0.025)))
                    q0975 = unname(quantile(D, probs = c(0.975)))
                    hdi_b = HDInterval::hdi(density(D), credMass = 0.95, allowSplit = T)
                    if (nrow(hdi_b) == 1){
                      HDI_L1 = unname(hdi_b[1, "begin"]); HDI_U1 = unname(hdi_b[1, "end"])
                      HDI_L2 = HDI_U2 = NA
                    } else {
                      HDI_L1 = unname(hdi_b[1, "begin"]); HDI_U1 = unname(hdi_b[1, "end"])
                      HDI_L2 = unname(hdi_b[2, "begin"]); HDI_U2 = unname(hdi_b[2, "end"])
                    }
                    
                    c(mean = mean_b, median = med_b, sd = sd_b, 
                      q025 = q025, q0975 = q0975,
                      HDI_L1 = HDI_L1, HDI_U1 = HDI_U1, HDI_L2 = HDI_L2, HDI_U2 = HDI_U2)
                  }
      ), t(out_optim))
      
    }
  
  res = simplify2array(res_tmp)
  results[[s]] = res
  cat("Update: ", format(Sys.time()), "\t", sep = "")
  cat(round(100*s/nrow(sim_settings), 2), "%\n", sep = "")
}
stopCluster(cl)
# saveRDS(results, "sim_study.RDS")

par(mfrow = c(3, 3))
for (j in seq(1, 27, by = 3)){
  min_x = -1 #min(c(results[[j]]["mean", 2, ], results[[j]]["BFGS", 2, ])) 
  max_x = 0 #max(c(results[[j]]["mean", 2, ], results[[j]]["BFGS", 2, ])) 
  plot(1, xlim = c(min_x, max_x), ylim = c(0, 10), type = "n", main = paste(sim_settings[j, ], collapse = " - "))
  lines(density(results[[j]]["mean", 2, ]))
  lines(density(results[[j]]["BFGS", 2, ]), col = "blue")
  lines(density(results[[j]]["CG", 2, ]), col = "green", lty = 2)
  abline(v = true_beta[2], col = "gold")
}

for (j in seq(2, 27, by = 3)){
  min_x =-1 #min(c(results[[j]]["mean", 2, ], results[[j]]["BFGS", 2, ])) 
  max_x = 0 #max(c(results[[j]]["mean", 2, ], results[[j]]["BFGS", 2, ])) 
  plot(1, xlim = c(min_x, max_x), ylim = c(0, 10), type = "n", main = paste(sim_settings[j, ], collapse = " - "))
  lines(density(results[[j]]["mean", 2, ]))
  lines(density(results[[j]]["BFGS", 2, ]), col = "blue")
  lines(density(results[[j]]["CG", 2, ]), col = "green", lty = 2)
  abline(v = true_beta[2], col = "gold")
}

for (j in seq(3, 27, by = 3)){
  min_x = -1 #min(c(results[[j]]["mean", 2, ], results[[j]]["BFGS", 2, ])) 
  max_x = 0 #max(c(results[[j]]["mean", 2, ], results[[j]]["BFGS", 2, ])) 
  plot(1, xlim = c(min_x, max_x), ylim = c(0, 10), type = "n", main = paste(sim_settings[j, ], collapse = " - "))
  lines(density(results[[j]]["mean", 2, ]))
  lines(density(results[[j]]["BFGS", 2, ]), col = "blue")
  lines(density(results[[j]]["CG", 2, ]), col = "green", lty = 2)
  abline(v = true_beta[2], col = "gold")
}


par(mfrow = c(3, 3))
for (j in seq(1, 27, by = 3)){
  min_x = -0.25 #min(c(results[[j]]["mean", 1, ], results[[j]]["BFGS", 1, ])) 
  max_x = 0.25 #max(c(results[[j]]["mean", 1, ], results[[j]]["BFGS", 1, ])) 
  plot(1, xlim = c(min_x, max_x), ylim = c(0, 15), type = "n", main = paste(sim_settings[j, ], collapse = " - "))
  lines(density(results[[j]]["mean", 1, ]))
  lines(density(results[[j]]["BFGS", 1, ]), col = "blue")
  lines(density(results[[j]]["CG", 1, ]), col = "green", lty = 2)
  abline(v = true_beta[1], col = "gold")
}

for (j in seq(2, 27, by = 3)){
  min_x = -0.25 #min(c(results[[j]]["mean", 1, ], results[[j]]["BFGS", 1, ])) 
  max_x = 0.25 #max(c(results[[j]]["mean", 1, ], results[[j]]["BFGS", 1, ])) 
  plot(1, xlim = c(min_x, max_x), ylim = c(0, 15), type = "n", main = paste(sim_settings[j, ], collapse = " - "))
  lines(density(results[[j]]["mean", 1, ]))
  lines(density(results[[j]]["BFGS", 1, ]), col = "blue")
  lines(density(results[[j]]["CG", 1, ]), col = "green", lty = 2)
  abline(v = true_beta[1], col = "gold")
}

for (j in seq(3, 27, by = 3)){
  min_x = -0.1 #min(c(results[[j]]["mean", 1, ], results[[j]]["BFGS", 1, ])) 
  max_x = 0.35 #max(c(results[[j]]["mean", 1, ], results[[j]]["BFGS", 1, ])) 
  plot(1, xlim = c(min_x, max_x), ylim = c(0, 15), type = "n", main = paste(sim_settings[j, ], collapse = " - "))
  lines(density(results[[j]]["mean", 1, ]))
  lines(density(results[[j]]["BFGS", 1, ]), col = "blue")
  lines(density(results[[j]]["CG", 1, ]), col = "green", lty = 2)
  abline(v = true_beta[1], col = "gold")
}


mean_mean = sapply(results, function(R) mean(R["mean", 2, ]))
sd_mean = sapply(results, function(R) sd(R["mean", 2, ]))
mean_sd = sapply(results, function(R) mean(R["sd", 2, ]))
coverage_SymmI = sapply(results, function(R) {
  I = R[c("q025", "q0975"), 2, ]
  mean(apply(I, 2, function(i) true_beta[2] >= i[1] & true_beta[2] <= i[2]))
})

coverage_HDI = sapply(results, function(R) {
  I1 = R[c("HDI_L1", "HDI_U1"), 2, ]
  I2 = R[c("HDI_L2", "HDI_U2"), 2, ]
  tmp1 = apply(I1, 2, function(i) true_beta[2] >= i[1] & true_beta[2] <= i[2])
  tmp2 = apply(I2, 2, function(i) true_beta[2] >= i[1] & true_beta[2] <= i[2])
  tmp2[is.na(tmp2)] = FALSE
  mean(tmp1 | tmp2)
})


library(tidyverse)
library(ggplot2)

df_res_beta1 = cbind(sim_settings, 
      mean_mean = round(mean_mean, 2), sd_mean = round(sd_mean, 2), 
      mean_sd = round(mean_sd, 2),
      coverage_SymmI,
      coverage_HDI) %>% 
  as_tibble() %>% 
  mutate(across(-(1:3), ~ as.numeric(as.character(.))))

df_res_beta1 %>% 
  ggplot() + 
  geom_point(aes(x = k, y = mean_mean)) +
  geom_errorbar(aes(x = k, ymin = mean_mean - sd_mean, ymax = mean_mean + sd_mean)) +
  geom_hline(aes(yintercept = true_beta[2]), col = "gold", lty = 2) +
  ggh4x::facet_nested(err_distr ~ unif_cov, scales = "free") +
  theme_bw()


df_res_beta1 %>% 
  select(err_distr:k, coverage_SymmI, coverage_HDI) %>%
  pivot_longer(coverage_SymmI:coverage_HDI, names_to = "Type_Int") %>% 
  mutate(Type_Int = str_split_i(Type_Int, "_", 2)) %>% 
  ggplot() + 
  geom_point(aes(x = k, y = value, col = Type_Int, shape = Type_Int)) +
  geom_hline(aes(yintercept = 0.95), col = "gold", lty = 2) +
  ggh4x::facet_nested(err_distr ~ unif_cov, scales = "free") +
  theme_bw()
