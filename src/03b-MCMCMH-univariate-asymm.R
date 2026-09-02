rm(list = ls())
source("src/01-utils.R")

# TRUE BETA AND SAMPLE SIZE ----
n = 1e3
true_beta = c(0, -0.5)

# SIMULATION SETTINGS ----
error_distr = c("Gaussian", "Gaussian_0.5", "Gamma_2_2")
covariate = c("0.5", "1", "asymm")
k = c("0.5", "init", "final")
sim_settings = as.matrix(expand.grid(error_distr, covariate, k))
colnames(sim_settings) = c("err_distr", "unif_cov", "k")

results = vector("list", nrow(sim_settings))
niter = 11e3
nrep = 100

# PARALLEL BACKEND ----
library(doParallel); library(foreach)
ncores = 8 #max(1, parallel::detectCores() - 1)
cl = makeForkCluster(ncores)
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
                    Gamma_2_2 = function(ndraws) rgamma(ndraws, 2, rate = 2),
                    Gamma_1.5_1.5 = function(ndraws) rgamma(ndraws, 1.5, rate = 1.5),
                    Gamma_1.5_3 = function(ndraws) rgamma(ndraws, 1.5, rate = 3))
      
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
      x = rcov(n)
      X = cbind(1, x)
      mu_x = mean(x); sd_x = sd(x)
      X_scaled = cbind(1, (x-mu_x)/sd_x)
      eta = X%*%true_beta
      
      # DRAW RESPONSE
      y = eta + ( rerr(n) - mode_shift)
      
      # OPTIMAL K
      const = (n^(5/12)) * (log(n)^(7/12)) / (n^(7/12))
      kinit = (n/const)^(1/7)
      res_qr = quantreg::rq(y ~ -1 + X_scaled, tau = 0.5)$residuals
      sd_quant = sd(res_qr)
      kfinal = kinit * (n^(1/7))
      k = switch(k_s,
                 "0.5" = 0.5,
                 "init" = kinit,
                 "final" = kfinal) 
      
      out_optim = sapply(c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"), 
                         function(m){
                           est = optim(par = c(0.05, -0.45), 
                                       fn = function(b) sum(loss_asymm2(y - X_scaled%*%b, 
                                                                        c(k/sd(res_qr[res_qr < 0]),
                                                                          k/sd(res_qr[res_qr > 0])), 
                                                                        1e-3)),
                                       method = m)$par
                           b0 = est[1]; b1 = est[2]
                           c(b0-b1*mu_x/sd_x, b1/sd_x)
                         }
                         
      )
      
      out_MCMC = MCMC_MH(niter = niter, X = X_scaled, y = y, 
                         k = c(kinit/sd(res_qr[res_qr < 0]), kinit/sd(res_qr[res_qr > 0])), 
                         c = 1e-3, asymm = T,
                         p = 2,  beta0 = c(0.1, -0.25), prior_sd = 1/3)
      
      draws_tmp = out_MCMC$beta[-(1:1000), ]
      draws = draws_tmp
      draws[ , 1] = draws_tmp[ , 1] - draws_tmp[ , 2]*mu_x/sd_x
      draws[ , 2] = draws_tmp[ , 2]/sd_x
      
      draws_thin5 = draws[seq(1, 1e4, by = 5), ]
      draws_thin10 = draws[seq(1, 1e4, by = 10), ]
      
      summ = rbind(apply(draws_thin10, 2,
                         function(D){
                           mean_b = mean(D)
                           med_b = median(D)
                           sd_b = sd(D)
                           q025 = unname(quantile(D, probs = c(0.025)))
                           q975 = unname(quantile(D, probs = c(0.975)))
                           hdi_b = HDInterval::hdi(density(D), credMass = 0.95, allowSplit = T)
                           if (nrow(hdi_b) == 1){
                             HDI_L1 = unname(hdi_b[1, "begin"]); HDI_U1 = unname(hdi_b[1, "end"])
                             HDI_L2 = HDI_U2 = NA
                           } else {
                             HDI_L1 = unname(hdi_b[1, "begin"]); HDI_U1 = unname(hdi_b[1, "end"])
                             HDI_L2 = unname(hdi_b[2, "begin"]); HDI_U2 = unname(hdi_b[2, "end"])
                           }
                           
                           c(mean = mean_b, median = med_b, sd = sd_b, 
                             q025 = q025, q975 = q975,
                             HDI_L1 = HDI_L1, HDI_U1 = HDI_U1, HDI_L2 = HDI_L2, HDI_U2 = HDI_U2)
                         }
      ), t(out_optim))
      
      diagn = rbind(apply(draws_thin10, 2, function(x) acf(x, plot = F)[[1]][2]),
                    apply(draws_thin10, 2, posterior::rhat),
                    apply(draws_thin10, 2, posterior::ess_basic),
                    apply(draws_thin5, 2, function(x) acf(x, plot = F)[[1]][2]),
                    apply(draws_thin5, 2, posterior::rhat),
                    apply(draws_thin5, 2, posterior::ess_basic),
                    apply(draws, 2, function(x) acf(x, plot = F)[[1]][2]),
                    apply(draws, 2, posterior::rhat),
                    apply(draws, 2, posterior::ess_basic))
      
      p = length(true_beta)
      rownames(diagn) = c(outer(c("acf_thin", "rhat_thin", "ess_thin"), c(10, 5, 1), paste0))
      
      rbind(summ, diagn)
      
    }
  
  res = simplify2array(res_tmp)
  results[[s]] = res
  cat("Update: ", format(Sys.time()), "\t", sep = "")
  cat(round(100*s/nrow(sim_settings), 2), "%\n", sep = "")
}
stopCluster(cl)
saveRDS(results, "sim_study_nonCalibrated_asymm.RDS")


# # Diagnostics ----
# # Rhat
# sapply(results, function(X) summary(X["rhat_thin10", 2, ]))
# 
# sapply(results, function(X) summary(X["rhat_thin1", 2, ]))
# sapply(results, function(X) summary(X["rhat_thin10", 1, ]))
# sapply(results, function(X) summary(X["rhat_thin1", 1, ]))
# # ESS
# cbind(sim_settings, round(t(sapply(results, function(X) summary(X["ess_thin10", 2, ]))), 2))
# cbind(sim_settings, round(t(sapply(results, function(X) summary(X["ess_thin1", 2, ]))), 2))
# cbind(sim_settings, round(t(sapply(results, function(X) summary(X["ess_thin10", 1, ]))), 2))
# cbind(sim_settings, round(t(sapply(results, function(X) summary(X["ess_thin1", 1, ]))), 2))
# # ACF
# cbind(sim_settings, round(t(sapply(results, function(X) summary(X["acf_thin10", 2, ]))), 2))
# cbind(sim_settings, round(t(sapply(results, function(X) summary(X["acf_thin1", 2, ]))), 2))
# cbind(sim_settings, round(t(sapply(results, function(X) summary(X["acf_thin10", 1, ]))), 2))
# cbind(sim_settings, round(t(sapply(results, function(X) summary(X["acf_thin1", 1, ]))), 2))