rm(list = ls())
source("src/01-utils.R")

# TRUE BETA AND SAMPLE SIZE ----
n = 1e3
true_beta = c(0, -0.5, +0.7)

# SIMULATION SETTINGS ----
error_distr = c("Gaussian", "Gamma_2_rad2", "Gaussian_rad2", "Gamma_2_1")
covariate = c("Unif_rad2", "Gamma_1.5_1.5", "BVN")
k = "final"
sim_settings = as.matrix(expand.grid(error_distr, covariate, k))
colnames(sim_settings) = c("err_distr", "unif_cov", "k")

results = vector("list", nrow(sim_settings))
niter = 11e3
nrep = 100

# PARALLEL BACKEND ----
library(doParallel); library(foreach)
ncores = 25 #max(1, parallel::detectCores() - 1)
registerDoParallel(ncores)

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
                    Gaussian_rad2 = function(ndraws) rnorm(ndraws, 0, sqrt(2)),
                    Gamma_2_2 = function(ndraws) rgamma(ndraws, 2, 2),
                    Gamma_2_rad2 = function(ndraws) rgamma(ndraws, 2, sqrt(2)),
                    Gamma_2_1 = function(ndraws) rgamma(ndraws, 2, 1))
      
      mode_shift = switch(ed_s,
                          Gaussian = 0,
                          Gaussian_0.5 = 0,
                          Gaussian_rad2 = 0,
                          Gamma_2_2 = (2-1)/2,
                          Gamma_1.5_1.5 = (1.5-1)/1.5,
                          Gamma_2_1 = (2-1)/1,
                          Gamma_2_rad2 = (2-1)/sqrt(2))
      
      rcov = switch(cd_s,
                    "0.5" = function(ndraws) runif(ndraws, -0.5, 0.5),
                    "1" = function(ndraws) runif(ndraws, -1, 1),
                    asymm = function(ndraws) runif(ndraws, 0, 2),
                    "Unif_rad2" = function(ndraws) runif(ndraws, -sqrt(2), sqrt(2)),
                    "Gamma_1.5_1.5" = function(ndraws) rgamma(ndraws, 1.5, 1.5),
                    "BVN" = function(ndraws) mvnfast::rmvn(ndraws, rep(0, 2), sigma = (26/12)*matrix(c(1, 0.75, 0.75, 1), nrow = 2), isChol = F))
      
      # DRAW COVARIATES
      if (cd_s != "BVN"){
        x1 = rcov(n)
        x2 = rcov(n)
        X = cbind(1, x1, x2)
      } else {
        X.tmp = rcov(n)
        x1 = X.tmp[ , 1]
        x2 = X.tmp[ , 2]
        X = cbind(1, X.tmp)
      }
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
      k = switch(k_s,
                 "0.5" = 0.5,
                 "init" = kinit,
                 "final" = kfinal)
      
      # I use 2*(..) so that if the err distribution is symmetric, then
      # factorp = factorm = sd(res) and we have the original k
      varp = var(res_qr[res_qr > 0])
      varn = var(res_qr[res_qr < 0])
      mup = mean(res_qr[res_qr > 0]) - mean(res_qr)
      mun = mean(res_qr[res_qr < 0]) - mean(res_qr)
      factorp = sqrt(varp+mup^2)
      factorn = sqrt(varn+mun^2)
      
      k_multi = c( 
        k/factorn,
        k/factorp
      )
      
      out_optim = sapply(c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"),
                         function(m){
                           est = optim(par = c(0.1, -0.25, +0.4),
                                       fn = function(b) sum(loss_asymm2(y - X_scaled%*%b, 
                                                                        k_multi, 
                                                                        1e-3)),
                                       method = m)$par
                           b0 = est[1]; b1 = est[2]; b2 = est[3]
                           c(b0-b1*mu_x1/sd_x1-b2*mu_x2/sd_x2, b1/sd_x1, b2/sd_x2)
                         }
                         
      )
      
      out_MCMC = MCMC_MH(niter = niter, X = X_scaled, y = y, 
                         k = k_multi, 
                         c = 1e-3, asymm = T,
                         p = 3,  beta0 = c(0.1, -0.25, +0.4), prior_sd = 1/3)
      
      draws_tmp = out_MCMC$beta[-(1:1000), ]
      draws = draws_tmp
      draws[ , 1] = draws_tmp[ , 1] - draws_tmp[ , 2]*mu_x1/sd_x1 - draws_tmp[ , 3]*mu_x2/sd_x2
      draws[ , 2] = draws_tmp[ , 2]/sd_x1
      draws[ , 3] = draws_tmp[ , 3]/sd_x2
      
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
      colnames(diagn) = paste0("beta", 1:p)
      
      list(draws = draws, summ = summ, diagn = diagn, w = out_MCMC$w, k = k_multi,
           X = X, X_scaled = X_scaled, y = y)
      
    }
  
  res = res_tmp
  results[[s]] = res
  cat("Update: ", format(Sys.time()), "\t", sep = "")
  cat(round(100*s/nrow(sim_settings), 2), "%\n", sep = "")
}
saveRDS(results, "sim_study_nonC-GBI/02-w_unit_loss/sim_study_nonCalibrated_multiv_asymm_TL.RDS")
