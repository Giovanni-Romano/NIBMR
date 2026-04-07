rm(list = ls())
source("src/01-utils.R")

# TRUE BETA AND SAMPLE SIZE ----
n = 1e3

# SIMULATION SETTINGS ----
error_distr = c("Gaussian", "Gamma")
cov_distr = c("Unif_rad2", "Beta_2_4") #"Gamma_1.5_1.5") #, "BVN")
transform_x = c("parabola", "cubic", "trigonometric")
SNR = c(2, 1)
sim_settings = as.matrix(expand.grid(error_distr, cov_distr, transform_x, SNR))
colnames(sim_settings) = c("err_distr", "cov_distr", "transf_x", "SNR")

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
  tr_s = unname(s_s[3])
  snr_s = as.numeric(unname(s_s[4]))
  k_s = "init"
  
  res_tmp = foreach(
    rep      = 1:nrep,
    .packages = c("quantreg", "HDInterval")) %dopar% {
      
      # DEFINE ERROR SAMPLER
      rerr = switch(ed_s,
                    Gaussian = function(ndraws, var_eta, SNR) rnorm(ndraws, 0, sd = sqrt(var_eta/SNR)),
                    Gamma = function(ndraws, var_eta, SNR) rgamma(ndraws, 2, sqrt(2*SNR/var_eta)))
      
      rcov = switch(cd_s,
                    "Unif_rad2" = function(ndraws) runif(ndraws, -sqrt(2), sqrt(2)),
                    "Beta_2_4" = function(ndraws) runif(ndraws, 2, 4))
      
      transf = switch(tr_s,
                      "parabola" = function(x) x^2,
                      "cubic" = function(x) x^3,
                      "trigonometric" = function(x) sin(2*(4*x - 2)) + 2*exp(-16^2*(x-0.5)^2))
      
      
      # Use tryCatch to avoid weird datasets with outliers on x and hence singular Z
      n_retry <- 0L
      
      stopTry = FALSE
      
      repeat {
        tryCatch(
          {
            # DRAW COVARIATES
            x1 = rcov(n)
            mu_x1 = mean(x1); sd_x1 = sd(x1)
            eta = transf(x1)
            var_eta = var(eta)
            
            # Mode shift
            mode_shift = switch(ed_s,
                                Gaussian = 0,
                                Gamma = (2 - 1) / sqrt(2*SNR/var_eta))
            
            # DRAW RESPONSE
            y = eta + rerr(n, var_eta, SNR = snr_s) - mode_shift
            
            # CREATE DESIGN MATRIX
            DR <- construct_DR_basis(x1, d = 10)
            Z <- DR$Z   # nonlinear DR basis
            X_scaled <- cbind(1, (x1 - mu_x1)/sd_x1, Z)
            
            # OPTIMAL K
            const = (n^(5/12)) * (log(n)^(7/12)) / (n^(7/12))
            kinit = (n/const)^(1/7)
            res_qr = quantreg::rq(y ~ -1 + X_scaled, tau = 0.5)$residuals
            
            stopTry = TRUE
          }, 
          error = function(e) {
            n_retry <<- n_retry + 1L
            NULL
          })
        
        if (stopTry) break
      }
      
      
      # kfinal = kinit * (n^(1/7))
      k = switch(k_s,
                 "0.5" = 0.5,
                 "init" = kinit,
                 "final" = kfinal)
      
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
      
      out_optim_scaled = sapply(c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"),
                                function(m){
                                  est = optim(par = rep(0, ncol(X_scaled)),
                                              fn = function(b) sum(loss_asymm2(y - X_scaled%*%b, 
                                                                               k_multi, 
                                                                               1e-3)),
                                              method = m)$par
                                  est
                                }
                                
      )
      
      p = 1; d = 10
      out_MCMC = MCMC_MH(niter = niter, X = X_scaled, y = y, 
                         k = k_multi, c = 1e-3, 
                         p = p, d = d,
                         asymm = T, verbose = F,
                         beta0 = out_optim_scaled[ , "SANN"])
      
      draws = out_MCMC$beta[-(1:1000), ]
      draws_thin5 = draws[seq(5, niter-1000, by = 5), ]
      draws_thin10 = draws[seq(10, niter-1000, by = 10), ]
      
      draws_tau = out_MCMC$tau[-(1:1000)]
      
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
      ), t(out_optim_scaled))
      
      diagn = rbind(apply(draws_thin10, 2, function(x) acf(x, plot = F)[[1]][2]),
                    apply(draws_thin10, 2, posterior::rhat),
                    apply(draws_thin10, 2, posterior::ess_basic),
                    apply(draws_thin5, 2, function(x) acf(x, plot = F)[[1]][2]),
                    apply(draws_thin5, 2, posterior::rhat),
                    apply(draws_thin5, 2, posterior::ess_basic),
                    apply(draws, 2, function(x) acf(x, plot = F)[[1]][2]),
                    apply(draws, 2, posterior::rhat),
                    apply(draws, 2, posterior::ess_basic))
      
      rownames(diagn) = c(outer(c("acf_thin", "rhat_thin", "ess_thin"), c(10, 5, 1), paste0))
      
      list(draws = draws, draws_tau = draws_tau,
           summ = summ, diagn = diagn, 
           w = out_MCMC$w, k = k_multi,
           x1 = x1, X_scaled = X_scaled, y = y, Z = Z,
           setting = s_s, n_retry = n_retry,
           sd_prop = out_MCMC$sd_prop, alfa = out_MCMC$alfa)
      
    }
  
  res = res_tmp
  results[[s]] = res
  cat("Update: ", format(Sys.time()), "\t", sep = "")
  cat(round(100*s/nrow(sim_settings), 2), "%\n", sep = "")
}
saveRDS(results, "sim_study_nonC-GBI/splines_univ/sim_study_nonCalibrated_asymm_TL_init.RDS")
