rm(list = ls())
source("src/01-utils.R")

# READ TYPE OF K ----
cmdline = commandArgs(trailingOnly = TRUE)

# SAMPLE SIZE ----
n = as.integer(cmdline[1])
# Number fo covariates
p = 1
# Number of splines
d = 20
# Type of k (fixed to "init")
k_s = "init"


# SIMULATION SETTINGS ----
error_distr = c("Gaussian", "Gamma")
cov_distr = c("Unif_rad2", "Beta_2_4") #"Gamma_1.5_1.5") #, "BVN")
transform_x = c("parabola", "cubic", "trigonometric")
# if (n == 1e3) {SNR = c(2, 1)} else {SNR = c(5, 2)}
SNR = c(5, 2)
sim_settings = as.matrix(expand.grid(error_distr, cov_distr, transform_x, SNR))
colnames(sim_settings) = c("err_distr", "cov_distr", "transf_x", "SNR")
ev_s = as.character(cmdline[2])  #c("hetero", "homo")

results = vector("list", nrow(sim_settings))
niter = 15e3
burnin = 5e3
nrep = 100

# PARALLEL BACKEND ----
library(doParallel); library(foreach); library(doRNG)
ncores = 25 #max(1, parallel::detectCores() - 1)
registerDoParallel(ncores)

cat("Starting time: ", format(Sys.time()), "\n", sep = "")

for (s in 1:nrow(sim_settings)){
  
  s_s = sim_settings[s, ]
  ed_s = unname(s_s[1])
  cd_s = unname(s_s[2])
  tr_s = unname(s_s[3])
  snr_s = as.numeric(unname(s_s[4]))
  
  res_tmp = foreach(
    rep      = 1:nrep,
    .packages = c("quantreg", "HDInterval"),
    .options.RNG = 76137) %dorng% {
      
      # DEFINE ERROR SAMPLER
      rerr = switch(ed_s,
                    Gaussian = function(ndraws, var_eta, SNR, hc) rnorm(ndraws, 0, sd = sqrt(var_eta/SNR)*sqrt(hc)),
                    Gamma = function(ndraws, var_eta, SNR, hc) rgamma(ndraws, 2, sqrt(2*SNR/var_eta)/sqrt(hc)))
      
      verr = switch(ev_s,
                    "homo" = function(x) 1,
                    "hetero" = function(x) 0.25 + x^2)
      
      rcov = switch(cd_s,
                    "Unif_rad2" = function(ndraws) runif(ndraws, -sqrt(2), sqrt(2)),
                    "Beta_2_4" = function(ndraws) 2*sqrt(2)*(rbeta(ndraws, 2, 4)-0.5))
      
      transf = switch(tr_s,
                      "parabola" = function(x) x^2,
                      "cubic" = function(x) x^3,
                      "trigonometric" = function(x) sin(2.5*x) + 3*exp(-(5^2)*((x)^2))
      )
      
      
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
            
            hetero_correction = verr(x1)/mean(verr(x1))
            
            # Mode shift
            mode_shift = switch(ed_s,
                                Gaussian = 0,
                                Gamma =  (2 - 1) / sqrt(2*snr_s/var_eta) * sqrt(hetero_correction))
            
            # DRAW RESPONSE
            y = eta + rerr(n, var_eta, SNR = snr_s, hc = hetero_correction) - mode_shift
            
            # CREATE DESIGN MATRIX
            DR <- construct_DR_basis(x1, d = d)
            Z <- DR$Z   # nonlinear DR basis
            X_scaled <- cbind(1, (x1 - mu_x1)/sd_x1, Z)
            
            # OPTIMAL K
            const = (n^(5/12)) * (log(n)^(7/12)) / (n^(7/12))
            kinit = (n/const)^(1/7)
            kfinal = kinit * (n^(1/7))
            res_qr = quantreg::rq(y ~ -1 + X_scaled, tau = 0.5)$residuals
            
            stopTry = TRUE
          }, 
          error = function(e) {
            n_retry <<- n_retry + 1L
            NULL
          })
        
        if (stopTry) break
      }
      
      
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
        sqrt(k)/factorn,
        sqrt(k)/factorp
      )
      
      out_optim_scaled = sapply(c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"),
                                function(m){
                                  est = optim(par = rep(0, ncol(X_scaled)),
                                              fn = function(b) sum(loss_asymm(y - X_scaled%*%b, 
                                                                              (k_multi)^(1/2), 
                                                                              1e-1)),
                                              method = m)$par
                                  est
                                }
      )
      
      loss_optim = apply(out_optim_scaled, 2, function(b){
        loss_asymm_pop(b, y, X_scaled, k_multi, 1e-1)
      })
      opt = which.min(loss_optim)
      b_init = out_optim_scaled[ , opt]
      
      tryCatch({
        out_MCMC = MCMC_IWLS(niter = niter, X = X_scaled, y = y, 
                             k = k_multi, k_deriv = k_multi,
                             c = 1e-1, c_deriv = 1e-1, 
                             p = p, d = d,
                             a_tau = c(2.1, 1e-3), b_tau = c(qgamma(0.001, shape = 2.1, rate = 1) * 5^2, 1e-3),
                             verbose = 0,
                             beta0 = b_init)
        
        draws = out_MCMC$beta[-(1:burnin), ]
        draws_thin5 = draws[seq(5, niter-burnin, by = 5), ]
        draws_thin10 = draws[seq(10, niter-burnin, by = 10), ]
        
        draws_tau = out_MCMC$tau[-(1:burnin), ]
        draws_tau_thin10 = draws_tau[seq(10, niter-burnin, by = 10), ]
        
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
        
        list(draws = draws_thin10, draws_tau = draws_tau_thin10,
             summ = summ, diagn = diagn, 
             w = out_MCMC$w, 
             k = kmulti, k_deriv = k_multi,
             c = 1e-1, c_deriv = 1e-1, 
             x1 = x1, X_scaled = X_scaled, y = y, Z = Z,
             DR = DR[c("Trafo", "smooth_object", "Gram", "Penalty")],
             setting = s_s, n_retry = n_retry,
             sd_prop = out_MCMC$sd_prop, alfa = out_MCMC$alfa, acc_prop = out_MCMC$acc_prop)
      }, error = function(e) {
        message("Task ", rep, " failed: ", e$message)
        return(list(draws = NULL, draws_tau = NULL,
                    summ = NULL, diagn = NULL, 
                    w =NULL, k = k_multi,
                    x1 = x1, X_scaled = X_scaled, y = y, Z = Z,
                    setting = s_s, n_retry = n_retry,
                    sd_prop = NULL, alfa = NULL))  # or NULL
      })
    }
  
  res = res_tmp
  results[[s]] = res
  cat("Update: ", format(Sys.time()), "\t", sep = "")
  cat(round(100*s/nrow(sim_settings), 2), "%\n", sep = "")
}

SAVE_PATH = paste0("sim_study_nonC-GBI/splines_univ/ss_p1_n", n, "_", ev_s, ".RDS")
saveRDS(list(res = results,
             k_and_c = list("k" = "sqrt(kmulti)", "k_deriv" = "sqrt(kmulti)", 
                            "c" = "1e-1", "c_deriv" = "1e-1",
                            "hyperpar tau" = c("lin" = "informative (Nico)",
                                               "nonlin" = "non=informative"))), SAVE_PATH)
