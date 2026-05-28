rm(list = ls())
source("src/01-utils.R")

# READ TYPE OF K ----
cmdline = commandArgs(trailingOnly = TRUE)

# SAMPLE SIZE ----
n = as.integer(cmdline[1])
# Number fo covariates
p = 2
# Number of splines
d = rep(20, p)
# Type of k (fixed to "init")
k_s = "init"


# SIMULATION SETTINGS ----
error_distr = c("Gaussian", "Gamma")
cov_distr = "Unif_rad2" #c("Unif_rad2", "Beta_2_4") #"Gamma_1.5_1.5") #, "BVN")
transform_x1 = c("parabola", "cubic", "trigonometric")
transform_x2 = c("parabola", "cubic", "trigonometric")
SNR = c(5, 2)
sim_settings = as.matrix(expand.grid(error_distr, cov_distr, transform_x1, transform_x2, SNR))
colnames(sim_settings) = c("err_distr", "cov_distr", "transf_x1", "transf_x2", "SNR")
sim_settings = sim_settings[sim_settings[ , "transf_x1"] != sim_settings[ , "transf_x2"], ]
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
  tr1_s = unname(s_s[3])
  tr2_s = unname(s_s[4])
  snr_s = as.numeric(unname(s_s[5]))
  
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
      
      transf1 = switch(tr1_s,
                       "parabola" = function(x) x^2,
                       "cubic" = function(x) x^3,
                       "trigonometric" = function(x) sin(2.5*x) + 3*exp(-(5^2)*((x)^2))
      )
      
      transf2 = switch(tr2_s,
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
            x2 = rcov(n)
            mu_x1 = mean(x1); sd_x1 = sd(x1)
            mu_x2 = mean(x2); sd_x2 = sd(x2)
            eta = transf1(x1) + transf2(x2)
            var_eta = var(eta)
            
            hetero_correction = 1 #verr(x1)/mean(verr(x1))
            
            # Mode shift
            mode_shift = switch(ed_s,
                                Gaussian = 0,
                                Gamma =  (2 - 1) / sqrt(2*snr_s/var_eta) * sqrt(hetero_correction))
            
            # DRAW RESPONSE
            y = eta + rerr(n, var_eta, SNR = snr_s, hc = hetero_correction) - mode_shift
            
            # CREATE DESIGN MATRIX
            DR1 <- construct_DR_basis(x1, d = d[1])
            Z1 <- DR1$Z
            DR2 <- construct_DR_basis(x2, d = d[2])
            Z2 <- DR2$Z
            X_scaled <- cbind(1, (x1 - mu_x1)/sd_x1, (x2 - mu_x2)/sd_x2,
                              Z1, Z2)
            
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
      
      out_optim_scaled = sapply(
        c("Nelder-Mead", "CG", "L-BFGS-B", "SANN"),
        # c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"),
        function(m){
          est = optim(par = rep(0, ncol(X_scaled)),
                      fn = function(b) sum(loss_asymm(y - X_scaled%*%b, 
                                                      k_multi, 
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
        a_tau1 = 2.1; a_tau2 = 1e-3
        b_tau1 = qgamma(0.001, shape = 2.1, rate = 1) * 5^2; b_tau2 = 1e-3
        
        k_post = k_multi; k_prop = (k_multi)^(4/5)
        c_post = 1e-1; c_prop = 1e-1
        
        out_MCMC = MCMC_IWLS(niter = niter, X = X_scaled, y = y, 
                             k = k_post, k_deriv = k_prop,
                             c = c_post, c_deriv = c_prop, 
                             p = p, d = d,
                             a_tau = c(a_tau1, a_tau1, a_tau2, a_tau2),
                             b_tau = c(b_tau1, b_tau1, b_tau2, b_tau2),
                             verbose = 0, K_fold = 50, print_step = 2.5e3,
                             beta0 = b_init)
        
        draws = out_MCMC$beta[-(1:burnin), ]
        draws_thin5 = draws[seq(5, niter-burnin, by = 5), ]
        draws_thin10 = draws[seq(10, niter-burnin, by = 10), ]
        
        draws_tau = out_MCMC$tau[-(1:burnin), ]
        draws_tau_thin10 = draws_tau[seq(10, niter-burnin, by = 10), ]
        
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
             diagn = diagn, 
             w = out_MCMC$w, 
             k = k_post, k_deriv = k_prop,
             c = c_post, c_deriv = c_prop, 
             x1 = x1, x2 = x2, y = y,
             Z1 = Z1, Z2 = Z2,
             DR1 = DR1[c("Trafo", "smooth_object")],
             DR2 = DR2[c("Trafo", "smooth_object")],
             setting = s_s, n_retry = n_retry,
             acc_prop = out_MCMC$acc_prop)
      }, error = function(e) {
        message("Task ", rep, " failed: ", e$message)
        return(list(draws = NULL, draws_tau = NULL,
                    summ = NULL, diagn = NULL, 
                    w =NULL, k = k_multi,
                    x1 = x1, X_scaled = X_scaled, y = y, Z1 = Z1, Z2 = Z2,
                    setting = s_s, n_retry = n_retry,
                    sd_prop = NULL, alfa = NULL))  # or NULL
      })
    }
  
  res = res_tmp
  results[[s]] = res
  cat("Update: ", format(Sys.time()), "\t", sep = "")
  cat(round(100*s/nrow(sim_settings), 2), "%\n", sep = "")
}

SAVE_PATH = paste0("sim_study_nonC-GBI/splines_biv/ss_p2_n", n, "_", ev_s, ".RDS")
saveRDS(list(res = results,
             k_and_c = list("k" = "sqrt(kinit)/factor", "k_deriv" = "sqrt(kinit)/factor", 
                            "c" = "1e-1", "c_deriv" = "1e-1",
                            "hyperpar tau" = c("lin" = "informative (Nico)",
                                               "nonlin" = "non=informative"))), SAVE_PATH)
