rm(list = ls())

# UTILS ----
# Loss function
loss = function(eps, k, c, g = 1){
  h = 2*g
  1 - exp(c^(1/h)-( (k*eps)^h + c )^(1/h) )
  # 1 - exp(c^(1/h)-( (k*eps)^h + c ) )
}

loss_pop = function(beta, y, x, k, c= 1e-3, g = 1){
  n = length(y)
  sum(sapply(1:n, function(i) loss(y[i]-x[i]*beta, k, c, g)))
}

loss_pop_V = Vectorize(loss_pop, "beta")

# Derivative of loss wrt residual
compute_t = function(r, k, c, g = 1){
  h = 2*g
  s = (k*r)^(h) + c
  t = exp(c^(1/h)) * exp(-s^(1/h)) * s^(1/h-1) * k^(h) * r^(h-1)
  return(t)
}



MCMC_MH = function(niter, x, y, k, c, beta0, prior_sd, calibrated = T, debug = F){
  
  p = 1
  n = length(x)
  
  beta_sample = rep(NA, niter)
  beta_init = beta0
  last = beta_init
  acc = 0
  alfa_sample = rep(NA, niter)
  loss_sample = mbar_ind_sample = array(NA, dim = c(niter, n, 2))
  dimnames(loss_sample)[[3]] = dimnames(mbar_ind_sample)[[3]] = c("last", "prop")
  sd_prop_sample = c()
  
  for (m in 1:niter){
    
    if (m %% (niter/10) == 0){
      cat(10 * m / (niter/10), "% \t")
    }
    
    if (m == 1){
      RM = 0
    } else { 
      RM = log(sd_prop) + (1/m)^(1/2) * (acc/(m-1) - 0.234)
    }
    
    sd_prop = exp(RM); sd_prop_sample = c(sd_prop_sample, sd_prop)
    prop = rnorm(1, mean = last, sd = sd_prop)
    
    if (debug){
      cat("\n\n beta last:", last, " |  beta prop:", prop, "\n\n")
    }
    
    
    if (debug){
      plot(x, y); abline(a = 0, b = true_beta, col = 2); points(x, last*x, col = 3); points(x, prop*x, col = 4) 
    }
    
    r_last = y - last*x
    r_prop = y - prop*x
    
    loss_last = loss(r_last, k, c)
    loss_prop = loss(r_prop, k, c)
    loss_sample[m, , ] = matrix(c(loss_last, loss_prop), nrow = n, ncol = 2)
    
    if (debug){
      plot(density(r_last)); lines(density(r_prop), col = 2)
      cat("\n\n Mean losses", colMeans(loss_sample[m, , ]), "\n\n")
    }
    
    if (calibrated){
      mbar_ind_last = -compute_t(r_last, k, c)*x
      mbar_ind_prop = -compute_t(r_prop, k, c)*x
      mbar_ind_sample[m , ,] = matrix(c(mbar_ind_last, mbar_ind_prop), nrow = n, ncol = 2)
      
      
      if (debug){
        check_prop = numDeriv::grad(function(x) loss(x, k, c), x = r_prop)*x
        check_last = numDeriv::grad(function(x) loss(x, k, c), x = r_last)*x
        cat("\n\n Check mbar", all(round(check_prop, 6) == round(mbar_ind_prop, 6)), " | ", all(round(check_last, 6) == round(mbar_ind_last, 6)), "\n\n")
      }
      mbar_last = sum(mbar_ind_last)/n
      mbar_prop = sum(mbar_ind_prop)/n
      
      W_last = var(mbar_ind_last)*(n-1)/n
      W_prop = var(mbar_ind_prop)*(n-1)/n
      
      Q_last = 0.5*mbar_last^2/W_last
      Q_prop = 0.5*mbar_prop^2/W_prop
      
      loglik_last = -0.5*log(W_last) - n*Q_last
      loglik_prop = -0.5*log(W_prop) - n*Q_prop
    } else {
      loglik_last = -sum(loss_last)
      loglik_prop = -sum(loss_prop)
    }
    
    if (debug){
      cat("loglik last:", loglik_last, "\t | \t loglik prop:", loglik_prop, "\n")
    }
    
    target_last = loglik_last - 0.5*last^2/(prior_sd^2)
    target_prop = loglik_prop - 0.5*prop^2/(prior_sd^2)
    
    alfa = alfa_sample[m] = exp(min(target_prop - target_last, 0))
    u = runif(1, 0, 1)
    
    if (u < alfa){
      beta_sample[m] = prop
      last = prop
      is_acc = 1
    } else {
      beta_sample[m] = last
      is_acc = 0
    }
    acc = acc + is_acc
  }
  
  return(list(beta = beta_sample, 
              acc_prop = acc/niter, alfa = alfa_sample, sd_prop = sd_prop_sample, 
              loss = loss_sample, mbar_ind = mbar_ind_sample))
}

# TRUE BETA AND SAMPLE SIZE ----
n = 1e3
true_beta = -0.5

# SIMULATION SETTINGS ----
error_distr = c("Gaussian", "Gaussian_0.5", "Gamma_2_2")
covariate = c("0.5", "1", "asymm")
k = c("init", "final")
calibration = c(T)
sim_settings = as.matrix(expand.grid(error_distr, covariate, k, calibration))

results = vector("list", nrow(sim_settings))
niter = 2e4
nrep = 100

# PARALLEL BACKEND ----
library(doParallel); library(foreach)
ncores = 5#max(1, parallel::detectCores() - 1)
cl = makeCluster(ncores)
registerDoParallel(cl)

cat("Starting time: ", format(Sys.time()), "\n", sep = "")

for (s in 1:nrow(sim_settings)){
  
  s_s = sim_settings[s, ]
  ed_s = unname(s_s[1])
  cd_s = unname(s_s[2])
  k_s = unname(s_s[3])
  cal_s = unname(s_s[4]); isCal = as.logical(cal_s)
  
  # column layout (you had ncol=14 but only stored 9 values)
  coln <- c("mean","median","sd","q025","q0975","HDI_L1","HDI_U1","HDI_L2","HDI_U2",
            "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN")
  
  res_tmp = foreach(
    rep      = 1:nrep,
    .combine = rbind,
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
      x = rcov(n)
      eta = x*true_beta
      
      # DRAW RESPONSE
      y = x*true_beta + rerr(n) - mode_shift
      
      
      # OPTIMAL K
      const = (n^(5/12)) * (log(n)^(7/12)) / (n^(7/12))
      kinit = (n/const)^(1/7)
      sd_quant = sd(quantreg::rq(y ~ x, tau = 0.5)$residuals)
      kfinal = kinit * (n^(1/7)) / sd_quant
      k = switch(k_s,
                 "0.5" = 0.5,
                 "init" = kinit,
                 "final" = kfinal) 
      
      out_optim = sapply(c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"), 
                         function(m)
                           optim(par = -0.25, 
                                 fn = function(b) sum(loss(y - x*b, k, 1e-3)),
                                 method = m)$par
      )
      
      out_MCMC = MCMC_MH(niter = niter, x = x, y = y, k = k, c = 1e-3, 
                         beta0 = -0.25, prior_sd = 1/3, calibrated = isCal)
      
      draws = out_MCMC$beta[-(1:1000)]
      mean_b = mean(draws)
      med_b = median(draws)
      sd_b = sd(draws)
      q025 = unname(quantile(draws, probs = c(0.025)))
      q0975 = unname(quantile(draws, probs = c(0.975)))
      hdi_b = HDInterval::hdi(density(draws), credMass = 0.95, allowSplit = T)
      if (nrow(hdi_b) == 1){
        HDI_L1 = unname(hdi_b[1, "begin"]); HDI_U1 = unname(hdi_b[1, "end"])
        HDI_L2 = HDI_U2 = NA
      } else {
        HDI_L1 = unname(hdi_b[1, "begin"]); HDI_U1 = unname(hdi_b[1, "end"])
        HDI_L2 = unname(hdi_b[2, "begin"]); HDI_U2 = unname(hdi_b[2, "end"])
      }
      
      c(mean = mean_b, median = med_b, sd = sd_b, 
        q025 = q025, q0975 = q0975,
        HDI_L1 = HDI_L1, HDI_U1 = HDI_U1, HDI_L2 = HDI_L2, HDI_U2 = HDI_U2,
        out_optim)
    }
  
  colnames(res_tmp) = coln
  results[[s]] <- res_tmp
  cat("Update: ", format(Sys.time()), "\t", sep = "")
  cat(round(100*s/nrow(sim_settings), 2), "%\n", sep = "")
}
stopCluster(cl)
# saveRDS(results, "sim_study.RDS")

par(mfrow = c(3, 3))
for (j in 1:9){
  min = -0.7 #min(c(results[[j]][, "mean"], results[[j]][ , "Nelder-Mead"])) 
  max = 0 #max(c(results[[j]][, "mean"], results[[j]][ , "Nelder-Mead"])) 
  plot(density(results[[j]][, "mean"]), main = paste(sim_settings[j, -4], collapse = " - "), xlim = c(min, max))
  lines(density(results[[j]][ , "Nelder-Mead"]), col = "blue")
  lines(density(results[[j]][ , "BFGS"]), col = "green", lty = 2)
  abline(v = true_beta, col = "gold")
}

for (j in 10:18){
  min = min(c(results[[j]][, "mean"], results[[j]][ , "Nelder-Mead"])) 
  max = max(c(results[[j]][, "mean"], results[[j]][ , "Nelder-Mead"])) 
  plot(density(results[[j]][, "mean"]), main = paste(sim_settings[j, -4], collapse = " - "), xlim = c(min, max))
  lines(density(results[[j]][ , "Nelder-Mead"]), col = "blue")
  lines(density(results[[j]][ , "BFGS"]), col = "green", lty = 2)
  abline(v = true_beta, col = "gold")
}

mean_mean = sapply(results, function(R) mean(R[, "mean"]))
sd_mean = sapply(results, function(R) sd(R[, "mean"]))
mean_sd = sapply(results, function(R) mean(R[, "sd"]))
coverage_SymmI = sapply(results, function(R) {
  I = R[ , c("q025", "q0975")]
  mean(apply(I, 1, function(i) true_beta >= i[1] & true_beta <= i[2]))
})
cbind(sim_settings[ , -4], 
      mean_mean = round(mean_mean, 2), sd_mean = round(sd_mean, 2), 
      mean_sd = round(mean_sd, 2),
      coverage_SymmI)

# beta_val = seq(-2, 1, by = 0.025)
# 
# loss_val = sapply(beta_val, function(b) 
#   mean(sapply(1:n, function(i) loss(y[i] - b*x[i], k, 1e-3))))
# mi_val = sapply(beta_val, function(b) 
#   sapply(1:n, function(i) -x[i]*compute_t(y[i]-b*x[i], k, 1e-3)))
# mbar_val = colMeans(mi_val); w_val = apply(mi_val, 2, var); Q_val = -0.5*mbar_val^2/w_val
# mbar_val_gr = numDeriv::grad(function(b) loss_pop(b, y, x, k), x = beta_val)/n
# nDenv = new.env(); nDenv$beta = beta_val; nDenv$y = y; nDenv$x = x; nDenv$k = k
# mbar_val_nD = drop(attr(numericDeriv(
#   expr = quote(loss_pop(beta, y, x, k)),
#   theta = "beta",
#   rho = nDenv
# ), "gradient"))/n
# mi_val_gr = matrix(NA, nrow = n, ncol = length(beta_val))
# for (j in seq_along(beta_val)){
#   for (i in 1:n){
#     res = y[i] - beta_val[j]*x[i]
#     mi_val_gr[i , j] = -x[i]*numDeriv::grad(function(eps) loss(eps, k = k, c = 1e-3), x = res)
#   }
# }
# mbar_val_grsum = colSums(mi_val_gr)/n
# 
# 
# par(mfrow = c(2, 2))
# plot(beta_val, loss_val, type = "b")
# plot(beta_val, mbar_val, type = "b", ylim = c(min(mbar_val)*1.25, max(mbar_val)*1.25), cex = 1.3, pch = 19)
# points(beta_val, mbar_val_nD, pch = 0, col = 3); points(beta_val, mbar_val_gr, pch = 4, col = 2)
# points(beta_val, mbar_val_grsum, pch = 17, col = 4, type = "b")
# plot(beta_val, -n*loss_val, type = "b")
# plot(beta_val, -Q_val, type = "b")
# 
# 
# 
# 
# sapply(c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"), 
#        function(m)
#          optim(par = 0, 
#                fn = function(b) sum(loss(y - b*x, k, 1e-3)),
#                method = m)$par
# )
# 
# mean(out$beta[-(1:1000)])
# hist(out$beta[-(1:1000)], breaks = 50)
# 
# 
# 
