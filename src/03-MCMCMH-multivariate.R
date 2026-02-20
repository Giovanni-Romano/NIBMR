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
  return(drop(t))
}


v1 = function(c){
  c^2 / 2 / 27 * ((27-2*c) + sqrt(27*(27-4*c)))
}

v2 = function(c){
  c^2 / 2 / 27 * ((27-2*c) - sqrt(27*(27-4*c)))
}

thresh = function(c){
  sqrt(v1(c)^(1/3) + v2(c)^(1/3) - c/3)
}


MCMC_MH = function(niter, X, y, k, c, p, beta0, prior_sd, calibrated = T, debug = F){
  
  if (p != length(beta0)){stop("Length of beta0 is not equal to p")}
  
  n = nrow(X)
  
  beta_sample = matrix(NA, nrow = niter, ncol = p)
  beta_init = beta0
  last = beta_init
  acc = 0
  alfa_sample = rep(NA, niter)
  loss_sample = array(NA, dim = c(niter, n, 2))
  mbar_ind_sample = array(NA, dim = c(niter, n, p, 2))
  dimnames(loss_sample)[[3]] = dimnames(mbar_ind_sample)[[4]] = c("last", "prop")
  sd_prop_sample = c()
  
  Sigma0 = diag(1, nrow = p)
  cholSigma0 = chol(Sigma0)
  
  for (m in 1:niter){
    
    if (m %% (niter/10) == 0){
      cat(10 * m / (niter/10), "% \t")
    }
    
    if (m == 1){
      RM = log(2.38/sqrt(p))
    } else {
      RM = log(sd_prop) + (1/m)^(1/2) * (acc/(m-1) - 0.234)
    }
    
    sd_prop = exp(RM); sd_prop_sample = c(sd_prop_sample, sd_prop)
    prop = t(mvnfast::rmvn(n = 1, mu = last, sigma = sd_prop*cholSigma0, isChol = T))
    
    r_last = y - X%*%last
    r_prop = y - X%*%prop
    
    loss_last = loss(r_last, k, c)
    loss_prop = loss(r_prop, k, c)
    loss_sample[m, , ] = matrix(c(loss_last, loss_prop), nrow = n, ncol = 2)
    
    if (debug){
      plot(density(r_last)); lines(density(r_prop), col = 2)
      cat("\n\n Mean losses", colMeans(loss_sample[m, , ]), "\n\n")
    }
    
    if (calibrated){
      mbar_ind_last = -X*compute_t(r_last, k, c)
      mbar_ind_prop = -X*compute_t(r_prop, k, c)
      mbar_ind_sample[m , , , "last"] = matrix(mbar_ind_last, nrow = n, ncol = p)
      mbar_ind_sample[m , , , "prop"] = matrix(mbar_ind_prop, nrow = n, ncol = p)
      
      if (debug){
        check_prop = -X*numDeriv::grad(function(x) loss(x, k, c), x = r_prop)
        check_last = -X*numDeriv::grad(function(x) loss(x, k, c), x = r_last)
        cat("\n\n Check mbar", all(round(check_prop, 5) == round(mbar_ind_prop, 5)), " | ", all(round(check_last, 5) == round(mbar_ind_last, 5)), "\n\n")
      }
      
      mbar_last = colSums(mbar_ind_last)/n
      mbar_prop = colSums(mbar_ind_prop)/n
      
      W_last = var(mbar_ind_last)*(n-1)/n
      W_prop = var(mbar_ind_prop)*(n-1)/n
      
      Q_last = drop(0.5*t(mbar_last) %*% solve(W_last) %*% mbar_last)
      Q_prop = drop(0.5*t(mbar_prop) %*% solve(W_prop) %*% mbar_prop)
      
      loglik_last = -0.5 * log(det(W_last)) - 0.5*log(det(W_last)) - n*Q_last
      loglik_prop = -0.5 * log(det(W_prop)) - 0.5*log(det(W_prop)) - n*Q_prop
    } else {
      loglik_last = -0.5 * sum(loss_last)
      loglik_prop = -0.5 * sum(loss_prop)
    }
    
    if (debug){
      cat("loglik last:", loglik_last, "\t | \t loglik prop:", loglik_prop, "\n")
    }
    
    target_last = loglik_last - 0.5*sum(last^2/(prior_sd^2))
    target_prop = loglik_prop - 0.5*sum(prop^2/(prior_sd^2))
    
    alfa = alfa_sample[m] = exp(min(target_prop - target_last, 0))
    u = runif(1, 0, 1)
    
    if (u < alfa){
      beta_sample[m, ] = prop
      last = prop
      is_acc = 1
    } else {
      beta_sample[m, ] = last
      is_acc = 0
    }
    acc = acc + is_acc
  }
  
  return(list(beta = beta_sample, 
              acc_prop = acc/niter, alfa = alfa_sample, sd_prop = sd_prop_sample, 
              loss = loss_sample, mbar_ind = mbar_ind_sample))
}


# Simulation setting
n = 1e3
true_beta = c(0, -0.5, +0.7)
X = cbind(1, runif(n, -1, 1), runif(n, -1, 1))
eta = X%*%true_beta
y = eta + rnorm(n, 0, sd = 1)
par(mfrow = c(1, 2)); plot(X[, 2], y); plot(X[, 3], y); par(mfrow = c(1, 1))

median_reg = quantreg::rq(y ~ -1 + ., data = data.frame(X), tau = 0.5)
sd_mr = sd(median_reg$residuals)
const = (n^(5/12)) * ((log(n))^(7/12)) / (n^(7/12))
kinit = (n/const)^(1/7)
k = kinit
thresh(1e-3)/k
round(quantile(abs(median_reg$residuals), probs = seq(0, 1, by = 0.1)), 2)

sapply(c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"), 
       function(m)
         optim(par = c(0.25, 0, 0), 
               fn = function(b) sum(loss(y - X%*%b, k, 1e-3)),
               method = m)$par
)

niter = 1e4
out = MCMC_MH(niter = niter, X = X, y = y, k = k, c = 1e-3, p = 3, beta0 = c(-0.5, 0, 0), prior_sd = 1)

par(mfrow = c(1, 3))
for (j in 1:3){
  plot(1:niter, out$beta[ , j], type = "l")
  lines(1001:niter, cumsum(out$beta[1001:niter, j])/(1:(niter-1000)), col = "red")
  abline(h = true_beta[j], col = "gold")
}

par(mfrow = c(1, 1))
out$acc_prop
plot(cumsum(out$alfa)/(1:10^4), type = "l"); abline(h = 0.234, col = "gold")
plot(out$sd_prop[-(1:1000)], type = "l")

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
