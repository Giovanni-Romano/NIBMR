# UTILS ----
# Loss function
loss_symm = function(eps, k, c, g = 1, theta = 0){
  h = 2*g
  1 - exp( c^(1/h)-(((k*eps)^h + c)^(1/h)) )
}

loss_asymm = function(eps, k, c, g = 1, theta = 0){
  h = 2*g
  m1 = tan(pi/4 - theta)
  1 - exp( c^(1/h)-(m1^sign(eps))*(((k*eps)^h + c)^(1/h)) )
}


loss_asymm2 = function(eps, k, c, g = 1, theta = pi/6){
  k1 = k[1]
  k2 = k[2]
  h = 2*g
  (eps < 0) * (1 - exp(c^(1/h)-( (k1*eps)^h + c )^(1/h) )) + 
    (eps > 0) * (1 - exp(c^(1/h)-( (k2*eps)^h + c )^(1/h) ))
}


loss_pop = function(beta, y, X, k, c= 1e-3, g = 1){
  n = length(y)
  sum(sapply(1:n, function(i) loss(y[i]-X[i, ]%*%beta, k, c, g)))
}


MCMC_MH = function(niter, X, y, k, c, p, beta0, prior_sd, asymm = F, debug = F){
  
  if (p != length(beta0)){stop("Length of beta0 is not equal to p")}
  
  n = nrow(X)
  
  loss = ifelse(asymm, loss_asymm2, loss_symm)
  
  beta_sample = matrix(NA, nrow = niter, ncol = p)
  beta_init = beta0
  last = beta_init
  acc = 0
  alfa_sample = rep(NA, niter)
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
    
    loss_last = sapply(r_last, loss, k = k, c = c)
    loss_prop = sapply(r_prop, loss, k = k, c = c)
    
    if (debug){
      plot(density(r_last)); lines(density(r_prop), col = 2)
      cat("\n\n Mean losses", colMeans(loss_sample[m, , ]), "\n\n")
    }
    
    loglik_last = -0.5 * sum(loss_last)
    loglik_prop = -0.5 * sum(loss_prop)
    
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
              acc_prop = acc/niter, alfa = alfa_sample, sd_prop = sd_prop_sample))
}
