library(dplyr)
library(ggplot2)

# Loss function
loss = function(y, x, beta, sigma2){
  # -dnorm(y, mean = x %*% beta, sd = sqrt(sigma2), log = T)
  -dnorm(y, mean = x * beta, sd = sqrt(sigma2), log = T)
}

compute_m_i = function(y_i, x_i, beta, sigma2){
  r_i = y_i - x_i * beta
  nabla_beta = (1/2/sigma2) * 2 * r_i * (-x_i)
  nabla_sigma2 = (1/2) / sigma2 - (r_i^2)/(2 * sigma2^2)
  c(nabla_beta, nabla_sigma2)
}

loss_CGBI = function(y, x, beta, sigma2){
  n = length(y)
  m_i = t(sapply(1:n, function(i) compute_m_i(y[i], x[i], beta, sigma2)))
  mbar = (1/n) * colSums(m_i)
  W = var(m_i) + diag(1e-5, 2)
  Q = drop((1/2) * t(mbar) %*% solve(W) %*% mbar)
  +0.5*log(det(W)) + n*Q
}

# Simulation setting
n = 1e3
true_beta = 2
x = runif(n, 0, 2)
eta = x*true_beta
gam = 0
true_var = 2/3 + 1/3 * abs(x)^gam
y = eta + rnorm(n, 0, sqrt(true_var))

# CHECK OPTIM
L = function(params) sum(loss(y = y, x = x, beta = params[1], sigma2 = exp(params[2])))
L_CGBI = function(params) loss_CGBI(y = y, x = x, beta = params[1], sigma2 = exp(params[2]))

optim_D = sapply(c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"), 
       function(m)
       {tmp = optim(par = c(2, 0), 
                    fn = L,
                    method = m)
       c(tmp$par[1], exp(tmp$par[2]), tmp$convergence)
       }
)

optim_CGBI = sapply(c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"), 
                    function(m)
                    {tmp = optim(par = c(1.8, 0), 
                                 fn = L_CGBI,
                                 method = m)
                    c(tmp$par[1], exp(tmp$par[2]), tmp$convergence)
                    }
)

beta_values = seq(0, 4, length = 51)
sigma2_values = seq(0.25, 2.75, length = 25) #quantile(true_var, probs = seq(0, 1, length = 25))


# setwd("Q_vs_D_NormLik")
# pdf("Q_vs_D_fixedsigma.pdf", width = 12, height = 24)

for (s in 1:nrow(settings)){
  
  loglik_GBI = loglik_D = matrix(NA, nrow = length(beta_values), ncol = length(sigma2_values))
  dimnames(loglik_GBI) = dimnames(loglik_D) = list(beta = beta_values, sigma2 = sigma2_values)
  
  for (j in seq_along(beta_values)){
    for (k in seq_along(sigma2_values)){
      
      be = beta_values[j]
      si = sigma2_values[k]
      
      m_i = t(sapply(1:n, function(i) compute_m_i(y[i], x[i], be, si)))
      mbar = (1/n) * colSums(m_i)
      W = diag(c(mean(x^2)/si, 1/2/si^2), 2) #var(m_i)
      Q = drop((1/2) * t(mbar) %*% solve(W) %*% mbar)
      loglik_GBI[j, k] = -0.5*log(det(W)) - n*Q
      
      loglik_D[j, k] = -sum(loss(y, x, be, si))
    }
  }
  
  df_D <- loglik_D %>% 
    reshape2::melt()
  D_breaks <- quantile(df_D$value,
                       probs = seq(0, 1, length.out = 51),
                       na.rm = TRUE)
  
  df_GBI <- loglik_GBI %>% 
    reshape2::melt()
  GBI_breaks <- quantile(df_GBI$value,
                         probs = seq(0, 1, length.out = 51),
                         na.rm = TRUE)
  
  ggplot(df_D) +
    geom_contour_filled(
      aes(x = beta, y = sigma2, z = value),
      col = "white", lwd = 0.1,
      breaks = q_breaks,
      show.legend = FALSE) + 
    geom_point(data = data.frame(x = true_beta, y = 1), 
               mapping = aes(x = x, y = 1),
               col = "red") +
    labs(title = sprintf("-D (n=%s)", n))
  
  ggplot(df_GBI) +
    geom_contour_filled(
      aes(x = beta, y = sigma2, z = value),
      col = "white", lwd = 0.1,
      breaks = GBI_breaks,
      show.legend = FALSE) +
    geom_point(data = data.frame(x = true_beta, y = 1), 
               mapping = aes(x = x, y = 1),
               col = "red") +
    labs(title = sprintf("CGBI loglik (n=%s)", n))
  
  # Plots likelihoods
  plot(beta_values, loglik_D, xlab = "beta", ylab = "-D",
       main = sprintf("n=%s", n), type = "l")
  plot(beta_values, loglik_GBI, xlab = "beta", ylab = "CGBI loglik",
       main = sprintf("n=%s", n), type = "l")
  plot(loglik_D, loglik_GBI, xlab = "-D", ylab = "Calibrated GBI", col = pal,
       main = sprintf("n=%s", n), type = "p")
  legend("bottomright", pch = 1, col = pal[c(1, 25, 75, 101)], 
         legend = paste0("beta = ", c(0, 1, 3, 4)))
}
dev.off()
