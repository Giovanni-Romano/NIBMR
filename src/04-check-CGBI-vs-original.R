rm(list = ls())

# UTILS ----
# Loss function
loss = function(eps, k, c, g){
  h = 2*g
  1 - exp(c^(1/h)-( (k*eps)^h + c )^(1/h))
}

# Derivative of loss wrt residual
compute_t_i = function(r, k, c, g){
  h = 2*g
  s = (k*r)^(h) + c
  t = exp(c^(1/h)) * exp(-s^(1/h)) * s^(1/h-1) * k^(h) * r^(h-1)
  return(t)
}

# numDeriv::grad(function(r) loss(r, 0.5, 1e-3, 1), seq(-1, 1, length = 10))
# compute_t_i(seq(-1, 1, length = 10), 0.5, 1e-3, 1)

# Individual contribution to 
compute_m_i = function(y_i, x_i, beta, k, c, g){
  t_i = compute_t_i(y_i - x_i*beta, k, c, g)
  -t_i * x_i
}

# DRAW COVARIATES ----
n = 1e3
true_beta = -0.5
x = runif(n, 0, 2)
eta = x*true_beta

# DEFINE ERROR SAMPLER ----
err_distr = "Gaussian_0.5"
rerr = switch(err_distr,
              Gaussian = function(ndraws) rnorm(ndraws, 0, 1),
              Gaussian_0.5 = function(ndraws) rnorm(ndraws, 0, 0.5),
              Gamma_2_2 = function(ndraws) rgamma(ndraws, 2, 2),
              Gamma_1.5_1.5 = function(ndraws) rgamma(ndraws, 1.5, 1.5),
              Gamma_1.5_3 = function(ndraws) rgamma(ndraws, 1.5, 3))

# DRAW RESPONSE ----
mode_shift = switch(err_distr,
                    Gaussian = 0,
                    Gaussian_0.5 = 0,
                    Gamma_2_2 = (2-1)/2,
                    Gamma_1.5_1.5 = (1.5-1)/1.5,
                    Gamma_1.5_3 = (1.5-1)/3)
y = x*true_beta + rerr(n) - mode_shift

# OPTIMAL K ----
# const = (n^(5/12)) * (log(n)^(7/12)) / (n^(7/12))
# kinit = (n/const)^(1/7)
sd_quant = sd(quantreg::rq(y ~ x, tau = 0.5)$residuals)
# kopt = kinit * (n^(1/7)) / sd
kopt = n^(1/7)/sd_quant

# BOOTSTRAP ESTIMATE OF FISHER INFORMATION ----
B = 1e3
beta_values = seq(-3, 2, length = 101)
y_boot = replicate(B, eta + rerr(n) - mode_shift)
x_boot = matrix(x, nrow = n, ncol = B) # just B replicates of x to be of same shape as y

# SETTINGS ANALYSIS ----
k_val = c(0.5, "opt/2", "opt", "2*opt")
g_val = c(5, 1)
settings = expand.grid(k_val, g_val)

# Create palette of 100 colors from blue to red through white
pal <- colorRampPalette(c("blue", "white", "red"))(101)

setwd("Q_vs_D_NIL")
if (!dir.exists(err_distr)){dir.create(err_distr)}
setwd(err_distr)

pdf("Q_vs_D_c1e3.pdf", width = 12, height = 24)
c = 1e-3
par(mfrow = c(8, 3))
par(cex.lab = 1.5)
for (i in 1:nrow(settings)){
  k_chr = settings[i, 1]
  if (k_chr == "0.5") {k = 0.5} 
  if (k_chr == "opt/2") {k = kopt/2}
  if (k_chr == "opt") {k = kopt}
  if (k_chr == "2*opt") {k = 2*kopt}
  
  g = as.numeric(settings[i, 2])
  
  mbar_boot = sapply(beta_values, function(b)
    mbar = colSums(compute_m_i(y_boot, x_boot, b, k, c, g))/n
  )
  W_estimate = apply(mbar_boot, 2, function(mb) var(sqrt(n)*mb))
  
  # Compare original loss L (or D in Frazier) and Q
  D = sapply(beta_values, function(b) sum(loss(y - b*x, k, c, g)))
  GBI = sapply(seq_along(beta_values), function(j) {
    b = beta_values[j]
    m_i = compute_m_i(y, x, b, k, c, g)
    mbar = (1/n) * sum(m_i)
    W = var(m_i)
    Q = (1/2) * mbar * (1/W) * mbar
    return(c(Q = Q, W = W))
  })
  
  GBI_boot = Q = sapply(seq_along(beta_values), function(j) {
    b = beta_values[j]
    W = W_estimate[j]
    mbar = (1/n) * sum(compute_m_i(y, x, b, k, c, g))
    Q = (1/2) * mbar * (1/W) * mbar
    return(c(Q = Q, W = W))
  })
  
  # Compute likelihoods
  loglik_D = -D
  loglik_GBI = -0.5*log(GBI["W", ]) - n*GBI["Q", ]
  loglik_GBI_boot = -0.5*log(GBI_boot["W", ]) - n*GBI_boot["Q", ]
  
  # Plots likelihoods
  plot(beta_values, loglik_D, xlab = "beta", ylab = "-D",
       main = sprintf("k=%s - c=%s - g=%s - n=%s", k_chr, c, g, n), type = "l")
  abline(v = true_beta, col = 2)
  plot(beta_values, loglik_GBI, xlab = "beta", ylab = "CGBI loglik",
       main = sprintf("k=%s - c=%s - g=%s - n=%s", k_chr, c, g, n), type = "l")
  lines(beta_values, loglik_GBI_boot, col = "green4", lty = 2); legend("bottomright", lty = c(1, 2), col = c(1, "green4"), legend = c("bootstrap", "sample var"))
  abline(v = true_beta, col = 2)
  plot(loglik_D, loglik_GBI_boot, xlab = "-D", ylab = "CGBI loglik", col = pal,
       main = sprintf("k=%s - c=%s - g=%s - n=%s", k_chr, c, g, n))
  legend("bottomright", pch = 1, col = pal[c(1, 25, 75, 101)],
         legend = paste0("beta = ", c(-2, -1, 0, 1)))
}
dev.off()


pdf("Q_vs_D_c1e5.pdf", width = 12, height = 24)
c = 1e-5
par(mfrow = c(8, 3))
par(cex.lab = 1.5)
for (i in 1:nrow(settings)){
  k_chr = settings[i, 1]
  if (k_chr == "0.5") {k = 0.5} 
  if (k_chr == "opt/2") {k = kopt/2}
  if (k_chr == "opt") {k = kopt}
  if (k_chr == "2*opt") {k = 2*kopt}
  
  g = as.numeric(settings[i, 2])
  
  mbar_boot = sapply(beta_values, function(b)
    mbar = colSums(compute_m_i(y_boot, x_boot, b, k, c, g))/n
  )
  W_estimate = apply(mbar_boot, 2, function(mb) var(sqrt(n)*mb))
  
  # Compare original loss L (or D in Frazier) and Q
  D = sapply(beta_values, function(b) sum(loss(y - b*x, k, c, g)))
  GBI = sapply(seq_along(beta_values), function(j) {
    b = beta_values[j]
    m_i = compute_m_i(y, x, b, k, c, g)
    mbar = (1/n) * sum(m_i)
    W = var(m_i)
    Q = (1/2) * mbar * (1/W) * mbar
    return(c(Q = Q, W = W))
  })
  
  GBI_boot = Q = sapply(seq_along(beta_values), function(j) {
    b = beta_values[j]
    W = W_estimate[j]
    mbar = (1/n) * sum(compute_m_i(y, x, b, k, c, g))
    Q = (1/2) * mbar * (1/W) * mbar
    return(c(Q = Q, W = W))
  })
  
  # Compute likelihoods
  loglik_D = -D
  loglik_GBI = -0.5*log(GBI["W", ]) - n*GBI["Q", ]
  loglik_GBI_boot = -0.5*log(GBI_boot["W", ]) - n*GBI_boot["Q", ]
  
  # Plots likelihoods
  plot(beta_values, loglik_D, xlab = "beta", ylab = "-D",
       main = sprintf("k=%s - c=%s - g=%s - n=%s", k_chr, c, g, n), type = "l")
  abline(v = true_beta, col = 2)
  plot(beta_values, loglik_GBI, xlab = "beta", ylab = "CGBI loglik",
       main = sprintf("k=%s - c=%s - g=%s - n=%s", k_chr, c, g, n), type = "l")
  lines(beta_values, loglik_GBI_boot, col = "green4", lty = 2); legend("bottomright", lty = c(1, 2), col = c(1, "green4"), legend = c("bootstrap", "sample var"))
  abline(v = true_beta, col = 2)
  plot(loglik_D, loglik_GBI_boot, xlab = "-D", ylab = "CGBI loglik", col = pal,
       main = sprintf("k=%s - c=%s - g=%s - n=%s", k_chr, c, g, n))
  legend("bottomright", pch = 1, col = pal[c(1, 25, 75, 101)], 
         legend = paste0("beta = ", c(-2, -1, 0, 1)))
}
dev.off()

setwd("../..")

