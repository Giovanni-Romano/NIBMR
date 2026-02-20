compute_t = function(r, k, c){
  s = sqrt(k^2 * r^2 + c)
  t = exp(-s)/s * k^2 * r
  return(t)
}


GBI_post = function(beta, data, k, c){
  y = data$y
  x = data$x
  n = length(x)
  
  r = y - beta*x
  
  mbar_ind = compute_t(r, k, c)*x
  mbar = -sum(mbar_ind)/n
  
  W = var(mbar_ind)*(n-1)/n
  
  Q = 0.5*mbar^2/W
  
  loglik = -0.5*log(W) - n*Q
  logprior = - 0.5*beta^2
  post = exp(loglik + logprior)
  
  return(post)
}

n = 10^4
x = runif(n, 0, 2)
true_beta = 2
y = x * true_beta + rnorm(n, 0, 1)
beta_val = seq(true_beta-2, true_beta+2, length = 1e4)
c = 1e-3
const = (n^(5/12)) * (log(n)^(7/12)) / (n^(7/12))
kopt = (n/const)^(1/7)

k_val = c(1, kopt, 10)
par(mfrow = c(1, 3))
for (j in seq_along(k_val)){
  k = k_val[j]
  GBI_post_val = sapply(beta_val, GBI_post, data = list(y = y, x = x), k = k, c = 1e-3)
  plot(beta_val, GBI_post_val, type = "l", col = j, main = sprintf("k = %.2f, n = %s", k, n))
  abline(v = true_beta, col = "gold")
}
