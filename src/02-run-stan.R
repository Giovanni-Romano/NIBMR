rm(list = ls())

n = 10^3
p = 1
x1 = runif(n, 0, 2)
# x2 = runif(n, 0, 2)
# X_tmp = matrix(c(x1, x2), nrow = n, ncol = p, byrow = F)
X_tmp = matrix(x1, nrow = n)
X = X_tmp
# X = apply(X_tmp, 2, scale)
# beta = c(+0.2, -2)
beta = array(c(-2), dim = c(1L))
y = c(X %*% beta) + rgamma(n, shape = 2, rate = 1) - 1#+ rnorm(n, mean = 0, sd = 1)

const = (n^(5/12)) * (log(n)^(7/12)) / (n^(7/12))
kinit = (n/const)^(1/7)
sd = sd(quantreg::rq(y ~ x1, tau = 0.5)$residuals)
kopt = kinit * (n^(1/7)) / sd

# Loss function
loss = function(eps, k, c, g = 1){
  h = 2*g
  1 - exp(c^(1/h)-( (k*eps)^h + c )^(1/h) )
}

sapply(c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"), 
       function(m)
         optim(par = 0, 
               fn = function(b) mean(loss(y - b*x1, kopt/2, 1e-3)),
               method = m)$par
)


library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

c = 1e-3 
k = kopt/4

dat_cgbi <- list(
  n = nrow(X),
  # J = 2L,
  # q = array(c(1L, 1L), dim = c(2L)),
  J = 1L,
  q = array(c(1L), dim = c(1L)),
  p = p,
  B = X,
  y = as.vector(y),
  k = k,
  c = c
)

sm_cgbi <- stan_model(file = "src/01-utils-stan.stan")

init_fun <- function(...) list("beta" = beta)

fit_cgbi <- sampling(
  object = sm_cgbi,
  data   = dat_cgbi,
  chains = 8,
  iter   = 2000,
  warmup = 1000,
  seed   = 1,
  init = init_fun
  # ,
  # control = list(adapt_delta = 0.99, max_treedepth = 12)
)

print(fit_cgbi, pars = "beta", probs = c(0.1, 0.5, 0.9))

rstan::check_hmc_diagnostics(fit_cgbi)

sum <- summary(fit_cgbi)$summary

rstan::stan_trace(fit_cgbi, pars = c("lp__", "beta[1]")) 

beta_post = extract(fit_cgbi)

plot(x1, y)
points(x1, colMeans(beta_post$beta) * x1, col = "2")
lines(x1, c(beta) * x1, col = 4)
