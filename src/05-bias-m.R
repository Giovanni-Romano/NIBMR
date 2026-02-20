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
  -x_i*t_i
}

n = 1e3
true_beta = 2
x = runif(n, 0, 2)
eta = x*true_beta

# DEFINE ERROR SAMPLER ----
err_distr = "Gaussian"
rerr = switch(err_distr,
              Gaussian = function(ndraws) rnorm(ndraws, 0, 1),
              Gamma_2_2 = function(ndraws) rgamma(ndraws, 2, 2),
              Gamma_2_1 = function(ndraws) rgamma(ndraws, 2, 1))

# DRAW RESPONSE ----
mode_shift = switch(err_distr,
                    Gaussian = 0,
                    Gamma_2_2 = (2-1)/2,
                    Gamma_2_1 = (2-1)/1)



m_n = replicate(1e3, {
  y = eta + rerr(n) - mode_shift
  sd_qr = sd(quantreg::rq(y~x)$residuals)
  mean(sapply(1:n, function(i) compute_m_i(y[i], x[i], true_beta, n^(1/7)*sd_qr, 1e-3, 1)))
})

summary(m_n)
hist(m_n, breaks = 50)