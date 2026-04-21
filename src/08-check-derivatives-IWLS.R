
source("src/01-utils.R")
rerr = switch("Gaussian",
              Gaussian = function(ndraws, var_eta, SNR) rnorm(ndraws, 0, sd = sqrt(var_eta/SNR)),
              Gamma = function(ndraws, var_eta, SNR) rgamma(ndraws, 2, sqrt(2*SNR/var_eta)))

rcov = switch("Beta_2_4",
              "Unif_rad2" = function(ndraws) runif(ndraws, -sqrt(2), sqrt(2)),
              "Beta_2_4" = function(ndraws) 2*sqrt(2)*(rbeta(ndraws, 2, 4)-0.5))

transf = switch("parabola",
                "parabola" = function(x) x^2,
                "cubic" = function(x) x^3,
                "trigonometric" = function(x) sin(2.5*x) + 3*exp(-(5^2)*((x)^2))
)

n = 1e3
x1 = rcov(n)
mu_x1 = mean(x1); sd_x1 = sd(x1)
eta = transf(x1)
var_eta = var(eta)

# Mode shift
mode_shift = switch("Gamma",
                    Gaussian = 0,
                    Gamma = (2 - 1) / sqrt(2*SNR/var_eta))

# DRAW RESPONSE
y = eta + rerr(n, var_eta, SNR = 5) - mode_shift

# CREATE DESIGN MATRIX
d = 20
DR <- construct_DR_basis(x1, d = d)
Z <- DR$Z   # nonlinear DR basis
X_scaled <- cbind(1, (x1 - mu_x1)/sd_x1, Z)
outerX = aperm(simplify2array(lapply(1:nrow(X_scaled), function(i) tcrossprod(X_scaled[i, ]))), c(3, 1, 2))

const = (n^(5/12)) * (log(n)^(7/12)) / (n^(7/12))
kinit = (n/const)^(1/7)
kfinal = kinit * (n^(1/7))
res_qr = quantreg::rq(y ~ -1 + X_scaled, tau = 0.5)$residuals
k = kinit

varp = var(res_qr[res_qr > 0])
varn = var(res_qr[res_qr < 0])
mup = mean(res_qr[res_qr > 0]) - mean(res_qr)
mun = mean(res_qr[res_qr < 0]) - mean(res_qr)
factorp = sqrt(varp+mup^2)
factorn = sqrt(varn+mun^2)

k_multi = c( 
  k/factorn,
  k/factorp
)

out_optim_scaled = sapply(c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"),
                          function(m){
                            est = optim(par = rep(0, ncol(X_scaled)),
                                        fn = function(b) sum(loss_asymm2(y - X_scaled%*%b, 
                                                                         k_multi, 
                                                                         1e-3)),
                                        method = m)$par
                            est
                          }
                          
)

beta_val = seq(0.35, 1, length = 100)
ll_val1 = sapply(beta_val, function(b) loss_asymm2_pop(c(out_optim_scaled[-22, "SANN"], b), y, X_scaled, c(kinit, kinit), 1e-3))
score_val1 = sapply(beta_val, function(b) score_pop(c(out_optim_scaled[-22, "SANN"], b), y, X_scaled, c(kinit, kinit), 1e-3))
H_val1 = lapply(beta_val, function(b) H_pop(c(out_optim_scaled[-22, "SANN"], b), y, X_scaled, outerX, c(kinit, kinit), 1e-3)) %>% simplify2array()
score_val1_ND = sapply(beta_val, function(b) numDeriv::grad(function(beta) -loss_asymm2_pop(beta, y, X_scaled, c(kinit, kinit), 1e-3), c(out_optim_scaled[-22, "SANN"], b)))
H_val1_ND = lapply(beta_val, function(b) numDeriv::hessian(function(beta) -loss_asymm2_pop(beta, y, X_scaled, c(kinit, kinit), 1e-3), c(out_optim_scaled[-22, "SANN"], b))) %>% simplify2array()

ll_val2 = sapply(beta_val, function(b) loss_asymm2_pop(c(out_optim_scaled[-22, "SANN"], b), y, X_scaled, c(1, 1), 1e-3))
score_val2 = sapply(beta_val, function(b) score_pop(c(out_optim_scaled[-22, "SANN"], b), y, X_scaled, c(1, 1), 1e-3))
H_val2 = lapply(beta_val, function(b) H_pop(c(out_optim_scaled[-22, "SANN"], b), y, X_scaled, outerX, c(1, 1), 1e-3)) %>% simplify2array()
score_val2_ND = sapply(beta_val, function(b) numDeriv::grad(function(beta) -loss_asymm2_pop(beta, y, X_scaled, c(1, 1), 1e-3), c(out_optim_scaled[-22, "SANN"], b)))
H_val2_ND = lapply(beta_val, function(b) numDeriv::hessian(function(beta) -loss_asymm2_pop(beta, y, X_scaled, c(1, 1), 1e-3), c(out_optim_scaled[-22, "SANN"], b))) %>% simplify2array()

ll_val3 = sapply(beta_val, function(b) loss_asymm2_pop(c(out_optim_scaled[-22, "SANN"], b), y, X_scaled, c(kinit, kinit), 1e-1))
score_val3 = sapply(beta_val, function(b) score_pop(c(out_optim_scaled[-22, "SANN"], b), y, X_scaled, c(kinit, kinit), 1e-1))
H_val3 = lapply(beta_val, function(b) H_pop(c(out_optim_scaled[-22, "SANN"], b), y, X_scaled, outerX, c(kinit, kinit), 1e-1)) %>% simplify2array()
score_val3_ND = sapply(beta_val, function(b) numDeriv::grad(function(beta) -loss_asymm2_pop(beta, y, X_scaled, c(kinit, kinit), 1e-1), c(out_optim_scaled[-22, "SANN"], b)))
H_val3_ND = lapply(beta_val, function(b) numDeriv::hessian(function(beta) -loss_asymm2_pop(beta, y, X_scaled, c(kinit, kinit), 1e-1), c(out_optim_scaled[-22, "SANN"], b))) %>% simplify2array()

ll_val4 = sapply(beta_val, function(b) loss_asymm2_pop(c(out_optim_scaled[-22, "SANN"], b), y, X_scaled, c(1, 1), 1e-1))
score_val4 = sapply(beta_val, function(b) score_pop(c(out_optim_scaled[-22, "SANN"], b), y, X_scaled, c(1, 1), 1e-1))
H_val4 = lapply(beta_val, function(b) H_pop(c(out_optim_scaled[-22, "SANN"], b), y, X_scaled, outerX, c(1, 1), 1e-1)) %>% simplify2array()
score_val4_ND = sapply(beta_val, function(b) numDeriv::grad(function(beta) -loss_asymm2_pop(beta, y, X_scaled, c(1, 1), 1e-1), c(out_optim_scaled[-22, "SANN"], b)))
H_val4_ND = lapply(beta_val, function(b) numDeriv::hessian(function(beta) -loss_asymm2_pop(beta, y, X_scaled, c(1, 1), 1e-1), c(out_optim_scaled[-22, "SANN"], b))) %>% simplify2array()

# ll_val5 = sapply(beta_val, function(b) loss_asymm2_pop(c(out_optim_scaled[-22, "SANN"], b), y, X_scaled, c(kinit, kinit), 1e-5))
# score_val5 = sapply(beta_val, function(b) score_pop(c(out_optim_scaled[-22, "SANN"], b), y, X_scaled, c(kinit, kinit), 1e-5))
# H_val5 = lapply(beta_val, function(b) H_pop(c(out_optim_scaled[-22, "SANN"], b), y, X_scaled, outerX, c(kinit, kinit), 1e-5)) %>% simplify2array()
# score_val5_ND = sapply(beta_val, function(b) numDeriv::grad(function(beta) -loss_asymm2_pop(beta, y, X_scaled, c(kinit, kinit), 1e-5), c(out_optim_scaled[-22, "SANN"], b)))
# H_val5_ND = lapply(beta_val, function(b) numDeriv::hessian(function(beta) -loss_asymm2_pop(beta, y, X_scaled, c(kinit, kinit), 1e-5), c(out_optim_scaled[-22, "SANN"], b))) %>% simplify2array()

pdf(file = "sim_study_nonC-GBI/splines_univ/checkDeriv_IWLS.pdf", width = 12, height = 6)
par(mfrow = c(1, 3))
plot(beta_val, -ll_val1 + min(ll_val1), type = "l",
     ylab = "-logloss"); abline(v = out_optim_scaled[22, "SANN"], col = "gold")
lines(beta_val, -ll_val2 + min(ll_val2), col = 2)
lines(beta_val, -ll_val3 + min(ll_val3), col = 3)
lines(beta_val, -ll_val4 + min(ll_val4), col = 4)
# lines(beta_val, -ll_val5 + min(ll_val5), col = 6)
legend("bottomleft", 
       legend = c("k=init(2.69), c=1e-3", "k=init(2.69), c=1e-1", "k=1, c=1e-3", "k=1, c=1e-1"),
       col = c(1, 3, 2, 4), lty = 1, pch = 4)

plot(beta_val, score_val1[22, ], type = "l",
     ylab = "gradient"); abline(v = out_optim_scaled[22, "SANN"], col = "gold"); abline(h = 0, col = "gold")
points(beta_val, score_val1_ND[22, ], col = 1, pch = 4)
lines(beta_val, score_val2[22, ], col = 2); points(beta_val, score_val2_ND[22, ], col = 2, pch = 4)
lines(beta_val, score_val3[22, ], col = 3); points(beta_val, score_val3_ND[22, ], col = 3, pch = 4)
lines(beta_val, score_val4[22, ], col = 4); points(beta_val, score_val4_ND[22, ], col = 4, pch = 4)
# lines(beta_val, score_val5[22, ], col = 6); points(beta_val, score_val5_ND[22, ], col = 6, pch = 4)
legend("bottomleft", 
       legend = c("k=init(2.69), c=1e-3", "k=init(2.69), c=1e-1", "k=1, c=1e-3", "k=1, c=1e-1"),
       col = c(1, 3, 2, 4), lty = 1, pch = 4)

plot(beta_val, H_val1[22, 22, ], type = "l",
     ylab = "Hessian"); abline(v = out_optim_scaled[22, "SANN"], col = "gold")
points(beta_val, H_val1_ND[22, 22, ], col = 1, pch = 4)
lines(beta_val, H_val2[22, 22, ], col = 2); points(beta_val, H_val2_ND[22, 22, ], col = 2, pch = 4)
lines(beta_val, H_val3[22, 22, ], col = 3); points(beta_val, H_val3_ND[22, 22, ], col = 3, pch = 4)
lines(beta_val, H_val4[22, 22, ], col = 4); points(beta_val, H_val4_ND[22, 22, ], col = 4, pch = 4)
# lines(beta_val, H_val5[22, 22, ], col = 6); points(beta_val, H_val5_ND[22, 22, ], col = 6, pch = 4)
legend("topleft", 
       legend = c("k=init(2.69), c=1e-3", "k=init(2.69), c=1e-1", "k=1, c=1e-3", "k=1, c=1e-1"),
       col = c(1, 3, 2, 4), lty = 1, pch = 4)
dev.off()

