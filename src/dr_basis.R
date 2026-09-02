## this function returns the nonlinear Demmler-Reinsch spline basis functions
# Construct Demmler Reinsch basis according to Bach and Klein 2024
construct_DR_basis <- function(x, d = 10, m = 3, penalty = 2, rescale = TRUE) {
  
  require(mgcv)
  
  n <- length(x)
  # set up P-spline design
  fit=smoothCon(s(x,bs="ps",k=d+2),data=data.frame(x),scale.penalty=FALSE)
  
  B=fit[[1]]$X
  K=fit[[1]]$S[[1]]
  
  # set up constraint matrix
  C <- t(cbind(1, scale(x))) %*% B
  V0 <- svd(C,nv=d+2)$v[, 3:(d+2), drop = FALSE]
  
  # transformed penalty and gram matrices
  K_tilde <- t(V0) %*% K %*% V0
  G_tilde <- t(V0) %*% (t(B) %*% B / n) %*% V0
  # solve generalized eigenvalue problem
  eig <- geigen::geigen(G_tilde, K_tilde)
  A <- eig$vectors
  Trf <- V0 %*% A
  
  # perform change of basis
  Z <- B %*% Trf
  # normalize to unit norm
  scaling_factor <- 1
  if (rescale) {
    scaling_factor <- sqrt(n / sum(Z^2))
    Z <- Z * scaling_factor
    Trf <- Trf * scaling_factor
  }
  # you don't really need the full list but I returned everything and kept what I needed
  list(Z = Z, Trafo = Trf, Gram = G_tilde, Penalty = K_tilde)
}

# Example
set.seed(123)
# sample size
n <- 400

# covariate
x <- runif(n, 0, 1)

# true nonlinear function
f_true <- sin(2*pi*x)
flin=(-6/pi)*(x-1/2)
fnonlin=f_true-flin 

# noise level
sigma <- 0.2

# response
y <- f_true + rnorm(n, 0, sigma)

dr <- construct_DR_basis(x, d = 10)

Z <- dr$Z   # nonlinear DR basis
X <- cbind(1, scale(x))

fit <- lm(y ~ X + Z - 1)

coef_all <- coef(fit)

print(coef_all)

beta_hat <- coef_all[1:2]
gamma_hat <- coef_all[-(1:2)]

# fitted values
y_hat <- X %*% beta_hat + Z %*% gamma_hat

# diagnostics
par(mfrow = c(2, 2))
plot(y, y_hat)
abline(0, 1, col = "red", lwd = 2)

ord <- order(x)
plot(x[ord], y[ord], col = "grey")
lines(x[ord], y_hat[ord], col = "red", lwd = 2)
lines(x[ord], f_true[ord], col = "blue", lwd = 2)

plot(x[ord], as.vector(Z%*%gamma_hat)[ord])
lines(x[ord], fnonlin[ord], col = "red", lwd = 2)

plot(x[ord], (X %*% beta_hat)[ord])
lines(x[ord], flin[ord], col = "red", lwd = 2)
