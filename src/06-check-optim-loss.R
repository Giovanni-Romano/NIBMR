n = 1e3
x = runif(n, 0, 2)
true_beta = 2
y = x * true_beta + rnorm(n, 0, 1)

meth = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN")
c = 1e-3
const = n^(5/12)*(log(n)^(7/12))/(n^(7/12))
koptim = (n/const)^(1/7)


kval = c(koptim, 1e0, 1e1, 1e2, 1e3, 1e4)
cval = c(1e-1, 1e-3, 1e-5)
par_val = expand.grid(kval, cval, meth)

optim_out = apply(par_val, 1, function(p){
  k = as.numeric(p[1])
  c = as.numeric(p[2])
  optim(par = 0, 
        fn = function(beta) mean(loss(y - beta*x, k, c)),
        method = p[3])$par
})


library(ggplot2)
library(dplyr)

as.data.frame(cbind(par_val, optim_out)) %>% 
  ggplot() + 
  geom_point(aes(y = optim_out, x = factor(round(Var1, 2)), col = Var3)) +
  geom_hline(aes(yintercept = true_beta), col = "gold") +
  labs(x = "k", y = "estimate optim") +
  facet_wrap(~ Var2)
