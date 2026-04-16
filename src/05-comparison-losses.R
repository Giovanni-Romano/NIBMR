library(plotly)
library(htmlwidgets)
source("src/01-utils.R")

# TRUE BETA AND SAMPLE SIZE ----
n = 1e3
true_beta = c(0, -1, +1.5)

# SIMULATION SETTINGS ----
error_distr = c("Gaussian", "Gamma_2_rad2", "Gaussian_rad2", "Gamma_2_1")
covariate = c("Unif_rad2", "Gamma_1.5_1.5", "BVN")
k = "final"
sim_settings = as.matrix(expand.grid(error_distr, covariate, k))
colnames(sim_settings) = c("err_distr", "unif_cov", "k")

# Pick a setting
s = 2
s_s = sim_settings[s, ]
ed_s = unname(s_s[1])
cd_s = unname(s_s[2])
k_s = unname(s_s[3])


# DEFINE ERROR SAMPLER
rerr = switch(ed_s,
              Gaussian = function(ndraws) rnorm(ndraws, 0, 1),
              Gaussian_0.5 = function(ndraws) rnorm(ndraws, 0, 0.5),
              Gaussian_rad2 = function(ndraws) rnorm(ndraws, 0, sqrt(2)),
              Gamma_2_2 = function(ndraws) rgamma(ndraws, 2, 2),
              Gamma_2_rad2 = function(ndraws) rgamma(ndraws, 2, sqrt(2)),
              Gamma_2_1 = function(ndraws) rgamma(ndraws, 2, 1))

mode_shift = switch(ed_s,
                    Gaussian = 0,
                    Gaussian_0.5 = 0,
                    Gaussian_rad2 = 0,
                    Gamma_2_2 = (2-1)/2,
                    Gamma_1.5_1.5 = (1.5-1)/1.5,
                    Gamma_2_1 = (2-1)/1,
                    Gamma_2_rad2 = (2-1)/sqrt(2))

rcov = switch(cd_s,
              "0.5" = function(ndraws) runif(ndraws, -0.5, 0.5),
              "1" = function(ndraws) runif(ndraws, -1, 1),
              asymm = function(ndraws) runif(ndraws, 0, 2),
              "Unif_rad2" = function(ndraws) runif(ndraws, -sqrt(2), sqrt(2)),
              "Gamma_1.5_1.5" = function(ndraws) rgamma(ndraws, 1.5, 1.5),
              "BVN" = function(ndraws) mvnfast::rmvn(ndraws, rep(0, 2), sigma = (26/12)*matrix(c(1, 0.75, 0.75, 1), nrow = 2), isChol = F))

# DRAW COVARIATES
if (cd_s != "BVN"){
  x1 = rcov(n)
  x2 = rcov(n)
  X = cbind(1, x1, x2)
} else {
  X.tmp = rcov(n)
  x1 = X.tmp[ , 1]
  x2 = X.tmp[ , 2]
  X = cbind(1, X.tmp)
}

mu_x1 = mean(x1); sd_x1 = sd(x1)
mu_x2 = mean(x2); sd_x2 = sd(x2)
X_scaled = cbind(1, (x1-mu_x1)/sd_x1, (x2-mu_x2)/sd_x2)
eta = X%*%true_beta

# DRAW RESPONSE
y = eta + rerr(n) - mode_shift


# OPTIMAL K
const = (n^(5/12)) * (log(n)^(7/12)) / (n^(7/12))
kinit = (n/const)^(1/7)
res_qr = quantreg::rq(y ~ -1 + X_scaled, tau = 0.5)$residuals
kfinal = kinit * (n^(1/7))


varp = var(res_qr[res_qr > 0])
varn = var(res_qr[res_qr < 0])
mup = mean(res_qr[res_qr > 0]) - mean(res_qr)
mun = mean(res_qr[res_qr < 0]) - mean(res_qr)

factorp_sd = sqrt(varp)
factorn_sd = sqrt(varn)

factorp_TL = sqrt(varp+mup^2)
factorn_TL = sqrt(varn+mun^2)

k_sd_init = c( 
  kinit/factorn_sd,
  kinit/factorp_sd
)

k_TL_init = c( 
  kinit/factorn_TL,
  kinit/factorp_TL
)

k_sd_final = c( 
  kfinal/factorn_sd,
  kfinal/factorp_sd
)

k_TL_final = c( 
  kfinal/factorn_TL,
  kfinal/factorp_TL
)

w_num = 1/2
LOO_sd_init <- sapply(1:n, function(i) {
  est <- optim(par = true_beta,
               fn = function(b) sum(loss_asymm2(y[-i] - X[-i, ] %*% b, k_sd_init, 1e-3)),
               method = "BFGS")$par
  loss_asymm2(y[i] - drop(X[i, ] %*% est), k_sd_init, 1e-3)
})
w_sd_init = w_num / (sum(LOO_sd_init)/n)

LOO_sd_final <- sapply(1:n, function(i) {
  est <- optim(par = true_beta,
               fn = function(b) sum(loss_asymm2(y[-i] - X[-i, ] %*% b, k_sd_final, 1e-3)),
               method = "BFGS")$par
  loss_asymm2(y[i] - drop(X[i, ] %*% est), k_sd_final, 1e-3)
})
w_sd_final = w_num / (sum(LOO_sd_final)/n)

LOO_TL_init <- sapply(1:n, function(i) {
  est <- optim(par = true_beta,
               fn = function(b) sum(loss_asymm2(y[-i] - X[-i, ] %*% b, k_TL_init, 1e-3)),
               method = "BFGS")$par
  loss_asymm2(y[i] - drop(X[i, ] %*% est), k_TL_init, 1e-3)
})
w_TL_init = w_num / (sum(LOO_TL_init)/n)

LOO_TL_final <- sapply(1:n, function(i) {
  est <- optim(par = true_beta,
               fn = function(b) sum(loss_asymm2(y[-i] - X[-i, ] %*% b, k_TL_final, 1e-3)),
               method = "BFGS")$par
  loss_asymm2(y[i] - drop(X[i, ] %*% est), k_TL_final, 1e-3)
})
w_TL_final = w_num / (sum(LOO_TL_final)/n)


# Optim ----
optim_sd_init = optim(par = c(0.05, -0.8, +1.2),
                      fn = function(b) sum(loss_asymm2(y - X%*%b, 
                                                       k_sd_init, 
                                                       1e-3)),
                      method = "SANN")$par

optim_TL_init = optim(par = c(0.05, -0.8, +1.2),
                      fn = function(b) sum(loss_asymm2(y - X%*%b, 
                                                       k_TL_init, 
                                                       1e-3)),
                      method = "SANN")$par

optim_sd_final = optim(par = c(0.05, -0.8, +1.2),
                       fn = function(b) sum(loss_asymm2(y - X%*%b, 
                                                        k_sd_final, 
                                                        1e-3)),
                       method = "SANN")$par


optim_TL_final = optim(par = c(0.05, -0.8, +1.2),
                       fn = function(b) sum(loss_asymm2(y - X%*%b, 
                                                        k_TL_final, 
                                                        1e-3)),
                       method = "SANN")$par

# Posterior contours ----
beta1_val = seq(-2.5, 0.5, length = 101)
beta2_val = seq(0, 3, length = 101)
beta_val = as.matrix(expand.grid(beta1_val, beta2_val))
post_sd_init_val = apply(beta_val, 1, function(B){
  res = y - X%*%c(optim_sd_init[1], B)
  -w_sd_init*sum(sapply(res, function(r) loss_asymm2(r, k_sd_init, 1e-3)))  + 0.5*sum(B^2)
})

post_TL_init_val = apply(beta_val, 1, function(B){
  res = y - X%*%c(optim_TL_init[1], B)
  -w_TL_init*sum(sapply(res, function(r) loss_asymm2(r, k_TL_init, 1e-3))) + 0.5*sum(B^2)
})

post_sd_final_val = apply(beta_val, 1, function(B){
  res = y - X%*%c(optim_sd_final[1], B)
  -w_sd_final*sum(sapply(res, function(r) loss_asymm2(r, k_sd_final, 1e-3))) + 0.5*sum(B^2)
})

post_TL_final_val = apply(beta_val, 1, function(B){
  res = y - X%*%c(optim_TL_final[1], B)
  -w_TL_final*sum(sapply(res, function(r) loss_asymm2(r, k_TL_final, 1e-3))) + 0.5*sum(B^2)
})

library(plotly)
plt_post_cont = subplot(
  plot_ly(x = beta_val[, 1], y = beta_val[, 2], z = post_sd_init_val,
          type = "contour", showscale = FALSE, showlegend = FALSE) %>%
    add_trace(x = true_beta[2], y = true_beta[3],
              type = "scatter", mode = "markers",
              name = "Truth", showlegend = TRUE,
              marker = list(color = "red", size = 5, symbol = "circle")) %>%
    add_trace(x = optim_sd_init[2], y = optim_sd_init[3],
              type = "scatter", mode = "markers",
              name = "Frequentist estimate",
              showlegend = TRUE,
              marker = list(color = "blue", size = 5, symbol = "x")),
  
  plot_ly(x = beta_val[, 1], y = beta_val[, 2], z = post_sd_final_val,
          type = "contour", showscale = FALSE, showlegend = FALSE) %>%
    add_trace(x = true_beta[2], y = true_beta[3],
              type = "scatter", mode = "markers",
              name = "Truth",
              showlegend = FALSE,
              marker = list(color = "red", size = 5, symbol = "circle")) %>%
    add_trace(x = optim_sd_final[2], y = optim_sd_final[3],
              type = "scatter", mode = "markers",
              name = "Frequentist estimate",
              showlegend = FALSE,
              marker = list(color = "blue", size = 5, symbol = "x")),
  
  plot_ly(x = beta_val[, 1], y = beta_val[, 2], z = post_TL_init_val,
          type = "contour", showscale = FALSE, showlegend = FALSE) %>%
    add_trace(x = true_beta[2], y = true_beta[3],
              type = "scatter", mode = "markers",
              name = "Truth",
              showlegend = FALSE,
              marker = list(color = "red", size = 5, symbol = "circle")) %>%
    add_trace(x = optim_TL_init[2], y = optim_TL_init[3],
              type = "scatter", mode = "markers",
              name = "Frequentist estimate",
              showlegend = FALSE,
              marker = list(color = "blue", size = 5, symbol = "x")),
  
  plot_ly(x = beta_val[, 1], y = beta_val[, 2], z = post_TL_final_val,
          type = "contour", showscale = FALSE, showlegend = FALSE) %>%
    add_trace(x = true_beta[2], y = true_beta[3],
              type = "scatter", mode = "markers",
              name = "Truth",
              showlegend = FALSE,
              marker = list(color = "red", size = 5, symbol = "circle")) %>%
    add_trace(x = optim_TL_final[2], y = optim_TL_final[3],
              type = "scatter", mode = "markers",
              name = "Frequentist estimate",
              showlegend = FALSE,
              marker = list(color = "blue", size = 5, symbol = "x")),
  nrows = 2, shareX = TRUE, shareY = TRUE
)

layout(plt_post_cont,
       title = paste0("Posterior Contours: [sd top, TL bottom] + [init left, final right] | ", sim_settings[s, 1], " - ", sim_settings[s, 2]),
       showlegend = TRUE)



# Posterior surfaces ----
plt_post_surf1 = plot_ly(x = unique(beta_val[, 1]), y = unique(beta_val[, 2]),
        z = matrix(post_TL_init_val,
                   nrow = length(unique(beta_val[, 2])),
                   ncol = length(unique(beta_val[, 1])),
                   byrow = TRUE),
        type = "surface", showscale = FALSE, showlegend = FALSE) %>%
  layout(title = paste0("Posterior surface [TL init] | ", sim_settings[s, 1], " - ", sim_settings[s, 2]), 
         scene = list(
           xaxis = list(title = "beta 1"),
           yaxis = list(title = "beta 2"),
           zaxis = list(title = "posterior")
         ))

# saveWidget(plt_post_surf1, "sim_study_nonC-GBI/02-w_unit_loss_SNR/post_surf_TL_init_Gauss_UnifRad2.html", selfcontained = TRUE)
saveWidget(plt_post_surf1, "sim_study_nonC-GBI/02-w_unit_loss_SNR/theoretical_comparisons/post_surf_TL_init_Gamma_UnifRad2.html", selfcontained = TRUE)

plt_post_surf2 = plot_ly(x = unique(beta_val[, 1]), y = unique(beta_val[, 2]),
        z = matrix(post_TL_final_val,
                   nrow = length(unique(beta_val[, 2])),
                   ncol = length(unique(beta_val[, 1])),
                   byrow = TRUE),
        type = "surface", showscale = FALSE, showlegend = FALSE) %>%
  layout(title = paste0("Posterior surface [TL final] |", sim_settings[s, 1], " - ", sim_settings[s, 2]),
         scene = list(
           xaxis = list(title = "beta 1"),
           yaxis = list(title = "beta 2"),
           zaxis = list(title = "posterior")
         ))

# saveWidget(plt_post_surf2, "sim_study_nonC-GBI/02-w_unit_loss_SNR/post_surf_TL_final_Gauss_UnifRad2.html", selfcontained = TRUE)
saveWidget(plt_post_surf2, "sim_study_nonC-GBI/02-w_unit_loss_SNR/theoretical_comparisons/post_surf_TL_final_Gamma_UnifRad2.html", selfcontained = TRUE)

plt_post_surf3 = plot_ly(x = unique(beta_val[, 1]), y = unique(beta_val[, 2]),
        z = matrix(post_sd_init_val,
                   nrow = length(unique(beta_val[, 2])),
                   ncol = length(unique(beta_val[, 1])),
                   byrow = TRUE),
        type = "surface", showscale = FALSE, showlegend = FALSE) %>%
  layout(title = paste0("Posterior surface [sd init] | ", sim_settings[s, 1], " - ", sim_settings[s, 2]),
         scene = list(
           xaxis = list(title = "beta 1"),
           yaxis = list(title = "beta 2"),
           zaxis = list(title = "posterior")
         ))

# saveWidget(plt_post_surf3, "sim_study_nonC-GBI/02-w_unit_loss_SNR/post_surf_sd_init_Gauss_UnifRad2.html", selfcontained = TRUE)
saveWidget(plt_post_surf3, "sim_study_nonC-GBI/02-w_unit_loss_SNR/theoretical_comparisons/post_surf_sd_init_Gamma_UnifRad2.html", selfcontained = TRUE)

plt_post_surf4 = plot_ly(x = unique(beta_val[, 1]), y = unique(beta_val[, 2]),
        z = matrix(post_sd_final_val,
                   nrow = length(unique(beta_val[, 2])),
                   ncol = length(unique(beta_val[, 1])),
                   byrow = TRUE),
        type = "surface", showscale = FALSE, showlegend = FALSE) %>%
  layout(title = paste0("Posterior surface [sd final] | ", sim_settings[s, 1], " - ", sim_settings[s, 2]),
         scene = list(
           xaxis = list(title = "beta 1"),
           yaxis = list(title = "beta 2"),
           zaxis = list(title = "posterior")
         ))

# saveWidget(plt_post_surf4, "sim_study_nonC-GBI/02-w_unit_loss_SNR/post_surf_sd_final_Gauss_UnifRad2.html", selfcontained = TRUE)
saveWidget(plt_post_surf4, "sim_study_nonC-GBI/02-w_unit_loss_SNR/theoretical_comparisons/post_surf_sd_final_Gamma_UnifRad2.html", selfcontained = TRUE)

# Curvature ----
## Hessian ----
library(parallel)
n_i <- 1000
n_b <- nrow(beta_val)
beta_aug <- cbind(0, beta_val)

## Final ----
H_sd_final.tmp <- mclapply(seq_len(n_i), function(i) {
  tmp <- vapply(
    seq_len(n_b),
    function(j) c(H_i(beta_aug[j, ], y[i], X[i, ], k_sd_final, 1e-3)),
    numeric(9)
  )
  array(tmp, dim = c(3, 3, n_b))
}, mc.cores = parallel::detectCores() - 1)

H_sd_final <- array(unlist(H_sd_final.tmp), dim = c(3, 3, n_b, n_i)) %>% aperm(c(4, 1, 2, 3)) %>% colMeans

H_TL_final.tmp <- mclapply(seq_len(n_i), function(i) {
  tmp <- vapply(
    seq_len(n_b),
    function(j) c(H_i(beta_aug[j, ], y[i], X[i, ], k_TL_final, 1e-3)),
    numeric(9)
  )
  array(tmp, dim = c(3, 3, n_b))
}, mc.cores = parallel::detectCores() - 1)

H_TL_final <- array(unlist(H_TL_final.tmp), dim = c(3, 3, n_b, n_i)) %>% aperm(c(4, 1, 2, 3)) %>% colMeans

## Init ----
H_sd_init.tmp <- mclapply(seq_len(n_i), function(i) {
  tmp <- vapply(
    seq_len(n_b),
    function(j) c(H_i(beta_aug[j, ], y[i], X[i, ], k_sd_init, 1e-3)),
    numeric(9)
  )
  array(tmp, dim = c(3, 3, n_b))
}, mc.cores = parallel::detectCores() - 1)

H_sd_init <- array(unlist(H_sd_init.tmp), dim = c(3, 3, n_b, n_i)) %>% aperm(c(4, 1, 2, 3)) %>% colMeans

H_TL_init.tmp <- mclapply(seq_len(n_i), function(i) {
  tmp <- vapply(
    seq_len(n_b),
    function(j) c(H_i(beta_aug[j, ], y[i], X[i, ], k_TL_init, 1e-3)),
    numeric(9)
  )
  array(tmp, dim = c(3, 3, n_b))
}, mc.cores = parallel::detectCores() - 1)

H_TL_init <- array(unlist(H_TL_init.tmp), dim = c(3, 3, n_b, n_i)) %>% aperm(c(4, 1, 2, 3)) %>% colMeans

## Plots ----
### H beta1 ----
s_H_beta1 = subplot(
  plot_ly(x = beta_val[ , 1], y = beta_val[, 2], z = H_sd_init[2, 2, ], type = "contour", showscale = F),
  plot_ly(x = beta_val[ , 1], y = beta_val[, 2], z = H_sd_final[2, 2, ], type = "contour", showscale = F),
  plot_ly(x = beta_val[ , 1], y = beta_val[, 2], z = H_TL_init[2, 2, ], type = "contour", showscale = F),
  plot_ly(x = beta_val[ , 1], y = beta_val[, 2], z = H_TL_final[2, 2, ], type = "contour", showscale = F),
  nrows=2, shareX=T, shareY=T
)

layout(s_H_beta1, title = "Hessian beta1-beta1: [sd top, TL bottom] + [init left, final right]")

### H beta2 ----
s_H_beta2 = subplot(
  plot_ly(x = beta_val[ , 1], y = beta_val[, 2], z = H_sd_init[3, 3, ], type = "contour", showscale = F),
  plot_ly(x = beta_val[ , 1], y = beta_val[, 2], z = H_sd_final[3, 3, ], type = "contour", showscale = F),
  plot_ly(x = beta_val[ , 1], y = beta_val[, 2], z = H_TL_init[3, 3, ], type = "contour", showscale = F),
  plot_ly(x = beta_val[ , 1], y = beta_val[, 2], z = H_TL_final[3, 3, ], type = "contour", showscale = F),
  nrows=2, shareX=T, shareY=T
)

layout(s_H_beta2, title = "Hessian beta2-beta2: [sd top, TL bottom] + [init left, final right]")

### Eigenvalues ----
ev_sd_final = apply(H_sd_final, 3, function(M) eigen(M)$values)
ev_TL_final = apply(H_TL_final, 3, function(M) eigen(M)$values)
ev_sd_init = apply(H_sd_init, 3, function(M) eigen(M)$values)
ev_TL_init = apply(H_TL_init, 3, function(M) eigen(M)$values)

s_ev = subplot(
  plot_ly(x = beta_val[ , 1], y = beta_val[, 2], z = ev_sd_init[1, ], type = "contour", showscale = F),
  plot_ly(x = beta_val[ , 1], y = beta_val[, 2], z = ev_sd_init[2, ], type = "contour", showscale = F),
  plot_ly(x = beta_val[ , 1], y = beta_val[, 2], z = ev_sd_init[3, ], type = "contour", showscale = F),
  plot_ly(x = beta_val[ , 1], y = beta_val[, 2], z = ev_TL_init[1, ], type = "contour", showscale = F),
  plot_ly(x = beta_val[ , 1], y = beta_val[, 2], z = ev_TL_init[2, ], type = "contour", showscale = F),
  plot_ly(x = beta_val[ , 1], y = beta_val[, 2], z = ev_TL_init[3, ], type = "contour", showscale = F),
  nrows=2, shareX=T, shareY=T
)

s_ev

### Convexity ----
convex_sd_init = apply(ev_sd_init, 2, function(vals) {
  if (all(vals>0)){
    return(+1)
  } else if (all(vals<0)){
    return(-1)
  } else {
    return(0)
  }
})
convex_TL_init = apply(ev_TL_init, 2, function(vals) {
  if (all(vals>0)){
    return(+1)
  } else if (all(vals<0)){
    return(-1)
  } else {
    return(0)
  }
})
convex_sd_final = apply(ev_sd_final, 2, function(vals) {
  if (all(vals>0)){
    return(+1)
  } else if (all(vals<0)){
    return(-1)
  } else {
    return(0)
  }
})
convex_TL_final = apply(ev_TL_final, 2, function(vals) {
  if (all(vals>0)){
    return(+1)
  } else if (all(vals<0)){
    return(-1)
  } else {
    return(0)
  }
})

disc_scale <- list(
  c(0/3, "blue"),  c(1/3, "blue"),
  c(1/3, "white"), c(2/3, "white"),
  c(2/3, "red"),   c(3/3, "red")
)

s_convex = subplot(
  plot_ly(x = beta_val[, 1], y = beta_val[, 2], z = convex_sd_init,
          type = "heatmap",
          zmin = -1.5, zmax = 1.5,
          colorscale = disc_scale,
          colorbar = list(title = "Value", tickvals = c(-1, 0, 1), ticktext = c("concave", "saddle", "convex")),
          showscale = TRUE),
  plot_ly(x = beta_val[, 1], y = beta_val[, 2], z = convex_sd_final,
          type = "heatmap",
          zmin = -1.5, zmax = 1.5,
          colorscale = disc_scale, showscale = FALSE),
  plot_ly(x = beta_val[, 1], y = beta_val[, 2], z = convex_TL_init,
          type = "heatmap",
          zmin = -1.5, zmax = 1.5,
          colorscale = disc_scale, showscale = FALSE),
  plot_ly(x = beta_val[, 1], y = beta_val[, 2], z = convex_TL_final,
          type = "heatmap",
          zmin = -1.5, zmax = 1.5,
          colorscale = disc_scale, showscale = FALSE),
  nrows = 2, shareX = TRUE, shareY = TRUE
)

s_convex

# # Fisher Info ----
# score = mclapply(seq_len(n_i), function(i) {
#   tmp <- vapply(
#     seq_len(n_b),
#     function(j) score_i(beta_aug[j, ], y[i], X[i, ], k, 1e-3),
#     numeric(3)
#   )
#   array(tmp, dim = c(3, n_b))
# }, mc.cores = parallel::detectCores() - 1)
# 
# score = res_score %>% abind::abind(along = 3) %>% aperm(c(3, 1, 2))
# 
# Fisher = apply(score, 3, function(M) {crossprod(M) / nrow(M)}, simplify = F) %>% abind::abind(along = 3)
# 
# 
