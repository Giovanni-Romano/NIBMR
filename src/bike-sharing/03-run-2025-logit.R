rm(list = ls())

# PACKAGES AND SOURCE  ----
suppressPackageStartupMessages(library(tidyverse))
source("src/01a-utils.R")
source("src/01b-gibbs.R")

# LOAD DATA AND FILTER ----
## Load data ----
data = readRDS("data/bike-sharing/milan_bikemi_data_for_regression.rds") |> 
  mutate(Weekend = as.numeric(weekdays(commit_at) %in% c("Saturday", "Sunday")),
         Rain = as.numeric(Precipitation > 0))
## Filter 2025 and no bikes ----
nobikes2025 = data |> filter(year(commit_at) == 2025 & unavail_type == "no_bikes")

# Exclude stations close to Pero: 334, 500 and 501. Actually 500 does not have any observation in 2025.
D = nobikes2025 |> filter(ID != 334, ID != 501) |> 
  select(-commit_at, -Precipitation, -Temperature, -Rain)


# CREATE OBJECTS FOR REGRESSION ----
## Response ----
y = D$time_numeric_shifted / 100 / 720
## Binary covariates ----
X_dummy = (D |> select(Weekend)  |> as.matrix())
## Continuous covariates ----
X_cont = D |> select(population_density, starts_with("density_"), starts_with("closest_")) |> as.matrix()
X_cont.scaled = scale(X_cont)
## ID for rand. eff. ----
id_re = as.integer(factor(D$station, levels = unique(D$station[order(D$ID)])))

# Dimensions ----
n = length(y)
p = ncol(X_cont)
q = ncol(X_dummy)
d = rep(5, p)

## Splines for continuous ----
splines_list = list()
for (j in 1:p){
  splines_list[[j]] = construct_DR_basis(X_cont.scaled[ , j], d = d[j], type_of_splines = "gps")
}
Z_cont = do.call(cbind, lapply(splines_list, function(x) x$Z))
colnames(Z_cont) = c(sapply(1:p, function(j) paste0(colnames(X_cont)[j], "_DR", 1:d[j])))
## Spatial component ----
X_spat = D |> select(station, longitude, latitude)
# I build the tensor-product splines on the unique values of the coordinates, 
# then I match the rows of the design matrix with the original data
X_spat_unique = X_spat |> distinct()
coord_utm = sf::st_as_sf(X_spat_unique, 
                         coords = c("longitude", "latitude"), crs = 4326) |> 
  sf::st_transform(crs = 32632) |> 
  sf::st_coordinates()

# # Option 1: Normalize with observed observed min and max x and y coordinates
coord_scaled = apply(coord_utm, 2, function(x) (x - min(x)) /  (max(x) - min(x)))

# # Option 2: Normalize with Milan min and max x and y coordinates
# library(sf)
# milano <- st_read("https://dati.comune.milano.it/dataset/e75d91fa-eca6-4ee5-b96e-08bcdbb8d6f0/resource/f56cb432-83e6-48de-ae30-d39b4be61e85/download/confine_comune_milano_layer_0_confine_comune_milano.geojson")
# milano_utm <- st_transform(milano, crs = 32632)
# milano_coords <- st_coordinates(milano_utm)
# min_x_milano = min(milano_coords[, 1]); max_x_milano = max(milano_coords[, 1])
# min_y_milano = min(milano_coords[, 2]); max_y_milano = max(milano_coords[, 2])
# coord_scaled = coord_utm
# coord_scaled[ , 1] = (coord_utm[ , 1] - min_x_milano) / (max_x_milano - min_x_milano)
# coord_scaled[ , 2] = (coord_utm[ , 2] - min_y_milano) / (max_y_milano - min_y_milano)

out_spatial_DR = construct_DR_Spatial(coord_scaled[ , 1], coord_scaled[ , 2], k1 = 5, k2 = 5)
Z_spat_unique = bind_cols(station = X_spat_unique$station, 
                          out_spatial_DR$Z)
colnames(Z_spat_unique) = c("station", paste0("spatial_DR", 1:ncol(out_spatial_DR$Z)))
s = ncol(Z_spat_unique)-1


data_for_regression = bind_cols(station = D$station, id_re = id_re, y = y, 
                                X_dummy, X_cont.scaled, Z_cont) |> 
  left_join(Z_spat_unique, by = "station")


# GROUPED_DATA
data_for_regression_grouped = data_for_regression |> 
  group_by_all() |> 
  summarise(freq = n()) |> 
  ungroup()

y_grouped = data_for_regression_grouped$y
id_re_grouped = data_for_regression_grouped$id_re
X_dummy_grouped = data_for_regression_grouped |> select(all_of(colnames(X_dummy))) |> as.matrix()
X_cont.scaled_grouped = data_for_regression_grouped |> select(all_of(colnames(X_cont.scaled))) |> as.matrix()
Z_cont_grouped = data_for_regression_grouped |> select(all_of(colnames(Z_cont))) |> as.matrix()
Z_spat_grouped = data_for_regression_grouped |> select(all_of(colnames(Z_spat_unique)[-1])) |> as.matrix()
freq_grouped = data_for_regression_grouped$freq


# Checks
sum(freq_grouped) == n
ncol(X_dummy_grouped) == q; ncol(X_cont.scaled_grouped) == p; ncol(Z_cont_grouped) == sum(d); ncol(Z_spat_grouped) == s


df = data_for_regression_grouped |> select(-station, -freq) |> as.data.frame()
df$y[df$y == 0] = 1e-3
df$y[df$y == 1] = 1-1e-3
X_complete_grouped = cbind("Intercept" = 1, df |> select(-y) |> as.matrix())

qr_reg = lqmm::lqmm(fixed = as.formula(paste("qlogis(y) ~ ", paste(colnames(df)[-(1:2)], collapse= "+"))), 
                    random = ~ 1, group = id_re, tau = 0.5, 
                    data = df, weights = freq_grouped)

# qr_reg = quantreg::rq(y ~ .,
#                       data = df,
#                       weights = freq_grouped,
#                       tau = 0.5)
# res_qr = qr_reg$residuals
# beta_qr = qr_reg$coefficients; names(beta_qr) = colnames(X_complete_grouped)
# unique_fit_qr = drop(X_complete_unique %*% beta_qr)


X_complete_unique = unique(X_complete_grouped)

res_qr = residuals(qr_reg)
beta_qr = coef(qr_reg); names(beta_qr) = colnames(X_complete_grouped)[-2] # Exclude id_re
unique_fit_qr = plogis(drop(X_complete_unique[ , names(beta_qr)] %*% beta_qr) + 
                         ranef(qr_reg)[["(Intercept)"]][X_complete_unique[ , "id_re"]])

labs_grouped = data_for_regression_grouped |> select(station, Weekend) |>  #, Rain) |> 
  mutate(Weekend = ifelse(Weekend == 1, "Weekend", "Work")) |> #,
  # Rain = ifelse(Rain == 1, "Rain", "Dry")) |>
  apply(1, paste, collapse = " | ")
labs = data_for_regression |> select(station, Weekend) |>  #, Rain) |> 
  mutate(Weekend = ifelse(Weekend == 1, "Weekend", "Work")) |> #,
  # Rain = ifelse(Rain == 1, "Rain", "Dry")) |>
  apply(1, paste, collapse = " | ")
labs_unique = unique(labs_grouped)

names(y) = labs
names(y_grouped) = labs_grouped
names(unique_fit_qr) = labs_unique

# par(mfrow = c(3, 4), ask = T)
# for (name in unique(labs_grouped)){
#   hist(y[names(y) == name], breaks = seq(-1/160, 1+1/160, by = 1/80),
#        main = name, lty = "blank", freq = F)
#   abline(v = unique_fit_qr[name], col = "purple", lwd = 2)
#   abline(v = median(y_grouped[names(y_grouped) == name]), col = "green3", lty = 2, lwd = 2)
# }


# OPTIMAL K
const = (n^(5/12)) * (log(n)^(7/12)) / (n^(7/12))
kinit = (n/const)^(1/7)
k = kinit
# kfinal = kinit * (n^(1/7))
varp = var(res_qr[res_qr > 0])
varn = var(res_qr[res_qr < 0])
mup = mean(res_qr[res_qr > 0]) - mean(res_qr)
mun = mean(res_qr[res_qr < 0]) - mean(res_qr)
factorp = sqrt(varp+mup^2)
factorn = sqrt(varn+mun^2)
k_multi = c( 
  sqrt(k)/factorn,
  sqrt(k)/factorp
)

# out_optim_scaled = sapply(
#   c("Nelder-Mead", "CG", "L-BFGS-B", "SANN"),
#   # c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"),
#   function(m){
#     est = optim(par = qr_reg$coefficients,
#                 fn = function(b){
#                   predictor_unique = logit_link(X_complete_grouped%*%b)
#                   sum(freq_grouped * 
#                         loss_asymm(y_grouped - predictor_unique, 
#                                    k_multi, 
#                                    1e-1))
#                 }, method = m)$par
#     est
#   }
# )
# 
# which.min(apply(out_optim_scaled, 2, function(b){
#   predictor_unique = logit_link(X_complete_grouped%*%b)
#   sum(freq_grouped * 
#         loss_asymm(y_grouped - predictor_unique, 
#                    k_multi, 
#                    1e-1))
# }))
# b_init = out_optim_scaled[ , "L-BFGS-B"]

b_init = optim(par = beta_qr,
               fn = function(b){
                 predictor_unique = logit_link(X_complete_grouped[ , -2]%*%b)
                 sum(freq_grouped * 
                       loss_asymm(y_grouped - predictor_unique, 
                                  k_multi, 
                                  1e-1))
               }, method = "L-BFGS-B")$par

unique_fit_init = logit_link(drop(X_complete_unique[ , -2] %*% b_init))
names(unique_fit_init) = labs_unique
# par(mfrow = c(3, 4), ask = T)
# for (name in unique(labs_grouped)){
#   hist(y[names(y) == name],
#        breaks = seq(-1/160, 1+1/160, by = 1/80),
#        main = name, lty = "blank")
#   abline(v = unique_fit_init[name], col = "red", lwd = 2)
#   abline(v = unique_fit_qr[name], col = "purple", lwd = 2)
#   abline(v = names(which.max(table(y[names(y) == name]))),
#          col = "green")
# }


# Find theta value for SDP for rand. eff.
library(sdPrior)
Z_th = diag(1, max(id_re))
Kinv_th = diag(1, ncol(Z_th))
theta_SDP = get_theta(Z = Z_th, Kinv = Kinv_th)


a_tau1 = 2.1; a_tau2 = 1e-3
b_tau1 = qgamma(0.001, shape = 2.1, rate = 1) * 5^2; b_tau2 = 1e-3

k_post = k_multi; k_prop = (k_multi)^(4/5)
c_post = 1e-1; c_prop = 1e-1

niter = 20e3
out_MCMC = MCMC_IWLS_bikes(niter = niter,
                           X_dummy = X_dummy_grouped, 
                           X_cont = X_cont.scaled_grouped, Z_cont = Z_cont_grouped,
                           Z_spat = Z_spat_grouped,
                           id = id_re_grouped,
                           y = y_grouped, 
                           freq = freq_grouped,
                           link = logit_link, link_prime = logit_prime, link_second = logit_second, 
                           k = k_post, k_deriv = k_prop,
                           c = c_post, c_deriv = c_prop,
                           d = d,
                           # 1 per ogni binary (q), 1 per ogni linear (p), 1 per ogni non-linear spline group (p), 1 per spatial (1)
                           a_tau = c(rep(a_tau1, q+p), rep(a_tau2, p+1)),
                           b_tau = c(rep(a_tau1, q+p), rep(a_tau2, p+1)),
                           theta_re = theta_SDP$root,
                           verbose = 2, K_fold = 5, print_step = 50,
                           step = 0.1,
                           beta0 = b_init,
                           seed = 1708)
saveRDS(out_MCMC, "application/bike-sharing/rds_mcmc/out_MCMC_SpatCont_RandEff.RDS")
# out_MCMC = readRDS("application/bike-sharing/rds_mcmc/out_MCMC_spatial_logit_norain_blockspatial.RDS")


burnin = 5000
update_idx <- c(1,
                if (q > 0) 2:(q+1) else integer(0),
                if (p > 0) (q + 2):(q + 1 + p) else integer(0),
                if (p > 0) q + 1 + p + rep(seq_len(p), d) else integer(0),
                if (s > 0) q+1+p+p + c(rep(1:(floor(s / 5)), each = 5), rep(floor(s / 5)+(s %% 5 > 0), s %% 5)) else integer(0)
)


IMG_FOLD = "application/bike-sharing/pics/results/SpatCont_RandEff/"

pdf(paste0(IMG_FOLD, "traceplot_spatial_logit.pdf"), width = 11, height = 7)
par(mfrow =c(5, 5), mar = c(3, 3, 1, 1))
for (j in 1:ncol(out_MCMC$beta)){
  plot(out_MCMC$beta[ , j], type = "l", 
       main = paste0(colnames(X_complete_grouped)[-2][j], " | AP: ", 
                     round(out_MCMC$acc_prop$beta[update_idx[j]], 2)),
       xlab = "Iteration", ylab = "", cex.main = 0.85)
  abline(v = burnin, col = 2, lty = 1)
}
dev.off()


pdf(paste0(IMG_FOLD, "traceplot_tau_spatial_logit.pdf"), width = 11, height = 7)
par(mfrow =c(3, 3), mar = c(3, 3, 1, 1))
for (j in 1:ncol(out_MCMC$tau)){
  plot(out_MCMC$tau[ , j], type = "l", main = paste("tau", j, sep = "_"),
       xlab = "Iteration", ylab = "Tau")
  abline(v = burnin, col = 2, lty = 1)
}
dev.off()


labs_re = data_for_regression_grouped |> select(station, id_re) |> 
  distinct() |> arrange(id_re)

pdf(paste0(IMG_FOLD, "traceplot_u_spatial_logit.pdf"), width = 11, height = 7)
par(mfrow =c(5, 5), mar = c(3, 3, 1, 1))
for (j in 1:ncol(out_MCMC$u)){
  st = labs_re$station[labs_re$id_re == j]
  plot(out_MCMC$u[ , j], type = "l", main = st, cex.main = 0.75,
       xlab = "Iteration", ylab = "Random effect")
  abline(v = burnin, col = 2, lty = 1)
}
dev.off()



pdf(paste0(IMG_FOLD, "traceplot_sigma2u_spatial_logit.pdf"), width = 11, height = 7)
par(mfrow =c(1, 1), mar = c(3, 3, 1, 1))
plot(out_MCMC$sigma2u, type = "l", main = "sigma2u",
     xlab = "Iteration", ylab = "Variance of random effects")
abline(v = burnin, col = 2, lty = 1)
dev.off()


# Caterpillar plot for random effect ----
u_mean = colMeans(out_MCMC$u[-(1:burnin), ])
u_median = apply(out_MCMC$u[-(1:burnin), ], 2, median)
u_lower = apply(out_MCMC$u[-(1:burnin), ], 2, quantile, probs = 0.025)
u_upper = apply(out_MCMC$u[-(1:burnin), ], 2, quantile, probs = 0.975)
u_df = labs_re |> 
  mutate(mean = u_mean, median = u_median, lower = u_lower, upper = u_upper)

# Create a data frame containing the group names, mean estimates, and lower/upper bounds
ggplot(u_df, aes(x = reorder(station, mean))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.5) +
  geom_point(aes(y = mean)) +
  geom_point(aes(y = median), color = "blue") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0) +
  geom_text(aes(y = lower, 
                label = if_else(lower > -1, station, "")), hjust = 1, size = 1.75) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(size = 6)) +
  labs(x = "Group", y = "Random Effect Estimate")

ggsave(path = IMG_FOLD,
       filename = "caterpillar_plot_random_effects_logit.pdf",
       width = 11, height = 22)


# Compute estimates ----
est = apply(out_MCMC$beta[-(1:burnin), ], 2, mean)
unique_fit_mode = logit_link(drop(X_complete_unique[ , -2] %*% est + colMeans(out_MCMC$u[-(1:burnin), X_complete_unique[ , 2]])))
q025 = apply(logit_link(X_complete_unique[ , -2] %*% t(out_MCMC$beta[-(1:burnin), ]) +
                          + t(out_MCMC$u[-(1:burnin), X_complete_unique[ , 2]])), 1, quantile, probs = 0.025)
q975 = apply(logit_link(X_complete_unique[ , -2] %*% t(out_MCMC$beta[-(1:burnin), ]) +
                          + t(out_MCMC$u[-(1:burnin), X_complete_unique[ , 2]])), 1, quantile, probs = 0.975)
names(q975) = names(q025) = names(unique_fit_mode) = labs_unique


# Fitted vs observed ----
## Hist ----
pdf(paste0(IMG_FOLD, "hist_vs_fitted_spatial_logit.pdf"), width = 11, height = 7)
par(mfrow = c(3, 4))
for (name in labs_unique){
  hist(y[names(y) == name], 
       xlim = c(-0.1, 1.1),
       breaks = seq(-1/160, 1+1/160, by = 1/80), 
       freq = F, cex.main = 0.85, lty = "blank",
       main = name, xlab = "Unavailability times (numeric and normalized)")
  ttt = table(y[names(y) == name])
  max_ttt = max(ttt)
  obs_mode = names(which(ttt == max_ttt))
  abline(v = obs_mode, 
         col = "green3", lty = 2, lwd = 1)
  abline(v = unique_fit_init[name], col = "blue", lwd = 1)
  abline(v = unique_fit_mode[name], col = "red", lwd = 1)
  abline(v = q025[name], col = "orange", lwd = 1)
  abline(v = q975[name], col = "orange", lwd = 1)
}
dev.off()

## Scatterplot ----
modal_times = D |>
  mutate(lab = paste(station, 
                     ifelse(Weekend == 1, "Weekend", "Work"), sep = " | ")) |>  #, 
  group_by(lab, time_numeric_shifted) |> summarise(n = n()) |>
  group_by(lab) |> summarise(mode_time_numeric_shifted = time_numeric_shifted[which.max(n)],
                             tot_events = sum(n))

predict_df = bind_cols(lab = labs_unique, id_re = X_complete_unique[ , 2], X_complete_unique[ , -2])
estimates = predict_df |>
  mutate(fitted = unique_fit_mode[predict_df$lab]*100*720) |>
  left_join(modal_times, by = "lab") |>
  select(lab, fitted, mode_time_numeric_shifted, tot_events)


plt_fit_vs_est = estimates |>
  mutate(out_of_bounds = factor(ifelse(fitted < 0 | fitted > 60*60*20, "OOF", "OK"))) |>
  # filter(out_of_bounds == "OK") |>
  ggplot() +
  geom_point(aes(x = mode_time_numeric_shifted, y = fitted, 
                 alpha = out_of_bounds, color = out_of_bounds, size = tot_events)) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  scale_x_continuous(breaks = seq(0, 60*60*20, by = 60*60*2),
                     labels = c("06:00", "08:00", "10:00", "12:00", "14:00", "16:00", "18:00",
                                "20:00", "22:00", "00:00", "02:00")) +
  scale_y_continuous(breaks = seq(0, 60*60*20, by = 60*60*2),
                     labels = c("06:00", "08:00", "10:00", "12:00", "14:00", "16:00", "18:00",
                                "20:00", "22:00", "00:00", "02:00"),
                     limits = c(0, 60*60*20)) +
  scale_color_manual(name = "Out of bounds", labels = c("In bounds", "Out of bounds"),
                     values = c("navy", "red3")) +
  scale_alpha_manual(name = "Out of bounds", labels = c("In bounds", "Out of bounds"),
                     values = c(1, 0.25)) +
  scale_size_continuous(name = "Total events", range = c(1, 5)) +
  labs(x = "Observed modal time", y = "Fitted modal time")

pdf(paste0(IMG_FOLD, "scatter_fit_vs_obs_logit.pdf"), width = 11, height = 7)
print(plt_fit_vs_est)
dev.off()

# Covariates effects ----
pdf(paste0(IMG_FOLD, "covariates_effects_logit.pdf"), width = 11, height = 7)
par(mfrow = c(3, 4))
for (j in 1:p){
  X_tmp = X_complete_unique[ , -2]
  x = X_tmp[ , 1+q+j]
  idx = 1 + q + p + sum(d[min(1, j-1):(j-1)]) + (1:d[j])
  S = X_tmp[ , idx]
  b = est[c(1, 1+q+j, idx)]
  X_pred = cbind(1, x, S)
  MCMC_fit = logit_link(X_pred[order(x), ] %*% t(out_MCMC$beta[-(1:burnin), c(1, 1+q+j, idx)]))
  mean_fit = logit_link(X_pred[order(x), ] %*% b)
  plot(sort(x), mean_fit, type = "l", ylim = c(0, 1),
       main = colnames(X_tmp)[1+q+j], xlab = colnames(X_tmp)[1+q+j], ylab = "")
  lines(sort(x), apply(MCMC_fit, 1, quantile, probs = 0.025), lty = 2)
  lines(sort(x), apply(MCMC_fit, 1, quantile, probs = 0.975), lty = 2)
  rug(x, col = "black")
}
dev.off()

# Correlation in posterior samples ----
plt_post_cor = cor(out_MCMC$beta[-(1:burnin), ]) |> reshape2::melt() |> 
  mutate(gr1 = update_idx[Var1], gr2 = update_idx[Var2]) |>
  ggplot() +
  geom_tile(aes(x = Var1, y = fct_rev(Var2), fill = value)) +
  # geom_text(aes(x = Var1, y = fct_rev(Var2), label = round(value, 2)), size = 3) +
  facet_grid(gr2 ~ gr1, scales = "free") +
  scale_fill_gradient2(low = "#4575b4", mid = "#f7f7f7", high = "#d73027", 
                       midpoint = 0, limits = c(-1, 1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

pdf(paste0(IMG_FOLD, "posterior_correlation_coefficients_logit.pdf"), width = 22, height = 14)
print(plt_post_cor)
dev.off()


# Spatial effect ----
## Paul's plot
require(sf)
milano <- st_read("https://dati.comune.milano.it/dataset/e75d91fa-eca6-4ee5-b96e-08bcdbb8d6f0/resource/f56cb432-83e6-48de-ae30-d39b4be61e85/download/confine_comune_milano_layer_0_confine_comune_milano.geojson")
milano_utm <- st_transform(milano, crs = 32632)
milano_coords <- st_coordinates(milano_utm)
xmin = min(coord_utm[ , 1]); xmax = max(coord_utm[ , 1])
ymin = min(coord_utm[ , 2]); ymax = max(coord_utm[ , 2])
milano_coords_scaled = milano_coords
milano_coords_scaled[, 1] <- (milano_coords[, 1] - xmin) / (xmax - xmin)
milano_coords_scaled[, 2] <- (milano_coords[, 2] - ymin) / (ymax - ymin)

pdf(paste0(IMG_FOLD, "spatial_effect_Paul_noconvex_logit.pdf"), width = 11, height = 7)
plot_effect_spat(coord_scaled, est[75:98], out_spatial_DR$Trafo,
                 knots1 = out_spatial_DR$knots1, knots2 = out_spatial_DR$knots2,
                 xlim = c(-1, 1.5), ylim = c(-0.75, 1.25), convex_hull = F)
lines(x = c(0, 1, 1, 0, 0), y = c(0, 0, 1, 1, 0), col = "grey", lwd = 2) # Add square [0, 1] x [0, 1]
lines(milano_coords_scaled)
dev.off()

pdf(paste0(IMG_FOLD, "spatial_effect_Paul_convex_logit.pdf"), width = 11, height = 7)
plot_effect_spat(coord_scaled, est[75:98], out_spatial_DR$Trafo, ngrid = 100,
                 knots1 = out_spatial_DR$knots1, knots2 = out_spatial_DR$knots2,
                 xlim = c(-1, 1.5), ylim = c(-0.75, 1.25), convex_hull = T)
lines(x = c(0, 1, 1, 0, 0), y = c(0, 0, 1, 1, 0), col = "grey", lwd = 2) # Add square [0, 1] x [0, 1]
lines(milano_coords_scaled)
dev.off()

pdf(paste0(IMG_FOLD, "spatial_effect_Paul_concave_logit.pdf"), width = 11, height = 7)
plot_effect_spat(coord_scaled, est[75:98], out_spatial_DR$Trafo, ngrid = 100,
                 knots1 = out_spatial_DR$knots1, knots2 = out_spatial_DR$knots2,
                 xlim = c(-1, 1.5), ylim = c(-0.75, 1.25), concave_hull = T)
lines(x = c(0, 1, 1, 0, 0), y = c(0, 0, 1, 1, 0), col = "grey", lwd = 2) # Add square [0, 1] x [0, 1]
lines(milano_coords_scaled)
dev.off()

## My plot ----
MCMC_fit_spatial = drop(as.matrix(Z_spat_unique[ , -1])) %*% t(out_MCMC$beta[-(1:burnin), 75:98])
pred_spat = data.frame(x = coord_scaled[ , 1], y = coord_scaled[ , 2],
                       spat_eff = drop(as.matrix(Z_spat_unique[ , -1]) %*% est[75:98]),
                       q025 = apply(MCMC_fit_spatial, 1, quantile, 0.025),
                       q975 = apply(MCMC_fit_spatial, 1, quantile, 0.975))
pred_spat |> 
  ggplot() +
  geom_point(aes(x = x, y = y, fill = spat_eff), size = 2, pch = 21, color = "grey") + 
  scale_fill_gradient2(low = "#4575b4", mid = "#f7f7f7", high = "#d73027", midpoint = 0) +
  geom_path(data = milano_coords_scaled, aes(x = X, y = Y), color = "black", lwd = 0.3)

ggsave(filename = "spatial_effect_my_logit.pdf",
       path = IMG_FOLD,
       width = 11, height = 7)

