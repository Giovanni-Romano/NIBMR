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

D = nobikes2025 |>
  select(-ID, -commit_at, -Precipitation, -Temperature)


# CREATE OBJECTS FOR REGRESSION ----
## Response ----
y = D$time_numeric_shifted / 100 / 720
## Binary covariates ----
X_dummy = (D |> select(Weekend, Rain)  |> as.matrix())
## Continuous covariates ----
X_cont = D |> select(population_density, starts_with("density_"), starts_with("closest_")) |> as.matrix()
X_cont.scaled = scale(X_cont)

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
coord_scaled = apply(coord_utm, 2, function(x) (x - min(x)) /  (max(x) - min(x)))
out_spatial_DR = construct_DR_Spatial(coord_scaled[ , 1], coord_scaled[ , 2], k1 = 5, k2 = 5)
Z_spat_unique = bind_cols(station = X_spat_unique$station, 
                          out_spatial_DR$Z)
colnames(Z_spat_unique) = c("station", paste0("spatial_DR", 1:ncol(out_spatial_DR$Z)))

s = ncol(Z_spat_unique)-1



data_for_regression = bind_cols(station = D$station, y = y, 
                            X_dummy, X_cont.scaled, Z_cont) |> 
  left_join(Z_spat_unique, by = "station")


# GROUPED_DATA
data_for_regression_grouped = data_for_regression |> 
  group_by_all() |> 
  summarise(freq = n()) |> 
  ungroup()

y_grouped = data_for_regression_grouped$y
X_dummy_grouped = data_for_regression_grouped |> select(all_of(colnames(X_dummy))) |> as.matrix()
X_cont.scaled_grouped = data_for_regression_grouped |> select(all_of(colnames(X_cont.scaled))) |> as.matrix()
Z_cont_grouped = data_for_regression_grouped |> select(all_of(colnames(Z_cont))) |> as.matrix()
Z_spat_grouped = data_for_regression_grouped |> select(all_of(colnames(Z_spat_unique)[-1])) |> as.matrix()
freq_grouped = data_for_regression_grouped$freq

# Checks
sum(freq_grouped) == n
ncol(X_dummy_grouped) == q; ncol(X_cont.scaled_grouped) == p; ncol(Z_cont_grouped) == sum(d); ncol(Z_spat_grouped) == s


df = data_for_regression_grouped |> select(-station, -freq) |> as.data.frame()
X_complete_grouped = cbind("Intercept" = 1, df |> select(-y) |> as.matrix())
qr_reg = quantreg::rq(y ~ ., 
                      data = df,
                      weights = freq_grouped,
                      tau = 0.5)

X_complete_unique = unique(X_complete_grouped)

res_qr = qr_reg$residuals
beta_qr = qr_reg$coefficients; names(beta_qr) = colnames(X_complete_grouped)
unique_fit_qr = drop(X_complete_unique %*% beta_qr)

labs_grouped = data_for_regression_grouped |> select(station, Weekend, Rain) |> 
  mutate(Weekend = ifelse(Weekend == 1, "Weekend", "Work"),
         Rain = ifelse(Rain == 1, "Rain", "Dry")) |>
  apply(1, paste, collapse = " | ")
labs = data_for_regression |> select(station, Weekend, Rain) |> 
  mutate(Weekend = ifelse(Weekend == 1, "Weekend", "Work"),
         Rain = ifelse(Rain == 1, "Rain", "Dry")) |>
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

out_optim_scaled = sapply(
  c("Nelder-Mead", "CG", "L-BFGS-B", "SANN"),
  # c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"),
  function(m){
    est = optim(par = qr_reg$coefficients,
                fn = function(b){
                  predictor_unique = X_complete_grouped%*%b
                  sum(freq_grouped * 
                        loss_asymm(y_grouped - predictor_unique, 
                                   k_multi, 
                                   1e-1))
                }, method = m)$par
    est
  }
)

which.min(apply(out_optim_scaled, 2, function(b){
  predictor_unique = X_complete_grouped%*%b
  sum(freq_grouped * 
        loss_asymm(y_grouped - predictor_unique, 
                   k_multi, 
                   1e-1))
  }))
b_init = out_optim_scaled[ , "L-BFGS-B"]


unique_fit_init = drop(X_complete_unique %*% b_init)
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



a_tau1 = 2.1; a_tau2 = 1e-3
b_tau1 = qgamma(0.001, shape = 2.1, rate = 1) * 5^2; b_tau2 = 1e-3

k_post = k_multi; k_prop = (k_multi)^(4/5)
c_post = 1e-1; c_prop = 1e-1

niter = 20e3
out_MCMC = MCMC_IWLS_bikes(niter = niter,
                           X_dummy = X_dummy_grouped, 
                           X_cont = X_cont.scaled_grouped, Z_cont = Z_cont_grouped,
                           Z_spat = Z_spat_grouped,
                           # Z_spat = Z_spat_grouped,
                           y = y_grouped, 
                           freq = freq_grouped,
                           k = k_post, k_deriv = k_prop,
                           c = c_post, c_deriv = c_prop,
                           d = d,
                           # 1 per ogni binary (q), 1 per ogni linear (p), 1 per ogni non-linear spline group (p), 1 per spatial (1)
                           a_tau = c(rep(a_tau1, q+p), rep(a_tau2, p+1)),
                           b_tau = c(rep(a_tau1, q+p), rep(a_tau2, p+1)),
                           verbose = 2, K_fold = 5, print_step = 250,
                           step = 0.1,
                           beta0 = b_init) #rep(0, 1+q+p+sum(d)+s))
saveRDS(out_MCMC, "application/bike-sharing/out_MCMC_spatial.RDS")
# out_MCMC = readRDS("application/bike-sharing/out_MCMC_spatial.RDS")


burnin = 5e3
update_idx <- c(1,
                if (q > 0) 2:(q+1) else integer(0),
                if (p > 0) (q + 2):(q + 1 + p) else integer(0),
                if (p > 0) q + 1 + p + rep(seq_len(p), d) else integer(0),
                if (s > 0) q+1+p+p + c(rep(1:(floor(s / 5)), each = 5), rep(floor(s / 5)+(s %% 5 > 0), s %% 5)) else integer(0)
)


pdf("application/bike-sharing/traceplot_spatial.pdf", width = 11, height = 7)
par(mfrow =c(5, 5), mar = c(3, 3, 1, 1))
for (j in 1:ncol(out_MCMC$beta)){
  plot(out_MCMC$beta[ , j], type = "l", 
       main = paste0(colnames(X_complete_grouped)[j], " | AP: ", 
                     round(out_MCMC$acc_prop[update_idx[j]], 2)),
       xlab = "Iteration", ylab = "", cex.main = 0.85)
  # abline(h = mean(out_MCMC$beta[-(1:burnin), j]), col = "red")
  abline(v = burnin, col = 2, lty = 1)
}
dev.off()

pdf("application/bike-sharing/traceplot_tau_spatial.pdf", width = 11, height = 7)
par(mfrow =c(3, 3), mar = c(3, 3, 1, 1))
for (j in 1:ncol(out_MCMC$tau)){
  plot(out_MCMC$tau[ , j], type = "l", main = paste("tau", j, sep = "_"),
       xlab = "Iteration", ylab = "Tau")
  # abline(h = mean(out_MCMC$beta[-(1:burnin), j]), col = "red")
  abline(v = burnin, col = 2, lty = 1)
}
dev.off()



est = apply(out_MCMC$beta[-(1:burnin), ], 2, mean)
unique_fit_mode = drop(X_complete_unique %*% est)
q025 = apply(X_complete_unique %*% t(out_MCMC$beta[-(1:burnin), ]), 1, quantile, probs = 0.025)
q975 = apply(X_complete_unique %*% t(out_MCMC$beta[-(1:burnin), ]), 1, quantile, probs = 0.975)
names(q975) = names(q025) = names(unique_fit_mode) = labs_unique

pdf("application/bike-sharing/hist_vs_fitted_spatial.pdf", width = 11, height = 7)
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

pdf("application/bike-sharing/covariates_effects.pdf", width = 11, height = 7)
par(mfrow = c(3, 4))
for (j in 1:p){
  x = X_complete_unique[ , 1+q+j]
  idx = 1 + q + p + sum(d[min(1, j-1):(j-1)]) + (1:d[j])
  S = X_complete_unique[ , idx]
  b = est[c(1, 1+q+j, idx)]
  X_pred = cbind(1, x, S)
  plot(sort(x), X_pred[order(x), ] %*% b, type = "l",
       main = colnames(X_complete_unique)[1+q+j], xlab = colnames(X_complete_unique)[1+q+j], ylab = "")
  # lines(sort(x), x[order(x)] * b[1], col = "blue")
  # lines(sort(x), S[order(x), ] %*% b[-1], col = "red")
  rug(x, col = "black")
}
dev.off()

# Fitted vs observed ----
modal_times = D |>
  mutate(lab = paste(station, 
                     ifelse(Weekend == 1, "Weekend", "Work"), 
                     ifelse(Rain == 1, "Rain", "Dry"), sep = " | ")) |>
  group_by(lab, time_numeric_shifted) |> summarise(n = n()) |>
  group_by(lab) |> summarise(mode_time_numeric_shifted = time_numeric_shifted[which.max(n)],
                             tot_events = sum(n))

predict_df = bind_cols(lab = labs_unique, X_complete_unique)
estimates = predict_df |>
  mutate(fitted = drop(as.matrix(predict_df[ , -1]) %*% est)*100*720) |>
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
                                "20:00", "22:00", "00:00", "02:00")) +
  scale_color_manual(name = "Out of bounds", labels = c("In bounds", "Out of bounds"),
                     values = c("navy", "red3")) +
  scale_alpha_manual(name = "Out of bounds", labels = c("In bounds", "Out of bounds"),
                     values = c(1, 0.25)) +
  scale_size_continuous(name = "Total events", range = c(1, 5)) +
  labs(x = "Observed modal time", y = "Fitted modal time")

pdf("application/bike-sharing/scatter_fit_vs_obs.pdf", width = 11, height = 7)
print(plt_fit_vs_est)
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

pdf("application/bike-sharing/posterior_correlation_coefficients.pdf", width = 22, height = 14)
print(plt_post_cor)
dev.off()
