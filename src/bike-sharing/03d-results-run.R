# === === === === === === === === === === === === === === === === === === === ===
# 03c-results-2025-logit.R
#
# Diagnostics, posterior summaries and figures. Reads only saved objects, so it
# can be re-run cheaply without touching the MCMC.
#
# Input : prep_<MODEL_TAG>.rds, out_MCMC_<MODEL_TAG>.rds
# Output: figures under application/bike-sharing/pics/results/<MODEL_TAG>/
# === === === === === === === === === === === === === === === === === === === ===

suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
})

source("src/bike-sharing/03a-config-run.R")
source("src/01a-utils.R")
source("src/01b-gibbs.R")

prep     <- readRDS(PATH_PREP)
out_MCMC <- readRDS(PATH_MCMC)

p <- prep$dims$p; q <- prep$dims$q; d <- prep$dims$d
idx        <- prep$idx
update_idx <- prep$update_idx

y                 <- prep$y
D                 <- prep$D
X_complete_unique <- prep$X_complete_unique
X_design_unique   <- drop_id_col(X_complete_unique)
id_unique         <- X_complete_unique[, "id_re"]

labs_unique <- prep$labels$unique
labs_re     <- prep$labels$re

unique_fit_init <- out_MCMC$init$unique_fit_init

# Post-burn-in draws, extracted once.
keep   <- -seq_len(BURNIN)
beta_s <- out_MCMC$beta[keep, , drop = FALSE]
u_s    <- out_MCMC$u[keep, , drop = FALSE]

# ---- Traceplots --------------------------------------------------------------
coef_names <- colnames(X_design_unique)

with_pdf(paste0(DIR_IMG, "traceplot_beta.pdf"), {
  par(mfrow = c(5, 5), mar = c(3, 3, 1, 1))
  for (j in seq_len(ncol(out_MCMC$beta))) {
    plot(out_MCMC$beta[, j], type = "l", cex.main = 0.85,
         main = paste0(coef_names[j], " | AP: ",
                       round(out_MCMC$acc_prop$beta[update_idx[j]], 2)),
         xlab = "Iteration", ylab = "")
    abline(v = BURNIN, col = 2)
  }
})

with_pdf(paste0(DIR_IMG, "traceplot_beta_noburnin_thinned.pdf"), {
  par(mfrow = c(5, 5), mar = c(3, 3, 1, 1))
  for (j in seq_len(ncol(out_MCMC$beta))) {
    plot(out_MCMC$beta[seq(BURNIN + 1, N_ITER, by = 10), j], type = "l", cex.main = 0.85,
         main = paste0(coef_names[j], " | AP: ",
                       round(out_MCMC$acc_prop$beta[update_idx[j]], 2)),
         xlab = "Iteration", ylab = "")
    abline(v = BURNIN, col = 2)
  }
})

with_pdf(paste0(DIR_IMG, "traceplot_tau.pdf"), {
  par(mfrow = c(3, 3), mar = c(3, 3, 1, 1))
  for (j in seq_len(ncol(out_MCMC$tau))) {
    plot(out_MCMC$tau[, j], type = "l", main = paste0("tau_", j),
         xlab = "Iteration", ylab = "Tau")
    abline(v = BURNIN, col = 2)
  }
})

with_pdf(paste0(DIR_IMG, "traceplot_tau_noburnin_thinned.pdf"), {
  par(mfrow = c(3, 3), mar = c(3, 3, 1, 1))
  for (j in seq_len(ncol(out_MCMC$tau))) {
    plot(out_MCMC$tau[seq(BURNIN + 1, N_ITER, by = 10), j], type = "l", main = paste0("tau_", j),
         xlab = "Iteration", ylab = "Tau")
    abline(v = BURNIN, col = 2)
  }
})


with_pdf(paste0(DIR_IMG, "traceplot_u.pdf"), {
  par(mfrow = c(5, 5), mar = c(3, 3, 1, 1))
  for (j in seq_len(ncol(out_MCMC$u))) {
    plot(out_MCMC$u[, j], type = "l", cex.main = 0.75,
         main = labs_re$station[labs_re$id_re == j],
         xlab = "Iteration", ylab = "Random effect")
    abline(v = BURNIN, col = 2)
  }
})

with_pdf(paste0(DIR_IMG, "traceplot_u_noburnin_thinned.pdf"), {
  par(mfrow = c(5, 5), mar = c(3, 3, 1, 1))
  for (j in seq_len(ncol(out_MCMC$u))) {
    plot(out_MCMC$u[seq(BURNIN + 1, N_ITER, by = 10), j], type = "l", cex.main = 0.75,
         main = labs_re$station[labs_re$id_re == j],
         xlab = "Iteration", ylab = "Random effect")
    abline(v = BURNIN, col = 2)
  }
})


with_pdf(paste0(DIR_IMG, "traceplot_sigma2u.pdf"), {
  par(mfrow = c(1, 1), mar = c(3, 3, 1, 1))
  plot(out_MCMC$sigma2u, type = "l", main = "sigma2u",
       xlab = "Iteration", ylab = "Variance of random effects")
  abline(v = BURNIN, col = 2)
})

with_pdf(paste0(DIR_IMG, "traceplot_sigma2u_noburnin_thinned.pdf"), {
  par(mfrow = c(1, 1), mar = c(3, 3, 1, 1))
  plot(out_MCMC$sigma2u[seq(BURNIN + 1, N_ITER, by = 10)], type = "l", main = "sigma2u",
       xlab = "Iteration", ylab = "Variance of random effects")
  abline(v = BURNIN, col = 2)
})


# ---- Caterpillar plot of the random effects ----------------------------------
u_df <- labs_re |>
  mutate(
    mean   = colMeans(u_s),
    median = apply(u_s, 2, median),
    lower  = apply(u_s, 2, quantile, probs = 0.025),
    upper  = apply(u_s, 2, quantile, probs = 0.975)
  )

plt_caterpillar <- ggplot(u_df, aes(x = reorder(station, mean))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.5) +
  geom_point(aes(y = mean)) +
  geom_point(aes(y = median), color = "blue") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0) +
  geom_text(aes(y = lower, label = if_else(lower > -1, station, "")),
            hjust = 1, size = 1.75) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(size = 6)) +
  labs(x = "Group", y = "Random effect estimate")

ggsave(plot = plt_caterpillar, path = DIR_IMG,
       filename = "caterpillar_random_effects.pdf", width = 11, height = 22)

# ---- Estimates ---------------------------------------------------------------
est <- colMeans(beta_s)
names(est) <- coef_names

# Draws of the linear predictor, computed once and reused for mean and quantiles.
eta_draws <- X_design_unique %*% t(beta_s) + t(u_s[, id_unique, drop = FALSE])
fit_draws <- logit_link(eta_draws)

unique_fit_mode <- logit_link(drop(X_design_unique %*% est+ colMeans(u_s)[id_unique]))
q025 <- apply(fit_draws, 1, quantile, probs = 0.025)
q975 <- apply(fit_draws, 1, quantile, probs = 0.975)
names(unique_fit_mode) <- names(q025) <- names(q975) <- labs_unique

# ---- Fitted vs observed ------------------------------------------------------
## Histograms per station x day-type ----
with_pdf(paste0(DIR_IMG, "hist_vs_fitted.pdf"), {
  par(mfrow = c(3, 4))
  for (name in labs_unique) {
    y_sub <- y[names(y) == name]
    hist(y_sub, xlim = c(-0.1, 1.1),
         breaks = seq(-1 / 160, 1 + 1 / 160, by = 1 / 80),
         freq = FALSE, cex.main = 0.85, lty = "blank", main = name,
         xlab = "Unavailability times (numeric and normalized)")
    
    tab <- table(y_sub)
    obs_mode <- as.numeric(names(which(tab == max(tab))))
    
    abline(v = obs_mode,                col = "green3", lty = 2)
    abline(v = unique_fit_init[name],   col = "blue")
    abline(v = unique_fit_mode[name],   col = "red")
    abline(v = c(q025[name], q975[name]), col = "orange")
  }
})

## Scatterplot Fitted vs observed modal time ----
modal_times <- D |>
  mutate(lab = paste(station, if_else(Weekend == 1, "Weekend", "Work"), sep = " | ")) |>
  count(lab, time_numeric_shifted, name = "n") |>
  group_by(lab) |>
  summarise(
    mode_time_numeric_shifted = time_numeric_shifted[which.max(n)],
    tot_events = sum(n),
    .groups = "drop"
  )

HOURS_LABS <- c("06:00", "08:00", "10:00", "12:00", "14:00", "16:00",
                "18:00", "20:00", "22:00", "00:00", "02:00")
T_MAX <- 60 * 60 * 20

estimates <- tibble(lab = labs_unique, fitted = unique_fit_mode * SCALE_TIME) |>
  left_join(modal_times, by = "lab")

plt_fit_vs_est <- estimates |>
  mutate(out_of_bounds = factor(if_else(fitted < 0 | fitted > T_MAX, "OOF", "OK"))) |>
  ggplot() +
  geom_point(aes(x = mode_time_numeric_shifted, y = fitted,
                 alpha = out_of_bounds, color = out_of_bounds, size = tot_events)) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  scale_x_continuous(breaks = seq(0, T_MAX, by = 60 * 60 * 2), labels = HOURS_LABS) +
  scale_y_continuous(breaks = seq(0, T_MAX, by = 60 * 60 * 2), labels = HOURS_LABS,
                     limits = c(0, T_MAX)) +
  scale_color_manual(name = "Out of bounds", labels = c("In bounds", "Out of bounds"),
                     values = c("navy", "red3")) +
  scale_alpha_manual(name = "Out of bounds", labels = c("In bounds", "Out of bounds"),
                     values = c(1, 0.25)) +
  scale_size_continuous(name = "Total events", range = c(1, 5)) +
  labs(x = "Observed modal time", y = "Fitted modal time") +
  theme_bw()

ggsave(plot = plt_fit_vs_est, path = DIR_IMG,
       filename = "scatter_fit_vs_obs.pdf", width = PDF_W, height = PDF_H)

# ---- Covariate effects (linear + non-linear part per covariate) --------------
with_pdf(paste0(DIR_IMG, "covariates_effects.pdf"), {
  par(mfrow = c(3, 4))
  for (j in seq_len(p)) {
    i_lin <- idx$linear[j]
    i_dr  <- idx$dr_block[[j]]
    i_all <- c(idx$intercept, i_lin, i_dr)
    
    x <- X_design_unique[, i_lin]
    ord <- order(x)
    X_pred <- X_design_unique[, i_all, drop = FALSE][ord, , drop = FALSE]
    
    x = x*scaling_factors$scale[j] + scaling_factors$center[j]
    
    switch(TRANSF_CONT_X,
           "identity" = {x <- x},
           "log" = {x <- exp(x) - 1})
    
    mean_fit <- X_pred %*% est[i_all]
    MCMC_fit <- X_pred %*% t(beta_s[, i_all, drop = FALSE])
    lower    <-  apply(MCMC_fit, 1, quantile, probs = 0.025)
    upper    <-  apply(MCMC_fit, 1, quantile, probs = 0.975)
    
    ml <- min(lower); mu <- max(upper); r <- mu - ml
    ylim <- c(ml*(1-sign(ml)*0.1*r), mu*(1+sign(mu)*0.1*r))
    
    plot(x[ord], mean_fit, type = "l", ylim = ylim,
         main = coef_names[i_lin], xlab = coef_names[i_lin], ylab = "")
    lines(x[ord], apply(MCMC_fit, 1, quantile, probs = 0.025), lty = 2)
    lines(x[ord], apply(MCMC_fit, 1, quantile, probs = 0.975), lty = 2)
    rug(x)
  }
})



# ---- Posterior correlation between coefficients ------------------------------
plt_post_cor <- cor(beta_s) |>
  reshape2::melt() |>
  mutate(gr1 = update_idx[Var1], gr2 = update_idx[Var2]) |>
  ggplot() +
  geom_tile(aes(x = Var1, y = fct_rev(Var2), fill = value)) +
  facet_grid(gr2 ~ gr1, scales = "free") +
  scale_fill_gradient2(low = "#4575b4", mid = "#f7f7f7", high = "#d73027",
                       midpoint = 0, limits = c(-1, 1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(plot = plt_post_cor, path = DIR_IMG,
       filename = "posterior_correlation_coefficients.pdf", width = 22, height = 14)

# ---- Spatial effect ----------------------------------------------------------
coord_utm      <- prep$spatial$coord_utm
coord_scaled   <- prep$spatial$coord_scaled
out_spatial_DR <- prep$spatial$out_spatial_DR
est_spat       <- est[idx$spatial]

milano_scaled <- rescale_coords(
  load_milano_utm(), range(coord_utm[, 1]), range(coord_utm[, 2])
) |>
  as.data.frame()

## Basis-surface plots, with and without hull clipping ----
hull_variants <- list(
  noconvex = list(convex_hull = FALSE),
  convex   = list(convex_hull = TRUE),
  concave  = list(concave_hull = TRUE)
)

for (nm in names(hull_variants)) {
  with_pdf(paste0(DIR_IMG, "spatial_effect_surface_", nm, ".pdf"), {
    do.call(plot_effect_spat, c(
      list(coord_scaled, est_spat, out_spatial_DR$Trafo, ngrid = 100,
           knots1 = out_spatial_DR$knots1, knots2 = out_spatial_DR$knots2,
           xlim = c(-1, 1.5), ylim = c(-0.75, 1.25)),
      hull_variants[[nm]]
    ))
    # basis domain [0, 1] x [0, 1]
    lines(x = c(0, 1, 1, 0, 0), y = c(0, 0, 1, 1, 0), col = "grey", lwd = 2)
    lines(milano_scaled$X, milano_scaled$Y)
  })
}

## Station-level spatial effect ----
Z_spat_mat <- as.matrix(prep$Z_spat_unique[, -1])
spat_draws <- Z_spat_mat %*% t(beta_s[, idx$spatial, drop = FALSE])

pred_spat <- data.frame(
  x        = coord_scaled[, 1],
  y        = coord_scaled[, 2],
  spat_eff = drop(Z_spat_mat %*% est_spat),
  q025     = apply(spat_draws, 1, quantile, probs = 0.025),
  q975     = apply(spat_draws, 1, quantile, probs = 0.975)
)

plt_spat <- ggplot(pred_spat) +
  geom_path(data = milano_scaled, aes(x = X, y = Y), color = "black", lwd = 0.3) +
  geom_point(aes(x = x, y = y, fill = spat_eff), size = 2, pch = 21, color = "grey") +
  scale_fill_gradient2(low = "#4575b4", mid = "#f7f7f7", high = "#d73027", midpoint = 0) +
  coord_equal() +
  theme_bw() +
  labs(x = NULL, y = NULL, fill = "Spatial effect")

ggsave(plot = plt_spat, path = DIR_IMG,
       filename = "spatial_effect_stations.pdf", width = PDF_W, height = PDF_H)

message("Figures written to ", DIR_IMG)