plot_grid_comparison <- function(truth, estimate, new_x,
                                 effect = c("fcont", "f"),
                                 region = NULL,
                                 model = "Estimate",
                                 truth_shift = 0,
                                 theta = 40, phi = 25) {
  
  effect <- match.arg(effect)
  
  z_true <- truth[[effect]]
  z_est  <- estimate$point_est_grid[[effect]]
  
  # For spatial f, select the same region from both lists.
  if (is.list(z_true) || is.list(z_est)) {
    if (is.null(region) || length(region) != 1L ||
        !is.numeric(region) || !is.finite(region) ||
        region < 1 || region != as.integer(region)) {
      stop("Supply region as a positive integer list position.")
    }
    
    if (is.list(z_true)) {
      if (region > length(z_true)) stop("Region absent from truth.")
      z_true <- z_true[[region]]
    }
    
    # A model without regional effects has one common surface.
    if (is.list(z_est)) {
      if (region > length(z_est)) stop("Region absent from estimate.")
      z_est <- z_est[[region]]
    }
  }
  
  G <- length(new_x)
  
  stopifnot(
    G >= 2L,
    all(is.finite(new_x)),
    all(diff(new_x) > 0),
    is.matrix(z_true), is.matrix(z_est),
    all(dim(z_true) == c(G, G)),
    all(dim(z_est) == c(G, G)),
    all(is.finite(z_true)), all(is.finite(z_est)),
    length(truth_shift) == 1L, is.finite(truth_shift)
  )
  
  # Align truth with the identifying convention of the estimate.
  z_true <- z_true - truth_shift
  
  zlim <- range(z_true, z_est)
  if (diff(zlim) == 0) zlim <- zlim + c(-0.5, 0.5)
  
  region_label <- if (is.null(region)) "" else {
    paste0(" — region ", region)
  }
  
  old_par <- par(c("mar", "cex", "mex"))
  on.exit(par(old_par), add = TRUE)
  par(mar = c(2, 0, 3, 0.5))
  
  # True surface.
  projection <- persp(
    x = new_x, y = new_x, z = z_true,
    zlim = zlim,
    theta = theta, phi = phi,
    expand = 0.7,
    col = "lightblue",
    border = NA,
    shade = 0.3,
    ticktype = "detailed",
    xlab = "x1", ylab = "x2", zlab = effect,
    main = paste0("Truth and ", model, region_label,
                  "\n(theta = ", theta, "°, phi = ", phi, "°)")
  )
  
  # Estimated surface: use fewer grid lines for readability.
  wire_idx <- unique(round(seq(
    1, length(new_x),
    length.out = min(20L, length(new_x))
  )))
  
  # Lines along x2, at fixed x1.
  for (i in wire_idx) {
    coords <- trans3d(
      x = rep(new_x[i], length(new_x)),
      y = new_x,
      z = z_est[i, ],
      pmat = projection
    )
    lines(coords, col = "firebrick", lwd = 0.8)
  }
  
  # Lines along x1, at fixed x2.
  for (j in wire_idx) {
    coords <- trans3d(
      x = new_x,
      y = rep(new_x[j], length(new_x)),
      z = z_est[, j],
      pmat = projection
    )
    lines(coords, col = "firebrick", lwd = 0.8)
  }
  
  legend(
    "topright",
    legend = c("Truth", model),
    col = c("lightblue3", "firebrick"),
    lty = c(1, 1),
    lwd = c(6, 1),
    bty = "n",
    cex = 0.85
  )
  
  invisible(list(
    truth = z_true,
    estimate = z_est,
    difference = z_est - z_true
  ))
}


plot_grid_heatmap <- function(truth, estimate, new_x,
                              effect = c("fcont", "f"),
                              region = NULL,
                              model = "Estimate",
                              truth_shift = 0,
                              n_colors = 80L,
                              n_contours = 10L) {
  
  effect <- match.arg(effect)
  
  z_true <- truth[[effect]]
  z_est  <- estimate$point_est_grid[[effect]]
  
  if (is.list(z_true) || is.list(z_est)) {
    if (is.null(region) || length(region) != 1L ||
        !is.numeric(region) || !is.finite(region) ||
        region < 1 || region != as.integer(region)) {
      stop("Supply region as a positive integer list position.")
    }
    
    if (is.list(z_true)) {
      if (region > length(z_true)) stop("Region absent from truth.")
      z_true <- z_true[[region]]
    }
    
    if (is.list(z_est)) {
      if (region > length(z_est)) stop("Region absent from estimate.")
      z_est <- z_est[[region]]
    }
  }
  
  G <- length(new_x)
  
  stopifnot(
    G >= 2L,
    all(is.finite(new_x)),
    all(diff(new_x) > 0),
    is.matrix(z_true), is.matrix(z_est),
    all(dim(z_true) == c(G, G)),
    all(dim(z_est) == c(G, G)),
    all(is.finite(z_true)), all(is.finite(z_est)),
    length(truth_shift) == 1L, is.finite(truth_shift),
    n_colors >= 2L, n_contours >= 1L
  )
  
  z_true <- z_true - truth_shift
  
  # Shared scales for both panels.
  zlim <- range(z_true, z_est)
  if (diff(zlim) == 0) zlim <- zlim + c(-0.5, 0.5)
  
  breaks <- seq(zlim[1], zlim[2], length.out = n_colors + 1L)
  colors <- grDevices::hcl.colors(n_colors, palette = "Viridis")
  contour_levels <- pretty(zlim, n = n_contours)
  contour_levels <- contour_levels[
    contour_levels > zlim[1] & contour_levels < zlim[2]
  ]
  
  region_label <- if (is.null(region)) "" else {
    paste0(" — region ", region)
  }
  
  difference <- z_est - z_true
  
  # Symmetric scale around zero for the difference.
  error_max <- max(abs(difference))
  if (error_max == 0) error_max <- 1e-8
  
  error_breaks <- seq(
    -error_max, error_max,
    length.out = n_colors + 1L
  )
  error_colors <- grDevices::colorRampPalette(
    c("#2166AC", "white", "#B2182B")
  )(n_colors)
  
  old_par <- par(c("mfrow", "mar", "oma", "cex", "mex"))
  on.exit({
    par(mfrow = old_par$mfrow)
    par(mar = old_par$mar, oma = old_par$oma,
        cex = old_par$cex, mex = old_par$mex)
  }, add = TRUE)
  
  par(mfrow = c(1, 1))
  
  draw_panel <- function(z, title, palette, cuts, levels) {
    image(
      x = new_x, y = new_x, z = z,
      col = palette, breaks = cuts,
      xlim = range(new_x), ylim = range(new_x),
      asp = 1,
      xlab = "x1", ylab = "x2",
      main = paste0(title, region_label),
      useRaster = FALSE
    )
    
    levels <- levels[levels > min(z) & levels < max(z)]
    
    if (length(levels)) {
      contour(
        x = new_x, y = new_x, z = z,
        levels = levels, add = TRUE,
        col = "grey20", lwd = 0.7, labcex = 0.65
      )
    }
    box()
  }
  
  draw_colorbar <- function(cuts, palette, title) {
    plot.new()
    plot.window(
      xlim = c(0, 1), ylim = range(cuts),
      xaxs = "i", yaxs = "i"
    )
    rect(
      xleft = 0, ybottom = head(cuts, -1L),
      xright = 1, ytop = tail(cuts, -1L),
      col = palette, border = NA
    )
    axis(4, las = 1, cex.axis = 0.75)
    box()
    title(main = title, cex.main = 0.85)
  }
  
  # Truth.
  par(fig = c(0, 0.28, 0, 1), mar = c(4, 4, 3, 1))
  draw_panel(z_true, "Truth", colors, breaks, contour_levels)
  
  # Estimate, sharing the truth's color scale.
  par(fig = c(0.28, 0.56, 0, 1), new = TRUE,
      mar = c(4, 4, 3, 1))
  draw_panel(z_est, model, colors, breaks, contour_levels)
  
  # Shared truth/estimate color bar.
  par(fig = c(0.56, 0.64, 0, 1), new = TRUE,
      mar = c(4, 0.5, 3, 3))
  draw_colorbar(breaks, colors, effect)
  
  # Difference.
  par(fig = c(0.64, 0.92, 0, 1), new = TRUE,
      mar = c(4, 4, 3, 1))
  draw_panel(
    difference, "Estimate − truth",
    error_colors, error_breaks,
    pretty(c(-error_max, error_max), n = n_contours)
  )
  
  # Separate difference color bar.
  par(fig = c(0.92, 1, 0, 1), new = TRUE,
      mar = c(4, 0.5, 3, 3))
  draw_colorbar(error_breaks, error_colors, "Difference")
  
  invisible(list(
    truth = z_true,
    estimate = z_est,
    difference = z_est - z_true
  ))
}


# Surface plots ----
truth <- grid_truth(
  s_X2$scenario$predictor,
  new_x,
  sp_eff = s$scenario$gamma
)

# Center using training observations, not grid points.
covariate_shift <- mean(s_X3$scenario$f0$f1) +
  mean(s_X3$scenario$f0$f2)

par(mfrow = c(2, 2))
plot_grid_comparison(truth, res_tr_X3$m3, new_x,
                     effect = "fcont",
                     model = "M3",
                     truth_shift = covariate_shift,
                     theta = 40)
plot_grid_comparison(truth, res_tr_X3$m3, new_x,
                     effect = "fcont",
                     model = "M3",
                     truth_shift = covariate_shift,
                     theta = 130)
plot_grid_comparison(truth, res_tr_X3$m3, new_x,
                     effect = "fcont",
                     model = "M3",
                     truth_shift = covariate_shift,
                     theta = 220)
plot_grid_comparison(truth, res_tr_X3$m3, new_x,
                     effect = "fcont",
                     model = "M3",
                     truth_shift = covariate_shift,
                     theta = 310)

par(mfrow = c(2, 2))
plot_grid_comparison(truth, res_tr_X3$m3, new_x,
                     effect = "f",
                     region = 1,
                     model = "M3",
                     truth_shift = covariate_shift + mean(s_X3$scenario$gamma),
                     theta = 40)
plot_grid_comparison(truth, res_tr_X3$m3, new_x,
                     effect = "f",
                     region = 1,
                     model = "M3",
                     truth_shift = covariate_shift + mean(s_X3$scenario$gamma),
                     theta = 130)
plot_grid_comparison(truth, res_tr_X3$m3, new_x,
                     effect = "f",
                     region = 1,
                     model = "M3",
                     truth_shift = covariate_shift + mean(s_X3$scenario$gamma),
                     theta = 220)
plot_grid_comparison(truth, res_tr_X3$m3, new_x,
                     effect = "f",
                     region = 1,
                     model = "M3",
                     truth_shift = covariate_shift + mean(s_X3$scenario$gamma),
                     theta = 310)

# Contour / pheatmap plots ----
plot_grid_heatmap(truth, res_tr_X3$m3, new_x,
                  effect = "fcont",
                  model = "M3",
                  truth_shift = covariate_shift)

plot_grid_heatmap(truth, res_tr_X3$m3, new_x,
                  effect = "f",
                  region = 1,
                  model = "M3",
                  truth_shift = covariate_shift + mean(s_X3$scenario$gamma))