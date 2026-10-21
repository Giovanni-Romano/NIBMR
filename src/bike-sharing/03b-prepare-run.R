# === === === === === === === === === === === === === === === === === === === ===
# 03a-prepare-2025-logit.R
#
# Builds the design matrices for the 2025 "no bikes" unavailability model and
# saves a single `prep` object consumed by 03c (run) and 03d (results).
#
# Input : data/bike-sharing/milan_bikemi_data_for_regression.rds
# Output: application/bike-sharing/rds_prep/prep_<MODEL_TAG>.rds
# === === === === === === === === === === === === === === === === === === === ===

suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
})

source("src/bike-sharing/03a-config-run.R")
source("src/01a-utils.R")
source("src/01b-gibbs.R")

# ---- Load and filter ---------------------------------------------------------
raw <- readRDS(PATH_RAW_DATA) |>
  mutate(
    Weekend = as.numeric(weekdays(commit_at) %in% c("Saturday", "Sunday")),
    Rain    = as.numeric(Precipitation > 0)
  )

D <- raw |>
  filter(
    year(commit_at) == YEAR_KEEP,
    unavail_type == UNAVAIL_TYPE,
    !ID %in% STATIONS_EXCLUDE
  ) |>
  select(-commit_at, -Precipitation, -Temperature, -Rain, -starts_with("count_"))

str(D)

# ---- Response, covariates, random-effect id ----------------------------------
y <- D$time_numeric_shifted / SCALE_TIME

X_dummy <- D |> select(Weekend) |> as.matrix()

X_cont <- D |>
  select(population_density, starts_with("density_"), starts_with("closest_")) |>
  as.matrix() |> 
  transf_cont_x()

X_cont_scaled <- scale(X_cont)
scaling_factors <- list(center = attr(X_cont_scaled, "scaled:center"),
                        scale = attr(X_cont_scaled, "scaled:scale"))

id_re <- as.integer(factor(D$station, levels = unique(D$station[order(D$ID)])))

n <- length(y)
q <- ncol(X_dummy)
p <- ncol(X_cont)
d <- rep(N_DR_PER_COV, p)

# ---- DR spline bases for the continuous covariates ---------------------------
splines_list <- lapply(seq_len(p), function(j) {
  construct_DR_basis(X_cont_scaled[, j], d = d[j], type_of_splines = "gps")
})

Z_cont <- do.call(cbind, lapply(splines_list, `[[`, "Z"))
colnames(Z_cont) <- unlist(lapply(seq_len(p), function(j) {
  paste0(colnames(X_cont)[j], "_DR", seq_len(d[j]))
}))

# ---- Spatial tensor-product basis --------------------------------------------
# Built on the unique station coordinates, then matched back to the rows of the
# data by station name.
X_spat_unique <- D |> select(station, longitude, latitude) |> distinct()

coord_utm <- X_spat_unique |>
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
  st_transform(crs = 32632) |>
  st_coordinates()

# Normalisation to [0, 1]^2.
#   "observed": min/max of the observed station coordinates  (current choice)
#   "milano"  : min/max of the Milan municipal boundary
COORD_NORM <- "observed"

coord_range <- switch(
  COORD_NORM,
  observed = list(x = range(coord_utm[, 1]), y = range(coord_utm[, 2])),
  milano   = {
    mc <- load_milano_utm()
    list(x = range(mc[, 1]), y = range(mc[, 2]))
  },
  stop("Unknown COORD_NORM: ", COORD_NORM)
)

coord_scaled <- rescale_coords(coord_utm, coord_range$x, coord_range$y)

out_spatial_DR <- construct_DR_Spatial(
  coord_scaled[, 1], coord_scaled[, 2], k1 = SPAT_K1, k2 = SPAT_K2
)

Z_spat_unique <- bind_cols(
  station = X_spat_unique$station,
  as.data.frame(out_spatial_DR$Z)
)
colnames(Z_spat_unique) <- c(
  "station", 
  paste0("spatial_DR", seq_len(ncol(out_spatial_DR$Z)))
)
s <- ncol(Z_spat_unique) - 1L

# ---- Assemble, then collapse to unique covariate patterns --------------------
data_for_regression <- bind_cols(
  station = D$station, id_re = id_re, y = y,
  X_dummy, X_cont_scaled, Z_cont
) |>
  left_join(Z_spat_unique, by = "station")

data_grouped <- data_for_regression |>
  group_by(across(everything())) |>
  summarise(freq = n(), .groups = "drop")

y_grouped             <- data_grouped$y
id_re_grouped         <- data_grouped$id_re
freq_grouped          <- data_grouped$freq
X_dummy_grouped       <- data_grouped |> select(all_of(colnames(X_dummy)))       |> as.matrix()
X_cont_scaled_grouped <- data_grouped |> select(all_of(colnames(X_cont_scaled))) |> as.matrix()
Z_cont_grouped        <- data_grouped |> select(all_of(colnames(Z_cont)))        |> as.matrix()
Z_spat_grouped        <- data_grouped |> select(all_of(colnames(Z_spat_unique)[-1])) |> as.matrix()

stopifnot(
  sum(freq_grouped) == n,
  ncol(X_dummy_grouped)       == q,
  ncol(X_cont_scaled_grouped) == p,
  ncol(Z_cont_grouped)        == sum(d),
  ncol(Z_spat_grouped)        == s
)

# ---- Full design matrix and labels -------------------------------------------
# `df` keeps id_re as a column because lqmm() needs it as a grouping variable;
# it is dropped by name (never by position) wherever a design matrix is needed.
df_lqmm <- data_grouped |> select(-station, -freq) |> as.data.frame()

# qlogis() is undefined at the boundary, so clamp the grouped response here.
# NOTE: `y` and `y_grouped` are deliberately left unclamped, they are the data
# the MCMC and the histograms use.
df_lqmm$y[df_lqmm$y == 0] <- EPS_BOUND
df_lqmm$y[df_lqmm$y == 1] <- 1 - EPS_BOUND

X_complete_grouped <- cbind("Intercept" = 1, df_lqmm |> select(-y) |> as.matrix())
X_complete_unique  <- unique(X_complete_grouped)

labs         <- make_labels(data_for_regression)
labs_grouped <- make_labels(data_grouped)
labs_unique  <- unique(labs_grouped)

names(y)         <- labs
names(y_grouped) <- labs_grouped

labs_re <- data_grouped |> select(station, id_re) |> distinct() |> arrange(id_re)

# ---- Save preparation object -------------------------------------------------
prep <- list(
  D            = D,
  dims         = list(n = n, p = p, q = q, s = s, d = d),
  idx          = coef_index(q, p, d, s),
  update_idx   = update_index(q, p, d, s),
  y            = y,
  grouped      = list(
    y       = y_grouped,
    id_re   = id_re_grouped,
    freq    = freq_grouped,
    X_dummy = X_dummy_grouped,
    X_cont  = X_cont_scaled_grouped,
    Z_cont  = Z_cont_grouped,
    Z_spat  = Z_spat_grouped
  ),
  df_lqmm            = df_lqmm,
  X_complete_grouped = X_complete_grouped,
  X_complete_unique  = X_complete_unique,
  Z_spat_unique      = Z_spat_unique,
  spatial = list(
    coord_utm      = coord_utm,
    coord_scaled   = coord_scaled,
    coord_range    = coord_range,
    out_spatial_DR = out_spatial_DR
  ),
  scaling_factors = scaling_factors,
  labels = list(
    obs     = labs,
    grouped = labs_grouped,
    unique  = labs_unique,
    re      = labs_re
  )
)

saveRDS(prep, PATH_PREP)
message("Saved preparation object to ", PATH_PREP,
        "\n(n = ", n, ", p = ", p, ", q = ", q, ", s = ", s,
        ", # stations = ", nrow(labs_re),
        ", # unique patterns = ", nrow(X_complete_unique), ")")