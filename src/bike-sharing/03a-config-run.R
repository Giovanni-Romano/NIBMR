# === === === === === === === === === === === === === === === === === === === ===
# 03a-config.R
# Shared paths, constants and small helpers for the 2025 "no bikes" logit model.
# Sourced by 03b / 03c / 03d. Contains no side effects other than creating
# output directories.
# === === === === === === === === === === === === === === === === === === === ===

# ---- Paths -------------------------------------------------------------------
PATH_RAW_DATA <- "data/bike-sharing/milan_bikemi_data_for_regression.rds"

MODEL_TAG <- "SpatCont_RandEff_linX"   # one knob controls where everything is written

DIR_PREP <- "application/bike-sharing/rds_prep/"
DIR_MCMC <- "application/bike-sharing/rds_mcmc/"
DIR_IMG  <- file.path("application/bike-sharing/pics/results", MODEL_TAG, "")

PATH_PREP <- file.path(DIR_PREP, paste0("prep_", MODEL_TAG, ".rds"))
PATH_MCMC <- file.path(DIR_MCMC, paste0("out_MCMC_", MODEL_TAG, ".rds"))

for (dd in c(DIR_PREP, DIR_MCMC, DIR_IMG)) {
  dir.create(dd, recursive = TRUE, showWarnings = TRUE)
}
rm(dd)

# Milan municipal boundary (downloaded once, then cached locally)
URL_MILANO_GEOJSON <- paste0(
  "https://dati.comune.milano.it/dataset/e75d91fa-eca6-4ee5-b96e-08bcdbb8d6f0/",
  "resource/f56cb432-83e6-48de-ae30-d39b4be61e85/download/",
  "confine_comune_milano_layer_0_confine_comune_milano.geojson"
)
PATH_MILANO_GEOJSON <- "data/bike-sharing/confine_comune_milano.geojson"

# ---- Data selection ----------------------------------------------------------
YEAR_KEEP        <- 2025
UNAVAIL_TYPE     <- "no_bikes"
# Stations close to Pero, excluded from the analysis.
# (500 has no observations in 2025 at all, so it needs no explicit filter.)
STATIONS_EXCLUDE <- c(334, 501)
TRANSF_CONT_X <- "identity"

# ---- Model settings ----------------------------------------------------------
SCALE_TIME    <- 100 * 720   # divisor mapping time_numeric_shifted into (0, 1)
N_DR_PER_COV  <- 5           # basis dimension per continuous covariate
SPAT_K1       <- 5
SPAT_K2       <- 5
EPS_BOUND     <- 1e-3        # clamp for y exactly at 0 or 1 before qlogis()

# ---- MCMC settings -----------------------------------------------------------
N_ITER  <- 20e3
BURNIN  <- 5000
SEED    <- 1708

# Inverse-gamma hyperparameters for the smoothing variances.
# Group 1: binary + linear continuous effects. Group 2: spline + spatial blocks.
A_TAU_LINEAR <- 2.1
B_TAU_LINEAR <- qgamma(0.001, shape = 2.1, rate = 1) * 5^2
A_TAU_SMOOTH <- 1e-3
B_TAU_SMOOTH <- 1e-3

C_POST <- 1e-1
C_PROP <- 1e-1
K_PROP_EXP <- 4 / 5          # k_prop = k_post^K_PROP_EXP

MCMC_STEP     <- 0.1
MCMC_K_FOLD   <- 5
MCMC_VERBOSE  <- 2
MCMC_PRINT_EVERY <- 250

# ---- Plotting ----------------------------------------------------------------
PDF_W <- 11
PDF_H <- 7

# === === === === === === === === === === === === === === === === === === === ===
# Helpers
# === === === === === === === === === === === === === === === === === === === ===

# Transformation of the covariates, based on TRANSF_CONT_X
transf_cont_x = switch(TRANSF_CONT_X,
                       "identity" = function(x) x,
                       "log" = function(x) log1p(x))

# Column blocks of the coefficient vector.
#
# The design matrix (and therefore `beta`) is laid out as
#   [intercept] [binary] [linear continuous] [DR spline blocks] [spatial]
# Every script derives its indices from here instead of hard-coding numbers.
coef_index <- function(q, p, d, s) {
  i_int  <- 1L
  i_bin  <- if (q > 0) seq.int(2L, 1L + q) else integer(0)
  i_lin  <- if (p > 0) seq.int(q + 2L, q + 1L + p) else integer(0)
  i_dr   <- if (p > 0) seq.int(q + p + 2L, q + p + 1L + sum(d)) else integer(0)
  i_spat <- if (s > 0) seq.int(q + p + sum(d) + 2L, length.out = s) else integer(0)
  
  # DR indices belonging to covariate j
  dr_block <- lapply(seq_len(p), function(j) {
    (q + p + 1L) + sum(d[seq_len(j - 1L)]) + seq_len(d[j])
  })
  
  list(intercept = i_int, binary = i_bin, linear = i_lin,
       dr = i_dr, dr_block = dr_block, spatial = i_spat)
}

# Block labels used to group coefficients in the traceplots / correlation plot.
# One id per coefficient: separate id per scalar effect, one id per DR block,
# and spatial coefficients chunked in groups of 5.
update_index <- function(q, p, d, s) {
  c(1L,
    if (q > 0) seq.int(2L, q + 1L) else integer(0),
    if (p > 0) seq.int(q + 2L, q + 1L + p) else integer(0),
    if (p > 0) q + 1L + p + rep(seq_len(p), d) else integer(0),
    if (s > 0) {
      q + 1L + 2L * p +
        c(rep(seq_len(floor(s / 5)), each = 5),
          rep(floor(s / 5) + (s %% 5 > 0), s %% 5))
    } else integer(0)
  )
}

# Drop the `id_re` bookkeeping column from a design matrix, by name.
drop_id_col <- function(X, id_col = "id_re") {
  X[, setdiff(colnames(X), id_col), drop = FALSE]
}

# Open a pdf device, evaluate `expr`, always close the device.
with_pdf <- function(file, expr, width = PDF_W, height = PDF_H) {
  grDevices::pdf(file, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(expr)
  invisible(file)
}

# Milan boundary as a matrix of UTM coordinates, cached on disk.
load_milano_utm <- function(url = URL_MILANO_GEOJSON,
                            path = PATH_MILANO_GEOJSON,
                            crs_utm = 32632) {
  if (!file.exists(path)) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    utils::download.file(url, path, quiet = TRUE)
  }
  sf::st_read(path, quiet = TRUE) |>
    sf::st_transform(crs = crs_utm) |>
    sf::st_coordinates()
}

# Rescale a 2-column coordinate matrix to [0, 1]^2 using supplied ranges.
rescale_coords <- function(coords, x_range, y_range) {
  out <- coords
  out[, 1] <- (coords[, 1] - x_range[1]) / diff(x_range)
  out[, 2] <- (coords[, 2] - y_range[1]) / diff(y_range)
  out
}

# Create labels for combination of station and workday/weekend
make_labels <- function(dat) {
  dat |>
    select(station, Weekend) |>
    mutate(Weekend = if_else(Weekend == 1, "Weekend", "Work")) |>
    apply(1, paste, collapse = " | ")
}