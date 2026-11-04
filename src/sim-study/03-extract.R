# === === === === === === === === === === === === === === === === === === ===
# 04-simulation/03-extract.R
# Load saved per-replicate `s`, apply effect_intervals (B1) and test_eval (B2),
# aggregate over replicates, write result tables.
# === === === === === === === === === === === === === === === === === === ===
source("04-simulation/02-worker.R")     # provides effect_intervals, test_eval, ...
OUT_DIR <- "sim-study/results/"; OUT_TAB <- "sim-study/results_tables/"

cell_path <- function(dgp, predictor, n, r) sprintf("%s_%s_n%d_rep%d", dgp, predictor, n, r)

rep_summary <- function(dgp, predictor, n, r, root = OUT_DIR) {
  s <- readRDS(file.path(root, cell_path(dgp, predictor, n, r), "cell.rds"))
  ci <- effect_intervals(s)                       # B1 point + HDI/quantile CIs
  te <- test_eval(s)                              # B2 predictions + scores
  list(b1 = ci, b2 = te)   # structured, aggregate-able
}

main <- function(argv = commandArgs(trailingOnly = TRUE),
                 R = 100, out = OUT_TAB) {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  grid <- expand.grid(dgp = paste0("D", 1:4), predictor = paste0("X", 1:3),
                      n = c(250, 1000), stringsAsFactors = FALSE)
  res <- lapply(seq_len(nrow(grid)), function(i)
    lapply(1:R, function(r)
      rep_summary(grid$dgp[i], grid$predictor[i], grid$n[i], r)))
  # aggregate B1 (mean point/coverage per cell) and B2 (mean RMSE/MAE, logp),
  # then write CSV/RDS tables.
}
# guarded main() so it only runs when executed as a script, not when sourced
if (sys.nframe() == 0L && !interactive()) main()