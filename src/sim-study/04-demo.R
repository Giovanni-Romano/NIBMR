# === === === === === === === === === === === === === === === === === === ===
# 04-simulation/04-demo.R
# One-cell end-to-end demo + diagnostics.
#   * run_cell on ONE disposable cell            (D1, X1, n = 250, rep = 1)
#   * print diagnostics needed to finalize:
#       (a) bamlss predict(type="term") term names    -> one_model_bamlss
#       (b) bamlss predict(type="parameter") params   -> M1 posterior-predictive
#       (c) loo::loo() success for all 3 GUD members  (log_lik exposure)
#       (d) HDInterval::hdi returns a single interval
#   * run effect_intervals(s) and test_eval(s), print concise summaries
#   * write_cell(s) and confirm the artifacts land
# === === === === === === === === === === === === === === === === === === ===
source("src/sim-study/02-worker.R")

cat("\n==== running run_cell(D1, X1, n=250, rep=1) ====\n")
s <- run_cell(dgp = "D3", predictor = "X1", n = 250, replicate = 1)
dat <- s$scenario$dat

cat("\n-- fit structures --\n")
cat("m0 draws:", paste(dim(s$fits$m0$draws), collapse = " x "), "\n")
cat("m1 draws:", paste(dim(s$fits$m1$draws), collapse = " x "), "\n")
cat("m2 draws:", paste(dim(s$fits$m2$draws), collapse = " x "), "\n")
cat("m3 beta :", length(s$fits$m3$beta), " coefs, conv =", s$fits$m3$convergence, "\n")
cat("m4 members:", names(s$fits$m4$draws), " | selected:", s$fits$m4$selected, "\n")

cat("\n======== (a) bamlss predict(type='term') term names ========\n")
trm <- predict(s$fits$m1$fit, newdata = dat, model = "mu", type = "term",
               FUN = function(x) x)
print(names(trm)); print(lapply(trm, dim))   # <-- we need these names to finalize one_model_bamlss

cat("\n======== (b) bamlss predict(type='parameter') param names ========\n")
prm <- predict(s$fits$m1$fit, newdata = dat, type = "parameter")
print(names(prm)); print(lapply(prm, dim))   # mu? sigma? -> for M1 posterior-predictive

cat("\n======== (c) loo() for the 3 GUD members ========\n")
gud_loo <- sapply(names(s$fits$m4$draws), function(m) {
  tryCatch(loo::loo(s$fits$m4$... )$estimates["elpd_loo", "Estimate"],
           error = function(e) paste("ERR:", conditionMessage(e)))
})
# NOTE: the loo call needs the stanfit of each member (m4 retains stanfits?).
# We print whether it succeeds for TPSC/FG/DTP.

cat("\n======== (d) HDInterval::hdi single-interval check ========\n")
h <- HDInterval::hdi(rnorm(5000), credMass = 0.95)
cat("length(h) =", length(h), "(1  => single HPD interval; >1  => union)\n")

cat("\n======== effect_intervals(s) ========\n")
ci <- effect_intervals(s)              # provisional (may need the term names)
cat("M0 effects:", paste(names(ci$M0), collapse = ", "), "\n")
cat("M0 whole point[1:3]:", round(ci$M0$whole$point[1:3], 3), "\n")

cat("\n======== test_eval(s) ========\n")
te <- test_eval(s)
print(te$pred)                 # per-model RMSE/MAE on the 101st test set
print(te$logp)                 # M1, M4 posterior-predictive log scores

cat("\n======== write_cell(s) ========\n")
write_cell(s, dgp = "D1", predictor = "X1", n = 250, replicate = 1)
cat("cell dir:", file.path(OUT_DIR, "D1_X1_n250_rep1"), "\n")
cat("artifacts:", paste(list.files(file.path(OUT_DIR, "D1_X1_n250_rep1")),
                        collapse = ", "), "\n")