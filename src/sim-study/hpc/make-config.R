# Run from R/: Rscript src/sim-study/hpc/make-config.R [output.csv]
source("src/sim-study/08-pipeline.R")
args <- commandArgs(trailingOnly = TRUE)
path <- if (length(args)) args[1L] else "src/sim-study/hpc/runs.csv"
write.csv(study_config(), path, row.names = FALSE)
cat("Wrote 2,400 runs to", path, "\n")
