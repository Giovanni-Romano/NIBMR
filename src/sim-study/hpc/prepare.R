# Rscript .../prepare.R config.csv settings.R study-directory
source("src/sim-study/08-pipeline.R")
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L)
settings <- new.env(parent = globalenv())
sys.source(args[2L], settings)
prepare_study(read.csv(args[1L], stringsAsFactors = FALSE), settings$study_settings, args[3L])
cat("Prepared", args[3L], "\n")
