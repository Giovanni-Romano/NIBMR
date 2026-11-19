# Rscript .../collect.R study-directory
source("src/sim-study/08-pipeline.R")
source("src/sim-study/06-evaluate.R")
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
collect_study(args[1L])
