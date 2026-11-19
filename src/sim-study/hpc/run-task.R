# Rscript .../run-task.R study-directory config.csv row-index
# Slurm selects a CSV row; the scenario values then go into run_cell().
source("src/sim-study/08-pipeline.R")
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L)
config <- read.csv(args[2L], stringsAsFactors = FALSE)
row <- as.integer(args[3L]); stopifnot(row >= 1L, row <= nrow(config))
manifest <- readRDS(file.path(args[1L], "manifest.rds"))
index <- match(study_key(config[row, ]), study_key(manifest$config))
stopifnot(!is.na(index))
run_study_task(args[1L], index)
