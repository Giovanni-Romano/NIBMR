# Shared console output for cells, extraction and the batch pipeline.
# verbose=FALSE hides our status lines; engine_output=TRUE restores engine prints.
run_status <- function(text, verbose = TRUE) {
  if (verbose) {
    cat(format(Sys.time(), "%H:%M:%S"), "  ", text, "\n", sep = "")
    flush(stdout())
  }
  invisible(NULL)
}

run_output <- function(expr, engine_output = FALSE) {
  if (engine_output) return(force(expr))
  # Discard routine stdout without accumulating a potentially large log in RAM.
  # Warnings/messages/errors use their normal channel and remain visible.
  depth <- sink.number()
  output <- file(nullfile(), "wt")
  sink(output)
  on.exit({
    while (sink.number() > depth) sink()
    close(output)
  }, add = TRUE)
  force(expr)
}

run_stage <- function(label, expr, verbose = TRUE, engine_output = FALSE) {
  run_status(paste0(label, ": starting"), verbose)
  start <- proc.time()[["elapsed"]]
  result <- run_output(expr, engine_output)
  run_status(sprintf("%s: completed in %.1f s", label,
                     proc.time()[["elapsed"]] - start), verbose)
  result
}
