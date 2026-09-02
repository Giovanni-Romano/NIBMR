#' Function immitating \code{\link[base]{try}} in parallelized setting
#'
#' This function immitates the behaviour of \code{\link[base]{try}} in parallelized settings using
#' \code{\link[pbmcapply]{pbmclapply}} / \code{\link[pbmcapply]{pbmcmapply}}.
#'
#' @param f expression to be wrapped
#' @export

try_parallel <- function(f) {function(...) {tryCatch({f(...)}, error = function(e) e)}}
