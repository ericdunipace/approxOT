#' Number of threads used by approxOT
#'
#' Get or set the number of OpenMP threads used by approxOT's parallel
#' code (currently the log-scale Sinkhorn algorithm).
#'
#' @param n Number of threads to use. `NULL` (the default) only returns the
#' current setting. `0` restores the default of half the available processors.
#'
#' @details
#' When the package loads, the number of threads is taken from
#' `getOption("approxOT.threads")` if that is set. Otherwise it is 2 during
#' `R CMD check` (when `_R_CHECK_LIMIT_CORES_` is set) and half the
#' available processors in all other cases.
#'
#' The number of threads never exceeds the number of available processors
#' or the `OMP_THREAD_LIMIT` environment variable. If approxOT was built
#' without OpenMP support, it is always 1.
#'
#' @return If `n` is `NULL`, the number of threads in use. Otherwise the
#' previous number of threads, invisibly, so that it can be restored later.
#' @export
#'
#' @examples
#' approxOT_threads()
#' old <- approxOT_threads(1)
#' approxOT_threads(old)
approxOT_threads <- function(n = NULL) {
  if (is.null(n)) return(get_threads_())
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 0) {
    stop("`n` must be a single non-negative number")
  }
  old <- get_threads_()
  set_threads_(as.integer(n))
  invisible(old)
}
