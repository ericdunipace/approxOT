.onLoad <- function(libname, pkgname) {
  # thread count for OpenMP loops: option "approxOT.threads" if set,
  # 2 under R CMD check core limits, otherwise half the processors (n = 0)
  n <- getOption("approxOT.threads", NA_integer_)
  if (is.null(n) || is.na(n)) {
    lim <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
    n <- if (nzchar(lim) && !identical(toupper(lim), "FALSE")) 2L else 0L
  }
  set_threads_(as.integer(n))
  invisible()
}
