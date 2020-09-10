transport_options <- function() {
  return(c("exact", "sinkhorn", "greenkhorn",
           "randkhorn", "gandkhorn",
           "hilbert", "rank", "univariate",
           "univariate.approximation.pwr",
           "swapping"))
}