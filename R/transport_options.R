transport_options <- function() {
  return(c("exact", "networkflow","shortsimplex",
           "sinkhorn", "greenkhorn",
           "randkhorn", "gandkhorn",
           "hilbert", "rank", "univariate",
           "univariate.approximation.pwr",
           "swapping", "sliced"))
}