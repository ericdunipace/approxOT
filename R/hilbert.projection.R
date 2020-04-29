hilbert.projection <- function(X, Sigma = NULL) {
  X <- as.matrix(t(X))
  if(!is.null(Sigma)) {
    X <- backsolve(chol(Sigma), X, transpose = TRUE)
  }
  idx <- hilbert_proj_(X) + 1 # +1 to adjust C to R indexing
  return(idx)
}