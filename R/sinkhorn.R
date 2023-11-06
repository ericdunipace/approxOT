# sinkhorn_distance <- function(x1, x2, p = 1, ground_p = 2, eps = 0.05, niter = 100){
#   w1 <- rep(1/ncol(x1), ncol(x1))
#   w2 <- rep(1/ncol(x2), ncol(x2))
#   C <- cost_calc(x1, x2, ground_p)
#   epsilon <- eps * median(C^p)
#   wass <- sinkhorn_(w1, w2, C^p, epsilon, niter)
#   ### CORRECTION OF THE MARGINALS
#   # explained in the appendix of Coupling of Particle Filters, Jacob Lindsten Schon  (arXiv v2 appendix E)
#   Phat <- wass$transportmatrix
#   u <- rowSums(Phat)
#   utilde <- colSums(Phat)
#   alpha <- min(pmin(w1/u, w2/utilde))
#   r <- (w1 - alpha * u) / (1 - alpha)
#   rtilde <- (w2 - alpha * utilde) / (1 - alpha)
#   P <- alpha * Phat + (1 - alpha) * matrix(r, ncol = 1) %*% matrix(rtilde, nrow = 1)
#   return(list(uncorrected = (sum(Phat * C^p))^(1/p), corrected = (sum(P * C^p))^(1/p)))
# }

#' Test sinkhorn distance
#'
#' @param mass_x empiric measure of first sample
#' @param mass_y empiric measure of second sample
#' @param cost cost matrix
#' @param p power to raise the cost matrix by
#' @param eps epsilon of cost matrix
#' @param niter number of iterations
#'
#' @return a numeric value
#'
#' @keywords internal
sinkhorn_distance <- function(mass_x, mass_y, cost = NULL, p = 1, eps = 0.05, niter = 100){
  costp <- cost^p
  epsilon <- eps * stats::median(costp)
  wass <- sinkhorn_(mass_x, mass_y, costp, epsilon, niter)
  ### CORRECTION OF THE MARGINALS
  # explained in the appendix of Coupling of Particle Filters, Jacob Lindsten Schon  (arXiv v2 appendix E)
  Phat <- wass$transportmatrix
  u <- rowSums(Phat)
  utilde <- colSums(Phat)
  alpha <- min(pmin(mass_x/u, mass_y/utilde))
  r <- (mass_x - alpha * u) / (1 - alpha)
  rtilde <- (mass_y - alpha * utilde) / (1 - alpha)
  P <- if ( alpha < 1 ) {
    alpha * Phat + (1 - alpha) * matrix(r, ncol = 1) %*% matrix(rtilde, nrow = 1)
  } else {
    Phat
  }
  return(list(uncorrected = (sum(Phat * costp))^(1/p), corrected = (sum(P * costp))^(1/p)))
}

#' Test sinkhorn transportation plan
#'
#' @param mass_x empiric measure of first sample
#' @param mass_y empiric measure of second sample
#' @param cost cost matrix
#' @param eps epsilon of cost matrix
#' @param niterations number of iterations
#'
#' @return transportation plan as list with slots "from","to", and "mass"
#'
#' @keywords internal
sinkhorn_transport <- function(mass_x, mass_y, cost = NULL, eps = 0.05, niterations = 100){
  n1 <- length(mass_x)
  n2 <- length(mass_y)
  # costp <- cost^p
  epsilon <- eps * stats::median(cost)
  transp <- sinkhorn_(mass_x, mass_y, cost, epsilon, niterations)
  ### CORRECTION OF THE MARGINALS
  # explained in the appendix of Coupling of Particle Filters, Jacob Lindsten Schon  (arXiv v2 appendix E)
  Phat <- transp$transportmatrix
  u <- rowSums(Phat)
  utilde <- colSums(Phat)
  alpha <- min(pmin(mass_x/u, mass_y/utilde))
  r <- (mass_x - alpha * u) / (1 - alpha)
  rtilde <- (mass_y - alpha * utilde) / (1 - alpha)
  P <- alpha * Phat + (1 - alpha) * matrix(r, ncol = 1) %*% matrix(rtilde, nrow = 1)
  return(list(from = rep(1:n1, n2), 
              to = rep(1:n1, each = n2), 
              mass = c(P)))
}


col_logsumexp <- function(mat) {
  maxes <- apply(mat,2,max)
  return(
    maxes + log(colSums(exp(sweep(mat, 2, maxes))))
  )
}

row_logsumexp <- function(mat) {
  maxes <- apply(mat,1,max)
  return(
    maxes + log(rowSums(exp(sweep(mat, 1, maxes))))
  )
}

log_sinkhorn_test <- function(mass_x, mass_y, cost = NULL, eps = 0.05, niterations = 100){
  update_g <- function(cost, f, lambda, log_b) {
    -lambda * col_logsumexp(-sweep(cost, 1, f)/lambda) - lambda * log_b
  }
  update_f <- function(cost, g, lambda, log_a) {
    -lambda * row_logsumexp(-sweep(cost, 2, g)/lambda) - lambda * log_a
  }
  converge_log <- function(pot, pot_old, tol) {
    return(isTRUE(sum(abs(pot - pot_old)/abs(pot_old)) < tol))
  }
  n1 <- length(mass_x)
  n2 <- length(mass_y)
  log_x <- log(mass_x)
  log_y <- log(mass_y)
  
  # costp <- cost^p
  lambda <- eps * stats::median(cost)
  f_old <- rep(0,n1)
  g_old <- rep(0,n2)
  
  f <- -lambda * row_logsumexp(-cost/lambda) + lambda * log_x
  g <- update_g(cost, f, lambda, log_y)
  for( i in 1:niterations) {
    f <- update_f(cost, g, lambda, log_x)
    g <- update_g(cost, f, lambda, log_y)
    
    if (converge_log(f, f_old, tol = 1e-8)) {
      break
    } else {
      f_old <- f
      g_old <- g
    }
  }
  
  return(list(f = f, g = g))
}


col_softmin <- function(mat) {
  -log(colSums(exp(mat)))
}

row_softmin <- function(mat) {
  -log(rowSums(exp(mat)))
}

log_sinkhorn_test_nomax <- function(mass_x, mass_y, cost = NULL, eps = 0.05, niterations = 100){
  generate_S <- function(cost, f, g, lambda) {
    S <- -sweep(sweep(cost, 1, f), 2, g)/lambda
    return(S)
  }
  update_g <- function(cost, f, g, lambda, log_a, log_b) {
    S <- generate_S(cost, f, g, lambda)
    return(
      lambda * col_softmin(S) + g + lambda * log_b
    )
  }
  update_f <- function(cost, f, g, lambda, log_a, log_b) {
    S <- generate_S(cost, f, g, lambda)
    return(
      lambda * row_softmin(S) + f + lambda * log_a
    )
  }
  converge_log <- function(pot, pot_old, tol) {
    return(isTRUE(sum(abs(pot - pot_old))/sum(abs(pot_old)) < tol))
  }
  n1 <- length(mass_x)
  n2 <- length(mass_y)
  log_x <- log(mass_x)
  log_y <- log(mass_y)
  
  # costp <- cost^p
  lambda <- eps * stats::median(cost)
  f <- -lambda * row_logsumexp(-cost/lambda) + lambda * log_x
  f_old <- rep(0,n1)
  g <- -lambda * col_logsumexp(-sweep(cost,1,f)/lambda) + lambda * log_y
  g_old <- rep(0,n2)
  
  for( i in 1:niterations) {
    f <- update_f(cost, f, g, lambda, log_x, log_y)
    g <- update_g(cost, f, g, lambda, log_x, log_y)
    
    if(any(is.nan(f)) || any(is.nan(g))) browser()
    if (converge_log(f, f_old, tol = 1e-8)) {
      break
    } else {
      f_old <- f
      g_old <- g
    }
  }
  
  return(list(f = f, g = g))
}

log_sinkhorn_test_nomax_KL <- function(mass_x, mass_y, cost = NULL, eps = 0.05, niterations = 100){
  generate_S <- function(cost, f, g, lambda) {
    S <- -sweep(sweep(cost, 1, f), 2, g)/lambda
    return(S)
  }
  update_g <- function(cost, f, g, lambda, log_a, log_b) {
    S <- generate_S(cost, f, g, lambda)
    return(
      lambda * col_softmin(sweep(S,1,log_a, "+")) + g 
    )
  }
  update_f <- function(cost, f, g, lambda, log_a, log_b) {
    S <- generate_S(cost, f, g, lambda)
    return(
      lambda * row_softmin(sweep(S,2,log_b, "+")) + f
    )
  }
  converge_log <- function(pot, pot_old, tol) {
    return(isTRUE(sum(abs(pot - pot_old))/sum(abs(pot_old)) < tol))
  }
  n1 <- length(mass_x)
  n2 <- length(mass_y)
  log_x <- log(mass_x)
  log_y <- log(mass_y)
  
  # costp <- cost^p
  lambda <- eps * stats::median(cost)
  f <- -lambda * row_logsumexp(sweep(-cost/lambda,2,log_y,"+")) 
  f_old <- rep(0,n1)
  g <- -lambda * col_logsumexp(sweep(-sweep(cost,1,f)/lambda,1,log_x,"+"))
  g_old <- rep(0,n2)
  
  for( i in 1:niterations) {
    f <- update_f(cost, f, g, lambda, log_x, log_y)
    g <- update_g(cost, f, g, lambda, log_x, log_y)
    
    if(any(is.nan(f)) || any(is.nan(g))) browser()
    if (converge_log(f, f_old, tol = 1e-8)) {
      break
    } else {
      f_old <- f
      g_old <- g
    }
  }
  
  return(list(f = f, g = g))
}


#' Round transportation matrix to feasible set
#'
#' @param transport_matrix A transportation matrix returned by an approximate method
#' @param mass_x The distribution of the first margin
#' @param mass_y The distribution of the second margin
#'
#' @return Returns a transportation matrix projected to the feasible set.
#' @keywords internal
round_transport_matrix <- function(transport_matrix, mass_x, mass_y) {
  # set.seed(32423)
  # n <- 100
  # d <- 10
  # x <- matrix(rnorm(d*n), nrow=d, ncol=n)
  # y <- matrix(rnorm(d*n), nrow=d, ncol=n)
  # mass_x <- rep(1/n,n)
  # mass_y <- rep(1/n,n)
  # cost <- approxOT::cost_calc(x,y, 2.0)
  # tpot <- sinkhorn_pot(mass_x, mass_y, p = 2,
  #            cost=cost)
  # tmat <- exp(sweep(sweep(cost^2,1,tpot$f ),2,tpot$g)/(0.05 * median(cost^2)))
  # tmat <- tmat/sum(tmat)
  # 
  # rounded <- approxOT:::round_transport_matrix(tmat, mass_x = mass_x,
  # mass_y = mass_y)
  # all.equal(rowSums(rounded), mass_x)
  transport_matrix <- as.matrix(transport_matrix)
  stopifnot(nrow(transport_matrix) == length(mass_x))
  stopifnot(ncol(transport_matrix) == length(mass_y))
  a <- as.double(mass_x)
  b <- as.double(mass_y)
  
  
  tmat <- round_2_feasible_(transport_matrix, a, b)
  return(tmat)
}