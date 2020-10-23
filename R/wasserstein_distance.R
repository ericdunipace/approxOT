wasserstein <- function (X, Y, p = 2, ground_p = 2, observation.orientation = c("rowwise","colwise"), 
                         method = c("exact", "sinkhorn", "greenkhorn",
                                    "randkhorn", "gandkhorn",
                                    "hilbert", "rank", "sinkhorn2",
                                    "univariate.approximation", 
                                    "univariate.approximation.pwr","univariate",
                                    "sliced"), ... ) {
  obs <- match.arg(observation.orientation,  c("colwise","rowwise"))
  method <- match.arg(method)
  
  if(missing(X)) stop("Must specify X")
  if(missing(Y)) stop("Must specify Y")
  
  if (!is.matrix(X)) {
    # warning("Attempting to coerce X to a matrix")
    X <- as.matrix(X)
    if(dim(X)[2] == 1 & obs == "colwise") X <- t(X)
  }
  if (!is.matrix(Y)) {
    # warning("Attempting to coerce Y to a matrix")
    Y <- as.matrix(Y)
    if(dim(Y)[2] == 1 & obs == "colwise") Y <- t(Y)
  }
  p <- as.double(p)
  ground_p <- as.double(ground_p)
  
  if(!(p >= 1)) stop("p must be >= 1")
  
  if(obs == "rowwise"){
    X <- t(X)
    Y <- t(Y)
    obs <- "colwise"
  }
  stopifnot(nrow(X) == nrow(Y))
  stopifnot(all(is.finite(X)))
  stopifnot(all(is.finite(Y)))
  
  if (method == "univariate.approximation" ) {
    loss <- wasserstein_p_iid_(X,Y, p)
  } else if ( method == "univariate.approximation.pwr") {
    loss <- wasserstein_p_iid_p_(X,Y, p)
  } else if (method == "univariate" | method == "hilbert" | method == "rank") {
    tp <- transport_plan(X = X, Y = Y, p = p, ground_p = ground_p,
                         observation.orientation = obs, method = method, ...)
    # loss <- c((((colSums(abs(X[, tp$tplan$from, drop = FALSE] - Y[, tp$tplan$to, drop=FALSE])^ground_p))^(1/ground_p))^p %*% tp$tplan$mass)^(1/p))
    loss <- tp$cost
  } else if (method == "sliced") {
    dots <- list(...)
    tplan <- NULL
    nboot <- dots$nsim
    d     <- nrow(X)
    theta <- matrix(rnorm(d * nboot), d, nboot)
    theta <- sweep(theta, 2, STAT=apply(theta,2,function(x) sqrt(sum(x^2))), FUN = "/")
    X_theta <- crossprod(x = X, y = theta)
    Y_theta <- crossprod(x = Y, y = theta)
    # u     <- sort(runif(nboot))
    costs <- sapply(1:nboot, function(i) {
      # x <- quantile(c(X_theta[,i]), probs = u)
      # y <- quantile(c(Y_theta[,i]), probs = u)
      x <- c(X_theta[,i])
      y <- c(Y_theta[,i])
      trans <- general_1d_transport(t(x),t(y),"univariate")
      cost <- ((sum(abs(x[trans$from] - y[trans$to])^ground_p))^(1/ground_p))^p %*% tplan$mass
      return(cost)
    }
    )
    loss <- mean(costs)^(1/p)
  } else {
    n1 <- ncol(X)
    n2 <- ncol(Y)
    mass_x <- as.double(rep(1/n1, n1))
    mass_y <- as.double(rep(1/n2, n2))
    
    cost <- cost_calc(X, Y, ground_p)
    
    # if (method == "exact") {
      
      tplan <- transport_plan_given_C(mass_x, mass_y, p, cost, method, ...)
      loss <- wasserstein_(tplan$mass, cost, p, from=tplan$from, to = tplan$to)
      
    # } else if (method == "sinkhorn") {
    #   dots <- list(...)
    #   eps <- dots$eps
    #   niter <- dots$niter
    #   if (is.null(eps)) eps <- 0.05
    #   if (is.null(niter)) niter <- 100
    #   
    #   sink_out <- sinkhorn_distance(mass_x, mass_y, cost, p, eps, niter)
    #   loss <- sink_out$corrected
    # }
    
  }
  
  return(loss)
}

# wasserstein_post <- function(mod = NULL, idx= NULL, coefs = NULL, post = NULL, p = 2, n.param=NULL){
#   if(!is.null(mod)){
#     ml <- minLambda(mod)
#     idx <- ml$nzero
#     coefs <- ml$coefs
#   }
#   if(!is.matrix(coefs)) coefs <- as.matrix(coefs)
#   p_post <- transport::pp(post)
#   if(p < 1) stop("p must be greater or equal to 1")
#   
#   sink("temp_asterix1234.txt") #transport uses annoying Rprintf which can't be diverted normally
#   loss <- apply(coefs, 2, function(gamma) 
#     transport::wasserstein(transport::pp(post %*% diag(gamma)), 
#                            p_post, p=p, method="shortsimplex"))
#   sink()
#   file.remove("temp_asterix1234.txt")
#   
#   out <- rep(NA, n.param)
#   out[idx] <- loss^p
#   return(out)
# }

wasserstein_individual <- function(X,Y, ground_p, observation.orientation = c("colwise","rowwise")) {
  if(!is.matrix(X)) X <- as.matrix(X)
  if(!is.matrix(Y)) Y <- as.matrix(Y)
  obs <- match.arg(observation.orientation)
  if(obs == "rowwise"){
    X <- t(X)
    Y <- t(Y)
  }
  
  Xs <- apply(X,2,sort)
  Ys <- apply(Y,2,sort)
  
  loss <- colMeans((Xs - Ys)^ground_p)
  
  return(loss^(1/p))
  
}

general_dist <- function(X, Y) {
  idx_x <- hilbert_proj_(X) + 1
  idx_y <- hilbert_proj_(Y) + 1
  n <- ncol(X)
  m <- ncol(Y)
  
  mass_a <- rep(1/n, n)
  mass_b <- rep(1/m, m)
  cum_a <- c(cumsum(mass_a))[-n]
  cum_b <- c(cumsum(mass_b))[-m]
  mass <- diff(c(0,sort(c(cum_a, cum_b)),1))
  cum_m <- cumsum(mass)
  arep <- table(cut(cum_m, c(-Inf, cum_a, Inf)))
  brep <- table(cut(cum_m, c(-Inf, cum_b, Inf)))
  a_idx <- rep(idx_x, times = arep)
  b_idx <- rep(idx_y, times = brep)
  
  transport <- list(from = a_idx[order(b_idx)], to = sort(b_idx), mass = mass[order(b_idx)])
  
  return(transport)
  # test <- data.frame(from = a_idx, to = b_idx, mass = mass)
}


wasserstein_multimarg <- function (..., p = 2, ground_p = 2, observation.orientation = c("rowwise","colwise"), 
                         method = c("hilbert", "univariate")) {
  
  if (method == "univariate" | method == "hilbert" ) {
    tp <- transport_plan_multimarg(..., p = p, ground_p = ground_p,
                         observation.orientation = obs, method = method)
    # loss <- c((((colSums(abs(X[, tp$tplan$from, drop = FALSE] - Y[, tp$tplan$to, drop=FALSE])^ground_p))^(1/ground_p))^p %*% tp$tplan$mass)^(1/p))
    loss <- tp$cost
  } else {
    
    stop("Transport method", method, "not currently supported for multimarginal problems.")
    
  }
  
  return(loss)
}