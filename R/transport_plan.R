transport_plan_given_C <- function(mass_x, mass_y, p = 2, 
                                   cost=NULL, method = "exact", ...) {
  method <- match.arg(method, c("exact","sinkhorn","greenkhorn", 
                                "randkhorn", "gandkhorn", "sinkhorn2"))
  
  dots <- list(...)
  epsilon <- as.double(dots$epsilon)
  niter <- as.integer(dots$niter)
  stopifnot(all(is.finite(cost)))
  
  if(length(epsilon) == 0) epsilon <- as.double(0.05)
  if(length(niter) == 0) niter <- as.integer(100)
  
  if (is.null(cost) ) stop("Cost matrix must be provided")
  tplan <- if (method == "exact" | method == "greenkhorn" | method == "sinkhorn" |
               method == "randkhorn" | method == "gandkhorn") {
    
    n1 <- length(mass_x)
    n2 <- length(mass_y)
    
    if(n1 > 1 & n2 > 1) {
      transport_C_(mass_a_ = mass_x, mass_b_ = mass_y, cost_matrix_ = cost^p, 
                   method_ = method, epsilon_ = epsilon, niter_ = niter)
    } else if (n2 == 1) {
      list(from = 1:n1, to = rep(1,n1), mass = mass_x)
    } else if (n1 == 1) {
      list(from = rep(1,n2), to = 1:n2, mass = mass_y)
    } else {
      stop("Some error found in mass_x or mass_y length. Check mass input.")
    }
    
  } else if (method == "sinkhorn2") {
    
    sinkhorn_transport(mass_x = mass_x, mass_y = mass_y, cost = cost^p, 
                       eps = epsilon, niter = niter)
    
  } else {
    stop( paste0( "Transport method ", method, " not supported" ) )
  }
  return( tplan )
  
}

transport_plan <- function(X, Y, p = 2, ground_p = 2,
                           observation.orientation = c("rowwise", "colwise"), 
                           method = transport_options(),... ) {
  
  obs <- match.arg(observation.orientation)
  method <- match.arg(method)
  
  if (!is.matrix(X)) {
    X <- as.matrix(X)
    if(dim(X)[2] == 1) X <- t(X)
  }
  if (!is.matrix(Y)) {
    Y <- as.matrix(Y)
    if(dim(Y)[2] == 1) Y <- t(Y)
  }
  
  p <- as.double(p)
  ground_p <- as.double(ground_p)
  
  if(obs == "rowwise"){
    X <- t(X)
    Y <- t(Y)
  }
  stopifnot(all(is.finite(X)))
  stopifnot(all(is.finite(Y)))
  
  cost <- tplan <- NULL
  if (method == "univariate.approximation") {
    tplan <-  list(from = apply(X, 1, order), to = apply(Y,1,order), mass = rep(1/ncol(X), ncol(X)))
    cost <- sapply(1:nrow(X), function(i) 
      sum((X[i, tplan$from[,i],drop=FALSE] - 
             Y[i, tplan$to[,i],drop = FALSE] )^ground_p * tplan$mass )^(1.0/ground_p))
  } else if (method == "univariate.approximation.pwr") {
    dots <- list(...)
    if(is.null(dots$is.X.sorted)) dots$is.X.sorted <- FALSE
    is.A.sorted <- as.logical(dots$is.X.sorted)
    tplan <- transport_(A_ = X, B_ = Y, p = p, ground_p = ground_p, 
                        method_ = method, a_sort = is.A.sorted)
    cost <- sum((X[tplan$from] - 
                   Y[tplan$to] )^p * tplan$mass*1/nrow(Y))
  } else if (method == "exact" | method == "sinkhorn" | method == "greenkhorn" | method == "randkhorn" | method == "gandkhorn" | method == "sinkhorn2") {
    # tplan <- transport_(X, Y, p, ground_p, "shortsimplex")
    n1 <- ncol(X)
    n2 <- ncol(Y)
    mass_x <- as.double(rep(1/n1, n1))
    mass_y <- as.double(rep(1/n2, n2))
    
    cost <- cost_calc(X, Y, ground_p)
    tplan <- transport_plan_given_C(mass_x, mass_y, p, cost, method, ...)
  } else if (method == "univariate" | method == "hilbert" | method == "rank") {
    dots <- list(...)
    if(is.null(dots$is.X.sorted)) dots$is.X.sorted <- FALSE
    is.A.sorted <- as.logical(dots$is.X.sorted)
    if(ncol(X) == ncol(Y) ) {
      tplan <- transport_(A_ = X, B_ = Y, p = p, ground_p = ground_p, 
                          method_ = method, a_sort = is.A.sorted, epsilon = 0.0, niter = 0)
    } else if(method == "hilbert" | method == "univariate") {
      tplan <- general_1d_transport(X, Y, method = method)
    } else {
      stop("only measures with same number of atoms supported for rank methods.")
    }
    cost <- c((((colSums(abs(X[, tplan$from, drop=FALSE] - Y[, tplan$to, drop=FALSE])^ground_p))^(1/ground_p))^p %*% tplan$mass)^(1/p))
  } else if (method == "sliced") {
    dots <- list(...)
    tplan <- NULL
    nboot <- as.double(dots$nsim)
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
      cost <- ((sum(abs(x[tplan$from] - y[tplan$to])^ground_p))^(1/ground_p))^p %*% tplan$mass
      return(cost)
      }
      )
    cost <- mean(costs)^(1/p)
    tplan <- NULL
  } else if (method == "swapping") {
    dots <- list(...)
    epsilon <- as.double(dots$epsilon)
    niter <- as.integer(dots$niter)
    
    if(length(epsilon) == 0) epsilon <- as.double(0.001)
    if(length(niter) == 0) niter <- as.integer(100)
    
    if(ncol(X) > 1e6 | ncol(Y) > 1e6) {
      if(ncol(X) == ncol(Y)) {
        idx_x <- hilbert_proj_(X) 
        idx_y <- hilbert_proj_(Y)
        idxs  <- cbind(idx_x[order(idx_y)], 0:(ncol(Y)-1)) # to make similar
        idxs  <- cbind(idx_x, idx_y)
        mass <- rep(1/ncol(X), ncol(X))
      } else {
        tplan <- general_hilbert_transport(X, Y)
        idxs <- cbind(tplan$from, tplan$to)-1
        mass <- tplan$mass
      }
      if(mode(idxs) != "integer") mode(idxs) <- "integer"
      
      tplan <- transport_swap_(X,
                               Y,
                               idxs,
                               as.double(mass),
                               p, ground_p,
                               epsilon, niter)
      
    } else {
      is.A.sorted <- FALSE
      
      tplan <- transport_(A_ = X, B_ = Y, p = p, ground_p = ground_p, 
                          method_ = method, a_sort = is.A.sorted, epsilon = epsilon, niter = niter)
    }
    cost <- c((((colSums(abs(X[, tplan$from, drop=FALSE] - Y[, tplan$to, drop=FALSE])^ground_p))^(1/ground_p))^p %*% tplan$mass)^(1/p))
    
  } else {
    stop( paste0( "Transport method ", method, " not supported" ) )
  }
  
  return(list(tplan = tplan, cost = cost ))
  
}

general_1d_transport <- function(X, Y, method = c("hilbert", "univariate")) {
  method <- match.arg(method)
  
  if(method == "hilbert") {
    idx_x <- hilbert_proj_(X) + 1L
    idx_y <- hilbert_proj_(Y) + 1L
    
  } else if (method == "univariate") {
    idx_x <- order(X) 
    idx_y <- order(Y) 
  }
  n <- ncol(X)
  m <- ncol(Y)
  
  mass_a <- rep(1/n, n)
  mass_b <- rep(1/m, m)
  cum_a  <- c(cumsum(mass_a))[-n]
  cum_b  <- c(cumsum(mass_b))[-m]
  mass   <- diff(c(0,unique(sort(c(cum_a, cum_b))),1))
  # mass   <- unique(mass)
  cum_m  <- cumsum(mass)
  arep   <- table(cut(cum_m, c(-Inf, cum_a, Inf)))
  brep   <- table(cut(cum_m, c(-Inf, cum_b, Inf)))
  a_idx  <- rep(idx_x, times = arep)
  b_idx  <- rep(idx_y, times = brep)
  
  transport <- list(from = a_idx[order(b_idx)], 
                    to = sort(b_idx), 
                    mass = mass[order(b_idx)])
  
  return(transport)
  # test <- data.frame(from = a_idx, to = b_idx, mass = mass)
}

general_hilbert_transport <- function(X, Y) {
  idx_x <-  hilbert_proj_(X) + 1L
  idx_y <-  hilbert_proj_(Y) + 1L
  n <- ncol(X)
  m <- ncol(Y)
  
  mass_a <- rep(1/n, n)
  mass_b <- rep(1/m, m)
  cum_a  <- c(cumsum(mass_a))[-n]
  cum_b  <- c(cumsum(mass_b))[-m]
  mass   <- diff(c(0,unique(sort(c(cum_a, cum_b))),1))
  # mass   <- unique(mass)
  cum_m  <- cumsum(mass)
  arep   <- table(cut(cum_m, c(-Inf, cum_a, Inf)))
  brep   <- table(cut(cum_m, c(-Inf, cum_b, Inf)))
  a_idx  <- rep(idx_x, times = arep)
  b_idx  <- rep(idx_y, times = brep)
  
  transport <- list(from = a_idx[order(b_idx)], to = sort(b_idx), mass = mass[order(b_idx)])
  
  return(transport)
  # test <- data.frame(from = a_idx, to = b_idx, mass = mass)
}

transport_plan_multimarg <- function(..., p = 2, ground_p = 2,
                               observation.orientation = c("rowwise", "colwise"), 
                               method = c("hilbert", "univariate", "sliced"),
                               nsim = 1000) {
  obs <- match.arg(observation.orientation)
  method <- match.arg(method)

  if(...length() > 1) {
    data <- list(...)
  } else {
    data <- (...)
  }
  
  data <- lapply(data, function(mm) {
    if(!is.matrix(mm)) {
      mm <- as.matrix(mm)
      if(dim(mm)[2] == 1) mm <- t(mm)
    }
    return(mm) 
    })
  if(obs == "rowwise"){
    data <- lapply(data, t)
  }
  lapply(data, function(X) stopifnot(all(is.finite(X))))
  p <- as.double(p)
  ground_p <- as.double(ground_p)
  ds <- sapply(data, nrow)
  if(all(ds != ds[1])) stop("Dimension of input data is not all the same. Data can have different numbers of observations but must have the same number of covariates.")
  d <- ds[1]
  
  cost <- tplan <- NULL
  
  if(method == "hilbert") {
    idx <- lapply(data, approxOT::hilbert_proj_)
    idx <- lapply(idx, "+", 1L)
  } else if (method == "univariate") {
    idx <- lapply(data, order)
  } else if (method == "sliced") {
    if(is.null(nsim)) nsim <- 1e3
    nboot <- nsim
    theta <- matrix(rnorm(d * nboot), d, nboot)
    theta <- sweep(theta, 2, STAT=apply(theta,2,function(x) sqrt(sum(x^2))), FUN = "/")
    data_theta <- lapply(data, crossprod, y = theta)
    cost <- (mean(sapply(1:nboot, function(i)  transport_plan_multimarg(lapply(data_theta, function(j) t(j[,i])), 
                                  p = p, ground_p = ground_p,
                                  observation.orientation = "colwise", 
                                  method = "univariate")$cost^p)))^(1.0/p)
    return(list(tplan = NULL, cost = cost))
  }
  
  n     <- sapply(data, ncol)
  cmass <- lapply(n, function(nn) seq(1/nn,(nn-1)/nn, by = 1/nn))
  mass  <- diff(c(0,unique(sort(unlist(cmass))),1))
  cum_m <- cumsum(mass)
  reps  <- lapply(cmass, function(m) table(cut(cum_m, c(-Inf, m, Inf))))
  repidx<- mapply(function(i,r){rep(i, times = r)}, i = idx, r = reps, SIMPLIFY = FALSE)
  names(repidx) <- names(data)
  tplan <- list(indexes = repidx, mass = mass)
  
  cost  <- multi_marg_final_cost_(idx_ = tplan$indexes, 
                                  data_ = data, 
                                  mass_ = tplan$mass,
                                  M = length( tplan$indexes[[1]]),
                                  D = d,
                                  p = p,
                                  ground_p = ground_p
                                  )
  
  return(list(tplan = tplan, cost = cost))
}

