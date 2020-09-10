testthat::test_that("multimarg cost, l2", {
  set.seed(23423)
  n <- 100
  d <- 10
  x <- matrix(rnorm(n*d), n , d)
  y <- matrix(rnorm(n*d), n , d)
  z <- matrix(rnorm(n*d), n , d)
  
  data <- list(x,y,z)
  datat <- lapply(data, t)
  indexes <- lapply(1:3, function(i) sample.int(n,n))
  mass <- rep(1/n,n)
  p <- ground_p <- 2
  
  cost <- multi_marg_final_cost_(indexes, datat, mass = mass,
                                   M = n, D = d, p = p, ground_p = ground_p)
  
  dists <- sqrt(weighted.mean(rowSums(abs(data[[1]][indexes[[1]],] - data[[2]][indexes[[2]],])^p) + 
                      rowSums(abs(data[[1]][indexes[[1]],] - data[[3]][indexes[[3]],])^p) + 
                      rowSums(abs(data[[2]][indexes[[2]],] - data[[3]][indexes[[3]],])^p),
                      w = mass))
  
  testthat::expect_equivalent(dists, cost)
  
})

testthat::test_that("multimarg cost, l1", {
  set.seed(23423)
  n <- 100
  d <- 10
  x <- matrix(rnorm(n*d), n , d)
  y <- matrix(rnorm(n*d), n , d)
  z <- matrix(rnorm(n*d), n , d)
  
  data <- list(x,y,z)
  datat <- lapply(data, t)
  indexes <- lapply(1:3, function(i) sample.int(n,n))
  mass <- rep(1/n,n)
  p <- ground_p <- 1
  
  cost <- multi_marg_final_cost_(indexes, datat, mass = mass,
                                 M = n, D = d, p = p, ground_p = ground_p)
  
  dists <- (weighted.mean(rowSums(abs(data[[1]][indexes[[1]],] - data[[2]][indexes[[2]],])^p) + 
                                rowSums(abs(data[[1]][indexes[[1]],] - data[[3]][indexes[[3]],])^p) + 
                                rowSums(abs(data[[2]][indexes[[2]],] - data[[3]][indexes[[3]],])^p),
                              w = mass))
  
  testthat::expect_equivalent(dists, cost)
  
})

testthat::test_that("multimarg cost, l3", {
  set.seed(23423)
  n <- 100
  d <- 10
  x <- matrix(rnorm(n*d), n , d)
  y <- matrix(rnorm(n*d), n , d)
  z <- matrix(rnorm(n*d), n , d)
  
  data <- list(x,y,z)
  datat <- lapply(data, t)
  indexes <- lapply(1:3, function(i) sample.int(n,n))
  mass <- rep(1/n,n)
  p <- ground_p <- 3
  
  cost <- multi_marg_final_cost_(indexes, datat, mass = mass,
                                 M = n, D = d, p = p, ground_p = ground_p)
  
  dists <- (weighted.mean(rowSums(abs(data[[1]][indexes[[1]],] - data[[2]][indexes[[2]],])^p) + 
                            rowSums(abs(data[[1]][indexes[[1]],] - data[[3]][indexes[[3]],])^p) + 
                            rowSums(abs(data[[2]][indexes[[2]],] - data[[3]][indexes[[3]],])^p),
                          w = mass))^(1/3)
  
  testthat::expect_equivalent(dists, cost)
  
})

testthat::test_that("multimarg transportation equal, hilbert", {
  set.seed(23423)
  n <- 100
  d <- 10
  x <- matrix(rnorm(n*d), n , d)
  y <- matrix(rnorm(n*d), n , d)
  z <- matrix(rnorm(n*d), n , d)
  
  data <- list(x,y,z)
  indexes <- lapply(data, approxOT::hilbert.projection)
  mass <- rep(1/n,n)
  p <- ground_p <- 2
  
  dists <- sqrt(weighted.mean(rowSums(abs(data[[1]][indexes[[1]],] - data[[2]][indexes[[2]],])^p) + 
                                rowSums(abs(data[[1]][indexes[[1]],] - data[[3]][indexes[[3]],])^p) + 
                                rowSums(abs(data[[2]][indexes[[2]],] - data[[3]][indexes[[3]],])^p),
                              w = mass))
  
  # debugonce(transport_plan_multimarg)
  tplan <- transport_plan_multimarg(data, p = p, ground_p = ground_p,
                           observation.orientation = "rowwise",
                           method = "hilbert")
  
  testthat::expect_equivalent(dists, tplan$cost)
  testthat::expect_equivalent(indexes, tplan$tplan$indexes)
  
})

testthat::test_that("multimarg transportation equal, univariate", {
  set.seed(01212)
  n <- 100
  d <- 1
  x <- matrix(rnorm(n*d), n , d)
  y <- matrix(rnorm(n*d), n , d)
  z <- matrix(rnorm(n*d), n , d)
  
  data <- list(x,y,z)
  indexes <- lapply(data, order)
  mass <- rep(1/n,n)
  p <- ground_p <- 2
  
  dists <- sqrt(weighted.mean((abs(data[[1]][indexes[[1]],] - data[[2]][indexes[[2]],])^p) + 
                                (abs(data[[1]][indexes[[1]],] - data[[3]][indexes[[3]],])^p) + 
                                (abs(data[[2]][indexes[[2]],] - data[[3]][indexes[[3]],])^p),
                              w = mass))
  
  # debugonce(transport_plan_multimarg)
  tplan <- transport_plan_multimarg(data, p = p, ground_p = ground_p,
                                    observation.orientation = "rowwise",
                                    method = "univariate")
  
  testthat::expect_equivalent(dists, tplan$cost)
  testthat::expect_equivalent(indexes, tplan$tplan$indexes)
  
})

