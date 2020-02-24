test_that("multiplication works", {
  set.seed(234234)
  n1 <- 100
  n2 <- 197
  d<- 10
  
  x <- matrix(rnorm(n1*d),n1,d)
  y <- matrix(rnorm(n2*d),n2,d)
  
  debugonce(general_hilbert_transport)
  general_hilbert_transport(t(x), t(y), 2, 2)
})
