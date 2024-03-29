% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transport_plan.R
\name{transport_plan_multimarg}
\alias{transport_plan_multimarg}
\title{Multimarginal optimal transport plans}
\usage{
transport_plan_multimarg(
  ...,
  p = 2,
  ground_p = 2,
  observation.orientation = c("rowwise", "colwise"),
  method = c("hilbert", "univariate", "sliced"),
  nsim = 1000
)
}
\arguments{
\item{...}{Either data matrices as separate arguments or a list of data matrices. Arguments after the data must be specified by name.}

\item{p}{The power of the Wasserstein distance to use}

\item{ground_p}{The power of the Euclidean distance to use}

\item{observation.orientation}{Are observations by rows or columns}

\item{method}{One of "hilbert", "univariate", or "sliced"}

\item{nsim}{Number of simulations to use for the sliced method}
}
\value{
transport plan
}
\description{
Multimarginal optimal transport plans
}
\examples{
set.seed(23423)
n <- 100
d <- 10
p <- ground_p <- 2 #euclidean cost, p = 2
x <- matrix(stats::rnorm((n + 11)*d), n + 11 , d)
y <- matrix(stats::rnorm(n*d), n, d)
z <- matrix(stats::rnorm((n +455)*d), n +455, d)

# make data a list
data <- list(x,y,z)

tplan <- transport_plan_multimarg(data, p = p, ground_p = ground_p,
observation.orientation = "rowwise", method = "hilbert")

#' #transpose data works too
datat <- lapply(data, t)

tplan2 <- transport_plan_multimarg(datat, p = p, ground_p = ground_p,
observation.orientation = "colwise",method = "hilbert")
}
