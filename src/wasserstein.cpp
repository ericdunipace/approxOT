#include "approxOT/wasserstein.h"

//[[Rcpp::export]]
double wasserstein_(const Rcpp::NumericVector & mass_,
             const Rcpp::NumericMatrix & cost_, const double p, 
             const Rcpp::IntegerVector & from_,  const Rcpp::IntegerVector & to_) {
  
  int N = from_.size();
  const vecMap mass(Rcpp::as<vecMap >(mass_));
  // const vecMap mass_b(Rcpp::as<vecMap >(mass_b_));
  const matMap cost(Rcpp::as<matMap >(cost_));
  vectorI from(N);
  vectorI to(N);
  // const vecMapI from(Rcpp::as<vecMapI >(from_));
  // const vecMapI to(Rcpp::as<vecMapI >(to_));
  
  for(int n = 0; n < N; n ++) {
    to(n) = to_(n) - 1;
    from(n) = from_(n) - 1;
  }
  
  return( wasserstein(mass, cost, p, from, to) );
}

//[[Rcpp::export]]
double wasserstein_p_iid_(const SEXP & X_, const SEXP & Y_, double p) {
  
  const matMap Xtemp(Rcpp::as<matMap >(X_));
  const matMap Ytemp(Rcpp::as<matMap >(Y_));
  
  matrix X = Xtemp;
  matrix Y = Ytemp;
  
  if ( p == 2.0 ) {
    // Rcout << "w2";
    return wasserstein_2_iid(X,Y);
  } else if ( p == 1.0 ) {
    // Rcout << "wp";
    return wasserstein_1_iid(X,Y);
  } else {
    return wasserstein_p_iid(X,Y,p);
  }
}

//[[Rcpp::export]]
double wasserstein_p_iid_p_(const SEXP & X_, const SEXP & Y_, double p) {
  
  const matMap Xtemp(Rcpp::as<matMap >(X_));
  const matMap Ytemp(Rcpp::as<matMap >(Y_));
  
  matrix X = Xtemp;
  matrix Y = Ytemp;
  
  if ( p == 2.0 ) {
    // Rcout << "w2";
    return wasserstein_2_iid_2(X,Y);
  } else if ( p == 1.0 ) {
    // Rcout << "wp";
    return wasserstein_1_iid(X,Y);
  } else {
    return wasserstein_p_iid_p(X,Y,p);
  }
}