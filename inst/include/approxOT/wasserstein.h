#ifndef WASSERSTEIN_H
#define WASSERSTEIN_H

#include "../approxOT_types.h"
#include "sort.h"

static inline double wasserstein_p(const refVec & mass,
                     const refMat & cost, const double p,
                     const refVecI & from,  const refVecI & to) {
  
  int N = from.size();
  
  double loss = 0.0;
  
  for (int n = 0; n < N; n ++){
    int a_idx = from(n);
    int b_idx = to(n);
    // double mass = mass_a(n);
    double cur_cost = cost(a_idx, b_idx);
    
    loss += std::pow(cur_cost, p) * mass(n);
  }
  
  return(std::pow(loss, 1.0/p));
}

static inline double wasserstein_2(const refVec & mass,
                     const refMat & cost,
                     const refVecI & from,  const refVecI & to) {
  
  int N = from.size();
  
  double loss = 0.0;
  
  for (int n = 0; n < N; n ++){
    // int a_idx = from(n);
    // int b_idx = to(n);
    // double mass = mass_a(n);
    double cur_cost = cost(from(n), to(n));
    loss += cur_cost * cur_cost * mass(n);
  }
  
  return( std::sqrt(loss) );
}

static inline double wasserstein_1(const refVec & mass,
                     const refMat & cost,
                     const refVecI & from,  const refVecI & to) {
  
  int N = from.size();
  
  double loss = 0.0;
  
  for (int n = 0; n < N; n ++){
    int a_idx = from(n);
    int b_idx = to(n);
    // double mass = mass_a(n);
    // double cur_cost = cost(a_idx, b_idx);
    loss += cost(a_idx, b_idx) * mass(n);
  }
  
  return(loss);
}

//' Calculates exact or approximate optimal transport distances
//'
//' @param mass_a A reference to an Eigen::VectorXd 
//' with the empirical weights from sample 1
//' @param mass_b A reference to an Eigen::VectorXd 
//' with the empirical weights from sample 2
//' @param cost The cost matrix
//' @param p The power to raise the cost by
//' @param from The indexes of the first sample
//' @param to he indexes of the second sample
//' @return returns a double value denoting the distance
//' @keywords internal
 static inline double wasserstein(const refVec & mass,
                   const refMat & cost, const double p, 
                   const refVecI & from,  const refVecI & to) {
  double loss;
  if( p == 2.0) {
    loss = wasserstein_2(mass, cost, from, to);
  } else if (p == 1.0) {
    loss = wasserstein_1(mass, cost, from, to);
  } else {
    loss = wasserstein_p(mass, cost, p, from, to);
  }
  
  return(loss);
}

static inline double wasserstein_p_iid_p(refMat X, refMat Y, double p){
  
  if(X.cols() != Y.cols()){
    Rcpp::stop("Number of columns of first matrix don't match number of columns of second matrix");
  }
  
  if(X.rows() != Y.rows()){
    Rcpp::stop("Number of rows of first matrix don't match number of rows of second matrix");
  }
  
  sort_matrix_by_row(X);
  sort_matrix_by_row(Y);
  
  return (X-Y).cwiseAbs().array().pow(p).mean();
}

static inline double wasserstein_p_iid(refMat X, refMat Y, double p){
  
  if(X.cols() != Y.cols()){
    Rcpp::stop("Number of columns of first matrix don't match number of columns of second matrix");
  }
  
  if(X.rows() != Y.rows()){
    Rcpp::stop("Number of rows of first matrix don't match number of rows of second matrix");
  }
  
  sort_matrix_by_row(X);
  sort_matrix_by_row(Y);
  
  Eigen::VectorXd loss = (X-Y).cwiseAbs().array().pow(p).colwise().mean();
  loss = loss.array().pow(1.0/p);
  return loss.mean();
}

static inline double wasserstein_2_iid(refMat X, refMat Y){
  
  if(X.cols() != Y.cols()){
    Rcpp::stop("Number of columns of first matrix don't match number of columns of second matrix");
  }
  
  if(X.rows() != Y.rows()){
    Rcpp::stop("Number of rows of first matrix don't match number of rows of second matrix");
  }
  
  sort_matrix_by_row(X);
  sort_matrix_by_row(Y);
  
  Eigen::VectorXd loss = (X-Y).cwiseAbs2().array().colwise().mean().sqrt();
  return loss.mean();
}

//' Calculates exact or approximate optimal transport distances
//' averaged by observation
//'
//' @param X a matrix of observations from sample 1, already sorted
//' @param Y a matrix of observations from sample 2, already sorted
//' @return returns a double value denoting the distance
//' @details Just calculates the univariate by observation
//' squared wasserstein distance
//' @keywords internal
 static inline double wasserstein_2_iid_2(refMat X, refMat Y){
   
   if(X.cols() != Y.cols()) {
     Rcpp::stop("Number of columns of first matrix don't match number of columns of second matrix");
   }
   
   if(X.rows() != Y.rows()) {
     Rcpp::stop("Number of rows of first matrix don't match number of rows of second matrix");
   }
   
   sort_matrix_by_row(X);
   sort_matrix_by_row(Y);
   
   return (X-Y).cwiseAbs2().mean();
 }

static inline double wasserstein_1_iid(refMat X, refMat Y){
  
  if(X.cols() != Y.cols()){
    Rcpp::stop("Number of columns of first matrix don't match number of columns of second matrix");
  }
  
  if(X.rows() != Y.rows()){
    Rcpp::stop("Number of rows of first matrix don't match number of rows of second matrix");
  }
  
  sort_matrix_by_row(X);
  sort_matrix_by_row(Y);
  
  return (X-Y).cwiseAbs().array().mean();
}

#endif //WASSERSTEIN_H