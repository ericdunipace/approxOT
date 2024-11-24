#ifndef COST_H
#define COST_H

#include "../approxOT_types.h"

static inline int dist_2d_to_1d_(int i, int j, int n) {
  if((i >= 0) && (j >= 0) && (i < n) && (j < n)) {
    int ii = i;
    int jj = j;
    if(ii < jj) {
      ii = j;
      jj = i;
    }
    int k = (2 * n - jj - 1) * (jj) / 2 + (ii - jj) - 1;
    if(k < 0) {
      Rcpp::Rcout << 
        "i: " << i <<
          ", j: " << j <<
            ", ii: " << ii <<
              ", jj: " << jj <<
                ", n: " << n << 
                  ", and k: " << k << ". ";
      Rcpp::stop("Non-valid result in dist_2d_to_1d_ function");
    }
    return(k);
  } else {
    Rcpp::Rcout << 
      "i: " << i <<
        ", j: " << j <<
          ", and n: " << n << ". ";
    Rcpp::stop("Non-valid indexes in dist_2d_to_1d_ function");
  }
}


//' Fills in the values of a cost matrix
//'
//' @param A The data matrix for sample 1 with obsrevations unique by column
//' @param B The data matrix for sample 2 with obsrevations unique by column
//' @param cost_matrix The empty cost matrix of class 
//' Eigen::MatrixXd with dimensions A.col() by B.cols()
//' @return void
//' @details Calculates the \eqn{L_2} cost raised to the power of 2
//' @keywords internal
static inline void cost_calculation_L2sq(const refMatConst & A, const refMatConst & B, matrix & cost_matrix) {
   for (int j = 0; j < B.cols(); j++) { 
     Eigen::VectorXd bvec = B.col(j);
     for (int i = 0; i < A.cols(); i++) {
       cost_matrix(i,j) = (A.col(i)-bvec).squaredNorm();
     }
   }
 }


//' Fills in the values of a cost matrix
//'
//' @param A The data matrix for sample 1 with obsrevations unique by column
//' @param B The data matrix for sample 2 with obsrevations unique by column
//' @param cost_matrix The empty cost matrix of class 
//' Eigen::MatrixXd with dimensions A.col() by B.cols()
//' @return void
//' @details Calculates the \eqn{L_2} cost
//' @keywords internal
static inline void cost_calculation_L2(const refMatConst & A, const refMatConst & B, matrix & cost_matrix) {
   for (int j = 0; j < B.cols(); j++) { 
     Eigen::VectorXd bvec = B.col(j);
     for (int i = 0; i < A.cols(); i++) {
       cost_matrix(i,j) = (A.col(i)-bvec).norm();
     }
   }
 }

//' Fills in the values of a cost matrix
//'
//' @param A The data matrix for sample 1 with obsrevations unique by column
//' @param B The data matrix for sample 2 with obsrevations unique by column
//' @param cost_matrix The empty cost matrix of class 
//' Eigen::MatrixXd with dimensions A.col() by B.cols()
//' @return void
//' @details Calculates the \eqn{L_1} cost
//' @keywords internal
static inline void cost_calculation_L1(const refMatConst & A, const refMatConst & B, matrix & cost_matrix) {
   for (int j = 0; j < B.cols(); j++) { 
     Eigen::VectorXd bvec = B.col(j);
     for (int i = 0; i < A.cols(); i++) {
       cost_matrix(i,j) = (A.col(i)-bvec).cwiseAbs().sum();
     }
   }
 }

//' Fills in the values of a cost matrix
//'
//' @param A The data matrix for sample 1 with obsrevations unique by column
//' @param B The data matrix for sample 2 with obsrevations unique by column
//' @param cost_matrix The empty cost matrix of class 
//' Eigen::MatrixXd with dimensions A.col() by B.cols()
//' @return void
//' @details Calculates the \eqn{L_p} cost
//' @keywords internal
static inline void cost_calculation_Lp(const refMatConst & A, const refMatConst & B, matrix & cost_matrix, double p) {
   double p_inv = 1.0/p;
   
   for (int j = 0; j < B.cols(); j++) { 
     Eigen::VectorXd bvec = B.col(j);
     for (int i = 0; i < A.cols(); i++) {
       double cost_p = (A.col(i)-bvec).array().pow(p).sum();
       cost_matrix(i,j) = std::pow(cost_p, p_inv);
     }
   }
 }

#endif //COST_H