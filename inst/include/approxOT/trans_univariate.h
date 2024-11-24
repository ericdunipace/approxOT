#ifndef TRANS_UNIVARIATE_H
#define TRANS_UNIVARIATE_H

#include "../approxOT_types.h"
#include "sort.h"

//' Transport plan based for univariate measures
//'
//' @param A An Eigen::MatrixXd of the data in sample A
//' @param B An Eigen::MatrixXd of the data in sample B
//' @param N The columns of A
//' @param M The columns of B
//' @param idx A two column Eigen::MatrixXi giving the paired
//' indexes between samples.
//' @param mass An Eigen::VectorXd giving the mass between pairs
//' of observations.
//' @param a_sort Is the data in A already sorted? (bool)
//' @return void
//' @keywords internal
 static inline void  trans_univariate(const Eigen::VectorXd & A, const Eigen::VectorXd & B, int N, int M,
                        matrixI & idx, Eigen::VectorXd & mass, bool & a_sort) {
   if(N != M) {
     Rcpp::stop("Number of atoms of A and B must match for univariate method!");
   }
   idx.resize(N, 2);
   mass.resize(N);
   mass.fill(1.0/double(N));
   std::vector<size_t> idx_A(N);
   std::iota (idx_A.begin(), idx_A.end(), 0);
   std::vector<size_t> idx_B = sort_indexes(B);
   // std::vector<size_t> idx_B(N);
   // std::iota (idx_B.begin(), idx_B.end(), 0);
   // sort_indexes(B, idx_B);
   if (!a_sort) {
     sort_indexes(A, idx_A);
     // a_sort = true;
   }
   idx.col(1) = vectorI::LinSpaced(N,0,N-1);
   for ( int n = 0; n < N; n++ ) {
     // Rcpp::Rcout << idx_B[n] << ", " << idx_A[n] << "\n";
     idx(idx_B[n],0) = idx_A[n];
   }
 }
#endif //TRANS_UNIVARIATE_H