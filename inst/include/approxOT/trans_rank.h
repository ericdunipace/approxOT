#ifndef TRANS_RANK_H
#define TRANS_RANK_H

#include "../approxOT_types.h"
#include "sort.h"

//' Transport plan based on average rank
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
 static inline void trans_rank(const matrix & A, const matrix & B, int N, int M,
                  matrixI & idx, Eigen::VectorXd & mass, bool & a_sort) {
   if (N != M) {
     Rcpp::stop("Number of atoms of A and B must match for ranks method!");
   }
   idx.resize(N, 2);
   mass.resize(N);
   mass.fill(1.0/double(N));
   matrixI ranksB(B.rows(), M);
   rank_mat(B, ranksB);
   Eigen::VectorXd meanB = ranksB.colwise().sum().cast<double>();
   // meanB /= double(M); // add in normalization if indices get too large
   std::vector<size_t> idx_B = sort_indexes(meanB);
   std::vector<size_t> idx_A(N);
   
   if (!a_sort) {
     // Rcpp::Rcout << "a not sorted\n";
     matrixI ranksA(A.rows(), N);
     rank_mat(A, ranksA);
     Eigen::VectorXd meanA = ranksA.colwise().sum().cast<double>();
     // meanA /= double(N); // add in normalization if indices get too large
     sort_indexes(meanA, idx_A);
   } else {
     // Rcpp::Rcout << "a sorted\n";
     std::iota (idx_A.begin(), idx_A.end(), 0);
   }
   idx.col(1) = vectorI::LinSpaced(N,0,N-1);
   for ( int n = 0; n < N; n++ ) {
     idx(idx_B[n],0) = idx_A[n];
     // Rcpp::Rcout << idx_A[n] << ", ";
   }
 }

#endif //TRANS_RANK_H