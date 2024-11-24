#ifndef TRANS_HILBERT_H
#define TRANS_HILBERT_H

#include "../approxOT_types.h"
#include "hilbert_cgal.h"

//' Transport plan based on Hilbert sorting
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
static inline void trans_hilbert(const matrix & A, const matrix & B, int N, int M,
                    matrixI & idx, Eigen::VectorXd &  mass, bool & a_sort)
 {
   int K = A.rows();
   if(N != M) {
     Rcpp::stop("Number of atoms of A and B must match for current implementation of Hilbert sort!");
   }
   idx.resize(N, 2);
   mass.resize(N);
   mass.fill(1.0/double(N));
   mass.fill(1.0/double(N));
   std::vector<int> idx_A(N);
   std::vector<int> idx_B(N);
   
   if ( a_sort ) {
     std::iota (std::begin(idx_A), std::end(idx_A), 0);
   } else {
     hilbert_sort_cgal_fun(A.data(), K, N, &idx_A[0] );
     
     // a_sort = true;
   }
   hilbert_sort_cgal_fun( B.data(), K, N, &idx_B[0] );
   idx.col(1) = vectorI::LinSpaced(N,0,N-1);
   for ( int n = 0; n < N; n++ ) {
     idx(idx_B[n],0) = idx_A[n];
   }
 }

#endif //TRANS_HILBERT_H