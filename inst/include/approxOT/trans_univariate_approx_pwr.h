#ifndef TRANS_UNIVARIATE_APPROXIMATION_PWR_H
#define TRANS_UNIVARIATE_APPROXIMATION_PWR_H

#include "../approxOT_types.h"
#include "sort.h"

//' Transport plan based on a univariate approximation by row
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
 static inline void  trans_univariate_approx_pwr(const matrix & A, const matrix & B, int N, int M,
                                   matrixI & idx, Eigen::VectorXd & mass, bool & a_sort) {
   if(N != M) {
     Rcpp::stop("Number of atoms of A and B must match for univariate approximation method!");
   }
   int S = A.rows();
   
   idx.resize(N*S,2);
   mass.resize(N*S);
   mass.fill(1.0/double(N));
   matrixI idx_A(S,N);
   matrixI idx_B(S,N);
   
   if (a_sort) { 
     idx_A = vectorI::LinSpaced(N*S,0,(N*S)-1);
   } else {
     sort_indexes_byrow_totalentry(A, idx_A);
     a_sort = true;
   }
   sort_indexes_byrow_totalentry(B, idx_B);
   
   idx.col(1) = vectorI::LinSpaced(N*S,0,(N*S)-1);
   vecMapI idx_avec(idx_A.data(), idx_A.size());
   vecMapI idx_bvec(idx_B.data(), idx_B.size());
   for ( int i = 0; i < N*S; i ++ ) {
     idx(idx_bvec(i), 0) = idx_avec(i);
   }
 }


#endif //TRANS_UNIVARIATE_APPROXIMATION_PWR_H