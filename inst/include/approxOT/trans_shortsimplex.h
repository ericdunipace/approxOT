#ifndef TRANS_SHORTSIMPLEX_H
#define TRANS_SHORTSIMPLEX_H

#include "../approxOT_types.h"

#include "shortsimplex.h"

//' Generates exact optimal transport plans 
//' using the shortlist algorithm
//'
//' @param mass_a An Eigen::VectorXi
//' giving the empirical mass in sample 1
//' @param mass_b An Eigen::VectorXi
//' giving the empirical mass in sample 2
//' @param cost_matrix A reference to an Eigen::MatrixXd giving 
//' the cost between samples A and B
//' @param assign_mat The assignment matrix as
//' @param basis_mat An Eigen::MatrixXi giving a variable used to
//' construct the basis functions
//' @return void
//' @keywords internal
static inline void trans_shortsimplex(vectorI & mass_a, vectorI & mass_b, refMat cost_matrix, 
                         matrixI & assign_mat, matrixI &  basis_mat) {
   
   int N = mass_a.size();
   // int * NN;
   // NN = &N;
   
   int M = mass_b.size();
   // int * MM;
   // MM = &M;
   
   int * a = mass_a.data();
   int * b = mass_b.data();
   
   double * costm = cost_matrix.data();
   
   int slength = std::min(M, 
                          15 + std::max(0, 
                                        int(std::floor(15.0 * 
                                            std::log(double(M)/400.0)/std::log(2.0)))));
   if (slength > M) {
     slength = M;
     Rcpp::warning("Shortlist parameter 'slength' too large...  decreased to maximal value.");
   }
   
   int kfound = slength; 
   double psearched = 0.05;
   int * assignment = assign_mat.data();
   int * basis = basis_mat.data();
   
   shortsimplex( &slength, &kfound, &psearched,
                 &N, &M, a, b, costm, assignment, basis);
   
 }

#endif //TRANS_SHORTSIMPLEX_H