#ifndef TRANSPORT_H
#define TRANSPORT_H

#include "../approxOT_types.h"
#include "cost.h"
#include <string>
// #include "utils.h"
#include "trans_hilbert.h"
#include "networkflow.h"
#include "trans_shortsimplex.h"
#include "trans_rank.h"
#include "trans_univariate.h"
#include "trans_univariate_approx_pwr.h"
#include "systematic_sample.h"
#include "trans_approxOT.h"
#include "trans_swap.h"

static inline void which(const matrixI & basis, int N, int M, matrixI & index) { //check which function is working
  if ( (N*M) != index.rows() ) Rcpp::stop("Index matrix rows don't match number of possible assignments");
  int count = 0;
  // Rcpp::Rcout << index(0,0) << std::endl;
  // Rcpp::Rcout << index(0,1) << std::endl;
  if (basis.rows() != (N) ) Rcpp::stop("Basis matrix rows don't match cost matrix rows");
  if (basis.cols() != (M) ) Rcpp::stop("Basis matrix columns don't match cost matrix cols");
  
  for (int j = 0; j < M; j++) {
    for (int i = 0; i < N; i ++) {
      if(basis(i,j) == 1 ) { 
        index(count, 0) = i;
        index(count, 1) = j;
        count++;
      }
    }
  }
  if(count == 0) Rcpp::stop("No matchings found!");
  index.conservativeResize( count,	Eigen::NoChange );
}

static inline void which_nonzero(const matrix & basis, int N, int M, matrixI & index) { //check which function is working
  if ( (N*M) != index.rows() ) Rcpp::stop("Index matrix rows don't match number of possible assignments");
  int count = 0;
  // Rcpp::Rcout << index(0,0) << std::endl;
  // Rcpp::Rcout << index(0,1) << std::endl;
  if (basis.rows() != (N) ) Rcpp::stop("Assignment matrix rows don't match cost matrix rows");
  if (basis.cols() != (M) ) Rcpp::stop("Assignment matrix columns don't match cost matrix cols");
  
  for (int j = 0; j < M; j++) {
    for (int i = 0; i < N; i ++) {
      if(basis(i,j) != 0.0 ) { 
        index(count, 0) = i;
        index(count, 1) = j;
        count++;
      }
    }
  }
  if(count == 0) Rcpp::stop("No matchings found!");
  index.conservativeResize( count,	Eigen::NoChange );
}

static inline void adjustVal(vectorI & x, int N) {
  int sum_x = x.sum();
  
  if (sum_x < N) {
    int n_samp = N - sum_x;
    vectorI select(n_samp);
    Eigen::VectorXd weights(n_samp);
    weights.fill( 1.0/double(n_samp) );
    
    sample_systematic(select, weights, n_samp);
    
    for( int i = 0; i < select.size(); i ++) {
      x(select(i)) += 1;
    }
  } else {
    while (sum_x > N) {
      int n_samp = N - sum_x;
      vectorI select(n_samp);
      Eigen::VectorXd weights(n_samp);
      weights.fill( 1.0/double(n_samp) );
      
      sample_systematic(select, weights, n_samp);
      for( int i = 0; i < select.size(); i ++) {
        x(select(i)) -= 1;
      }
      sum_x = x.sum();
    }
  }
}

static inline void intNormalizeMass(const refVecConst & a, const refVecConst & b, 
                      vectorI & mass_a, vectorI & mass_b, double & max_val) {
  //based on "fudge" function from transport pacakge
  double sum_a = a.sum();
  double sum_b = b.sum();
  // double max_val = 1e09;
  
  for(int n = 0; n < a.size(); n++) mass_a(n) = int( a(n)/sum_a * max_val );
  for(int m = 0; m < b.size(); m++) mass_b(m) = int( b(m)/sum_b * max_val );
  
  adjustVal(mass_a, int(max_val));
  adjustVal(mass_b, int(max_val));
  
}

//' Calculate transportation plans given a cost matrix
//'
//' @param mass_a A reference to an Eigen::VectorXd
//' giving the empirical mass in sample A
//' @param mass_b A reference to an Eigen::VectorXi 
//' giving the empirical mass in sample B
//' @param cost_matrix A reference to an Eigen::MatrixXd giving 
//' the cost between samples A and B
//' @param idx A two column Eigen::MatrixXi giving the paired
//' indexes between samples.
//' @param mass An Eigen::VectorXd giving the mass between pairs
//' of observations.
//' @param method One of "sinkhorn", "greenkhorn", "randkhorn", 
//' grandkhorn", "hilbert", "shortsimplex", "exact",
//' "networkflow", or "univariate"
//' @param cost_matrix_A The cost matrix between observations in sample 1.
//' For use method is "sinkhorn" and unbiased is true.
//' @param cost_matrix_B The cost matrix between observations in sample 2.
//' For use method is "sinkhorn" and unbiased is true.
//' @param epsilon The value to multiple the median of the
//' cost_matrix by to give a value of lambda for the entropic
//' penalty
//' @param niter The iterations to use for the methods
//' @param unbiased Should the Sinkhorn distance be de-biased? (bool)
//' Note if its unbiased, the transport plan won't be rounded to
//' the feasible set
//' @param threads If method is "exact" or "networkflow", how many
//' threads to use in the optimization problem?
//' @return void
//' @keywords internal

// void (*shortsimplex_trans)(int*, int*, double*, int*, int*,
//       int*, int*, double*, int*,int*);
// 
// extern "C" void R_shortsimplex(DllInfo *dll) {
// shortsimplex_trans = (void (*) (int*, int*, double*, int*, int*,
//                       int*, int*, double*, int*,int*)) R_GetCCallable("transport", "shortsimplex");
// }

static inline void transport_C(const refVecConst & mass_a, const refVecConst & mass_b, 
                 refMat cost_matrix, matrixI & idx, Eigen::VectorXd & mass, const std::string & method,
                 refMat cost_matrix_A, refMat cost_matrix_B,
                 double epsilon, int niter, bool unbiased, int threads) {
  
  int N = mass_a.size();
  int M = mass_b.size();
  matrix assignment = matrix::Zero(N,M);
  
  // Rcpp::Rcout << mass.size() << std::endl;
  // if (mass.size() != (N*M)) {
  //   // Rcpp::Rcout << "resize mass\n";
  //   // Rcpp::Rcout << N*M<<std::endl;
  //   mass.resize(N*M);
  // } //Rcpp::stop("mass Eigen::VectorXd size not equal to nobs(A) * nobs(B)");
  if (idx.rows() != (N*M)) {
    idx.resize( (N*M), 2);
  } //Rcpp::stop("Index nrows not equal to size of nobs(A) * nobs(B)");
  // Rcpp::Rcout << mass.size() << std::endl;
  
  
  if (method == "shortsimplex" ) {
    double max_val = double(N);
    matrixI assign_mat = matrixI::Zero(N,M);
    matrixI basis_mat = matrixI::Zero(N,M);
    vectorI int_mass_a = vectorI::Ones(N);
    vectorI int_mass_b = vectorI::Ones(M);
    if (N != M || (mass_a.array() != 1.0/double(N)).all() ||  (mass_a.array() != 1.0/double(M)).all() ) {
      max_val = double(1e09);
      intNormalizeMass(mass_a, mass_b, int_mass_a, int_mass_b, max_val);
      // if(method == "univariate") {
      //   method = "shortsimplex";
      // }
    }
    
    trans_shortsimplex(int_mass_a, int_mass_b, cost_matrix, assign_mat, basis_mat);
    if (N != M ) {
      which(basis_mat, N, M, idx); // usage: which(const matrixI & basis, int N, int M, matrixI index )
    } else {
      which(assign_mat, N, M, idx);
    }
    double renorm = double(mass_a.sum())/max_val;
    assignment = assign_mat.cast<double>() * renorm;
  } else if (method == "networkflow" || method == "exact") {
    bool accuracy = false;
    trans_networkflow(mass_a, mass_b, cost_matrix, assignment, threads, accuracy, niter);
    which_nonzero(assignment, N, M, idx);
  } else if (method == "sinkhorn" || method == "sinkhorn_log" ||
    method == "greenkhorn" 
               // || method == "randkhorn" || method == "gandkhorn"
  ) {
    
    // void trans_approxOT(const refVecConst & mass_a, const refVecConst & mass_b, 
    //                     refMat cost_matrix, 
    //                     matrix & assign_mat,
    //                     double epsilon, int niterations,
    //                     const std::string & method);
    trans_approxOT(mass_a, mass_b, cost_matrix, assignment, 
                   epsilon, niter, unbiased, method,
                   cost_matrix_A, cost_matrix_B);
    which_nonzero(assignment, N, M, idx);
    // } else if (method == "randkhorn") {
    //   Rcpp::stop("transport method not found!");
    // } else if (method == "gandkhorn") {
    //   Rcpp::stop("transport method not found!");
  } else if (method == "hilbert") {
    Rcpp::stop("Hilbert method shouldn't rely on this function");
  } else if (method == "univariate") {
    Rcpp::stop("Univariate method shouldn't rely ont his function");
  } else {
    Rcpp::stop("transport method not found!");
  }
  
  if(mass.size() != idx.rows()){
    mass.resize(idx.rows());
  }
  
  for(int i = 0; i < idx.rows(); i ++) {
    mass(i) = assignment( idx(i,0), idx(i,1) );
  }
  
}

//' Calculate transportation plans where a cost matrix
//' needs to be calculated
//'
//' @param A An Eigen::MatrixXd of the data in sample A
//' @param B An Eigen::MatrixXd of the data in sample B
//' @param p A double greater than or equal to 1 giving the power
//' to raise the cost matrix by
//' @param ground_p A double greater than or equal to 1 
//' giving the power of the L_p norm
//' @param idx A two column Eigen::MatrixXi giving the paired
//' indexes between samples.
//' @param mass An Eigen::VectorXd giving the mass between pairs
//' of observations.
//' @param method One of "sinkhorn", "greenkhorn", "randkhorn", 
//' grandkhorn", "hilbert", "shortsimplex", "exact",
//' "networkflow", or "univariate"
//' @param a_sort Is the data in A already sorted? (bool)
//' @param epsilon The value to multiple the median of the
//' cost_matrix by to give a value of lambda for the entropic
//' penalty
//' @param niter The iterations to use for the methods
//' @param unbiased Should the Sinkhorn distance be de-biased? (bool)
//' Note if its unbiased, the transport plan won't be rounded to
//' the feasible set
//' @param threads If method is "exact" or "networkflow", how many
//' threads to use in the optimization problem?
//' @return void
//' @keywords internal
 static inline void transport(const matrix & A, const matrix & B, const double p, const double ground_p,
                matrixI & idx, Eigen::VectorXd & mass, const std::string & method, bool & a_sort,
                double epsilon, int niter, bool unbiased, int threads ) {
   
   int N = A.cols();
   int M = B.cols();
   bool univ = false;
   
   if ( method == "univariate" || ((A.rows() == 1) && (B.rows() == 1) && (N == M) ) ){
     univ = true;
   }
   if ( univ ) {
     if(A.rows() != 1 || B.rows() != 1) Rcpp::stop("rows of A and B must be 1 for univariate");
     
     vecMapConst avec(A.data(),N);
     vecMapConst bvec(B.data(),M);
     trans_univariate(avec, bvec, N, M, idx, mass, a_sort);
     
   } else if ( method == "hilbert" ) {
     trans_hilbert(A, B, N, M, idx, mass, a_sort);
   } else if ( method == "rank") {
     trans_rank(A, B, N, M, idx, mass, a_sort);
   } else if ( method == "univariate.approximation.pwr") {
     // void  trans_univariate_approx_pwr(const matrix & A, const matrix & B, int N, int M,
     //                                   matrixI & idx, Eigen::VectorXd & mass, bool & a_sort)
     trans_univariate_approx_pwr(A, B, N, M, idx, mass, a_sort);
   } else if (method == "swapping") {
     trans_hilbert(A, B, N, M, idx, mass, a_sort);
     trans_swap(A, B, N, M,
                idx, mass, ground_p,
                p, epsilon, niter);
   } else {
     
     matrix cost_matrix(N,M);
     
     Eigen::VectorXd mass_a(N);
     Eigen::VectorXd mass_b(M);
     
     mass_a.fill(1.0/double(N));
     mass_b.fill(1.0/double(M));
     
     //cost matrix calculation for non-univariate measures
     if (ground_p == 2.0) {
       cost_calculation_L2(A, B, cost_matrix);
     } else if (ground_p == 1.0){
       cost_calculation_L1(A, B, cost_matrix);
     } else {
       cost_calculation_Lp(A, B, cost_matrix, ground_p);
     }
     
     matrix cost_matrix_A;
     matrix cost_matrix_B;
     
     if (unbiased) {
       cost_matrix_A = matrix::Zero(N,N);
       cost_matrix_B = matrix::Zero(M,M);
       if (ground_p == 2.0) {
         cost_calculation_L2(B, B, cost_matrix_A);
         cost_calculation_L2(B, B, cost_matrix_A);
       } else if (ground_p == 1.0){
         cost_calculation_L1(B, B, cost_matrix_A);
         cost_calculation_L1(B, B, cost_matrix_A);
       } else {
         cost_calculation_Lp(B, B, cost_matrix_A, ground_p);
         cost_calculation_Lp(B, B, cost_matrix_A, ground_p);
       }
     } else {
       cost_matrix_A = matrix::Zero(0,0);
       cost_matrix_B = matrix::Zero(0,0);
     }
     
     cost_matrix.array() = cost_matrix.array().pow(p).eval();
     
     transport_C(mass_a, mass_b, 
                 cost_matrix, idx, mass, 
                 method, cost_matrix_A, cost_matrix_B,
                 epsilon, niter, 
                 unbiased, threads
     );
   }
   
 }

#endif //TRANSPORT_H