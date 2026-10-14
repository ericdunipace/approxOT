#ifndef TRANS_GREENKHORN_H
#define TRANS_GREENKHORN_H

#include "../approxOT_types.h"
#include <algorithm>

static inline Eigen::VectorXd rho_vec(const Eigen::VectorXd & a, const Eigen::VectorXd & b) {
  return ( b.array()  - a.array() +  a.array() * ( a.array().log() - b.array().log() ) );
}

static inline int argmax_rho (const Eigen::VectorXd & r) {
  // Eigen::VectorXd r = rho_vec(a,b);
  Eigen::Index maxIndex;
  
  r.maxCoeff(&maxIndex);
  
  int which_max = maxIndex;
  return(which_max);
}

static inline double dist_approx_ot(const refVecConst & mass_a, const refVecConst & mass_b,
                      const Eigen::VectorXd & r, const Eigen::VectorXd & c, int p) {
  double out = 0.0;
  if(p == 2) {
    out = (r - mass_a).norm() + (c - mass_b).norm();
  } else if ( p == 1) {
    out = (r - mass_a).lpNorm<1>() + (c - mass_b).lpNorm<1>();
  } else {
    Rcpp::stop("Other norms not supported");
  }
  return (out);
}

// Algorithm 4 of Altschuler, J., Weed, J., & Rigollet, P. (2017). Near-linear time approximation 
// algorithms for optimal transport via Sinkhorn iteration. 31st Conference on 
// Neural Information Processing Systems, (1), 1–11. Long Beach, CA.
// void trans_greenkhorn(const refVecConst & mass_a, const refVecConst & mass_b,
//                     const matrix & exp_cost,
//                     matrix & A,
//                     double eta, double epsilon, int niterations) {
//   int N = mass_a.size();
//   int M = mass_b.size();
// 
//   matrix A_0 = exp_cost.array()/exp_cost.lpNorm<1>();
//   A = A_0;
//   // matrix A_0 = exp_cost;
//   // A = exp_cost.array()/exp_cost.lpNorm<1>();
//   Eigen::VectorXd log_a = mass_a.array().log();
//   Eigen::VectorXd log_b = mass_b.array().log();
//   Eigen::VectorXd x = vector::Zero(N);
//   Eigen::VectorXd y = vector::Zero(M);
//   Eigen::VectorXd r = A.rowwise().sum();
//   Eigen::VectorXd c = A.colwise().sum();
//   Eigen::VectorXd rho_vals_r = rho_vec(mass_a, r);
//   Eigen::VectorXd rho_vals_c = rho_vec(mass_b, c);
// 
//   for( int i = 0; i < niterations; i ++){
//     int I = argmax_rho(rho_vals_r);
//     int J = argmax_rho(rho_vals_c);
// 
//     if (rho_vals_r(I) > rho_vals_c(J)) {
//       x(I) += log_a(I);
//       x(I) -= std::log(r(I));
//     } else {
//       y(J) += log_b(J);
//       y(J) -= std::log(c(J));
//     }
//     A = (x.array().exp()).matrix().asDiagonal() * A_0 * (y.array().exp()).matrix().asDiagonal();
//     // Rcpp::Rcout <<"dist: "<< dist_approx_ot(mass_a, mass_b, r, c, 2) << "\n";
//     // Rcpp::Rcout <<"A: " << A(0,0) << ", ";
//     // Rcpp::Rcout <<"A0: " << A_0(0,0) << ", ";
//     // Rcpp::Rcout <<"x: " << x(0) << ", ";
//     // Rcpp::Rcout <<"y: " << y(0) << "\n";
//     r = A.rowwise().sum();
//     c = A.colwise().sum();
//     if(dist_approx_ot(mass_a, mass_b, r, c, 2) <= epsilon) {
//       break;
//     }
//     
//     if (rho_vals_r(I) > rho_vals_c(J)) {
//       rho_vals_r(I) = 0.0;
//       rho_vals_c = rho_vec(mass_b, c);
//     } else {
//       rho_vals_r = rho_vec(mass_a, r);
//       rho_vals_c(J) = 0.0;
//     }
//   }
// }

//' Generates approximate optimal transport plans 
//' using the greenkhorn algorithm
//'
//' @param mass_a A reference to an Eigen::VectorXd
//' giving the empirical mass in sample 1
//' @param mass_b A reference to an Eigen::VectorXd
//' giving the empirical mass in sample 2
//' @param exp_cost A reference to an Eigen::MatrixXd giving 
//' the cost between samples A and B
//' @param A The assignment matrix
//' @param eta The inverse of lambda value used for 
//' the entropy penalty
//' @param epsilon The desired error bound
//' @param niterations The iterations to use for the methods
//' @return void
//' @keywords internal
// code adapted from https://github.com/JasonAltschuler/OptimalTransportNIPS17/blob/master/algorithms/greenkhorn.m
static inline void trans_greenkhorn(const refVecConst & mass_a, const refVecConst & mass_b,
                      const matrix & exp_cost,
                      matrix & A,
                      double eta, double epsilon, int niterations) {
  // int N = mass_a.size();
  // int M = mass_b.size();
  
  A = exp_cost.array()/exp_cost.lpNorm<1>(); //can probably just use sum
  // A = A_0;
  // matrix A_0 = exp_cost;
  // Eigen::VectorXd log_a = mass_a.array().log(); //not working on log-scale but maybe should?
  // Eigen::VectorXd log_b = mass_b.array().log();
  Eigen::VectorXd r = A.rowwise().sum();
  Eigen::VectorXd c = A.colwise().sum();
  // Eigen::VectorXd rho_vals_r = rho_vec(mass_a, log_a, r);
  // Eigen::VectorXd rho_vals_c = rho_vec(mass_b, log_b, c);
  
  for( int i = 0; i < niterations; i ++){
    // int I = argmax_rho(rho_vals_r);
    // int J = argmax_rho(rho_vals_c);
    Eigen::Index I, J;
    double viol_r = (r - mass_a).cwiseAbs().maxCoeff(&I);
    double viol_c = (c - mass_b).cwiseAbs().maxCoeff(&J);
    if (std::max(viol_r, viol_c) <= epsilon)  break;
    if (viol_r > viol_c) {
      double scale_x = mass_a(I)/r(I);
      
      c += (scale_x - 1.0) * A.row(I).transpose();
      A.row(I) *= scale_x;
      r(I) = mass_a(I);
      
    } else {
      double scale_y = mass_b(J)/c(J);
      
      r += (scale_y - 1.0) * A.col(J);
      A.col(J) *= scale_y;
      c(J) = mass_b(J);
    }
    
    if(i % 10 == 0) Rcpp::checkUserInterrupt();
      // if (dist_approx_ot(mass_a, mass_b, r, c, 1.0) <= epsilon) {
      //       // Rcpp::Rcout << " " << dist_approx_ot(mass_a, mass_b, r, c, 2.0) << "\n";
      //       // Rcpp::Rcout << dist_approx_ot(mass_a, mass_b, r, c, 1.0) << "\n";
      //       // Rcpp::Rcout << epsilon << "\n";
      //   break;
      // }
    
  }
  A /= A.sum();
}
#endif //TRANS_GREENKHORN_H