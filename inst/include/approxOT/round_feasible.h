#ifndef TRANS_ROUND_FEASIBLE_H
#define TRANS_ROUND_FEASIBLE_H

#include "../approxOT_types.h"

//' Rounds approximate transport plans to the feasible set
//'
//' @param F the transport plan, a matrix of class Eigen::MatrixXd 
//' @param mass_a The sample weights of the first sample of data
//' @param mass_b The sample weights of the second sample of data
//' @return void
//' @keywords internal
static inline void round_feasible(matrix & F, const refVecConst & mass_a, const refVecConst & mass_b) {
   Eigen::VectorXd a_f = F.rowwise().sum();
   Eigen::VectorXd b_f = F.colwise().sum();
   
   // Rcpp::Rcout << a_f << "\n";
   // Rcpp::Rcout << b_f << "\n";
   // 
   Eigen::VectorXd x = mass_a.cwiseQuotient(a_f).cwiseMin(1.0);
   Eigen::VectorXd y = mass_b.cwiseQuotient(b_f).cwiseMin(1.0);
   
   matrix F_prime = x.asDiagonal() * F * y.asDiagonal();
   
   Eigen::VectorXd err_a = mass_a;
   err_a -= F_prime.rowwise().sum();
   Eigen::VectorXd err_b = mass_b;
   err_b -= F_prime.colwise().sum();
   
   // Rcpp::Rcout << F_prime.colwise().sum() <<"\n";
   // Rcpp::Rcout << (mass_b - F_prime.colwise().sum().transpose()) <<"\n";
   // // Rcpp::Rcout << err_a.lpNorm<1>() <<"\n";
   // if((F_prime.array() < 0).any()) Rcpp::stop("F' < 0");
   // if((err_a.array() < 0).any()) Rcpp::stop("erra' < 0");
   // if((err_b.array() < 0).any()) Rcpp::stop("errb' < 0");
   // 
   F.noalias() = F_prime;
   F.noalias() += err_a * err_b.transpose() / err_a.lpNorm<1>();
   // if((F.array() < 0).any()) Rcpp::stop("F < 0");
 }

#endif //TRANS_ROUND_FEASIBLE_H