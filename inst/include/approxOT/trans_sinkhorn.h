#ifndef TRANS_SINKHORN_H
#define TRANS_SINKHORN_H

#include "../approxOT_types.h"
#include "utils.h"

static inline double LogSumExp(matrix Mat) {
  double max = Mat.maxCoeff();
  return std::log((Mat.array() - max).array().exp().sum()) + max;
}

static inline Eigen::VectorXd rowLogSumExp_(matrix Mat) {
  
  Eigen::VectorXd max = Mat.rowwise().maxCoeff();
  // Rcpp::Rcout << "max: " << max << "\n";
  // Rcpp::Rcout << "sweep:\n" << (Mat.colwise() - max) <<"\n";
  Eigen::VectorXd sum = (Mat.colwise() - max).array().exp().rowwise().sum().log();
  // Rcpp::Rcout << "sum: " << sum << "\n";
  // Rcpp::Rcout << "out: " << max + sum << "\n";
  
  return max + sum;
  // return (Mat.array().exp().rowwise().sum().log());
}

static inline Eigen::VectorXd colLogSumExp_(matrix Mat) {
  
  Eigen::RowVectorXd max = Mat.colwise().maxCoeff();
  // Rcpp::Rcout << "max: " << max << "\n";
  // Rcpp::Rcout << "sweep:\n" << (Mat.rowwise() - max) <<"\n";
  Eigen::VectorXd sum = (Mat.rowwise() - max).array().exp().colwise().sum().log();
  // Rcpp::Rcout << "out: " << max.transpose() + sum << "\n";
  return max.transpose() + sum;
  // return (Mat.array().exp().colwise().sum().log());
}

static inline double sinkhorn_converge(const Eigen::VectorXd & u, const Eigen::VectorXd & u_old) {
  Eigen::VectorXd f = u.array().log();
  Eigen::VectorXd f_old = u_old.array().log();
  Eigen::VectorXd diff = (f - f_old).array().abs();
  double out = (diff.array()/ f_old.array().abs()).sum();
  
  return(out);
}

static inline double sinkhorn_converge_log(const Eigen::VectorXd & f, const Eigen::VectorXd & f_old) {
  Eigen::VectorXd diff = (f - f_old).array().abs();
  double out = (diff.array()/ f_old.array().abs()).sum();
  
  return(out);
}

static inline matrix generate_S_(const matrix & cost, Eigen::VectorXd & f, Eigen::VectorXd & g, double eta) {
  matrix S = ((cost.colwise() - f).rowwise() - g.transpose()) * -eta;
  return(S);
}

static inline matrix generate_S_star_(const matrix & K, Eigen::VectorXd & f, Eigen::VectorXd & g) {
  matrix S = ((K.colwise() + f).rowwise() + g.transpose()) ;
  return(S);
}

static inline Eigen::VectorXd rowMin_eps_(const matrix & cost, Eigen::VectorXd & f, Eigen::VectorXd & g, double eta) {
  // return (rowLogSumExp((cost.rowwise() - g.transpose()) * - eta).array()/(-eta));
  matrix S = generate_S_(cost, f, g, eta);
  return(
    -(S.array().exp().rowwise().sum().log())/eta 
  );
}

static inline Eigen::VectorXd colMin_eps_(const matrix & cost, Eigen::VectorXd & f, Eigen::VectorXd & g, double eta) {
  
  // return (colLogSumExp((cost.colwise() - f) * - eta).array()/(-eta));
  matrix S = generate_S_(cost, f, g, eta);
  return(
    -(S.array().exp().colwise().sum().log())/eta 
  );
}

static inline Eigen::VectorXd rowMin_eps_star_(const matrix & K, Eigen::VectorXd & f, Eigen::VectorXd & g) {
  // return (rowLogSumExp((cost.rowwise() - g.transpose()) * - eta).array()/(-eta));
  matrix S = generate_S_star_(K, f, g);
  return(
    -(S.array().exp().rowwise().sum().log())
  );
}

static inline Eigen::VectorXd colMin_eps_star_(const matrix & K, Eigen::VectorXd & f, Eigen::VectorXd & g) {
  
  // return (colLogSumExp((cost.colwise() - f) * - eta).array()/(-eta));
  matrix S = generate_S_star_(K, f, g);
  return(
    -(S.array().exp().colwise().sum().log()) 
  );
}

static inline Eigen::VectorXd rowMin_eps_KL_(const matrix & cost, Eigen::VectorXd & f, Eigen::VectorXd & g, double eta,
                     Eigen::VectorXd & log_a, Eigen::VectorXd & log_b) {
  // return (rowLogSumExp((cost.rowwise() - g.transpose()) * - eta).array()/(-eta));
  matrix S = generate_S_(cost, f, g, eta);
  return(
    -((S.rowwise() + log_b.transpose()).array().exp().rowwise().sum().log())/eta 
  );
}

static inline Eigen::VectorXd colMin_eps_KL_(const matrix & cost, Eigen::VectorXd & f, Eigen::VectorXd & g, double eta,
                     Eigen::VectorXd & log_a, Eigen::VectorXd & log_b) {
  
  // return (colLogSumExp((cost.colwise() - f) * - eta).array()/(-eta));
  matrix S = generate_S_(cost, f, g, eta);
  return(
    -((S.colwise() + log_a).array().exp().colwise().sum().log())/eta 
  );
}


//' Generates approximate optimal transport plans 
//' using the sinkhorn algorithm
//'
//' @param mass_a A reference to an Eigen::VectorXd
//' giving the empirical mass in sample 1
//' @param mass_b A reference to an Eigen::VectorXd
//' giving the empirical mass in sample 2
//' @param exp_cost A reference to an Eigen::MatrixXd giving 
//' the exponentiated cost between samples A and B
//' @param A The assignment matrix
//' @param eta The inverse of lambda value used for 
//' the entropy penalty
//' @param epsilon The desired error bound
//' @param niterations The iterations to use for the methods
//' @param f The potentials of the dual problem
//' @param g The potentials of the dual problem
//' @return void
//' @keywords internal
// avoid log scaling
static inline void trans_sinkhorn(const refVecConst & mass_a, const refVecConst & mass_b,
                     const matrix & exp_cost,
                     matrix & A,
                     double eta, double epsilon, int niterations,
                     Eigen::VectorXd & f, Eigen::VectorXd & g) {
   int N = mass_a.size();
   int M = mass_b.size();
   
   // Eigen::VectorXd ones_n = Eigen::VectorXd::Ones(N);
   // Eigen::VectorXd ones_m = Eigen::VectorXd::Ones(M);
   
   Eigen::VectorXd u = Eigen::VectorXd::Ones(N); // first margins
   Eigen::VectorXd v = Eigen::VectorXd::Ones(M); // second margins
   
   Eigen::VectorXd u_old = u; // first margins
   // Eigen::VectorXd v_old = v; // second margins
   
   
   // matrix scaling
   for ( int i = 0; i < niterations; i ++){
     // column margins
     v = mass_b.cwiseQuotient(exp_cost.transpose() * u);
     
     // row margins
     u = mass_a.cwiseQuotient(exp_cost * v);
     
     // calc relative change in scaling vectors to see if approx converged
     if (i % 10) {
       if (sinkhorn_converge(u, u_old) <= epsilon) {
         break;
       } else {
         Rcpp::checkUserInterrupt();
       }
     }
     u_old = u;
     // v_old = v;
   }
   
   // get approximate assignment matrix
   A = u.asDiagonal() * exp_cost * v.asDiagonal();
   
   f = u.array().log() / eta;
   g = v.array().log() / eta;
 }

//' Self-sinkhorn distance
//'
//' @param f The potentials of the dual problem
//' @param mass_a A reference to an Eigen::VectorXd
//' giving the empirical mass in sample 1
//' @param eta The inverse of lambda value used for 
//' the entropy penalty
//' @param exp_cost A reference to an Eigen::MatrixXd giving 
//' the exponentiated cost between samples A and B
//' @param A The assignment matrix
//' @param epsilon The desired error bound
//' @param niterations The iterations to use for the methods
//' @return void
//' @keywords internal
static inline void trans_sinkhorn_self(Eigen::VectorXd & f, const refVecConst & mass_a,
                          double eta,
                          const matrix & exp_cost,
                          double epsilon, int niterations) {
   int N = mass_a.size();
   
   Eigen::VectorXd ones_n = Eigen::VectorXd::Ones(N);
   Eigen::VectorXd u = ones_n;
   Eigen::VectorXd u_old = ones_n; // first margins
   
   
   // matrix scaling
   for ( int i = 0; i < niterations; i ++){
     // row margins
     u = mass_a.cwiseQuotient(exp_cost * u);
     
     // calc relative change in scaling vectors to see if approx converged
     if (i % 10) {
       if (sinkhorn_converge(u, u_old) <= epsilon) {
         break;
       }
     }
     u_old = u;
   }
   
   f = u.array().log()/eta;
   
 }

//' Generates approximate optimal transport plans 
//' using the sinkhorn algorithm on the log scale
//'
//' @param mass_a A reference to an Eigen::VectorXd
//' giving the empirical mass in sample 1
//' @param mass_b A reference to an Eigen::VectorXd
//' giving the empirical mass in sample 2
//' @param cost A reference to an Eigen::MatrixXd giving 
//' the cost between samples A and B
//' @param Assign The assignment matrix
//' @param eta The inverse of lambda value used for 
//' the entropy penalty
//' @param epsilon The desired error bound
//' @param niterations The iterations to use for the methods
//' @param f_pot The potentials of the dual problem
//' @param g_pot The potentials of the dual problem
//' @return void
//' @keywords internal
 static inline void trans_sinkhorn_log(const refVecConst & mass_a, const refVecConst & mass_b,
                         const matrix & cost,
                         matrix & Assign,
                         double eta, double epsilon, int niterations,
                         Eigen::VectorXd & f_pot, Eigen::VectorXd & g_pot) {
   int N = mass_a.size();
   // int M = mass_b.size();
   
   Eigen::VectorXd log_a = mass_a.array().log();
   Eigen::VectorXd log_b = mass_b.array().log();
   
   // matrix K = -eta * cost.array();
   // Eigen::VectorXd A = -rowLogSumExp_(K) + log_a; // first margins
   // Eigen::VectorXd B = -colLogSumExp_(K.colwise() + A) + log_b; // second margins
   
   Eigen::VectorXd A = -rowLogSumExp_(-eta * cost.array()).array()/eta + log_a.array()/eta; // first margins
   Eigen::VectorXd B = -colLogSumExp_((-eta * (cost.colwise() - A).array())).array()/eta + log_b.array()/eta; // second margins
   
   // Eigen::VectorXd A = Eigen::VectorXd::Zero(N);
   // Eigen::VectorXd B = Eigen::VectorXd::Zero(M);
   
   Eigen::VectorXd A_old = Eigen::VectorXd::Zero(N); // first margins
   // Eigen::VectorXd B_old = Eigen::VectorXd::Zero(M); // second margins
   
   
   
   
   // Rcpp::Rcout << A << "\n";
   // Rcpp::Rcout << B << "\n";
   // matrix scaling
   for ( int i = 0; i < niterations; i ++) {
     
     // // row margins
     // A = -rowLogSumExp_(K.rowwise() + B.transpose()) + log_a;
     // // A += (colMin_eps_star(K, A, B) + log_a).eval();
     //   
     // // column margins
     // B = -rowLogSumExp_(K.colwise() + A) + log_b;
     // // B += (rowMin_eps_star(K, A, B) + log_b).eval();
     
     
     
     // // row margins
     // A += (rowMin_eps(cost, A, B, eta)  + log_a/eta).eval();
     // 
     // // column margins
     // B += (colMin_eps(cost, A, B, eta)  + log_b/eta).eval();
     
     // KL version of sinkhorn
     // row margins
     A += rowMin_eps_KL_(cost, A, B, eta, log_a, log_b).eval();
     
     // column margins
     B += colMin_eps_KL_(cost, A, B, eta, log_a, log_b).eval();
     
     
     // Rcpp::Rcout << i << ":\n ";
     // Rcpp::Rcout << "A: \n" << A.transpose() << "\n";
     // Rcpp::Rcout << "B: \n" << B.transpose() << "\n";
     // calc relative change in scaling vectors to see if approx converged
     // if (i % 10) {
     if (sinkhorn_converge_log(A, A_old) <= epsilon) {
       break;
     } else {
       Rcpp::checkUserInterrupt();
     }
     // }
     A_old = A;
     // B_old = B;
   }
   // f_pot = A / eta;
   // g_pot = B / eta;
   f_pot = A;
   g_pot = B;
   // Rcpp::Rcout << "epsilon: " << epsilon <<"\n";
   // Rcpp::Rcout << (mass_a.array() * f.array()).sum() +
   //   (mass_b.array() * g.array()).sum() << "\n";
   // Rcpp::Rcout << (mass_a.array() * (u.array().log() - u_pot.array().log()) / eta).sum() +
   //   (mass_b.array() * (v.array().log() - v_pot.array().log()) / eta).sum();
   // f = f - f_pot;
   // g = g - g_pot;
   
   // Rcpp::Rcout << "using biased potentials " << ((A.asDiagonal() * K.exp().matrix() * B.asDiagonal()).array() * cost.array()).sum() << "\n";
   
   // Rcpp::Rcout << (mass_a.array() * f.array()).sum() +
   //   (mass_b.array() * g.array()).sum() << "\n";
   
   // Rcpp::Rcout << "using biased potentials " << ((A.array().exp().matrix().asDiagonal() * K.array().exp().matrix() * B.array().exp().matrix().asDiagonal()).array() * cost.array()).sum() << "\n";
   
   // Rcpp::Rcout << "using biased potentials from paper " << ((K + ( A * Eigen::VectorXd::Ones(M).transpose() + Eigen::VectorXd::Ones(N) * B.transpose())).exp().array() * (mass_a * mass_b.transpose()).array() * cost.array()).sum() << "\n";
   
   // unbiased potentials if available
   // A = f * eta;
   // B = g * eta;
   
   // get approximate assignment matrix
   // Assign = ((K.colwise() + A).rowwise() + B.transpose()).array().exp();
   matrix S = generate_S_(cost, f_pot, g_pot, eta);
   Assign = (S.array() - LogSumExp(S)).array().exp();
   // Rcpp::Rcout << "using unbiased potentials " << ((A.array().exp().matrix().asDiagonal() * K.array().exp().matrix() * B.array().exp().matrix().asDiagonal()).array() * cost.array()).sum() << "\n";
   // Rcpp::Rcout << "using Assign from paper" << (Assign.array() * cost.array()).sum() << "\n";
   
 }

//' Autocorr sinkhorn distance
//'
//' @param f The potentials of the dual problem
//' @param mass_a A reference to an Eigen::VectorXd
//' giving the empirical mass in sample 1
//' @param eta The inverse of lambda value used for 
//' the entropy penalty
//' @param exp_cost A reference to an Eigen::MatrixXd giving 
//' the exponentiated cost between samples A and B
//' @param A The assignment matrix
//' @param epsilon The desired error bound
//' @param niterations The iterations to use for the methods
//' @return void
//' @keywords internal
 static inline void trans_sinkhorn_autocorr_log(Eigen::VectorXd & f, const refVecConst & mass_a,
                                  const matrix & cost,
                                  double eta, double epsilon, int niterations) {
   int N = mass_a.size();
   
   Eigen::VectorXd log_a = mass_a.array().log();
   
   // matrix K = -eta * cost.array();
   // Eigen::VectorXd A = -rowLogSumExp_(K) + log_a; // first margins
   // Eigen::VectorXd B = -colLogSumExp_(K.colwise() + A) + log_b; // second margins
   
   f = -rowLogSumExp_(-eta * cost.array()).array()/eta + log_a.array()/eta; // first margins
   
   // Eigen::VectorXd A = Eigen::VectorXd::Zero(N);
   // Eigen::VectorXd B = Eigen::VectorXd::Zero(M);
   
   Eigen::VectorXd f_old = Eigen::VectorXd::Zero(N); // second margins
   
   
   
   
   // Rcpp::Rcout << A << "\n";
   // Rcpp::Rcout << B << "\n";
   // matrix scaling
   for ( int i = 0; i < niterations; i ++) {
     
     
     // column margins
     f += colMin_eps_KL_(cost, f, f, eta, log_a, log_a).eval();
     
     
     if (sinkhorn_converge_log(f, f_old) <= epsilon) {
       break;
     }
     // }
     f_old = f;
   }
   // f_pot = A / eta;
   // g_pot = B / eta;
 }
#endif //TRANS_SINKHORN_H