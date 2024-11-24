#include "approxOT/trans_sinkhorn.h"

//[[Rcpp::export]]
Eigen::VectorXd rowLogSumExp(matrix Mat) {
  
  return rowLogSumExp_(Mat);
}

//[[Rcpp::export]]
Eigen::VectorXd colLogSumExp(matrix Mat) {
  
  return colLogSumExp_(Mat);
}

//[[Rcpp::export]]
matrix generate_S(const matrix & cost, Eigen::VectorXd & f, Eigen::VectorXd & g, double eta) {
  return generate_S_(cost, f, g, eta);
}

//[[Rcpp::export]]
Eigen::VectorXd rowMin_eps(const matrix & cost, Eigen::VectorXd & f, Eigen::VectorXd & g, double eta) {
  return rowMin_eps_(cost,f,g,eta);
}

//[[Rcpp::export]]
Eigen::VectorXd colMin_eps(const matrix & cost, Eigen::VectorXd & f, Eigen::VectorXd & g, double eta) {
  return colMin_eps_(cost,f,g,eta);
}


//[[Rcpp::export]]
Eigen::VectorXd rowMin_eps_KL(const matrix & cost, Eigen::VectorXd & f, Eigen::VectorXd & g, double eta,
                     Eigen::VectorXd & log_a, Eigen::VectorXd & log_b) {
  return rowMin_eps_KL_(cost,f,g,eta, log_a, log_b);
}

//[[Rcpp::export]]
Eigen::VectorXd colMin_eps_KL(const matrix & cost, Eigen::VectorXd & f, Eigen::VectorXd & g, double eta,
                     Eigen::VectorXd & log_a, Eigen::VectorXd & log_b) {
  
  return colMin_eps_KL_(cost, f, g, eta, log_a, log_b);
}

//[[Rcpp::export]]
Rcpp::List sinkhorn_pot_(const Eigen::VectorXd & mass_a, const Eigen::VectorXd & mass_b, 
                         const matrix & cost_matrix, 
                         double epsilon, int niterations,
                         bool unbiased,
                         const matrix & cost_matrix_A, 
                         const matrix & cost_matrix_B) {
  double med_cost = median(cost_matrix);
  double eta = 1.0 / (epsilon * med_cost); //avoid underflow
  // double eta = 4 * log(double(mass_a.size())) / epsilon;
  const matrix exp_cost = (-eta * cost_matrix.array() ).exp();
  double epsilon_prime = epsilon / (8 * cost_matrix.maxCoeff());
  
  
  int N = mass_a.size();
  int M = mass_b.size();
  
  matrix assign_mat = matrix::Zero(N,M);
  
  Eigen::VectorXd f = Eigen::VectorXd::Ones(N);
  Eigen::VectorXd g = Eigen::VectorXd::Ones(M);
  
  
  if (unbiased) {
    const matrix exp_cost_a = (-eta * cost_matrix_A.array() ).exp();
    const matrix exp_cost_b = (-eta * cost_matrix_B.array() ).exp();
    Eigen::VectorXd p = Eigen::VectorXd::Ones(N);
    Eigen::VectorXd p_unused = Eigen::VectorXd::Ones(N);
    Eigen::VectorXd q = Eigen::VectorXd::Ones(M);
    Eigen::VectorXd q_unused = Eigen::VectorXd::Ones(M);
    matrix assign_mat_a = matrix::Zero(N,N);
    matrix assign_mat_b = matrix::Zero(M,M);
    
    trans_sinkhorn(mass_a, mass_a, exp_cost_a, assign_mat_a, eta, epsilon_prime/2.0, niterations,
                   p, p_unused);
    trans_sinkhorn(mass_b, mass_b, exp_cost_b, assign_mat_b, eta, epsilon_prime/2.0, niterations,
                   q_unused, q);
    trans_sinkhorn(mass_a, mass_b, exp_cost, assign_mat, eta, epsilon_prime/2.0, niterations,
                       f, g);
    // Rcpp::Rcout << f << "\n";
    // Rcpp::Rcout << g << "\n";
    // Rcpp::Rcout << p << "\n";
    // Rcpp::Rcout << q << "\n";
    f = f - p;
    g = g - q;
    
  } else {
    trans_sinkhorn(mass_a, mass_b, exp_cost, assign_mat, eta, epsilon_prime/2.0, niterations,
                   f, g);
  }
  
  return(Rcpp::List::create(Rcpp::Named("f") = Rcpp::wrap(f),
                            Rcpp::Named("g") = Rcpp::wrap(g)));
}


//[[Rcpp::export]]
Rcpp::List sinkhorn_pot_log_(const Eigen::VectorXd & mass_a, const Eigen::VectorXd & mass_b, 
                         const matrix & cost_matrix, 
                         double epsilon, int niterations,
                         bool unbiased,
                         const matrix & cost_matrix_A, 
                         const matrix & cost_matrix_B) {
  double med_cost = median(cost_matrix);
  double eta = 1.0 / (epsilon * med_cost); //avoid underflow
  // double eta = 4 * log(double(mass_a.size())) / epsilon;
  // const matrix exp_cost = (-eta * cost_matrix.array() ).exp();
  double epsilon_prime = 1e-8;
  
  
  int N = mass_a.size();
  int M = mass_b.size();
  
  matrix assign_mat = matrix::Zero(N,M);
  
  Eigen::VectorXd f = Eigen::VectorXd::Ones(N);
  Eigen::VectorXd g = Eigen::VectorXd::Ones(M);
  
  
  if (unbiased) {
    // const matrix exp_cost_a = (-eta * cost_matrix_A.array() ).exp();
    // const matrix exp_cost_b = (-eta * cost_matrix_B.array() ).exp();
    Eigen::VectorXd p = Eigen::VectorXd::Ones(N);
    Eigen::VectorXd p_unused = Eigen::VectorXd::Ones(N);
    Eigen::VectorXd q = Eigen::VectorXd::Ones(M);
    Eigen::VectorXd q_unused = Eigen::VectorXd::Ones(M);
    // matrix assign_mat_a = matrix::Zero(N,N);
    // matrix assign_mat_b = matrix::Zero(M,M);
    
    // Eigen::VectorXd & f, const refVecConst & mass_a,
    // const matrix & cost,
    // double eta, double epsilon, int niterations
    trans_sinkhorn_autocorr_log(p, mass_a, cost_matrix_A, eta, epsilon_prime, niterations);
    trans_sinkhorn_autocorr_log(q, mass_b, cost_matrix_B, eta, epsilon_prime, niterations);
    
    
    // const refVecConst & mass_a, const refVecConst & mass_b,
    // const matrix & cost,
    // matrix & Assign,
    // double eta, double epsilon, int niterations,
    // const refVecConst & f_pot, const refVecConst & g_pot
    trans_sinkhorn_log(mass_a, mass_b, cost_matrix, assign_mat, eta, epsilon_prime, niterations,
                   f, g);
    // Rcpp::Rcout << f << "\n";
    // Rcpp::Rcout << g << "\n";
    // Rcpp::Rcout << p << "\n";
    // Rcpp::Rcout << q << "\n";
    f = f - p;
    g = g - q;
    
  } else {
    trans_sinkhorn_log(mass_a, mass_b, cost_matrix, assign_mat, eta, epsilon_prime, niterations,
                       f, g);
  }
  
  return(Rcpp::List::create(Rcpp::Named("f") = Rcpp::wrap(f),
                            Rcpp::Named("g") = Rcpp::wrap(g)));
}
