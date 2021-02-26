#include "../inst/include/trans_sinkhorn.h"

// algorithm 3 in Altschuler, J., Weed, J., & Rigollet, P. (2017). 
// Near-linear time approximation algorithms for optimal transport via Sinkhorn iteration. 
// 31st Conference on Neural Information Processing Systems, (1), 1–11. Long Beach, CA.
// void trans_sinkhorn(const refVecConst & mass_a, const refVecConst & mass_b, 
//                     const matrix & exp_cost, 
//                     matrix & A,
//                     double eta, double epsilon, int niterations) {
//   int N = mass_a.size();
//   int M = mass_b.size();
//   
//   
//   matrix A_0 = exp_cost;
//   A = exp_cost.array()/exp_cost.lpNorm<1>();
//   vector log_a = mass_a.array().log();
//   vector log_b = mass_b.array().log();
//   vector x = vector::Zero(N);
//   vector y = vector::Zero(M);
//   vector r = vector::Zero(N);
//   vector c = vector::Zero(M);
//   
//   for( int i = 0; i < niterations; i ++){
//     r = A.rowwise().sum();
//     c = A.colwise().sum();
//     
//     if ((i % 2) == 0) {
//       x.noalias() += log_a;
//       x.array() -= r.array().log();
//       
//     } else {
//       y.noalias() += log_b;
//       y.array() -= c.array().log();
//       
//     }
//     A = x.array().exp().matrix().asDiagonal() * A_0 * y.array().exp().matrix().asDiagonal();
//     // A /= A.sum();
//     if(dist_approx_ot(mass_a, mass_b, r, c, 2.0) <= epsilon) {
//       break;
//     }
//   }
// }


// void trans_sinkhorn(const refVecConst & mass_a, const refVecConst & mass_b,
//                     const matrix & exp_cost,
//                     matrix & A,
//                     double eta, double epsilon, int niterations) {
//   int N = mass_a.size();
//   int M = mass_b.size();
// 
// 
//   // matrix A_0 = exp_cost;
//   A = exp_cost.array()/exp_cost.lpNorm<1>();
//   vector log_a = mass_a.array().log();
//   vector log_b = mass_b.array().log();
//   // vector x = vector::Zero(N);
//   // vector y = vector::Zero(M);
//   vector r = A.rowwise().sum();
//   vector c = A.colwise().sum();
// 
//   for( int i = 0; i < niterations; i ++){
// 
//     if ((i % 2) == 0) {
//       vector x = log_a;
//       x.array() -= r.array().log();
// 
//       A = x.array().exp().matrix().asDiagonal() * A;
//       r = mass_a;
//       c = A.colwise().sum();
//       // Rcpp::Rcout << "Even: ";
//       // // Rcpp::Rcout << "x: " << x.transpose() << " ";
//       // Rcpp::Rcout << "Rsum: " << A.rowwise().sum().sum() << " ";
//       // Rcpp::Rcout << "Csum: " << A.colwise().sum().sum() << " ";
//     } else {
//       vector y = log_b;
//       y.array() -= c.array().log();
// 
//       A = A * y.array().exp().matrix().asDiagonal();
//       r = A.rowwise().sum();
//       c = mass_b;
//       // Rcpp::Rcout << "Odd: ";
//       // // Rcpp::Rcout << "y: " << y.transpose() << " ";
//       // Rcpp::Rcout << "Rsum: " << A.rowwise().sum().sum() << " ";
//       // Rcpp::Rcout << "Csum: " << A.colwise().sum().sum() << " ";
//     }
//     // A /= A.sum();
//     // r = A.rowwise().sum();
//     // c = A.colwise().sum();
//     // Rcpp::Rcout << "R: " << r.transpose() << " ";
//     // Rcpp::Rcout << "Csum: " << A.colwise().sum().sum() << " ";
//     // Rcpp::Rcout << "Asum " << A.sum() << "\n";
//     if(dist_approx_ot(mass_a, mass_b, r, c, 1.0) <= epsilon) {
//       // Rcpp::Rcout << epsilon << "\n";
//       // Rcpp::Rcout << dist_approx_ot(mass_a, mass_b, r, c, 2.0) << "\n";
//       break;
//     }
//   }
//   A /= A.sum();
// }


// avoid log scaling
void trans_sinkhorn(const refVecConst & mass_a, const refVecConst & mass_b,
                    const matrix & exp_cost,
                    matrix & A,
                    double eta, double epsilon, int niterations) {
  int N = mass_a.size();
  int M = mass_b.size();

  vector ones_n = vector::Ones(N);
  vector ones_m = vector::Ones(M);

  vector u = ones_n; // first margins
  vector v = ones_m; // second margins
  
  vector u_old = ones_n; // first margins
  vector v_old = ones_m; // second margins
  
  
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
      }
    }
    u_old = u;
    v_old = v;
  }
  
  // get approximate assignment matrix
  A = u.asDiagonal() * exp_cost * v.asDiagonal();
}

vector rowLogSumExp(matrix  Mat) {
  
  vector max = Mat.rowwise().maxCoeff();
  
  vector sum = (Mat.colwise() - max).array().exp().rowwise().sum().log();
  return max + sum;
}

vector colLogSumExp(matrix Mat) {
  
  Eigen::RowVectorXd max = Mat.colwise().maxCoeff();
  vector sum = (Mat.rowwise() - max).array().exp().colwise().sum().log();
  
  return max.transpose() + sum;
}

void trans_sinkhorn_log(const refVecConst & mass_a, const refVecConst & mass_b,
                    const matrix & cost,
                    matrix & Assign,
                    double eta, double epsilon, int niterations,
                    const refVecConst & f_pot, const refVecConst & g_pot) {
  int N = mass_a.size();
  int M = mass_b.size();
  
  vector A = vector::Zero(N); // first margins
  vector B = vector::Zero(M); // second margins
  
  vector A_old = vector::Zero(N); // first margins
  vector B_old = vector::Zero(M); // second margins
  
  vector log_a = mass_a.array().log();
  vector log_b = mass_b.array().log();
  
  matrix K = -eta * cost.array();
  
  // matrix scaling
  for ( int i = 0; i < niterations; i ++) {
    // column margins
    B = -colLogSumExp(K.colwise() + (A + log_a));
    
    // row margins
    A = -rowLogSumExp(K.rowwise() +  (B + log_b).transpose());
    
    // calc relative change in scaling vectors to see if approx converged
    // if (i % 10) {
      if (sinkhorn_converge_log(A, A_old) <= epsilon) {
        break;
      }
    // }
    A_old = A;
    B_old = B;
  }
  vector f = A / eta;
  vector g = B / eta;
  
  // Rcpp::Rcout << (mass_a.array() * f.array()).sum() +
  //   (mass_b.array() * g.array()).sum() << "\n";
  // Rcpp::Rcout << (mass_a.array() * (u.array().log() - u_pot.array().log()) / eta).sum() +
  //   (mass_b.array() * (v.array().log() - v_pot.array().log()) / eta).sum();
  f = f - f_pot;
  g = g - g_pot;
  
  // Rcpp::Rcout << "using biased potentials " << ((A.asDiagonal() * K.exp().matrix() * B.asDiagonal()).array() * cost.array()).sum() << "\n";
  
  // Rcpp::Rcout << (mass_a.array() * f.array()).sum() +
  //   (mass_b.array() * g.array()).sum() << "\n";
  
  // Rcpp::Rcout << "using biased potentials " << ((A.array().exp().matrix().asDiagonal() * K.array().exp().matrix() * B.array().exp().matrix().asDiagonal()).array() * cost.array()).sum() << "\n";
  
  // Rcpp::Rcout << "using biased potentials from paper " << ((K + ( A * vector::Ones(M).transpose() + vector::Ones(N) * B.transpose())).exp().array() * (mass_a * mass_b.transpose()).array() * cost.array()).sum() << "\n";
  
  // unbiased potentials if available
  A = f * eta;
  B = g * eta;
  
  // get approximate assignment matrix
  Assign = (K + ( A * vector::Ones(M).transpose() + vector::Ones(N) * B.transpose())).array().exp() * (mass_a * mass_b.transpose()).array();
  
  // Rcpp::Rcout << "using unbiased potentials " << ((A.array().exp().matrix().asDiagonal() * K.array().exp().matrix() * B.array().exp().matrix().asDiagonal()).array() * cost.array()).sum() << "\n";
  // Rcpp::Rcout << "using Assign from paper" << (Assign.array() * cost.array()).sum() << "\n";
  
}

void trans_sinkhorn_autocorr(vector & f, const refVecConst & mass_a,
                             const matrix & cost,
                             double eta, double epsilon, int niterations) {
  int N = mass_a.size();
  
  vector log_a = mass_a.array().log();

  vector A = vector::Zero(N); // first margins
  
  vector A_old = vector::Zero(N); // first margins save old iteration
  
  matrix K = (-cost * eta).colwise() + log_a;

  // matrix scaling
  for ( int i = 0; i < niterations; i ++) {
    
    // row margins
    A = 0.5 * (A - colLogSumExp(K.colwise() + A));
    // calc relative change in scaling vectors to see if approx converged
    // if (i % 3) { // calc every 10 iterations to avoid log calculation every step
    if (sinkhorn_converge_log(A, A_old) <= epsilon) {
      break;
    }
    // }
    
    // save current iteration
    A_old = A;
  }
  f = -colLogSumExp(K.colwise() + A).array() / eta;
  // get potential
  // f = u.array().log();
}

// void trans_sinkhorn_autocorr(vector & u_pot, const refVecConst & mass_a,
//                     const matrix & exp_cost,
//                     double eta, double epsilon, int niterations) {
//   int N = mass_a.size();
// 
//   vector ones_n = vector::Ones(N);
// 
//   vector u = ones_n; // first margins
// 
//   vector u_old = ones_n; // first margins save old iteration
// 
//   matrix K = exp_cost * mass_a.asDiagonal();
// 
//   // matrix scaling
//   for ( int i = 0; i < niterations; i ++){
// 
//     // row margins
//     u = (u.array() / (K * u).array()).sqrt();
// 
//     // calc relative change in scaling vectors to see if approx converged
//     // if (i % 3) { // calc every 10 iterations to avoid log calculation every step
//       if (sinkhorn_converge(u, u_old) <= epsilon) {
//         break;
//       }
//     // }
// 
//     // save current iteration
//     u_old = u;
//   }
// 
//   // get potential
//   // f = u.array().log() / eta;
//   u_pot = u;
// }
