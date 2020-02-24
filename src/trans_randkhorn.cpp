#include "../inst/include/trans_randkhorn.h"

// algorithm 3 in Altschuler, J., Weed, J., & Rigollet, P. (2017). 
// Near-linear time approximation algorithms for optimal transport via Sinkhorn iteration. 
// 31st Congerence on Neural Information Processing Systems, (1), 1–11. Long Beach, CA.
void trans_randkhorn(const refVecConst & mass_r, const refVecConst & mass_l, 
                    const matrix & exp_cost, 
                    matrix & B,
                    double eta, double epsilon, int niterations) {
  int N = mass_r.size();
  int M = mass_l.size();
  // if(N != M) Rcpp::stop("")
  double theta = 2.0;//R::runif(0.0, 2.0);
  matrix A = exp_cost/exp_cost.sum();
  
  // A = exp_cost.array()/exp_cost.lpNorm<1>();
  vector log_r = mass_r.array().log();
  vector log_l = mass_l.array().log();
  vector r = vector::Zero(N);
  vector r_bar = vector::Zero(N);
  vector l = vector::Zero(M);
  vector l_bar = vector::Zero(M);
  vector y_u = vector::Zero(N);
  vector y_v = vector::Zero(M);
  vector u = vector::Zero(N);
  vector u_hat = vector::Zero(N);
  vector u_tilde = vector::Zero(N);
  vector u_bar = vector::Zero(N);
  vector v = vector::Zero(M);
  vector v_hat = vector::Zero(M);
  vector v_tilde = vector::Zero(M);
  vector v_bar = vector::Zero(M);
  matrix B_bar = matrix::Zero(N, M);
  
  for( int i = 0; i < niterations; i ++){
    //step 1
    theta = theta/2.0 * (std::sqrt(theta * theta + 4.0) - theta);
    
    //step 2
    u_bar = (1 - theta) * y_u + theta * u_tilde;
    v_bar = (1 - theta) * y_v + theta * v_tilde;
    B_bar = (u_bar.array().exp()).matrix().asDiagonal() * A * v_bar.array().exp().matrix().asDiagonal();
    // B_bar /= B_bar.lpNorm<1>();
    
    //step 3
    r_bar = B_bar.rowwise().sum();
    l_bar = B_bar.colwise().sum();
    if ( rho_ot(mass_r, r_bar) >= rho_ot(mass_l, l_bar) ) {
      u_hat = u_bar + log_r;
      u_hat.array() -= r_bar.array().log();
      v_hat = v_bar;
    } else {
      u_hat = u_bar;
      v_hat = v_bar + log_l;
      v_hat.array() -= l_bar.array().log();
    }
    
    //step 4
    double xi = R::runif(0.0,1.0);
    // Rcpp::Rcout << xi << "\n";
    if(xi <= 0.5) {
      u_tilde -= eta/(8 * theta ) * (r_bar - mass_r); //their eta is 1.0/my_eta
      //v_tilde = v_tilde
    } else {
      //u_tilde = u_tilde
      v_tilde -= eta/(8 * theta) * (l_bar - mass_l);
    }
    
    //step 5
    argmin_f(mass_r, mass_l, exp_cost,
                  u, v, y_u, y_v, u_hat, v_hat);
    
    //step 6
    B = u.array().exp().matrix().asDiagonal() * A * v.array().exp().matrix().asDiagonal();
    // B /= B.lpNorm<1>();
    r = B.rowwise().sum();
    l = B.colwise().sum();
    
    // if((r.array() < 0).any()) Rcpp::stop("r < 0");
    // if((l.array() < 0).any()) Rcpp::stop("l < 0");
    
    if ( dist_approx_ot(mass_r, mass_l,r,l, 1.0) <= epsilon ) {
      // Rcpp::Rcout << "iter: " << i << "\n";
      break;
    }
    
    // Rcpp::Rcout <<"dist: "<< dist_approx_ot(mass_r, mass_l, r, l, 1) << "\n";
    // Rcpp::Rcout <<"B: " << B(0,0) << ", ";
    // Rcpp::Rcout <<"u: " << u(0) << ", ";
    // Rcpp::Rcout <<"v: " << v(0) << "\n";
    
    if( rho_ot(mass_r, r) >= rho_ot(mass_l, l)) {
      y_u = u + log_r;
      y_u.array() -= r.array().log();
      y_v = v;
    } else {
      y_u = u;
      y_v = v + log_l;
      y_v.array() -= l.array().log();
    }
  } //end for
  B /= B.lpNorm<1>();
  
}