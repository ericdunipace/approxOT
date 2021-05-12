#include "../inst/include/trans_gandkhorn.h"


//' Generates approximate optimal transport plans
//'
//' @param a An Eigen::VectorXd
//' @param b An Eigen::VectorXd
//' @return returns an integer giving the index of the
//' maximum difference.
//' @details Gives the index of the largest difference between
//' vectors a and b.
//' @keywords internal
int argmax_abs (const vector & a, const vector & b) {
  vector r = (a - b).cwiseAbs();
  Eigen::Index maxIndex;
  
  r.maxCoeff(&maxIndex);
  
  int which_max = maxIndex;
  return(which_max);
}

//' Generates approximate optimal transport plans 
//' using the gandkhorn algorithm
//'
//' @param mass_r A reference to an Eigen::VectorXi 
//' giving the sampled indicators
//' @param mass_b A reference to an Eigen::VectorXi 
//' giving the sampled indicators
//' @param exp_cost A reference to an Eigen::MatrixXd giving 
//' the cost between samples A and B
//' @param B The assignment matrix
//' @param eta The inverse of lambda value used for 
//' the entropy penalty
//' @param epsilon The desired error bound
//' @param niterations The iterations to use for the methods
//' @return void
//' @keywords internal
void trans_gandkhorn(const refVecConst & mass_r, 
                     const refVecConst & mass_l, 
                      const matrix & exp_cost, 
                      matrix & B,
                      double eta, double epsilon, int niterations) {
  int N = mass_r.size();
  int M = mass_l.size();
  double Nd = double(N);
  double theta = 2.0;//R::runif(0.0, 2.0);
  vector samp_weights = vector::Constant(N, 1.0/ Nd);

  matrix A = exp_cost.array()/exp_cost.lpNorm<1>();
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
    B_bar = u_bar.array().exp().matrix().asDiagonal() * A * v_bar.array().exp().matrix().asDiagonal();
    // B_bar /= B_bar.lpNorm<1>();
    
    //step 3
    r_bar = B_bar.rowwise().sum();
    l_bar = B_bar.colwise().sum();
    
    int I = argmax_abs(mass_r,r_bar);
    int J = argmax_abs(mass_l,l_bar);
    
    if (rho(mass_r(I), r_bar(I)) >= rho(mass_l(J), l_bar(J)) ) {
      u_hat = u_bar;
      u_hat(I) += log_r(I) - std::log(r_bar(I));
      v_hat = v_bar;
    } else {
      u_hat = u_bar;
      v_hat = v_bar;
      v_hat(J) += log_l(J) - std::log(l_bar(J));
    }
    
    //step 4
    double xi = R::runif(0.0,1.0);
    
    // void sample_systematic(vectorI & samps, const vector & weight, const int nsamp )
    vectorI pi(1);
    sample_systematic(pi, samp_weights, 1);

    if(xi <= 0.5) {
      u_tilde(pi(0)) -= eta/(8 * theta * Nd) * (r_bar(pi(0)) - mass_r(pi(0))); //note, eta here is reverse of paper my_eta = 1/their_eta
      //v_tilde = v_tilde
    } else {
      //u_tilde = u_tilde
      v_tilde(pi(0)) -= eta/(8 * theta * Nd) * (l_bar(pi(0)) - mass_l(pi(0)));
    }
    
    //step 5
    argmin_f(mass_r, mass_l, A,
                  u, v, y_u, y_v, u_hat, v_hat);
    
    //step 6
    B = u.array().exp().matrix().asDiagonal() * A * v.array().exp().matrix().asDiagonal();
    // B /= B.lpNorm<1>();
    r = B.rowwise().sum();
    l = B.colwise().sum();
    // if((r.array() < 0).any()) Rcpp::stop("r < 0");
    // if((l.array() < 0).any()) Rcpp::stop("l < 0");
    if ( dist_approx_ot(mass_r, mass_l, r, l, 1.0) <= epsilon ) {
      break;
    }
    // Rcpp::Rcout <<"dist: "<< dist_approx_ot(mass_r, mass_l, r, l, 1) << "\n";
    // Rcpp::Rcout <<"B: " << B(0,0) << ", ";
    // Rcpp::Rcout <<"u: " << u(0) << ", ";
    // Rcpp::Rcout <<"v: " << v(0) << "\n";
    
    I = argmax_abs(mass_r, r);
    J = argmax_abs(mass_l, l);
    
    if( rho(mass_r(I), r(I)) >= rho(mass_l(J), l(J))) {
      y_u = u;
      y_u(I) += log_r(I) - std::log(r(I));
      y_v = v;
    } else {
      y_u = u;
      y_v = v;
      y_v(J) += log_l(J) - std::log(l(J));
    }
  } //end for
  B /= B.lpNorm<1>();
  
}
