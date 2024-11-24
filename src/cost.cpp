#include "approxOT/cost.h"


//[[Rcpp::export]]
Rcpp::NumericMatrix cost_calculation_(const Rcpp::NumericMatrix & A_, const Rcpp::NumericMatrix & B_, const double p) {
  int N = A_.cols();
  int M = B_.cols();
  
  const matMap A(Rcpp::as<matMap >(A_));
  const matMap B(Rcpp::as<matMap >(B_));
  
  matrix cost_matrix(N,M);
  
  if(p == 2.0) {
    cost_calculation_L2(A, B, cost_matrix);
  } else if (p == 1.0){
    cost_calculation_L1(A, B, cost_matrix);
  } else {
    cost_calculation_Lp(A, B, cost_matrix, p);
  }
  
  return Rcpp::wrap(cost_matrix);
}


double multi_marg_final_cost_L2(const Rcpp::List & idx_, 
                                const Rcpp::List & data_, 
                                const Rcpp::NumericVector & mass_,
                                int M, int D,
                                double p, double ground_p) {
  int N = idx_.size();
  double cost = 0.0;
  double scale = double(N);
  
  for( int m = 0; m < M; m ++) {
    matrix temp(D,N);
    for(int n = 0; n < N; n ++) {
      int cur_idx = Rcpp::as<Rcpp::IntegerVector>(idx_[n])[m] - 1; // -1 to adjust R to C indexing
      
      for(int d = 0; d < D; d ++) temp(d,n) = Rcpp::as<Rcpp::NumericMatrix>(data_[n])(d, cur_idx);
    }
    // Rcpp::Rcout << temp.sum()<<", ";
    Eigen::VectorXd mean = temp.rowwise().mean();
    // Rcpp::Rcout << mean.sum()<<", ";
    temp.noalias() = temp.colwise() - mean;
    // Rcpp::Rcout << temp.sum()<<", ";
    // Rcpp::Rcout << temp.array().abs().pow(ground_p).sum() * scale << "\n";
    cost += std::pow(std::pow(temp.array().abs().pow(ground_p).sum() * scale, 1.0/ground_p), p) * mass_(m);
  }
  
  return std::pow(cost, 1.0/p);
}

//[[Rcpp::export]]
double multi_marg_final_cost_(const Rcpp::List & idx_, 
                              const Rcpp::List & data_, 
                              const Rcpp::NumericVector & mass_,
                              int M, int D,
                              double p, double ground_p) {
  
  if(ground_p == 2.0) {
    // Rcpp::Rcout <<"using l2";
    return(multi_marg_final_cost_L2(idx_, data_, mass_, M, D, p, ground_p));
  } else {
    // Rcpp::Rcout <<"using general";
    int N = idx_.size();
    double cost = 0.0;
    for( int m = 0; m < M; m ++) {
      matrix temp(D,N);
      Eigen::VectorXd cost_vec = Eigen::VectorXd::Zero(D);
      
      for(int n = 0; n < N; n ++) {
        int cur_idx = Rcpp::as<Rcpp::IntegerVector>(idx_[n])[m] - 1; // -1 to adjust R to C indexing
        
        for(int d = 0; d < D; d ++) temp(d,n) = Rcpp::as<Rcpp::NumericMatrix>(data_[n])(d, cur_idx);
      }
      
      for(int n = 0; n < (N-1); n ++) {
        Eigen::Ref<Eigen::VectorXd> temp_col = temp.col(n);
        for(int nn = n+1; nn < N; nn++) {
          cost_vec += (temp_col - temp.col(nn)).array().abs().pow(ground_p).matrix();
        }
      }
      cost += std::pow(std::pow(cost_vec.sum(), 1.0/ground_p), p) * mass_(m);
    }
    
    return std::pow(cost, 1.0/p);
  }
  
}


//[[Rcpp::export]]
double multi_marg_given_dist_(const Rcpp::List & idx_, 
                              const Rcpp::NumericVector & mass_,
                              const Rcpp::NumericVector & cost_,
                              int M, 
                              int N_cost,
                              double p) 
{
  
  int N = idx_.size();
  double cost = 0.0;
  for (int m = 0; m < M; m ++) {
    Rcpp::IntegerVector cur_idxs(N);
    double cur_mass = mass_(m);
    for (int n = 0 ; n < N; n++) {
      // Rcpp::Rcout << Rcpp::as<Rcpp::IntegerVector>(idx_[n])[m] << ", ";
      cur_idxs(n) = Rcpp::as<Rcpp::IntegerVector>(idx_[n])[m] - 1; //adapt to R index
    }
    for (int n = 0; n < (N-1); ++n) {
      for (int nn = n+1; nn < N; ++nn) {
        if(cur_idxs(n) != cur_idxs(nn)) {
          cost += std::pow(cost_(dist_2d_to_1d_(cur_idxs(n), cur_idxs(nn), N_cost)), p) * cur_mass;
        }
        // Rcpp::Rcout << "M: " << m <<", " <<
        //   cur_idxs(n) << ", " << cur_idxs(nn) << ", " <<
        //   dist_2d_to_1d_(cur_idxs(n), cur_idxs(nn), N_cost) << "\n";
      }
    }
  }
  
  return std::pow(cost, 1.0/p);
  
}