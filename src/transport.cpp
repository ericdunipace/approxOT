#include "approxOT/transport.h"


// .C("shortsimplex",
// as.integer(control$para$slength), 
// as.integer(control$para$kfound), 
// as.double(control$para$psearched),
// as.integer(N), 
// as.integer(N), 
// as.integer(rep.int(1,N)), 
// as.integer(rep.int(1,N)),
// as.double(dd), 
// assignment = as.integer(initassig), 
//   basis = as.integer(initbasis),
// DUP=TRUE, PACKAGE="transport")


//[[Rcpp::export]]
Rcpp::List transport_C_(const Rcpp::NumericVector & mass_a_, const Rcpp::NumericVector & mass_b_, 
                        const Rcpp::NumericMatrix & cost_matrix_,
                        const Rcpp::CharacterVector & method_,
                        double epsilon_,
                        int niter_,
                        bool unbiased_,
                        int threads_,
                        const Rcpp::NumericMatrix &  cost_matrix_A_, 
                        const Rcpp::NumericMatrix &  cost_matrix_B_) {
  
  const vecMap mass_a( Rcpp::as< vecMap >(mass_a_) );
  const vecMap mass_b( Rcpp::as< vecMap >(mass_b_) );
  const matMap cost_matrix(Rcpp::as< matMap> (cost_matrix_));
  const matMap cost_matrix_A(Rcpp::as< matMap> (cost_matrix_A_));
  const matMap cost_matrix_B(Rcpp::as< matMap> (cost_matrix_B_));
  
  std::string method(Rcpp::as<std::string>(method_(0)));
  
  int N = mass_a.size();
  int M = mass_b.size();
  
  matrixI idx(N*M,2);
  Eigen::VectorXd mass(N);
  
  transport_C(mass_a, mass_b, cost_matrix, idx, mass, method, 
              cost_matrix_A, cost_matrix_B,
              epsilon_, niter_, unbiased_, threads_);
  
  for(int i = 0; i < idx.size(); i++) idx(i) += 1;
  
  return  Rcpp::List::create(Rcpp::Named("from")  = Rcpp::wrap(idx.col(0)),
                             Rcpp::Named("to")    = Rcpp::wrap(idx.col(1)),
                             Rcpp::Named("mass")  = Rcpp::wrap(mass));
}


//[[Rcpp::export]]
Rcpp::List transport_(const Rcpp::NumericMatrix & A_, 
                      const Rcpp::NumericMatrix & B_, double p, double ground_p,
                      const Rcpp::CharacterVector & method_,
                      bool a_sort, double epsilon_ = 0.0, int niter_ = 0,
                      bool unbiased_ = false,
                      int threads_ = 1) {
  int N = A_.cols();
  int M = B_.cols();
  
  const matMap A(Rcpp::as<matMap>(A_));
  const matMap B(Rcpp::as<matMap>(B_));
  
  const std::string method(Rcpp::as<std::string>(method_(0)));
  
  matrixI idx(N * M, 2);
  Eigen::VectorXd mass(N * M);
  
  transport(A, B, p, ground_p,
            idx, mass, method, a_sort, epsilon_, niter_, 
            unbiased_,
            threads_);
  
  
  for(int i = 0; i < idx.size(); i++) idx(i) += 1;
  
  return  Rcpp::List::create(Rcpp::Named("from")  = Rcpp::wrap(idx.col(0)),
                             Rcpp::Named("to")    = Rcpp::wrap(idx.col(1)),
                             Rcpp::Named("mass")  = Rcpp::wrap(mass));
  
}

//[[Rcpp::export]]
Rcpp::List transport_swap_(const Rcpp::NumericMatrix & A_,
                      const Rcpp::NumericMatrix & B_,
                      Eigen::MatrixXi & idx_,
                      Eigen::VectorXd & mass_,
                      double p, double ground_p,
                      double tolerance_, int niter_ = 0) {
  int N = A_.cols();
  int M = B_.cols();

  const matMap A(Rcpp::as<matMap>(A_));
  const matMap B(Rcpp::as<matMap>(B_));


  matrixI idx = idx_;
  Eigen::VectorXd mass = mass_;
  trans_swap(A, B, N, M,
             idx, mass, ground_p,
             p, tolerance_, niter_);

  for(int i = 0; i < idx.size(); i++) idx(i) += 1;

  return  Rcpp::List::create(Rcpp::Named("from")  = Rcpp::wrap(idx.col(0)),
                             Rcpp::Named("to")    = Rcpp::wrap(idx.col(1)),
                             Rcpp::Named("mass")  = Rcpp::wrap(mass));

}

