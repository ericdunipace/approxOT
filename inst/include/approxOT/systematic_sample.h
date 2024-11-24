#ifndef SYSTEMATIC_SAMPLE_H
#define SYSTEMATIC_SAMPLE_H

#include "../approxOT_types.h"

//' Samples from a multinomial systematically
//'
//' @param samps An Eigen::VectorXi giving the sampled indicators
//' @param v An Eigen::VectorXd of sample weights
//' @param nsamp An int denoting the number of samples to take
//' @return void
//' @keywords internal
 static inline void sample_systematic(vectorI & samps, const Eigen::VectorXd & weight, const int nsamp ) {
   Rcpp::RNGScope scope;
   
   Rcpp::NumericVector draw = Rcpp::runif(1);
   double u = draw(0)/double(nsamp);
   double sampWeight = weight(0);
   int i = 0;
   
   for ( int j = 0; j < nsamp; j++ ){
     while(sampWeight < u) {
       i++;
       sampWeight += weight(i);
     }
     samps(j) = i;
     u += 1.0/double(nsamp);
   }
   
 }
  
#endif //SYSTEMATIC_SAMPLE_H