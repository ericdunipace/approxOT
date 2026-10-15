#include "approxOT/threads.h"

//[[Rcpp::export]]
int set_threads_(int n) {
  return approxOT_set_threads(n);
}

//[[Rcpp::export]]
int get_threads_() {
  return approxOT_nthreads();
}
