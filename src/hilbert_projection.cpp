#include "../inst/include/approxOT_types.h"
#include "../inst/include/hilbert_cgal.h"

// [[Rcpp::export]]
Rcpp::IntegerVector hilbert_proj_(const matrix & A)
{
  int K = A.rows();
  int N = A.cols();
  std::vector<int> idx(N);
  
  hilbert_sort_cgal_fun(A.data(), K, N, &idx[0] );
  return(Rcpp::wrap(idx));
}