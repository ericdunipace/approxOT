#include "approxOT/round_feasible.h"

//[[Rcpp::export]]
matrix round_2_feasible_(matrix & F, const Eigen::VectorXd & mass_a, const Eigen::VectorXd & mass_b) {
  matrix F_safe = F;
  round_feasible(F_safe, mass_a, mass_b);
  return(F_safe);
}
