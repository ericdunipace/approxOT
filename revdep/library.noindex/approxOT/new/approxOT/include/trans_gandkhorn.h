#ifndef TRANS_GRANDKHORN_H
#define TRANS_GRANDKHORN_H

#include "approxOT_types.h"
#include "utils.h"
#include "systematic_sample.h"

//' Generates approximate optimal transport plans 
//' using the gandkhorn algorithm
//'
//' @param mass_r A reference to an Eigen::VectorXd
//' giving the empirical mass in sample 1
//' @param mass_l A reference to an Eigen::VectorXd
//' giving the empirical mass in sample 2
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
                     double eta, double epsilon, int niterations);
#endif //TRANS_GRANDKHORN_H
