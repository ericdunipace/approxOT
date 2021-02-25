#ifndef TRANSPORT_H
#define TRANSPORT_H

#include "approxOT_types.h"
#include "cost.h"
#include <string>
#include "utils.h"
#include "trans_hilbert.h"
#include "networkflow.h"
#include "trans_shortsimplex.h"
#include "trans_rank.h"
#include "trans_univariate.h"
#include "trans_univariate_approx_pwr.h"
#include "systematic_sample.h"
#include "trans_approxOT.h"
#include "trans_swap.h"

void transport_C(const refVecConst & mass_a, const refVecConst & mass_b, 
                 refMat cost_matrix, matrixI & idx, vector & mass, const std::string & method,
                 refMat cost_matrix_A, refMat cost_matrix_B,
                 double epsilon = 0.0, int niter = 0,
                 bool unbiased = false, int threads = 1);
void transport(const matrix & A, const matrix & B, const double p, const double ground_p,
               matrixI & idx, vector & mass, const std::string & method, bool & a_sort,
               double epsilon = 0.0, int niter = 0,
               bool unbiased = false, int threads = 1);
  
#endif //TRANSPORT_H