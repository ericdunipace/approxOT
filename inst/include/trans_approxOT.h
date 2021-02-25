#ifndef TRANS_APPROX_OT_H
#define TRANS_APPROX_OT_H

#include "approxOT_types.h"
#include "round_feasible.h"
#include "trans_sinkhorn.h"
#include "trans_greenkhorn.h"
#include "trans_randkhorn.h"
#include "trans_gandkhorn.h"
#include "utils.h"

void trans_approxOT(const refVecConst & mass_a, const refVecConst & mass_b, 
                    refMat cost_matrix, 
                    matrix & assign_mat,
                    double epsilon, int niterations,
                    bool unbiased,
                    const std::string & method,
                    refMat cost_matrix_A, refMat cost_matrix_B);
#endif //TRANS_SINKHORN_H