#ifndef TRANS_RANDKHORN_H
#define TRANS_RANDKHORN_H

#include "approxOT_types.h"
#include "utils.h"

void trans_randkhorn(const refVecConst & mass_r, const refVecConst & mass_l, 
                     const matrix & exp_cost, 
                     matrix & B,
                     double eta, double epsilon, int niterations);
#endif //TRANS_RANDKHORN_H