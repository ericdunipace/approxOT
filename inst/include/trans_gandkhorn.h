#ifndef TRANS_GRANDKHORN_H
#define TRANS_GRANDKHORN_H

#include "approxOT_types.h"
#include "utils.h"
#include "systematic_sample.h"

void trans_gandkhorn(const refVecConst & mass_r, const refVecConst & mass_l, 
                     const matrix & exp_cost, 
                     matrix & B,
                     double eta, double epsilon, int niterations);
#endif //TRANS_GRANDKHORN_H