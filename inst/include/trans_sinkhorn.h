#ifndef TRANS_SINKHORN_H
#define TRANS_SINKHORN_H

#include "approxOT_types.h"
#include "utils.h"

void trans_sinkhorn(const refVecConst & mass_a, const refVecConst & mass_b, 
                    const matrix & exp_cost, 
                    matrix & A,
                    double eta, double epsilon, int niterations);
void trans_sinkhorn_self(vector & u, const refVecConst & mass_a,
                         const matrix & exp_cost,
                         double epsilon, int niterations);
void trans_sinkhorn_log(const refVecConst & mass_a, const refVecConst & mass_b,
                        const matrix & cost,
                        matrix & Assign,
                        double eta, double epsilon, int niterations,
                        const refVecConst & f_pot, const refVecConst & g_pot);
void trans_sinkhorn_autocorr(vector & f, const refVecConst & mass_a,
                             const matrix & exp_cost,
                             double eta, double epsilon, int niterations);
#endif //TRANS_SINKHORN_H