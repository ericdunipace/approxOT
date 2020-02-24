#ifndef TRANS_SWAP_H
#define TRANS_SWAP_H

#include "approxOT_types.h"
#include "cost.h"

void trans_swap(const matrix & A, const matrix & B, int N, int M,
                   matrixI & idx, vector &  mass, double ground_p,
                   double p, double tol, int niter);
#endif //TRANS_SWAP_H