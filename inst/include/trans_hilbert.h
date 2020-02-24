#ifndef TRANS_HILBERT_H
#define TRANS_HILBERT_H

#include "approxOT_types.h"
#include "hilbert_cgal.h"

void trans_hilbert(const matrix & A, const matrix & B, int N, int M,
                   matrixI & idx, vector &  mass, bool & a_sort) ;
#endif //TRANS_HILBERT_H