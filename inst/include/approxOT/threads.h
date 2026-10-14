#ifndef APPROXOT_THREADS_H
#define APPROXOT_THREADS_H

#include <algorithm>
#include <RcppEigen.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// Number of OpenMP threads used by approxOT's parallel loops.
// Not `static`: an inline function's local static is shared by every
// translation unit, so setting it once affects all parallel regions.
inline int & approxOT_nthreads() {
  static int n = 1;
  return n;
}

// Smallest problem (number of cost matrix entries) split across threads.
// Below this, starting the threads costs about as much as the work saved.
// Small problems then never start libgomp's thread pool, whose idle
// threads valgrind reports as "possibly lost" when R exits.
constexpr long long approxOT_min_par_size = 1LL << 14;

// n <= 0 uses the default of half the available processors.
// Never exceeds the processor count or OMP_THREAD_LIMIT, which avoids
// OpenMP runtime warnings (e.g., CRAN machines set OMP_THREAD_LIMIT=2).
// Not capped by omp_get_max_threads(): omp_set_num_threads() calls
// (e.g., in trans_networkflow) lower it for the whole session.
inline int approxOT_set_threads(int n) {
#ifdef _OPENMP
  if (n <= 0) n = omp_get_num_procs() / 2;
  int lim = std::min(omp_get_num_procs(), omp_get_thread_limit());
  approxOT_nthreads() = std::max(1, std::min(n, lim));
#else
  approxOT_nthreads() = 1;
#endif
  Eigen::setNbThreads(approxOT_nthreads()); // Eigen's own matrix-product threads
  return approxOT_nthreads();
}

#endif //APPROXOT_THREADS_H
