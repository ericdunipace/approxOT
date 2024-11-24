#ifndef _UTILS_H
#define _UTILS_H

#include <Rcpp.h>
#include <RcppEigen.h>
#include "../approxOT_types.h"

static inline double median(refMat A) {
  if ( A.size() == 0) {
    Rcpp::stop("Can't take the median of an empty matrix.");
  }
  int size = A.size();
  int middleIdx = size/2;
  matrix A_copy = A;
  double * begin = A_copy.data();
  double * end = begin + size;
  double * target = begin + middleIdx;
  std::nth_element( begin, target, end);
  double a = *target;
  //
  if (size % 2 != 0) { //Odd number of elements
    return (a);
  } else {            //Even number of elements
    double * targetNeighbor = target-1;
    std::nth_element(begin, targetNeighbor, end);
    double an = *targetNeighbor;
    return ( (a + an)/2.0);
  }
  // return(0.0);
}

static inline double median(const matrix & A) {
  if ( A.size() == 0) {
    Rcpp::stop("Can't take the median of an empty matrix.");
  }
  int size = A.size();
  int middleIdx = size/2;
  matrix A_copy = A;
  double * begin = A_copy.data();
  double * end = begin + size;
  double * target = begin + middleIdx;
  std::nth_element( begin, target, end);
  double a = *target;
  //
  if (size % 2 != 0) { //Odd number of elements
    return (a);
  } else {            //Even number of elements
    double * targetNeighbor = target-1;
    std::nth_element(begin, targetNeighbor, end);
    double an = *targetNeighbor;
    return ( (a + an)/2.0);
  }
  // return(0.0);
}

static inline double median( matrix & A) {
  if ( A.size() == 0) {
    Rcpp::stop("Can't take the median of an empty matrix.");
  }
  int size = A.size();
  int middleIdx = size/2;
  matrix A_copy = A;
  double * begin = A_copy.data();
  double * end = begin + size;
  double * target = begin + middleIdx;
  std::nth_element( begin, target, end);
  double a = *target;
  //
  if (size % 2 != 0) { //Odd number of elements
    return (a);
  } else {            //Even number of elements
    double * targetNeighbor = target-1;
    std::nth_element(begin, targetNeighbor, end);
    double an = *targetNeighbor;
    return ( (a + an)/2.0);
  }
  // return(0.0);
}

#endif