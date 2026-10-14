# approxOT

## Changes in Version 1.3
### Major changes
* Added OpenMP support for the Sinkhorn algorithm. The log-scale Sinkhorn updates were also rewritten to avoid N×M temporary matrices, giving roughly a 1.4× speedup even on a single thread.
* By default, approxOT uses half of the available processors (never more than `OMP_THREAD_LIMIT`, and at most 2 during `R CMD check`). Use `approxOT_threads(n)` to change this, `approxOT_threads()` to see the current value, or set `options(approxOT.threads = n)` before loading the package.
* Packages linking to the approxOT headers run single-threaded unless they call `approxOT_set_threads(n)` from C++.

### Minor Improvements and Bug fixes
* The `shortsimplex.h` header caused a false positive warning on gcc-asan due to the posibility of feeding negative integer to the function; however, this would never happen due to higher level code always passing a positive integer since it is set by the size of the corresponding mass (i.e., an empirical measure passed to the function always has length greater than 0). A guard has been placed to prevent this gcc-asan warning.

## Changes in Version 1.2
### Major changes
* C++ code is now header-only; headers have moved to `inst/include/approxOT/`. Packages linking to approxOT should update their include paths.
* Updated for compatibility with CGAL 6.0.
* Now requires C++17.

### Bug fixes
* `transport_plan()` now correctly handles `method = "exact"` (mapped to `"networkflow"`).
* Fixed argument matching for `observation.orientation` in `wasserstein_calc_cost()`.

### Documentation
* Clarified that `epsilon` is the entropy penalty itself (must be > 0), not a multiple of the median cost.

## Changes in Version 1.1.1

### Minor Improvements and Bug fixes
* Minor bug fix to update a few Eigen C++ functions (Thank you Andrew Johnson!) 
* Removing older functions that are no longer used to speed up build
* No changes to user interface

## Changes in Version 1.1

* Added a `NEWS.md` file to track changes to the package.
* Added a `as.matrix` function to convert transportation plans to matrices.
* Added a `as.transport.plan` to convert from transportation matrix to transportation plan.
* Added classes to support these two functions.
* Added support for log domain calculations for the sinkhorn method. This can be accessed with "sinkhorn_log" as the method.
* Fixed typo in C files incorrectly referencing Rcpp::stop
* updated DLL registrations in Namespace
