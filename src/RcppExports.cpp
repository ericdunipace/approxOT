// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/approxOT_types.h"
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// cost_calculation_
Rcpp::NumericMatrix cost_calculation_(const Rcpp::NumericMatrix& A_, const Rcpp::NumericMatrix& B_, const double p);
RcppExport SEXP _approxOT_cost_calculation_(SEXP A_SEXP, SEXP B_SEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type A_(A_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type B_(B_SEXP);
    Rcpp::traits::input_parameter< const double >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(cost_calculation_(A_, B_, p));
    return rcpp_result_gen;
END_RCPP
}
// multi_marg_final_cost_
double multi_marg_final_cost_(const Rcpp::List& idx_, const Rcpp::List& data_, const Rcpp::NumericVector& mass_, int M, int D, double p, double ground_p);
RcppExport SEXP _approxOT_multi_marg_final_cost_(SEXP idx_SEXP, SEXP data_SEXP, SEXP mass_SEXP, SEXP MSEXP, SEXP DSEXP, SEXP pSEXP, SEXP ground_pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type idx_(idx_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type data_(data_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type mass_(mass_SEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type D(DSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type ground_p(ground_pSEXP);
    rcpp_result_gen = Rcpp::wrap(multi_marg_final_cost_(idx_, data_, mass_, M, D, p, ground_p));
    return rcpp_result_gen;
END_RCPP
}
// multi_marg_given_dist_
double multi_marg_given_dist_(const Rcpp::List& idx_, const Rcpp::NumericVector& mass_, const Rcpp::NumericVector& cost_, int M, int N_cost, double p);
RcppExport SEXP _approxOT_multi_marg_given_dist_(SEXP idx_SEXP, SEXP mass_SEXP, SEXP cost_SEXP, SEXP MSEXP, SEXP N_costSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type idx_(idx_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type mass_(mass_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type cost_(cost_SEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type N_cost(N_costSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(multi_marg_given_dist_(idx_, mass_, cost_, M, N_cost, p));
    return rcpp_result_gen;
END_RCPP
}
// hilbert_proj_
Rcpp::IntegerVector hilbert_proj_(const matrix& A);
RcppExport SEXP _approxOT_hilbert_proj_(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const matrix& >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(hilbert_proj_(A));
    return rcpp_result_gen;
END_RCPP
}
// sinkhorn_
Rcpp::List sinkhorn_(Rcpp::NumericVector p_, Rcpp::NumericVector q_, Rcpp::NumericMatrix cost_matrix_, double epsilon, int niterations);
RcppExport SEXP _approxOT_sinkhorn_(SEXP p_SEXP, SEXP q_SEXP, SEXP cost_matrix_SEXP, SEXP epsilonSEXP, SEXP niterationsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type p_(p_SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type q_(q_SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type cost_matrix_(cost_matrix_SEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< int >::type niterations(niterationsSEXP);
    rcpp_result_gen = Rcpp::wrap(sinkhorn_(p_, q_, cost_matrix_, epsilon, niterations));
    return rcpp_result_gen;
END_RCPP
}
// transport_C_
Rcpp::List transport_C_(const Rcpp::NumericVector& mass_a_, const Rcpp::NumericVector& mass_b_, const Rcpp::NumericMatrix& cost_matrix_, const Rcpp::CharacterVector& method_, double epsilon_, int niter_, int threads_);
RcppExport SEXP _approxOT_transport_C_(SEXP mass_a_SEXP, SEXP mass_b_SEXP, SEXP cost_matrix_SEXP, SEXP method_SEXP, SEXP epsilon_SEXP, SEXP niter_SEXP, SEXP threads_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type mass_a_(mass_a_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type mass_b_(mass_b_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type cost_matrix_(cost_matrix_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type method_(method_SEXP);
    Rcpp::traits::input_parameter< double >::type epsilon_(epsilon_SEXP);
    Rcpp::traits::input_parameter< int >::type niter_(niter_SEXP);
    Rcpp::traits::input_parameter< int >::type threads_(threads_SEXP);
    rcpp_result_gen = Rcpp::wrap(transport_C_(mass_a_, mass_b_, cost_matrix_, method_, epsilon_, niter_, threads_));
    return rcpp_result_gen;
END_RCPP
}
// transport_
Rcpp::List transport_(const Rcpp::NumericMatrix& A_, const Rcpp::NumericMatrix& B_, double p, double ground_p, const Rcpp::CharacterVector& method_, bool a_sort, double epsilon_, int niter_, int threads_);
RcppExport SEXP _approxOT_transport_(SEXP A_SEXP, SEXP B_SEXP, SEXP pSEXP, SEXP ground_pSEXP, SEXP method_SEXP, SEXP a_sortSEXP, SEXP epsilon_SEXP, SEXP niter_SEXP, SEXP threads_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type A_(A_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type B_(B_SEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type ground_p(ground_pSEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type method_(method_SEXP);
    Rcpp::traits::input_parameter< bool >::type a_sort(a_sortSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon_(epsilon_SEXP);
    Rcpp::traits::input_parameter< int >::type niter_(niter_SEXP);
    Rcpp::traits::input_parameter< int >::type threads_(threads_SEXP);
    rcpp_result_gen = Rcpp::wrap(transport_(A_, B_, p, ground_p, method_, a_sort, epsilon_, niter_, threads_));
    return rcpp_result_gen;
END_RCPP
}
// transport_swap_
Rcpp::List transport_swap_(const Rcpp::NumericMatrix& A_, const Rcpp::NumericMatrix& B_, matrixI& idx_, vector& mass_, double p, double ground_p, double tolerance_, int niter_);
RcppExport SEXP _approxOT_transport_swap_(SEXP A_SEXP, SEXP B_SEXP, SEXP idx_SEXP, SEXP mass_SEXP, SEXP pSEXP, SEXP ground_pSEXP, SEXP tolerance_SEXP, SEXP niter_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type A_(A_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type B_(B_SEXP);
    Rcpp::traits::input_parameter< matrixI& >::type idx_(idx_SEXP);
    Rcpp::traits::input_parameter< vector& >::type mass_(mass_SEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type ground_p(ground_pSEXP);
    Rcpp::traits::input_parameter< double >::type tolerance_(tolerance_SEXP);
    Rcpp::traits::input_parameter< int >::type niter_(niter_SEXP);
    rcpp_result_gen = Rcpp::wrap(transport_swap_(A_, B_, idx_, mass_, p, ground_p, tolerance_, niter_));
    return rcpp_result_gen;
END_RCPP
}
// wasserstein_
double wasserstein_(const Rcpp::NumericVector& mass_, const Rcpp::NumericMatrix& cost_, const double p, const Rcpp::IntegerVector& from_, const Rcpp::IntegerVector& to_);
RcppExport SEXP _approxOT_wasserstein_(SEXP mass_SEXP, SEXP cost_SEXP, SEXP pSEXP, SEXP from_SEXP, SEXP to_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type mass_(mass_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type cost_(cost_SEXP);
    Rcpp::traits::input_parameter< const double >::type p(pSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type from_(from_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type to_(to_SEXP);
    rcpp_result_gen = Rcpp::wrap(wasserstein_(mass_, cost_, p, from_, to_));
    return rcpp_result_gen;
END_RCPP
}
// wasserstein_p_iid_
double wasserstein_p_iid_(const SEXP& X_, const SEXP& Y_, double p);
RcppExport SEXP _approxOT_wasserstein_p_iid_(SEXP X_SEXP, SEXP Y_SEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const SEXP& >::type X_(X_SEXP);
    Rcpp::traits::input_parameter< const SEXP& >::type Y_(Y_SEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(wasserstein_p_iid_(X_, Y_, p));
    return rcpp_result_gen;
END_RCPP
}
// wasserstein_p_iid_p_
double wasserstein_p_iid_p_(const SEXP& X_, const SEXP& Y_, double p);
RcppExport SEXP _approxOT_wasserstein_p_iid_p_(SEXP X_SEXP, SEXP Y_SEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const SEXP& >::type X_(X_SEXP);
    Rcpp::traits::input_parameter< const SEXP& >::type Y_(Y_SEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(wasserstein_p_iid_p_(X_, Y_, p));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_approxOT_cost_calculation_", (DL_FUNC) &_approxOT_cost_calculation_, 3},
    {"_approxOT_multi_marg_final_cost_", (DL_FUNC) &_approxOT_multi_marg_final_cost_, 7},
    {"_approxOT_multi_marg_given_dist_", (DL_FUNC) &_approxOT_multi_marg_given_dist_, 6},
    {"_approxOT_hilbert_proj_", (DL_FUNC) &_approxOT_hilbert_proj_, 1},
    {"_approxOT_sinkhorn_", (DL_FUNC) &_approxOT_sinkhorn_, 5},
    {"_approxOT_transport_C_", (DL_FUNC) &_approxOT_transport_C_, 7},
    {"_approxOT_transport_", (DL_FUNC) &_approxOT_transport_, 9},
    {"_approxOT_transport_swap_", (DL_FUNC) &_approxOT_transport_swap_, 8},
    {"_approxOT_wasserstein_", (DL_FUNC) &_approxOT_wasserstein_, 5},
    {"_approxOT_wasserstein_p_iid_", (DL_FUNC) &_approxOT_wasserstein_p_iid_, 3},
    {"_approxOT_wasserstein_p_iid_p_", (DL_FUNC) &_approxOT_wasserstein_p_iid_p_, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_approxOT(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
