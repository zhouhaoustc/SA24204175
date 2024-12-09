// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// mse_fun
List mse_fun(NumericVector beta, NumericVector est, Nullable<NumericMatrix> X_test);
RcppExport SEXP _SA24204175_mse_fun(SEXP betaSEXP, SEXP estSEXP, SEXP X_testSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type est(estSEXP);
    Rcpp::traits::input_parameter< Nullable<NumericMatrix> >::type X_test(X_testSEXP);
    rcpp_result_gen = Rcpp::wrap(mse_fun(beta, est, X_test));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SA24204175_mse_fun", (DL_FUNC) &_SA24204175_mse_fun, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_SA24204175(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
