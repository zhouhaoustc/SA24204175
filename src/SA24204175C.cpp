#include <Rcpp.h>
using namespace Rcpp;

//' @title computing the MSE
//' @description computing the MSE
//' @param beta coefficient vector
//' @param est estimator
//' @param X_test test dataset
//' @export
// [[Rcpp::export]]
List mse_fun(NumericVector beta, NumericVector est, Nullable<NumericMatrix> X_test = R_NilValue) {
  // Calculate estimation error
  double est_err = 0.0;
  for (int i = 0; i < beta.size(); ++i) {
    double diff = beta[i] - est[i];
    est_err += diff * diff;
  }

  // Initialize prediction error to NA
  double pred_err = NA_REAL;

  // If X.test is provided, calculate prediction error
  if (X_test.isNotNull()) {
    NumericMatrix x_test = as<NumericMatrix>(X_test);

    // Calculate (X.test %*% (beta - est))
    NumericVector diff_vec(beta.size());
    for (int i = 0; i < beta.size(); ++i) {
      diff_vec[i] = beta[i] - est[i];
    }

    NumericVector pred_vec = x_test * diff_vec;

    // Calculate mean squared prediction error
    pred_err = mean(pow(pred_vec, 2));
  }

  // Return results
  return List::create(
    Named("est.err") = est_err,
    Named("pred.err") = pred_err
  );
}
