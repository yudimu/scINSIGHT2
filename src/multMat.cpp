// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
NumericMatrix multMat(arma::mat m1, arma::mat m2) {
  arma::mat m = m1 * m2;
  return wrap(m);
}
