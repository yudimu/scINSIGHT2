// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <iostream>
#include <armadillo>
#include <tuple>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat update_U_arma(arma::mat U, arma::mat V, arma::mat ddiff,
                  arma::mat ddratio,
                  double alpha,
                  double sigma) {
  double th = 1e2;
  int n = U.n_rows;
  int m = U.n_cols;


  arma::uvec keep_rows = arma::sum(arma::abs(V), 1) < th;
  int nfull = arma::accu(keep_rows);

  ddiff.elem(arma::find_nonfinite(ddiff)).zeros();

  // Calculate count as the row sums of the binary matrix (abs(ddiff) < th)
  arma::vec count(n);
  arma::vec count1(n);
  for (int i = 0; i < n; i++) {
    count(i) = arma::accu(arma::abs(ddiff.row(i)) < th);
    count1(i) = arma::accu(arma::abs(ddratio.row(i)) < th);
  }

  // Calculate corrections as nfull divided by count
  arma::vec corrections = static_cast<double>(nfull) / count;

  // Calculate corrections1 as nfull divided by count1
  arma::vec corrections1 = static_cast<double>(nfull) / count1;


  arma::mat dU(n, m);
  arma::mat ddU(n, m);

  ddiff.elem(arma::find(arma::abs(ddiff) > th)).zeros();
  ddratio.elem(arma::find(arma::abs(ddratio) > th)).zeros();

  arma::mat ddiff_new(n, nfull);
  arma::mat ddratio_new(n, nfull);
  arma::mat V_new(nfull, m);

  // Initialize matrices
  int col_idx = 0;
  for (int i = 0; i < nfull; i++) {
    if (keep_rows(i)) {
      ddiff_new.col(col_idx) = ddiff.col(i);
      ddratio_new.col(col_idx) = ddratio.col(i);
      V_new.row(col_idx) = V.row(i);
      col_idx++;
    }
  }

  dU = - ddiff_new * V_new;
  dU = dU % repmat(corrections, 1, dU.n_cols);
  dU +=  U / (sigma * sigma);

  ddU = ddratio_new * (V_new % V_new) % repmat(corrections1, 1, ddU.n_cols) +  1 / (sigma * sigma);

  arma::mat Unew(n,m);
  Unew = U - alpha * (dU / ddU);

  return Unew;
}



// [[Rcpp::export]]
arma::mat update_par_arma(arma::mat U, arma::mat V, arma::mat ddiff,
                    arma::mat ddratio,
                    double alpha) {
  double th = 1e2;
  int n = U.n_rows;
  int m = U.n_cols;


  arma::uvec keep_rows = arma::sum(arma::abs(V), 1) < th;
  int nfull = arma::accu(keep_rows);
  ddiff.elem(arma::find_nonfinite(ddiff)).zeros();


  // Calculate count as the row sums of the binary matrix (abs(ddiff) < th)
  arma::vec count(n);
  arma::vec count1(n);
  for (int i = 0; i < n; i++) {
    count(i) = arma::accu(arma::abs(ddiff.row(i)) < th);
    count1(i) = arma::accu(arma::abs(ddratio.row(i)) < th);
  }

  // Calculate corrections as nfull divided by count
  arma::vec corrections = static_cast<double>(nfull) / count;

  // Calculate corrections1 as nfull divided by count1
  arma::vec corrections1 = static_cast<double>(nfull) / count1;


  arma::mat dU(n, m);
  arma::mat ddU(n, m);

  ddiff.elem(arma::find(arma::abs(ddiff) > th)).zeros();
  ddratio.elem(arma::find(arma::abs(ddratio) > th)).zeros();

  arma::mat ddiff_new(n, nfull);
  arma::mat ddratio_new(n, nfull);
  arma::mat V_new(nfull, m);

  // Initialize matrices
  int col_idx = 0;
  for (int i = 0; i < nfull; i++) {
    if (keep_rows(i)) {
      ddiff_new.col(col_idx) = ddiff.col(i);
      ddratio_new.col(col_idx) = ddratio.col(i);
      V_new.row(col_idx) = V.row(i);
      col_idx++;
    }
  }

  dU = - ddiff_new * V_new;
  dU = dU % repmat(corrections, 1, dU.n_cols);

  ddU = ddratio_new * (V_new % V_new) % repmat(corrections1, 1, ddU.n_cols) ;

  arma::mat Unew(n,m);
  Unew = U - alpha * (dU / ddU);

  return Unew;
}


// [[Rcpp::export]]
double matrix_deviance(arma::mat pred, arma::mat obs){


  if (pred.n_elem == 1) {

    // Create a predMatrix with the same dimensions as obs
    arma::mat predMatrix = obs;

    // Fill predMatrix with the same value as pred
    predMatrix.fill(pred(0));

    // Update pred to be predMatrix
    pred = predMatrix;
  }

  // Create a boolean matrix to identify NA values in obs
  arma::umat isna = (obs == arma::datum::nan);

  // Filter pred and obs to exclude NA values
  pred = pred.elem(arma::find(isna == 0)); // Remove rows with NA in obs
  obs = obs.elem(arma::find(isna == 0));   // Remove rows with NA in obs

  arma::vec predVector = arma::vectorise(pred);
  arma::vec obsVector = arma::vectorise(obs);

  arma::vec r = predVector;

  // Find indices where y > 0
  arma::uvec p = arma::find(obsVector > 0);

  // Calculate the right-hand side of the assignment for positive y values
  arma::vec rhs = (obsVector % arma::log(obsVector / predVector) - (obsVector - predVector));

  // Update r with the calculated values for positive y values
  r.elem(p) = rhs.elem(p);

  // Calculate 2 * r
  r = 2.0 * r;

  double a = arma::mean(arma::nonzeros(r));
  return a;

}







// [[Rcpp::export]]
List iteration(arma::mat U, arma::mat V, arma::mat Y, arma::mat beta, arma::mat X, double alpha, double sigma,
               arma::mat latentpart, arma::mat regpart, arma::vec logs, int maxIter, double d,
               int n, int m, double llast, double tol1, double tol2, double p){
  double l;
  int i;
  arma::vec dev(maxIter);
  arma::vec squaredU(maxIter);
  arma::vec meansquaredU(maxIter);
  arma::vec likelihood(maxIter);
  for (i = 1; i < maxIter; ++i) {
    if (i > 1) {
      // Update Y based on the condition
      arma::umat isna = (Y != Y);
      arma::uvec isnaIndices = arma::find(isna == 1); // Find indices where isna is true
      arma::mat expValues = arma::exp(latentpart + regpart + repmat(logs, 1, latentpart.n_cols));

      // Update elements of Y using the expValues where isna is true
      Y.elem(isnaIndices) = expValues.elem(isnaIndices);
    }

    arma::mat eta = latentpart + regpart + repmat(logs, 1, latentpart.n_cols);
    arma::mat M = arma::exp(eta);
    arma::mat A = arma::exp(eta);
    arma::mat dratio = A / M;
    dratio.elem(arma::find_nonfinite(dratio)).zeros();

    arma::mat ddiff = (Y - M) % dratio;
    ddiff.elem(arma::find_nonfinite(ddiff)).zeros();

    arma::mat ddratio = dratio % A;
    ddratio.elem(arma::find_nonfinite(ddratio)).zeros();


    ///// update ///
    arma::mat Uold = U;
    arma::mat Unew = update_U_arma(U,V,ddiff,ddratio,alpha,sigma);
    V = update_par_arma(V,U,ddiff.t(),ddratio.t(),alpha);

    U = Unew;
    U.elem(arma::find_nonfinite(U)).zeros();
    V.elem(arma::find_nonfinite(V)).zeros();

    if (d) {
      //Update 'beta' using your custom function
      beta = update_par_arma(beta.t(), X, ddiff.t(), ddratio.t(), alpha);
      beta = beta.t();
    }

    // Replace NA values in 'beta' with 0
    beta.elem(arma::find_nonfinite(beta)).zeros();


    double dsigma = (- trace(U * U.t())/ (sigma * sigma * sigma)) + (n*p / sigma);
    // Calculate ddsigma = 3*crossprod(as.vector(U), as.vector(U))/sigma^4 - n/sigma^2
    double ddsigma = (3 * trace(U * U.t()) / (sigma * sigma * sigma * sigma)) - (n*p/ (sigma * sigma));
    // Update sigma = sigma - alpha*dsigma/ddsigma
    sigma = sigma - (alpha * dsigma / ddsigma);


    if (d) {
      regpart = X * beta;
    }
    // Calculate latentpart = U * V.t()
    latentpart = U * V.t();
    arma::mat predicted = arma::exp(latentpart + regpart + repmat(logs, 1, latentpart.n_cols));

    double penalties = arma::norm(U, "fro") / (sigma * sigma);

    /// matrix deviance
    double matrix_dev = matrix_deviance(predicted, Y);

    // Calculate the product of the dimensions
    l = (matrix_dev + penalties) / (n*m);

    dev(i) = l;


    squaredU(i) = accu(square(U - Uold));
    meansquaredU(i) = mean(mean(square(U - Uold)));

    arma::mat eta_predicted = latentpart + regpart + repmat(logs, 1, latentpart.n_cols);
    likelihood(i) = accu(Y % eta_predicted - predicted) - trace(trans(U) * U) / 2 / pow(sigma, 2) - n * p * log(sigma);
    //
    //
    //if (i % 100 == 0) {
      //Rcpp::Rcout << "iteration " << i << ", PQL = " << l << std::endl;
    //}

    if (i > 100) {
      if (sqrt(meansquaredU(i)) < tol1 && mean(meansquaredU.subvec(i - 100, i)) < tol2){
        break;
      }
    }

    // if (std::abs(llast - l) / l < tol) {
    //   break;
    // } else {
    //   llast = l;
    // }

  }

  List results;
  results["U"] = U;
  results["V"] = V;
  results["squaredU"]  = squaredU;
  results["meansquaredU"] = meansquaredU;
  results["likelihood"] = likelihood;
  results["beta"] = beta;
  results["sigma"] = sigma;
  results["latentpart"] = latentpart;
  results["regpart"] = regpart;
  results["logs"] = logs;
  results["i"] = i;
  results["l"] = l;
  results["dev"] = dev;

  return results;
}



