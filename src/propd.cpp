#include <Rcpp.h>

using namespace Rcpp;

// Function for Box-Cox psuedo-vlr
// [[Rcpp::export]]
NumericVector boxRcpp(NumericMatrix & X,
                      const double a){

  // Raise all of X to the a power
  for(int i = 0; i < X.nrow(); i++){
    for(int j = 0; j < X.ncol(); j++){
      X(i, j) = pow(X(i, j), a);
    }
  }

  // Sweep out column means
  for(int j = 0; j < X.ncol(); j++){
    X(_, j) = X(_, j) / mean(X(_, j));
  }

  // Output a half-matrix
  int nfeats = X.ncol();
  int llt = nfeats * (nfeats - 1) / 2;
  NumericVector result(llt);
  int counter = 0;

  // Calculate sum([i - j]^2)
  for(int i = 1; i < nfeats; i++){
    for(int j = 0; j < i; j++){
      result(counter) = sum(pow(X(_, i) - X(_, j), 2));
      counter += 1;
    }
  }

  result = result / (pow(a, 2) * (X.nrow() - 1));
  return result;
}

// Function to count joint zero frequency
// [[Rcpp::export]]
NumericVector ctzRcpp(NumericMatrix & X){

  int nfeats = X.ncol();
  int nsubjs = X.nrow();
  int llt = nfeats * (nfeats - 1) / 2;

  // Count zero frequency per feature
  Rcpp::NumericVector zeroes(nfeats);
  for(int i = 0; i < nfeats; i++){
    for(int j = 0; j < nsubjs; j++){
      if(X(j, i) == 0){
        zeroes(i) += 1;
      }
    }
  }

  // Count joint zero frequency
  Rcpp::NumericVector result(llt);
  int counter = 0;
  for(int i = 1; i < nfeats; i++){
    for(int j = 0; j < i; j++){
      result(counter) = zeroes(i) + zeroes(j);
      counter += 1;
    }
  }

  return result;
}

// Function to calculate log-ratio mean
// [[Rcpp::export]]
NumericVector lrmRcpp(NumericMatrix & X){

  int nfeats = X.ncol();
  int llt = nfeats * (nfeats - 1) / 2;
  Rcpp::NumericVector result(llt);
  int counter = 0;

  for(int i = 1; i < nfeats; i++){
    for(int j = 0; j < i; j++){
      result(counter) = mean(log(X(_, i) / X(_, j)));
      counter += 1;
    }
  }

  return result;
}
