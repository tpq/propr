#include <Rcpp.h>
using namespace Rcpp;

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
