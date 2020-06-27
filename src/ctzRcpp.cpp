#include <Rcpp.h>
#include <stdint.h>
using namespace Rcpp;

// Function to count joint zero frequency
// [[Rcpp::export]]
NumericVector ctzRcpp(NumericMatrix & X){

  int32_t nfeats = X.ncol();
  int32_t nsubjs = X.nrow();
  int32_t llt = nfeats * (nfeats - 1) / 2;

  // Count zero frequency per feature
  Rcpp::NumericVector zeroes(nfeats);
  for(int32_t i = 0; i < nfeats; i++){
    for(int32_t j = 0; j < nsubjs; j++){
      if(X(j, i) == 0){
        zeroes(i) += 1;
      }
    }
  }

  // Count joint zero frequency
  Rcpp::NumericVector result(llt);
  int32_t counter = 0;
  for(int32_t i = 1; i < nfeats; i++){
    for(int32_t j = 0; j < i; j++){
      result(counter) = zeroes(i) + zeroes(j);
      counter += 1;
    }
  }

  return result;
}
