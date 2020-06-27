#include <Rcpp.h>
#include <stdint.h>
using namespace Rcpp;

// Calculate lrv weight modifier
// [[Rcpp::export]]
NumericVector omega(NumericMatrix & W){

  int32_t nfeats = W.ncol();
  int32_t llt = nfeats * (nfeats - 1) / 2;
  Rcpp::NumericVector result(llt);
  Rcpp::NumericVector Wij(nfeats);
  int32_t counter = 0;
  double n = 0;

  for(int32_t i = 1; i < nfeats; i++){
    for(int32_t j = 0; j < i; j++){
      Wij = W(_, i) * W(_, j);
      n = sum(Wij);
      result(counter) = n - sum(pow(Wij, 2)) / n;
      counter += 1;
    }
  }

  return result;
}

// Calculate lrv weight modifier (population-level for F-stat and F-mod)
// [[Rcpp::export]]
NumericVector Omega(NumericMatrix & W){

  int32_t nfeats = W.ncol();
  int32_t llt = nfeats * (nfeats - 1) / 2;
  Rcpp::NumericVector result(llt);
  Rcpp::NumericVector Wij(nfeats);
  int32_t counter = 0;

  for(int32_t i = 1; i < nfeats; i++){
    for(int32_t j = 0; j < i; j++){
      Wij = W(_, i) * W(_, j);
      result(counter) = sum(Wij);
      counter += 1;
    }
  }

  return result;
}
