#include <Rcpp.h>
using namespace Rcpp;

// Calculate lrv weight modifier
// [[Rcpp::export]]
NumericVector omega(NumericMatrix & X,
                    NumericMatrix & W){

  int nfeats = X.ncol();
  int llt = nfeats * (nfeats - 1) / 2;
  Rcpp::NumericVector result(llt);
  Rcpp::NumericVector Wij(nfeats);
  int counter = 0;
  double n = 0;

  for(int i = 1; i < nfeats; i++){
    for(int j = 0; j < i; j++){
      Wij = W(_, i) * W(_, j);
      n = sum(Wij);
      result(counter) = n - sum(pow(Wij, 2)) / n;
      counter += 1;
    }
  }

  return result;
}
