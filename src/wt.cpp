#include <Rcpp.h>

using namespace Rcpp;

// Calculate weighted mean
// [[Rcpp::export]]
double wtmRcpp(NumericVector x,
               NumericVector w){

  return sum(x * w) / sum(w);
}

// Calculate weighted var
// [[Rcpp::export]]
double wtvRcpp(NumericVector x,
               NumericVector w){

  double xbar = wtmRcpp(x, w);
  return sum(w * pow(x - xbar, 2)) * (sum(w) / (pow(sum(w), 2) - sum(pow(w, 2))));
}

// Calculate lrm with weights
// [[Rcpp::export]]
NumericVector lrm(NumericMatrix & X,
                  NumericMatrix & W,
                  bool weighted = false){

  int nfeats = X.ncol();
  int llt = nfeats * (nfeats - 1) / 2;
  Rcpp::NumericVector result(llt);
  Rcpp::NumericVector Wij(nfeats);
  int counter = 0;

  if(weighted){
    for(int i = 1; i < nfeats; i++){
      for(int j = 0; j < i; j++){
        Wij = W(_, i) * W(_, j);
        result(counter) = wtmRcpp(log(X(_, i) / X(_, j)), Wij);
        counter += 1;
      }
    }
  }else{
    for(int i = 1; i < nfeats; i++){
      for(int j = 0; j < i; j++){
        result(counter) = mean(log(X(_, i) / X(_, j)));
        counter += 1;
      }
    }
  }

  return result;
}

// Calculate lrv with weights
// [[Rcpp::export]]
NumericVector lrv(NumericMatrix & X,
                  NumericMatrix & W,
                  bool weighted = false){

  int nfeats = X.ncol();
  int llt = nfeats * (nfeats - 1) / 2;
  Rcpp::NumericVector result(llt);
  Rcpp::NumericVector Wij(nfeats);
  int counter = 0;

  if(weighted){
    for(int i = 1; i < nfeats; i++){
      for(int j = 0; j < i; j++){
        Wij = W(_, i) * W(_, j);
        result(counter) = wtvRcpp(log(X(_, i) / X(_, j)), Wij);
        counter += 1;
      }
    }
  }else{
    for(int i = 1; i < nfeats; i++){
      for(int j = 0; j < i; j++){
        result(counter) = var(log(X(_, i) / X(_, j)));
        counter += 1;
      }
    }
  }

  return result;
}

// Calculate lrv modifier
// [[Rcpp::export]]
NumericVector lrvMod(NumericMatrix & X,
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
