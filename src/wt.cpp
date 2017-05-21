#include <Rcpp.h>

using namespace Rcpp;

// Calculate weighted mean
// [[Rcpp::export]]
double wtmRcpp(NumericVector x,
               NumericVector w){

  return sum(x * w) / sum(w);
}

// Calculate weighted variance
// [[Rcpp::export]]
double wtvRcpp(NumericVector x,
               NumericVector w){

  double xbar = wtmRcpp(x, w);
  return sum(w * pow(x - xbar, 2)) * (sum(w) / (pow(sum(w), 2) - sum(pow(w, 2))));
}

// Calculate lrm with weights
// [[Rcpp::export]]
NumericMatrix lrmRcppWt(NumericMatrix & X,
                        NumericMatrix & W){
  return X;
}

// Calculate vlr with weights
// [[Rcpp::export]]
NumericMatrix vlrRcppWt(NumericMatrix & X,
                        NumericMatrix & W){
  return X;
}

// Calculate lrm with or without weights
// [[Rcpp::export]]
NumericMatrix lrm(NumericMatrix & X,
                  NumericMatrix & W){
  return X;
}

// Calculate vlr with or without weights
// [[Rcpp::export]]
NumericMatrix vlr(NumericMatrix & X,
                  NumericMatrix & W){
  return X;
}
