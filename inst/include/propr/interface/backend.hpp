#pragma once

#include <Rcpp.h>

namespace propr {
    Rcpp::NumericMatrix covRcpp(Rcpp::NumericMatrix & X, const int norm_type);
    Rcpp::NumericMatrix corRcpp(Rcpp::NumericMatrix X);
    double              wtmRcpp(Rcpp::NumericVector x, Rcpp::NumericVector w);
    double              wtvRcpp(Rcpp::NumericVector x, Rcpp::NumericVector w);
    Rcpp::NumericMatrix covRcpp(Rcpp::NumericMatrix X, int norm_type);
    Rcpp::NumericMatrix vlrRcpp(Rcpp::NumericMatrix X);
    Rcpp::NumericMatrix clrRcpp(Rcpp::NumericMatrix X);
    Rcpp::NumericMatrix alrRcpp(Rcpp::NumericMatrix X, int ivar);
    Rcpp::NumericMatrix symRcpp(Rcpp::NumericMatrix X);
    Rcpp::NumericMatrix phiRcpp(Rcpp::NumericMatrix X, bool sym);
    Rcpp::NumericMatrix rhoRcpp(Rcpp::NumericMatrix X, Rcpp::NumericMatrix lr, int ivar);
    std::vector<int> indexPairs(Rcpp::NumericMatrix X, Rcpp::String op, double ref);
    Rcpp::List indexToCoord(Rcpp::IntegerVector V, int N);
    Rcpp::IntegerVector coordToIndex(Rcpp::IntegerVector row, Rcpp::IntegerVector col, int N);
    Rcpp::NumericMatrix linRcpp(Rcpp::NumericMatrix rho, Rcpp::NumericMatrix lr);
    Rcpp::NumericVector lltRcpp(Rcpp::NumericMatrix X);
    Rcpp::NumericVector urtRcpp(Rcpp::NumericMatrix X);
    Rcpp::List labRcpp(int nfeats);
    Rcpp::NumericMatrix half2mat(Rcpp::NumericVector X);
    Rcpp::NumericMatrix vector2mat(Rcpp::NumericVector X, Rcpp::IntegerVector i, Rcpp::IntegerVector j, int nfeats);
    Rcpp::NumericMatrix ratiosRcpp(Rcpp::NumericMatrix X);
    Rcpp::NumericMatrix results2matRcpp(Rcpp::DataFrame results, int n, double diagonal);
}
