#pragma once

#include <Rcpp.h>

namespace propr {
    Rcpp::NumericMatrix covRcpp(Rcpp::NumericMatrix & X, const int norm_type);
    Rcpp::NumericMatrix propr::corRcpp(Rcpp::NumericMatrix X);
    double              wtmRcpp(Rcpp::NumericVector x, Rcpp::NumericVector w);
    double              wtvRcpp(Rcpp::NumericVector x, Rcpp::NumericVector w);
    Rcpp::NumericMatrix propr::covRcpp(Rcpp::NumericMatrix X, int norm_type);
    Rcpp::NumericMatrix propr::vlrRcpp(Rcpp::NumericMatrix X);
    Rcpp::NumericMatrix propr::clrRcpp(Rcpp::NumericMatrix X);
    Rcpp::NumericMatrix propr::alrRcpp(Rcpp::NumericMatrix X, int ivar);
    Rcpp::NumericMatrix propr::symRcpp(Rcpp::NumericMatrix X);
    Rcpp::NumericMatrix propr::phiRcpp(Rcpp::NumericMatrix X, bool sym);
    Rcpp::NumericMatrix propr::rhoRcpp(Rcpp::NumericMatrix X, Rcpp::NumericMatrix lr, int ivar);
    std::vector<int> propr::indexPairs(Rcpp::NumericMatrix X, Rcpp::String op, double ref);
    Rcpp::List propr::indexToCoord(Rcpp::IntegerVector V, int N);
    Rcpp::IntegerVector propr::coordToIndex(Rcpp::IntegerVector row, Rcpp::IntegerVector col, int N);
    Rcpp::NumericMatrix propr::linRcpp(Rcpp::NumericMatrix rho, Rcpp::NumericMatrix lr);
    Rcpp::NumericVector propr::lltRcpp(Rcpp::NumericMatrix X);
    Rcpp::NumericVector propr::urtRcpp(Rcpp::NumericMatrix X);
    Rcpp::List propr::labRcpp(int nfeats);
    Rcpp::NumericMatrix propr::half2mat(Rcpp::NumericVector X);
    Rcpp::NumericMatrix propr::vector2mat(Rcpp::NumericVector X, Rcpp::IntegerVector i, Rcpp::IntegerVector j, int nfeats);
    Rcpp::NumericMatrix propr::ratiosRcpp(Rcpp::NumericMatrix X);
    Rcpp::NumericMatrix propr::results2matRcpp(Rcpp::DataFrame results, int n, double diagonal);
}
