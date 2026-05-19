#pragma once

#include <Rcpp.h>
#include <vector>

namespace propr {
    Rcpp::NumericMatrix covRcpp(Rcpp::NumericMatrix& X, const int norm_type, Rcpp::String backend = "auto");
    Rcpp::NumericMatrix corRcpp(Rcpp::NumericMatrix X, Rcpp::String backend = "auto");
    double wtmRcpp(Rcpp::NumericVector x, Rcpp::NumericVector w, Rcpp::String backend = "auto");
    double wtvRcpp(Rcpp::NumericVector x, Rcpp::NumericVector w, Rcpp::String backend = "auto");
    Rcpp::NumericMatrix covRcpp(Rcpp::NumericMatrix X, int norm_type, Rcpp::String backend = "auto");
    Rcpp::NumericMatrix vlrRcpp(Rcpp::NumericMatrix X, Rcpp::String backend = "auto");
    Rcpp::NumericMatrix clrRcpp(Rcpp::NumericMatrix X, Rcpp::String backend = "auto");
    Rcpp::NumericMatrix alrRcpp(Rcpp::NumericMatrix X, int ivar, Rcpp::String backend = "auto");
    Rcpp::NumericMatrix symRcpp(Rcpp::NumericMatrix X, Rcpp::String backend = "auto");
    Rcpp::NumericMatrix phiRcpp(Rcpp::NumericMatrix X, bool sym, Rcpp::String backend = "auto");
    Rcpp::NumericMatrix rhoRcpp(Rcpp::NumericMatrix X, Rcpp::NumericMatrix lr, int ivar, Rcpp::String backend = "auto");
    std::vector<int> indexPairs(Rcpp::NumericMatrix X, Rcpp::String op, double ref, Rcpp::String backend = "auto");
    Rcpp::List indexToCoord(Rcpp::IntegerVector V, int N, Rcpp::String backend = "auto");
    Rcpp::IntegerVector coordToIndex(Rcpp::IntegerVector row, Rcpp::IntegerVector col, int N, Rcpp::String backend = "auto");
    Rcpp::NumericMatrix linRcpp(Rcpp::NumericMatrix rho, Rcpp::NumericMatrix lr, Rcpp::String backend = "auto");
    Rcpp::NumericVector lltRcpp(Rcpp::NumericMatrix X, Rcpp::String backend = "auto");
    Rcpp::NumericVector urtRcpp(Rcpp::NumericMatrix X, Rcpp::String backend = "auto");
    Rcpp::List labRcpp(int nfeats, Rcpp::String backend = "auto");
    Rcpp::NumericMatrix half2mat(Rcpp::NumericVector X, Rcpp::String backend = "auto");
    Rcpp::NumericMatrix vector2mat(Rcpp::NumericVector X, Rcpp::IntegerVector i, Rcpp::IntegerVector j, int nfeatsbool, Rcpp::String backend = "auto");
    Rcpp::NumericMatrix ratiosRcpp(Rcpp::NumericMatrix X, Rcpp::String backend = "auto");
    Rcpp::NumericMatrix results2matRcpp(Rcpp::DataFrame results, int n, double diagonal, Rcpp::String backend = "auto");
}
