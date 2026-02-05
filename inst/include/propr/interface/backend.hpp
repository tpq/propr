#pragma once

#include <Rcpp.h>

namespace propr {
    Rcpp::NumericMatrix covRcpp(Rcpp::NumericMatrix & X, const int norm_type,bool use_gpu = false);
    Rcpp::NumericMatrix corRcpp(Rcpp::NumericMatrix X,bool use_gpu = false);
    double              wtmRcpp(Rcpp::NumericVector x, Rcpp::NumericVector w,bool use_gpu = false);
    double              wtvRcpp(Rcpp::NumericVector x, Rcpp::NumericVector w,bool use_gpu = false);
    Rcpp::NumericMatrix covRcpp(Rcpp::NumericMatrix X, int norm_type,bool use_gpu = false);
    Rcpp::NumericMatrix vlrRcpp(Rcpp::NumericMatrix X,bool use_gpu = false);
    Rcpp::NumericMatrix clrRcpp(Rcpp::NumericMatrix X,bool use_gpu = false);
    Rcpp::NumericMatrix alrRcpp(Rcpp::NumericMatrix X, int ivar,bool use_gpu = false);
    Rcpp::NumericMatrix symRcpp(Rcpp::NumericMatrix X,bool use_gpu = false);
    Rcpp::NumericMatrix phiRcpp(Rcpp::NumericMatrix X, bool sym,bool use_gpu = false);
    Rcpp::NumericMatrix rhoRcpp(Rcpp::NumericMatrix X, Rcpp::NumericMatrix lr, int ivar,bool use_gpu = false);
    std::vector<int> indexPairs(Rcpp::NumericMatrix X, Rcpp::String op, double ref,bool use_gpu = false);
    Rcpp::List indexToCoord(Rcpp::IntegerVector V, int N,bool use_gpu = false);
    Rcpp::IntegerVector coordToIndex(Rcpp::IntegerVector row, Rcpp::IntegerVector col, int N,bool use_gpu = false);
    Rcpp::NumericMatrix linRcpp(Rcpp::NumericMatrix rho, Rcpp::NumericMatrix lr, bool use_gpu = false);
    Rcpp::NumericVector lltRcpp(Rcpp::NumericMatrix X,bool use_gpu = false);
    Rcpp::NumericVector urtRcpp(Rcpp::NumericMatrix X,bool use_gpu = false);
    Rcpp::List labRcpp(int nfeats,bool use_gpu = false);
    Rcpp::NumericMatrix half2mat(Rcpp::NumericVector X,bool use_gpu = false);
    Rcpp::NumericMatrix vector2mat(Rcpp::NumericVector X, Rcpp::IntegerVector i, Rcpp::IntegerVector j, int nfeatsbool, bool use_gpu = false);
    Rcpp::NumericMatrix ratiosRcpp(Rcpp::NumericMatrix X, bool use_gpu = false);
    Rcpp::NumericMatrix results2matRcpp(Rcpp::DataFrame results, int n, double diagonal, bool use_gpu = false);
}
