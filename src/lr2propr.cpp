#include <Rcpp.h>

#include <propr/interface/lr2propr.hpp>
#include <propr/kernels/cpu/dispatch/lr2propr.hpp>

using namespace propr;

// CUDA doesn't support lr2* ops, always use CPU.

// [[Rcpp::export]]
Rcpp::NumericMatrix lr2vlr(Rcpp::NumericMatrix lr, Rcpp::String backend = "auto") {
    const int nfeats = lr.ncol();
    Rcpp::NumericMatrix result(nfeats, nfeats);
    dispatch::cpu::lr2vlr(result, lr);
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix lr2phi(Rcpp::NumericMatrix lr, Rcpp::String backend = "auto") {
    const int nfeats = lr.ncol();
    Rcpp::NumericMatrix result(nfeats, nfeats);
    dispatch::cpu::lr2phi(result, lr);
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix lr2rho(Rcpp::NumericMatrix lr, Rcpp::String backend = "auto") {
    const int nfeats = lr.ncol();
    Rcpp::NumericMatrix result(nfeats, nfeats);
    dispatch::cpu::lr2rho(result, lr);
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix lr2phs(Rcpp::NumericMatrix lr, Rcpp::String backend = "auto") {
    const int nfeats = lr.ncol();
    Rcpp::NumericMatrix result(nfeats, nfeats);
    dispatch::cpu::lr2phs(result, lr);
    return result;
}
