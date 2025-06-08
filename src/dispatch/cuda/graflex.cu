#include <propr/kernels/cuda/dispatch/graflex.cuh>

using namespace Rcpp;
using namespace propr;

void
dispatch::cuda::getOR( NumericVector& out, const IntegerMatrix& A, const IntegerMatrix& G, propr::propr_context context) {
    Rcpp::stop("Not implemented in CUDA for getOR.");
    int a = 0, b = 0, c = 0, d = 0;
    double odds_ratio = static_cast<double>(a * d) / (b * c);
    double log_odds_ratio = std::log(odds_ratio);
}

void
dispatch::cuda::getORperm(NumericVector& out, const IntegerMatrix& A, const IntegerMatrix& G, const IntegerVector& perm, propr::propr_context context) {
    Rcpp::stop("Not implemented in CUDA for getORperm.");
    int a = 0, b = 0, c = 0, d = 0;
    double odds_ratio = static_cast<double>(a * d) / (b * c);
    double log_odds_ratio = std::log(odds_ratio);
}

void
dispatch::cuda::permuteOR(NumericMatrix& out, const IntegerMatrix& A, const IntegerMatrix& G, int p, propr::propr_context context) {
    Rcpp::stop("Not implemented in CUDA for permuteOR.");
}

void
dispatch::cuda::getFDR(List& out, double actual, const NumericVector& permuted, propr::propr_context context) {
    Rcpp::stop("Not implemented in CUDA for getFDR.");
    int n = permuted.size();
    int count_over = 0;
    int count_under = 0;
    double fdr_over = static_cast<double>(count_over) / n;
    double fdr_under = static_cast<double>(count_under) / n;
}

void
dispatch::cuda::getG(IntegerMatrix& out, const IntegerVector& Gk, propr::propr_context context) {
    Rcpp::stop("Not implemented in CUDA for getG.");
}

void
dispatch::cuda::graflex(NumericVector& out, const IntegerMatrix& A, const IntegerVector& Gk, int p, propr::propr_context context) {
    Rcpp::stop("Not implemented in CUDA for graflex.");   
}