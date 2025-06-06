#include <Rcpp.h>
#include <math.h>
#include <propr/kernels/cuda/dispatch/backend.cuh>
#include <propr/utils.hpp>

using namespace Rcpp;
using namespace propr;


void dispatch::cuda::wtmRcpp(const NumericVector& x, const NumericVector& w, double& out, propr_context context){
  Rcpp::stop("wtmRcpp is not implemented in CUDA. Falling back to CPU via dispatcher.");
}

void dispatch::cuda::wtvRcpp(const NumericVector& x, const NumericVector& w, double& out, propr_context context) {
  Rcpp::stop("wtvRcpp is not implemented in CUDA. Falling back to CPU via dispatcher.");
}

void centerNumericMatrix(const NumericMatrix & X, NumericMatrix& out, propr_context context){
  CHECK_MATRIX_DIMS(out, X.nrow(), X.ncol());
  Rcpp::stop("centerNumericMatrix is not implemented in CUDA. Falling back to CPU via dispatcher.");
}

void dispatch::cuda::corRcpp(const NumericMatrix & X, NumericMatrix& out, propr_context context) {
  CHECK_MATRIX_DIMS(out, X.ncol(), X.ncol());
  Rcpp::stop("corRcpp is not implemented in CUDA. Falling back to CPU via dispatcher.");
}

void dispatch::cuda::covRcpp(const NumericMatrix & X, const int norm_type, NumericMatrix& out, propr_context context) {
  CHECK_MATRIX_DIMS(out, X.ncol(), X.ncol());
  Rcpp::stop("covRcpp is not implemented in CUDA. Falling back to CPU via dispatcher.");
}

void dispatch::cuda::vlrRcpp(const NumericMatrix & X, NumericMatrix& out, propr_context context){
  CHECK_MATRIX_DIMS(out, X.ncol(), X.ncol());
  Rcpp::stop("vlrRcpp is not implemented in CUDA. Falling back to CPU via dispatcher.");
}

void dispatch::cuda::clrRcpp(const NumericMatrix & X, NumericMatrix& out, propr_context context){
  CHECK_MATRIX_DIMS(out, X.nrow(), X.ncol());
  Rcpp::stop("clrRcpp is not implemented in CUDA. Falling back to CPU via dispatcher.");
}

void dispatch::cuda::alrRcpp(const NumericMatrix & X, const int ivar, NumericMatrix& out, propr_context context){
  if(ivar == 0) {
    Rcpp::stop("Select non-zero ivar for alrRcpp.");
  }
  CHECK_MATRIX_DIMS(out, X.nrow(), X.ncol());
  Rcpp::stop("alrRcpp is not implemented in CUDA. Falling back to CPU via dispatcher.");
}

void dispatch::cuda::symRcpp(const NumericMatrix & X, NumericMatrix& out, propr_context context) {
  CHECK_MATRIX_DIMS(out, X.nrow(), X.ncol());
  Rcpp::stop("symRcpp is not implemented in CUDA. Falling back to CPU via dispatcher.");
}

void dispatch::cuda::phiRcpp(const NumericMatrix &X, const bool sym, NumericMatrix& out, propr_context context) {
  CHECK_MATRIX_DIMS(out, X.ncol(), X.ncol());
  Rcpp::stop("phiRcpp is not implemented in CUDA. Falling back to CPU via dispatcher.");
}

void dispatch::cuda::rhoRcpp(const NumericMatrix &X, const NumericMatrix &lr, const int ivar, NumericMatrix& out, propr_context context){
  CHECK_MATRIX_DIMS(out, X.ncol(), X.ncol());
  Rcpp::stop("rhoRcpp is not implemented in CUDA. Falling back to CPU via dispatcher.");
}


void dispatch::cuda::indexPairs(const NumericMatrix & X, const String op, const double ref, std::vector<int>& out, propr_context context){
    out.clear();
    Rcpp::stop("indexPairs is not implemented in CUDA. Falling back to CPU via dispatcher.");
}

void dispatch::cuda::indexToCoord(IntegerVector V, int N, List& out, propr_context context){
    Rcpp::stop("indexToCoord is not implemented in CUDA. Falling back to CPU via dispatcher.");
}

void dispatch::cuda::coordToIndex(IntegerVector row, IntegerVector col, int N, IntegerVector& out, propr_context context){
    CHECK_VECTOR_SIZE(out, row.length());
    Rcpp::stop("coordToIndex is not implemented in CUDA. Falling back to CPU via dispatcher.");
}

void dispatch::cuda::linRcpp(const NumericMatrix & rho, const NumericMatrix &lr, NumericMatrix& out, propr_context context){
    CHECK_MATRIX_DIMS(out, rho.nrow(), rho.ncol());
    Rcpp::stop("linRcpp is not implemented in CUDA. Falling back to CPU via dispatcher.");
}

void dispatch::cuda::lltRcpp(const NumericMatrix & X, NumericVector& out, propr_context context){
    int nfeats = X.nrow();
    int llt = nfeats * (nfeats - 1) / 2;
    CHECK_VECTOR_SIZE(out, llt);
    Rcpp::stop("lltRcpp is not implemented in CUDA. Falling back to CPU via dispatcher.");
}

void dispatch::cuda::urtRcpp(const NumericMatrix & X, NumericVector& out, propr_context context){
    int nfeats = X.nrow();
    int llt = nfeats * (nfeats - 1) / 2;
    CHECK_VECTOR_SIZE(out, llt);
    Rcpp::stop("urtRcpp is not implemented in CUDA. Falling back to CPU via dispatcher.");
}

void dispatch::cuda::labRcpp(int nfeats, List & out, propr_context context){
    Rcpp::stop("labRcpp is not implemented in CUDA. Falling back to CPU via dispatcher.");
}

void dispatch::cuda::half2mat(const NumericVector & X, NumericMatrix& out, propr_context context){
    int nfeats = round(sqrt(2 * X.length() + 0.25) + 0.5); // Re-calculate nfeats from X.length()
    CHECK_MATRIX_DIMS(out, nfeats, nfeats);
    Rcpp::stop("half2mat is not implemented in CUDA. Falling back to CPU via dispatcher.");
}

void dispatch::cuda::vector2mat(const NumericVector & X, const IntegerVector & i, const IntegerVector & j, int nfeats, NumericMatrix& out, propr_context context){
    CHECK_MATRIX_DIMS(out, nfeats, nfeats);
    Rcpp::stop("vector2mat is not implemented in CUDA. Falling back to CPU via dispatcher.");
}

void dispatch::cuda::ratiosRcpp(const NumericMatrix & X, NumericMatrix & out, propr_context context){
    int nfeats = X.ncol();
    int nsamps = X.nrow();
    int llt = nfeats * (nfeats - 1) / 2;
    CHECK_MATRIX_DIMS(out, nsamps, llt);
    Rcpp::stop("ratiosRcpp is not implemented in CUDA. Falling back to CPU via dispatcher.");
}

void dispatch::cuda::results2matRcpp(const DataFrame& results, int n, double diagonal, NumericMatrix & out, propr_context context){
    CHECK_MATRIX_DIMS(out, n, n);
    Rcpp::stop("results2matRcpp is not implemented in CUDA. Falling back to CPU via dispatcher.");
}