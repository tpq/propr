#include <Rcpp.h>
#include <vector>

#include <propr/interface/backend.hpp>
#include <propr/interface/device_selector.hpp>
#include <propr/context.h>

#include <propr/kernels/cpu/dispatch/backend.hpp>
#include <propr/kernels/cuda/dispatch/backend.cuh>


using namespace propr;

// [[Rcpp::export]]
double wtmRcpp(Rcpp::NumericVector x, Rcpp::NumericVector w) {
    double result;
    dispatch::cpu::wtmRcpp(result, x, w);
    return result;
}

// [[Rcpp::export]]
double wtvRcpp(Rcpp::NumericVector x, Rcpp::NumericVector w) {
    double result;
    dispatch::cpu::wtvRcpp(result, x, w);
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix corRcpp(Rcpp::NumericMatrix X) {
    int nfeats = X.ncol();
    Rcpp::NumericMatrix result(nfeats, nfeats);

    if (is_gpu_backend()) {
        dispatch::cuda::corRcpp(result, X);
    } else {
        dispatch::cpu::corRcpp(result, X);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix covRcpp(Rcpp::NumericMatrix X, int norm_type) {
    int nfeats = X.ncol();
    Rcpp::NumericMatrix result(nfeats, nfeats);

    if (is_gpu_backend()) {
        dispatch::cuda::covRcpp(result, X, norm_type);
    } else {
        dispatch::cpu::covRcpp(result, X, norm_type);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix vlrRcpp(Rcpp::NumericMatrix X) {
    int nfeats = X.ncol();
    Rcpp::NumericMatrix result(nfeats, nfeats);
    if (is_gpu_backend()) {
        dispatch::cuda::vlrRcpp(result, X);
    } else {
        dispatch::cpu::vlrRcpp(result, X);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix clrRcpp(Rcpp::NumericMatrix X) {
    int n_rows = X.nrow();
    int n_cols = X.ncol();
    Rcpp::NumericMatrix result(n_rows, n_cols);

    if (is_gpu_backend()) {
        dispatch::cuda::clrRcpp(result, X);
    } else {
        dispatch::cpu::clrRcpp(result, X);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix alrRcpp(Rcpp::NumericMatrix X, int ivar) {
    int n_rows = X.nrow();
    int n_cols = X.ncol();
    Rcpp::NumericMatrix result(n_rows, n_cols);

    if (is_gpu_backend()) {
        dispatch::cuda::alrRcpp(result, X, ivar);
    } else {
        dispatch::cpu::alrRcpp(result, X, ivar);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix symRcpp(Rcpp::NumericMatrix X) {
    int n_rows = X.nrow();
    int n_cols = X.ncol();
    Rcpp::NumericMatrix result(n_rows, n_cols);

    if (is_gpu_backend()) {
        dispatch::cuda::symRcpp(result, X);
    } else {
        dispatch::cpu::symRcpp(result, X);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix phiRcpp(Rcpp::NumericMatrix X, bool sym) {
    int nfeats = X.ncol();
    Rcpp::NumericMatrix result(nfeats, nfeats);

    if (is_gpu_backend()) {
        dispatch::cuda::phiRcpp(result, X, sym);
    } else {
        dispatch::cpu::phiRcpp(result, X, sym);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix rhoRcpp(Rcpp::NumericMatrix X, Rcpp::NumericMatrix lr, int ivar) {
    int nfeats = X.ncol();
    Rcpp::NumericMatrix result(nfeats, nfeats);

    if (is_gpu_backend()) {
        dispatch::cuda::rhoRcpp(result, X, lr, ivar);
    } else {
        dispatch::cpu::rhoRcpp(result, X, lr, ivar);
    }
    return result;
}

// [[Rcpp::export]]
std::vector<int> indexPairs(Rcpp::NumericMatrix X, Rcpp::String op, double ref) {
    std::vector<int> result;
    if (is_gpu_backend()) {
        dispatch::cuda::indexPairs( result, X, op, ref);
    } else {
        dispatch::cpu::indexPairs(result, X, op, ref);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::List indexToCoord(Rcpp::IntegerVector V, int N) {
    Rcpp::List result;
    if (is_gpu_backend()) {
        dispatch::cuda::indexToCoord(result, V, N);
    } else {
        dispatch::cpu::indexToCoord(result, V, N);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::IntegerVector coordToIndex(Rcpp::IntegerVector row, Rcpp::IntegerVector col, int N) {
    int size = row.length();
    Rcpp::IntegerVector result(size);
    if (is_gpu_backend()) {
        dispatch::cuda::coordToIndex(result, row, col, N);
    } else {
        dispatch::cpu::coordToIndex(result, row, col, N);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix linRcpp(Rcpp::NumericMatrix rho, Rcpp::NumericMatrix lr) {
    int n_rows = rho.nrow();
    int n_cols = rho.ncol();
    Rcpp::NumericMatrix result(n_rows, n_cols);
    if (is_gpu_backend()) {
        dispatch::cuda::linRcpp(result, rho, lr);
    } else {
        dispatch::cpu::linRcpp(result, rho, lr);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericVector lltRcpp(Rcpp::NumericMatrix X) {
    int nfeats = X.nrow();
    int llt = nfeats * (nfeats - 1) / 2;
    Rcpp::NumericVector result(llt);
     if (is_gpu_backend()) {
        dispatch::cuda::lltRcpp(result, X);
    } else {
        dispatch::cpu::lltRcpp(result, X);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericVector urtRcpp(Rcpp::NumericMatrix X) {
    int nfeats = X.nrow();
    int llt = nfeats * (nfeats - 1) / 2;
    Rcpp::NumericVector result(llt);
    if (is_gpu_backend()) {
        dispatch::cuda::urtRcpp(result, X);
    } else {
        dispatch::cpu::urtRcpp(result, X);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::List labRcpp(int nfeats) {
    Rcpp::List result;
     if (is_gpu_backend()) {
        dispatch::cuda::labRcpp(result, nfeats);
    } else {
        dispatch::cpu::labRcpp(result, nfeats);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix half2mat(Rcpp::NumericVector X) {
    int nfeats = round(sqrt(2 * X.length() + 0.25) + 0.5);
    Rcpp::NumericMatrix result(nfeats, nfeats);
    if (is_gpu_backend()) {
        dispatch::cuda::half2mat(result, X);
    } else {
        dispatch::cpu::half2mat(result, X);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix vector2mat(Rcpp::NumericVector X, Rcpp::IntegerVector i, Rcpp::IntegerVector j, int nfeats) {
    Rcpp::NumericMatrix result(nfeats, nfeats);
    if (is_gpu_backend()) {
        dispatch::cuda::vector2mat(result, X, i, j, nfeats);
    } else {
        dispatch::cpu::vector2mat(result, X, i, j, nfeats);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix ratiosRcpp(Rcpp::NumericMatrix X) {
    int nfeats = X.ncol();
    int nsamps = X.nrow();
    int llt = nfeats * (nfeats - 1) / 2;
    Rcpp::NumericMatrix result(nsamps, llt);

    if (is_gpu_backend()) {
        dispatch::cuda::ratiosRcpp(result, X);
    } else {
        dispatch::cpu::ratiosRcpp(result, X);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix results2matRcpp(Rcpp::DataFrame results, int n, double diagonal) {
    Rcpp::NumericMatrix result(n, n);
    if (is_gpu_backend()) {
        dispatch::cuda::results2matRcpp(result, results, n, diagonal);
    } else {
        dispatch::cpu::results2matRcpp(result, results, n, diagonal);
    }
    return result;
}
