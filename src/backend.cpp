#include <Rcpp.h>

#include <propr/interface/backend.hpp>
#include <propr/kernels/cpu/dispatch/backend.hpp>
#include <propr/runtime/cuda_executor.hpp>
#include <propr/runtime/dispatch.hpp>

#include <cmath>

using namespace propr;

// [[Rcpp::export]]
double wtmRcpp(Rcpp::NumericVector x, Rcpp::NumericVector w, Rcpp::String backend = "auto") {
    double result = 0.0;
    if (runtime::resolve_backend(backend) == runtime::Backend::CUDA) {
        runtime::cuda_executor::wtmRcpp(result, x, w);
    } else {
        dispatch::cpu::wtmRcpp(result, x, w);
    }
    return result;
}

// [[Rcpp::export]]
double wtvRcpp(Rcpp::NumericVector x, Rcpp::NumericVector w, Rcpp::String backend = "auto") {
    double result = 0.0;
    if (runtime::resolve_backend(backend) == runtime::Backend::CUDA) {
        runtime::cuda_executor::wtvRcpp(result, x, w);
    } else {
        dispatch::cpu::wtvRcpp(result, x, w);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix corRcpp(Rcpp::NumericMatrix X, Rcpp::String backend = "auto") {
    const int nfeats = X.ncol();
    Rcpp::NumericMatrix result(nfeats, nfeats);
    if (runtime::resolve_backend(backend) == runtime::Backend::CUDA) {
        runtime::cuda_executor::corRcpp(result, X);
    } else {
        dispatch::cpu::corRcpp(result, X);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix covRcpp(Rcpp::NumericMatrix X, int norm_type, Rcpp::String backend = "auto") {
    const int nfeats = X.ncol();
    Rcpp::NumericMatrix result(nfeats, nfeats);
    if (runtime::resolve_backend(backend) == runtime::Backend::CUDA) {
        runtime::cuda_executor::covRcpp(result, X, norm_type);
    } else {
        dispatch::cpu::covRcpp(result, X, norm_type);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix vlrRcpp(Rcpp::NumericMatrix X, Rcpp::String backend = "auto") {
    const int nfeats = X.ncol();
    Rcpp::NumericMatrix result(nfeats, nfeats);
    if (runtime::resolve_backend(backend) == runtime::Backend::CUDA) {
        runtime::cuda_executor::vlrRcpp(result, X);
    } else {
        dispatch::cpu::vlrRcpp(result, X);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix clrRcpp(Rcpp::NumericMatrix X, Rcpp::String backend = "auto") {
    const int n_rows = X.nrow();
    const int n_cols = X.ncol();
    Rcpp::NumericMatrix result(n_rows, n_cols);
    if (runtime::resolve_backend(backend) == runtime::Backend::CUDA) {
        runtime::cuda_executor::clrRcpp(result, X);
    } else {
        dispatch::cpu::clrRcpp(result, X);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix alrRcpp(Rcpp::NumericMatrix X, int ivar, Rcpp::String backend = "auto") {
    const int n_rows = X.nrow();
    const int n_cols = X.ncol();
    Rcpp::NumericMatrix result(n_rows, n_cols);
    if (runtime::resolve_backend(backend) == runtime::Backend::CUDA) {
        runtime::cuda_executor::alrRcpp(result, X, ivar);
    } else {
        dispatch::cpu::alrRcpp(result, X, ivar);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix symRcpp(Rcpp::NumericMatrix X, Rcpp::String backend = "auto") {
    const int n_rows = X.nrow();
    const int n_cols = X.ncol();
    Rcpp::NumericMatrix result(n_rows, n_cols);
    if (runtime::resolve_backend(backend) == runtime::Backend::CUDA) {
        runtime::cuda_executor::symRcpp(result, X);
    } else {
        dispatch::cpu::symRcpp(result, X);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix phiRcpp(Rcpp::NumericMatrix X, bool sym, Rcpp::String backend = "auto") {
    const int nfeats = X.ncol();
    Rcpp::NumericMatrix result(nfeats, nfeats);
    if (runtime::resolve_backend(backend) == runtime::Backend::CUDA) {
        runtime::cuda_executor::phiRcpp(result, X, sym);
    } else {
        dispatch::cpu::phiRcpp(result, X, sym);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix rhoRcpp(Rcpp::NumericMatrix X, Rcpp::NumericMatrix lr, int ivar, Rcpp::String backend = "auto") {
    const int nfeats = X.ncol();
    Rcpp::NumericMatrix result(nfeats, nfeats);
    if (runtime::resolve_backend(backend) == runtime::Backend::CUDA) {
        runtime::cuda_executor::rhoRcpp(result, X, lr, ivar);
    } else {
        dispatch::cpu::rhoRcpp(result, X, lr, ivar);
    }
    return result;
}

// [[Rcpp::export]]
std::vector<int> indexPairs(Rcpp::NumericMatrix X, Rcpp::String op, double ref, Rcpp::String backend = "auto") {
    std::vector<int> result;
    if (runtime::resolve_backend(backend) == runtime::Backend::CUDA) {
        runtime::cuda_executor::indexPairs(result, X, op, ref);
    } else {
        dispatch::cpu::indexPairs(result, X, op, ref);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::List indexToCoord(Rcpp::IntegerVector V, int N, Rcpp::String backend = "auto") {
    Rcpp::List result;
    if (runtime::resolve_backend(backend) == runtime::Backend::CUDA) {
        runtime::cuda_executor::indexToCoord(result, V, N);
    } else {
        dispatch::cpu::indexToCoord(result, V, N);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::IntegerVector coordToIndex(Rcpp::IntegerVector row, Rcpp::IntegerVector col, int N, Rcpp::String backend = "auto") {
    const int size = row.length();
    Rcpp::IntegerVector result(size);
    if (runtime::resolve_backend(backend) == runtime::Backend::CUDA) {
        runtime::cuda_executor::coordToIndex(result, row, col, N);
    } else {
        dispatch::cpu::coordToIndex(result, row, col, N);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix linRcpp(Rcpp::NumericMatrix rho, Rcpp::NumericMatrix lr, Rcpp::String backend = "auto") {
    const int n_cols = rho.ncol();
    Rcpp::NumericMatrix result(n_cols, n_cols);
    if (runtime::resolve_backend(backend) == runtime::Backend::CUDA) {
        runtime::cuda_executor::linRcpp(result, rho, lr);
    } else {
        dispatch::cpu::linRcpp(result, rho, lr);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericVector lltRcpp(Rcpp::NumericMatrix X, Rcpp::String backend = "auto") {
    const int nfeats = X.nrow();
    const int llt = nfeats * (nfeats - 1) / 2;
    Rcpp::NumericVector result(llt);
    if (runtime::resolve_backend(backend) == runtime::Backend::CUDA) {
        runtime::cuda_executor::lltRcpp(result, X);
    } else {
        dispatch::cpu::lltRcpp(result, X);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericVector urtRcpp(Rcpp::NumericMatrix X, Rcpp::String backend = "auto") {
    const int nfeats = X.nrow();
    const int llt = nfeats * (nfeats - 1) / 2;
    Rcpp::NumericVector result(llt);
    if (runtime::resolve_backend(backend) == runtime::Backend::CUDA) {
        runtime::cuda_executor::urtRcpp(result, X);
    } else {
        dispatch::cpu::urtRcpp(result, X);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::List labRcpp(int nfeats, Rcpp::String backend = "auto") {
    Rcpp::List result;
    if (runtime::resolve_backend(backend) == runtime::Backend::CUDA) {
        runtime::cuda_executor::labRcpp(result, nfeats);
    } else {
        dispatch::cpu::labRcpp(result, nfeats);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix half2mat(Rcpp::NumericVector X, Rcpp::String backend = "auto") {
    const int nfeats = static_cast<int>(std::round(std::sqrt(2.0 * X.length() + 0.25) + 0.5));
    Rcpp::NumericMatrix result(nfeats, nfeats);
    if (runtime::resolve_backend(backend) == runtime::Backend::CUDA) {
        runtime::cuda_executor::half2mat(result, X);
    } else {
        dispatch::cpu::half2mat(result, X);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix vector2mat(Rcpp::NumericVector X, Rcpp::IntegerVector i, Rcpp::IntegerVector j, int nfeats, Rcpp::String backend = "auto") {
    Rcpp::NumericMatrix result(nfeats, nfeats);
    if (runtime::resolve_backend(backend) == runtime::Backend::CUDA) {
        runtime::cuda_executor::vector2mat(result, X, i, j, nfeats);
    } else {
        dispatch::cpu::vector2mat(result, X, i, j, nfeats);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix ratiosRcpp(Rcpp::NumericMatrix X, Rcpp::String backend = "auto") {
    const int nfeats = X.ncol();
    const int nsamps = X.nrow();
    const int llt = nfeats * (nfeats - 1) / 2;
    Rcpp::NumericMatrix result(nsamps, llt);
    if (runtime::resolve_backend(backend) == runtime::Backend::CUDA) {
        runtime::cuda_executor::ratiosRcpp(result, X);
    } else {
        dispatch::cpu::ratiosRcpp(result, X);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix results2matRcpp(Rcpp::DataFrame results, int n, double diagonal, Rcpp::String backend = "auto") {
    Rcpp::NumericMatrix out(n, n);
    if (runtime::resolve_backend(backend) == runtime::Backend::CUDA) {
        runtime::cuda_executor::results2matRcpp(out, results, n, diagonal);
    } else {
        dispatch::cpu::results2matRcpp(out, results, n, diagonal);
    }
    return out;
}
