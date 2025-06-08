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
    dispatch::cpu::wtmRcpp(x, w, result);
    return result;
}

// [[Rcpp::export]]
double wtvRcpp(Rcpp::NumericVector x, Rcpp::NumericVector w) {
    double result;
    dispatch::cpu::wtvRcpp(x, w, result);
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix corRcpp(Rcpp::NumericMatrix X) {
    int nfeats = X.ncol();
    Rcpp::NumericMatrix result(nfeats, nfeats);

    if (is_gpu_backend()) {
        propr_context context;
        cudaStream_t stream;
        cudaError_t err = cudaStreamCreate(&stream);
        if (err != cudaSuccess) {
            Rcpp::warning("CUDA stream creation failed for corRcpp: %s. Falling back to CPU.", cudaGetErrorString(err));
            dispatch::cpu::corRcpp(X, result);
        } else {
            context.stream = stream;
            dispatch::cuda::corRcpp(X, result, context);
            cudaStreamDestroy(stream);
        }
    } else {
        dispatch::cpu::corRcpp(X, result);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix covRcpp(Rcpp::NumericMatrix X, int norm_type) {
    int nfeats = X.ncol();
    Rcpp::NumericMatrix result(nfeats, nfeats);

    if (is_gpu_backend()) {
        propr_context context;
        cudaStream_t stream;
        cudaError_t err = cudaStreamCreate(&stream);
        if (err != cudaSuccess) {
            Rcpp::warning("CUDA stream creation failed for covRcpp: %s. Falling back to CPU.", cudaGetErrorString(err));
            dispatch::cpu::covRcpp(X, norm_type, result);
        } else {
            context.stream = stream;
            dispatch::cuda::covRcpp(X, norm_type, result, context);
            cudaStreamDestroy(stream);
        }
    } else {
        dispatch::cpu::covRcpp(X, norm_type, result);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix vlrRcpp(Rcpp::NumericMatrix X) {
    int nfeats = X.ncol();
    Rcpp::NumericMatrix result(nfeats, nfeats);

    if (is_gpu_backend()) {
        propr_context context;
        cudaStream_t stream;
        cudaError_t err = cudaStreamCreate(&stream);
        if (err != cudaSuccess) {
            Rcpp::warning("CUDA stream creation failed for vlrRcpp: %s. Falling back to CPU.", cudaGetErrorString(err));
            dispatch::cpu::vlrRcpp(X, result);
        } else {
            context.stream = stream;
            dispatch::cuda::vlrRcpp(X, result, context);
            cudaStreamDestroy(stream);
        }
    } else {
        dispatch::cpu::vlrRcpp(X, result);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix clrRcpp(Rcpp::NumericMatrix X) {
    int n_rows = X.nrow();
    int n_cols = X.ncol();
    Rcpp::NumericMatrix result(n_rows, n_cols);

    if (is_gpu_backend()) {
        propr_context context;
        cudaStream_t stream;
        cudaError_t err = cudaStreamCreate(&stream);
        if (err != cudaSuccess) {
            Rcpp::warning("CUDA stream creation failed for clrRcpp: %s. Falling back to CPU.", cudaGetErrorString(err));
            dispatch::cpu::clrRcpp(X, result);
        } else {
            context.stream = stream;
            dispatch::cuda::clrRcpp(X, result, context);
            cudaStreamDestroy(stream);
        }
    } else {
        dispatch::cpu::clrRcpp(X, result);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix alrRcpp(Rcpp::NumericMatrix X, int ivar) {
    int n_rows = X.nrow();
    int n_cols = X.ncol();
    Rcpp::NumericMatrix result(n_rows, n_cols);

    if (is_gpu_backend()) {
        propr_context context;
        cudaStream_t stream;
        cudaError_t err = cudaStreamCreate(&stream);
        if (err != cudaSuccess) {
            Rcpp::warning("CUDA stream creation failed for alrRcpp: %s. Falling back to CPU.", cudaGetErrorString(err));
            dispatch::cpu::alrRcpp(X, ivar, result);
        } else {
            context.stream = stream;
            dispatch::cuda::alrRcpp(X, ivar, result, context);
            cudaStreamDestroy(stream);
        }
    } else {
        dispatch::cpu::alrRcpp(X, ivar, result);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix symRcpp(Rcpp::NumericMatrix X) {
    int n_rows = X.nrow();
    int n_cols = X.ncol();
    Rcpp::NumericMatrix result(n_rows, n_cols);

    if (is_gpu_backend()) {
        propr_context context;
        cudaStream_t stream;
        cudaError_t err = cudaStreamCreate(&stream);
        if (err != cudaSuccess) {
            Rcpp::warning("CUDA stream creation failed for symRcpp: %s. Falling back to CPU.", cudaGetErrorString(err));
            dispatch::cpu::symRcpp(X, result);
        } else {
            context.stream = stream;
            dispatch::cuda::symRcpp(X, result, context);
            cudaStreamDestroy(stream);
        }
    } else {
        dispatch::cpu::symRcpp(X, result);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix phiRcpp(Rcpp::NumericMatrix X, bool sym) {
    int nfeats = X.ncol();
    Rcpp::NumericMatrix result(nfeats, nfeats);

    if (is_gpu_backend()) {
        propr_context context;
        cudaStream_t stream;
        cudaError_t err = cudaStreamCreate(&stream);
        if (err != cudaSuccess) {
            Rcpp::warning("CUDA stream creation failed for phiRcpp: %s. Falling back to CPU.", cudaGetErrorString(err));
            dispatch::cpu::phiRcpp(X, sym, result);
        } else {
            context.stream = stream;
            dispatch::cuda::phiRcpp(X, sym, result, context);
            cudaStreamDestroy(stream);
        }
    } else {
        dispatch::cpu::phiRcpp(X, sym, result);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix rhoRcpp(Rcpp::NumericMatrix X, Rcpp::NumericMatrix lr, int ivar) {
    int nfeats = X.ncol();
    Rcpp::NumericMatrix result(nfeats, nfeats);

    if (is_gpu_backend()) {
        propr_context context;
        cudaStream_t stream;
        cudaError_t err = cudaStreamCreate(&stream);
        if (err != cudaSuccess) {
            Rcpp::warning("CUDA stream creation failed for rhoRcpp: %s. Falling back to CPU.", cudaGetErrorString(err));
            dispatch::cpu::rhoRcpp(X, lr, ivar, result);
        } else {
            context.stream = stream;
            dispatch::cuda::rhoRcpp(X, lr, ivar, result, context);
            cudaStreamDestroy(stream);
        }
    } else {
        dispatch::cpu::rhoRcpp(X, lr, ivar, result);
    }
    return result;
}

// [[Rcpp::export]]
std::vector<int> indexPairs(Rcpp::NumericMatrix X, Rcpp::String op, double ref) {
    std::vector<int> result;
    if (is_gpu_backend()) {
         propr_context context;
        cudaStream_t stream;
        cudaError_t err = cudaStreamCreate(&stream);
        if (err != cudaSuccess) {
             Rcpp::warning("CUDA stream creation failed for indexPairs: %s. Falling back to CPU.", cudaGetErrorString(err));
             dispatch::cpu::indexPairs(X, op, ref, result);
        } else {
             context.stream = stream;
             dispatch::cuda::indexPairs(X, op, ref, result, context);
             cudaStreamDestroy(stream);
        }
    } else {
        dispatch::cpu::indexPairs(X, op, ref, result);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::List indexToCoord(Rcpp::IntegerVector V, int N) {
    Rcpp::List result;
    if (is_gpu_backend()) {
         propr_context context;
        cudaStream_t stream;
        cudaError_t err = cudaStreamCreate(&stream);
        if (err != cudaSuccess) {
             Rcpp::warning("CUDA stream creation failed for indexToCoord: %s. Falling back to CPU.", cudaGetErrorString(err));
             dispatch::cpu::indexToCoord(V, N, result);
        } else {
             context.stream = stream;
             dispatch::cuda::indexToCoord(V, N, result, context);
             cudaStreamDestroy(stream);
        }
    } else {
        dispatch::cpu::indexToCoord(V, N, result);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::IntegerVector coordToIndex(Rcpp::IntegerVector row, Rcpp::IntegerVector col, int N) {
    int size = row.length();
    Rcpp::IntegerVector result(size);
    if (is_gpu_backend()) {
        propr_context context;
        cudaStream_t stream;
        cudaError_t err = cudaStreamCreate(&stream);
        if (err != cudaSuccess) {
            Rcpp::warning("CUDA stream creation failed for coordToIndex: %s. Falling back to CPU.", cudaGetErrorString(err));
            dispatch::cpu::coordToIndex(row, col, N, result);
        } else {
            context.stream = stream;
            dispatch::cuda::coordToIndex(row, col, N, result, context);
            cudaStreamDestroy(stream);
        }
    } else {
        dispatch::cpu::coordToIndex(row, col, N, result);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix linRcpp(Rcpp::NumericMatrix rho, Rcpp::NumericMatrix lr) {
    int n_rows = rho.nrow();
    int n_cols = rho.ncol();
    Rcpp::NumericMatrix result(n_rows, n_cols);
    if (is_gpu_backend()) {
        propr_context context;
        cudaStream_t stream;
        cudaError_t err = cudaStreamCreate(&stream);
        if (err != cudaSuccess) {
            Rcpp::warning("CUDA stream creation failed for linRcpp: %s. Falling back to CPU.", cudaGetErrorString(err));
            dispatch::cpu::linRcpp(rho, lr, result);
        } else {
            context.stream = stream;
            dispatch::cuda::linRcpp(rho, lr, result, context);
            cudaStreamDestroy(stream);
        }
    } else {
        dispatch::cpu::linRcpp(rho, lr, result);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericVector lltRcpp(Rcpp::NumericMatrix X) {
    int nfeats = X.nrow();
    int llt = nfeats * (nfeats - 1) / 2;
    Rcpp::NumericVector result(llt);
     if (is_gpu_backend()) {
        propr_context context;
        cudaStream_t stream;
        cudaError_t err = cudaStreamCreate(&stream);
        if (err != cudaSuccess) {
            Rcpp::warning("CUDA stream creation failed for lltRcpp: %s. Falling back to CPU.", cudaGetErrorString(err));
            dispatch::cpu::lltRcpp(X, result);
        } else {
            context.stream = stream;
            dispatch::cuda::lltRcpp(X, result, context);
            cudaStreamDestroy(stream);
        }
    } else {
        dispatch::cpu::lltRcpp(X, result);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericVector urtRcpp(Rcpp::NumericMatrix X) {
    int nfeats = X.nrow();
    int llt = nfeats * (nfeats - 1) / 2;
    Rcpp::NumericVector result(llt);
    if (is_gpu_backend()) {
        propr_context context;
        cudaStream_t stream;
        cudaError_t err = cudaStreamCreate(&stream);
        if (err != cudaSuccess) {
            Rcpp::warning("CUDA stream creation failed for urtRcpp: %s. Falling back to CPU.", cudaGetErrorString(err));
            dispatch::cpu::urtRcpp(X, result);
        } else {
            context.stream = stream;
            dispatch::cuda::urtRcpp(X, result, context);
            cudaStreamDestroy(stream);
        }
    } else {
        dispatch::cpu::urtRcpp(X, result);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::List labRcpp(int nfeats) {
    Rcpp::List result;
     if (is_gpu_backend()) {
        propr_context context;
        cudaStream_t stream;
        cudaError_t err = cudaStreamCreate(&stream);
        if (err != cudaSuccess) {
            Rcpp::warning("CUDA stream creation failed for labRcpp: %s. Falling back to CPU.", cudaGetErrorString(err));
            dispatch::cpu::labRcpp(nfeats, result);
        } else {
            context.stream = stream;
            dispatch::cuda::labRcpp(nfeats, result, context);
            cudaStreamDestroy(stream);
        }
    } else {
        dispatch::cpu::labRcpp(nfeats, result);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix half2mat(Rcpp::NumericVector X) {
    int nfeats = round(sqrt(2 * X.length() + 0.25) + 0.5);
    Rcpp::NumericMatrix result(nfeats, nfeats);
    if (is_gpu_backend()) {
        propr_context context;
        cudaStream_t stream;
        cudaError_t err = cudaStreamCreate(&stream);
        if (err != cudaSuccess) {
            Rcpp::warning("CUDA stream creation failed for half2mat: %s. Falling back to CPU.", cudaGetErrorString(err));
            dispatch::cpu::half2mat(X, result);
        } else {
            context.stream = stream;
            dispatch::cuda::half2mat(X, result, context);
            cudaStreamDestroy(stream);
        }
    } else {
        dispatch::cpu::half2mat(X, result);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix vector2mat(Rcpp::NumericVector X, Rcpp::IntegerVector i, Rcpp::IntegerVector j, int nfeats) {
    Rcpp::NumericMatrix result(nfeats, nfeats);
    if (is_gpu_backend()) {
        propr_context context;
        cudaStream_t stream;
        cudaError_t err = cudaStreamCreate(&stream);
        if (err != cudaSuccess) {
            Rcpp::warning("CUDA stream creation failed for vector2mat: %s. Falling back to CPU.", cudaGetErrorString(err));
            dispatch::cpu::vector2mat(X, i, j, nfeats, result);
        } else {
            context.stream = stream;
            dispatch::cuda::vector2mat(X, i, j, nfeats, result, context);
            cudaStreamDestroy(stream);
        }
    } else {
        dispatch::cpu::vector2mat(X, i, j, nfeats, result);
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
        propr_context context;
        cudaStream_t stream;
        cudaError_t err = cudaStreamCreate(&stream);
        if (err != cudaSuccess) {
             Rcpp::warning("CUDA stream creation failed for ratiosRcpp: %s. Falling back to CPU.", cudaGetErrorString(err));
             dispatch::cpu::ratiosRcpp(X, result);
        } else {
             context.stream = stream;
             dispatch::cuda::ratiosRcpp(X, result, context);
             cudaStreamDestroy(stream);
        }
    } else {
        dispatch::cpu::ratiosRcpp(X, result);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix results2matRcpp(Rcpp::DataFrame results, int n, double diagonal) {
    Rcpp::NumericMatrix result(n, n);
    if (is_gpu_backend()) {
        propr_context context;
        cudaStream_t stream;
        cudaError_t err = cudaStreamCreate(&stream);
        if (err != cudaSuccess) {
            Rcpp::warning("CUDA stream creation failed for results2matRcpp: %s. Falling back to CPU.", cudaGetErrorString(err));
            dispatch::cpu::results2matRcpp(results, n, diagonal, result);
        } else {
            context.stream = stream;
            dispatch::cuda::results2matRcpp(results, n, diagonal, result, context);
            cudaStreamDestroy(stream);
        }
    } else {
        dispatch::cpu::results2matRcpp(results, n, diagonal, result);
    }
    return result;
}
