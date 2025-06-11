#include <Rcpp.h>
#include <catch2/catch.hpp>

#include <propr/context.h>
#include <propr/kernels/cpu/dispatch/lrv.hpp>
#include <propr/kernels/cuda/dispatch/lrv.cuh>

#include "test_utils.hpp"

TEST_CASE("CUDA Kernels - lrv_weighted", "[cuda][lrv]") {
    Rcpp::NumericMatrix raw_counts = generate_counts(50, 20, 1, 5000, 0.01);
    Rcpp::NumericMatrix Y          = apply_zero_replacement(raw_counts);
    Rcpp::NumericMatrix W          = generate_weights_matrix(Y.nrow(), Y.ncol(), 0.01, 1.0);

    int nfeats = Y.ncol();
    int llt = nfeats * (nfeats - 1) / 2;
    Rcpp::NumericVector result_cuda(llt);

    propr::dispatch::cuda::lrv_weighted(result_cuda, Y, W);

    Rcpp::NumericVector cpu_result(llt);
    propr::dispatch::cpu::lrv(cpu_result, Y, W, true);

    REQUIRE(result_cuda.length() == llt);
    for (int i = 0; i < llt; ++i) {
        REQUIRE(cpu_result[i] == result_cuda[i]);
    }
}

TEST_CASE("CUDA Kernels - lrv_alpha", "[cuda][lrv]") {
    Rcpp::NumericMatrix Y     = generate_counts(50, 20, 1, 5000, 0.0);
    Rcpp::NumericMatrix Yfull = generate_counts(100, 20, 1, 5000, 0.0);
    double alpha = 0.5;

    int nfeats = Y.ncol();
    int llt = nfeats * (nfeats - 1) / 2;
    Rcpp::NumericVector result_cuda(llt);

    propr::propr_context context = propr::DEFAULT_GLOBAL_CONTEXT;
    propr::dispatch::cuda::lrv_alpha(result_cuda, Y, alpha, Yfull, context);

    Rcpp::NumericVector cpu_result(llt);
    Rcpp::NumericMatrix W_dummy(Y.nrow(), Y.ncol());
    propr::dispatch::cpu::lrv(cpu_result, Y, W_dummy, false, alpha, Yfull);

    REQUIRE(result_cuda.length() == llt);
    for (int i = 0; i < llt; ++i) {
        REQUIRE(result_cuda[i] == cpu_result[i]);
    }
}

TEST_CASE("CUDA Kernels - lrv_alpha_weighted matches CPU results with realistic data", "[cuda][lrv]") {
    Rcpp::NumericMatrix Y     = generate_counts(50, 20, 1, 5000, 0.0);
    Rcpp::NumericMatrix W     = generate_weights_matrix(Y.nrow(), Y.ncol(), 0.01, 1.0);
    Rcpp::NumericMatrix Yfull = generate_counts(100, 20, 1, 5000, 0.0);
    Rcpp::NumericMatrix Wfull = generate_weights_matrix(Yfull.nrow(), Yfull.ncol(), 0.01, 1.0);
    double alpha = 0.5;

    int nfeats = Y.ncol();
    int llt = nfeats * (nfeats - 1) / 2;
    Rcpp::NumericVector result_cuda(llt);

    propr::dispatch::cuda::lrv_alpha_weighted(result_cuda, Y, W, alpha, Yfull, Wfull);

    Rcpp::NumericVector cpu_result(llt);
    propr::dispatch::cpu::lrv(cpu_result, Y, W, true, alpha, Yfull, Wfull);

    REQUIRE(result_cuda.length() == llt);
    for (int i = 0; i < llt; ++i) {
        REQUIRE(result_cuda[i] == cpu_result[i]);
    }
}
