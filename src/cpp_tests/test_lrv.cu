#include <Rcpp.h>

#include <propr/context.h>
#include <propr/kernels/cpu/dispatch/lrv.hpp>
#include <propr/kernels/cuda/dispatch/lrv.cuh>
#include <testthat/testthat.h>


context("lrm_test") {

    test_that("Stupid case") {
        expect_true(1 == 1);
    }

    test_that("CUDA Kernels - lrv_alpha [cuda][lrv]") {
       /**  
        int np = 0;

        // protect Y and Yfull
        SEXP Y_sxp     = PROTECT(Rcpp::wrap(generate_counts(50, 20, 1, 5000, 0.0)));  np++;
        SEXP Yfull_sxp = PROTECT(Rcpp::wrap(generate_counts(100, 20, 1, 5000, 0.0))); np++;
        Rcpp::NumericMatrix Y     = Rcpp::as<Rcpp::NumericMatrix>(Y_sxp);
        Rcpp::NumericMatrix Yfull = Rcpp::as<Rcpp::NumericMatrix>(Yfull_sxp);

        double alpha = 0.5;
        int nfeats = Y.ncol();
        int llt    = nfeats * (nfeats - 1) / 2;

        // protect result vectors
        SEXP result_cuda_sxp = PROTECT(Rcpp::wrap(Rcpp::NumericVector(llt))); np++;
        Rcpp::NumericVector result_cuda = Rcpp::as<Rcpp::NumericVector>(result_cuda_sxp);

        propr::propr_context context = propr::DEFAULT_GLOBAL_CONTEXT;
        propr::dispatch::cuda::lrv_alpha(result_cuda, Y, alpha, Yfull, context);

        SEXP cpu_result_sxp = PROTECT(Rcpp::wrap(Rcpp::NumericVector(llt))); np++;
        Rcpp::NumericVector cpu_result = Rcpp::as<Rcpp::NumericVector>(cpu_result_sxp);

        // dummy weights matrix, also protect
        SEXP W_dummy_sxp = PROTECT(Rcpp::wrap(Rcpp::NumericMatrix(Y.nrow(), Y.ncol()))); np++;
        Rcpp::NumericMatrix W_dummy = Rcpp::as<Rcpp::NumericMatrix>(W_dummy_sxp);

        propr::dispatch::cpu::lrv(cpu_result, Y, W_dummy, false, alpha, Yfull);

        expect_true(result_cuda.length() == llt);
        for (int i = 0; i < llt; ++i) {
            expect_true(result_cuda[i] == cpu_result[i]);
        }

        // unprotect all
        UNPROTECT(np);
    }

    test_that("CUDA Kernels - lrv_alpha_weighted matches CPU results with realistic data [cuda][lrv]") {
        int np = 0;

        SEXP Y_sxp     = PROTECT(Rcpp::wrap(generate_counts(50, 20, 1, 5000, 0.0))); np++;
        SEXP W_sxp     = PROTECT(Rcpp::wrap(generate_weights_matrix(50, 20, 0.01, 1.0))); np++;
        SEXP Yfull_sxp = PROTECT(Rcpp::wrap(generate_counts(100, 20, 1, 5000, 0.0))); np++;
        SEXP Wfull_sxp = PROTECT(Rcpp::wrap(generate_weights_matrix(100, 20, 0.01, 1.0))); np++;

        Rcpp::NumericMatrix Y     = Rcpp::as<Rcpp::NumericMatrix>(Y_sxp);
        Rcpp::NumericMatrix W     = Rcpp::as<Rcpp::NumericMatrix>(W_sxp);
        Rcpp::NumericMatrix Yfull = Rcpp::as<Rcpp::NumericMatrix>(Yfull_sxp);
        Rcpp::NumericMatrix Wfull = Rcpp::as<Rcpp::NumericMatrix>(Wfull_sxp);

        double alpha = 0.5;
        int nfeats = Y.ncol();
        int llt    = nfeats * (nfeats - 1) / 2;

        SEXP result_cuda_sxp = PROTECT(Rcpp::wrap(Rcpp::NumericVector(llt))); np++;
        Rcpp::NumericVector result_cuda = Rcpp::as<Rcpp::NumericVector>(result_cuda_sxp);

        propr::dispatch::cuda::lrv_alpha_weighted(result_cuda, Y, W, alpha, Yfull, Wfull);

        SEXP cpu_result_sxp = PROTECT(Rcpp::wrap(Rcpp::NumericVector(llt))); np++;
        Rcpp::NumericVector cpu_result = Rcpp::as<Rcpp::NumericVector>(cpu_result_sxp);

        propr::dispatch::cpu::lrv(cpu_result, Y, W, true, alpha, Yfull, Wfull);

        expect_true(result_cuda.length() == llt);
        for (int i = 0; i < llt; ++i) {
            expect_true(result_cuda[i] == cpu_result[i]);
        }

        UNPROTECT(np);
        **/
    }
}