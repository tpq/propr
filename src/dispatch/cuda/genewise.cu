#include <Rcpp.h>

#include <propr/context.h>

#include <propr/kernels/cuda/dispatch/genewise.cuh>
#include <propr/kernels/cuda/detail/genewise.cuh>
#include <propr/kernels/cuda/traits/genewise.cuh>

#include <propr/utils/common/constants.h>
#include <propr/utils/rcpp/rcpp_checks.h>
#include <propr/utils/rcpp/rcpp_cuda.cuh>
#include <propr/utils/cuda/cuda_checks.h>
#include <propr/utils/profilers/cuda_profiler.cuh>

using namespace Rcpp;
using namespace propr;

namespace {
int sort_end_bit_for_keys(const int num_genes) {
    if (num_genes <= 1) return 1;

    unsigned int v = static_cast<unsigned int>(num_genes - 1);
    int bits = 0;
    while (v > 0) {
        ++bits;
        v >>= 1;
    }
    return bits;
}
}  // namespace

void propr::dispatch::cuda::genewise_connectivity(
    IntegerVector& per_gene_count,
    IntegerVector& per_gene_conn,
    NumericVector& per_gene_wconn,
    NumericVector& per_gene_fdr_sum,
    const IntegerVector& partner,
    const IntegerVector& pair,
    const NumericVector& theta,
    const NumericVector& fdr,
    int num_genes,
    double fdr_thresh,
    propr_context context) {
    using Config = propr::cuda::traits::genewise_connectivity_stats_config;

    if (num_genes < 0) {
        Rcpp::stop("num_genes must be non-negative.");
    }

    const int num_edges = partner.size();
    if (pair.size() != num_edges || theta.size() != num_edges || fdr.size() != num_edges) {
        Rcpp::stop("partner, pair, theta, and fdr must have identical lengths.");
    }

    PROPR_CHECK_VECTOR_SIZE(per_gene_count, num_genes);
    PROPR_CHECK_VECTOR_SIZE(per_gene_conn, num_genes);
    PROPR_CHECK_VECTOR_SIZE(per_gene_wconn, num_genes);
    PROPR_CHECK_VECTOR_SIZE(per_gene_fdr_sum, num_genes);

    if (num_genes == 0) return;

    IntegerVector partner_zero(num_edges);
    IntegerVector pair_zero(num_edges);

    for (int i = 0; i < num_edges; ++i) {
        const int p = partner[i];
        const int q = pair[i];

        if (p == NA_INTEGER || q == NA_INTEGER) {
            Rcpp::stop("partner/pair cannot contain NA values.");
        }
        if (p < 1 || p > num_genes || q < 1 || q > num_genes) {
            Rcpp::stop("partner/pair indices must be in [1, num_genes].");
        }

        partner_zero[i] = p - 1;
        pair_zero[i] = q - 1;
    }

    int* d_partner = nullptr;
    int* d_pair = nullptr;
    float* d_theta = nullptr;
    float* d_fdr = nullptr;

    if (num_edges > 0) {
        d_partner = RcppVectorToDevice<int>(partner_zero, num_edges);
        d_pair = RcppVectorToDevice<int>(pair_zero, num_edges);
        d_theta = RcppVectorToDevice<float>(theta, num_edges);
        d_fdr = RcppVectorToDevice<float>(fdr, num_edges);
    }

    int* d_count = nullptr;
    int* d_conn = nullptr;
    float* d_wconn = nullptr;
    float* d_fdr_sum = nullptr;

    PROPR_CUDA_CHECK(cudaMalloc(&d_count, static_cast<size_t>(num_genes) * sizeof(int)));
    PROPR_CUDA_CHECK(cudaMalloc(&d_conn, static_cast<size_t>(num_genes) * sizeof(int)));
    PROPR_CUDA_CHECK(cudaMalloc(&d_wconn, static_cast<size_t>(num_genes) * sizeof(float)));
    PROPR_CUDA_CHECK(cudaMalloc(&d_fdr_sum, static_cast<size_t>(num_genes) * sizeof(float)));

    PROPR_CUDA_CHECK(cudaMemsetAsync(d_count, 0, static_cast<size_t>(num_genes) * sizeof(int), context.stream));
    PROPR_CUDA_CHECK(cudaMemsetAsync(d_conn, 0, static_cast<size_t>(num_genes) * sizeof(int), context.stream));
    PROPR_CUDA_CHECK(cudaMemsetAsync(d_wconn, 0, static_cast<size_t>(num_genes) * sizeof(float), context.stream));
    PROPR_CUDA_CHECK(cudaMemsetAsync(d_fdr_sum, 0, static_cast<size_t>(num_genes) * sizeof(float), context.stream));

    if (num_edges > 0) {
        const int edges_per_block = Config::THREADS_PER_BLOCK * Config::PAIRS_PER_THREAD;
        const int grid = propr::ceil_div(num_edges, edges_per_block);
        const int sort_end_bit = sort_end_bit_for_keys(num_genes);

        {
            PROPR_PROFILE_CUDA("kernel", context.stream);
            propr::dispatch::cuda::genewise_connectivity_stats<
                Config::THREADS_PER_BLOCK,
                Config::PAIRS_PER_THREAD><<<grid, Config::THREADS_PER_BLOCK, 0, context.stream>>>(
                d_partner,
                d_pair,
                d_theta,
                d_fdr,
                num_edges,
                static_cast<float>(fdr_thresh),
                sort_end_bit,
                d_count,
                d_conn,
                d_wconn,
                d_fdr_sum);
            PROPR_CUDA_CHECK(cudaGetLastError());
            PROPR_STREAM_SYNCHRONIZE(context);
        }
    }

    copyToNumericVector(d_count, per_gene_count, num_genes);
    copyToNumericVector(d_conn, per_gene_conn, num_genes);
    copyToNumericVector(d_wconn, per_gene_wconn, num_genes);
    copyToNumericVector(d_fdr_sum, per_gene_fdr_sum, num_genes);

    PROPR_CUDA_CHECK(cudaFree(d_partner));
    PROPR_CUDA_CHECK(cudaFree(d_pair));
    PROPR_CUDA_CHECK(cudaFree(d_theta));
    PROPR_CUDA_CHECK(cudaFree(d_fdr));

    PROPR_CUDA_CHECK(cudaFree(d_count));
    PROPR_CUDA_CHECK(cudaFree(d_conn));
    PROPR_CUDA_CHECK(cudaFree(d_wconn));
    PROPR_CUDA_CHECK(cudaFree(d_fdr_sum));
}

void propr::dispatch::cuda::genewise_theta_stats(
    NumericVector& out_mean,
    NumericVector& out_median,
    const NumericVector& theta_edges,
    int num_genes,
    propr_context context) {
    using Config = propr::cuda::traits::genewise_theta_stats_config_for<float>;

    if (num_genes < 0) {
        Rcpp::stop("num_genes must be non-negative.");
    }

    PROPR_CHECK_VECTOR_SIZE(out_mean, num_genes);
    PROPR_CHECK_VECTOR_SIZE(out_median, num_genes);

    const size_t expected_edges =
        (num_genes <= 1) ? 0 : (static_cast<size_t>(num_genes) * static_cast<size_t>(num_genes - 1) / 2);
    if (theta_edges.size() != static_cast<R_xlen_t>(expected_edges)) {
        Rcpp::stop("theta_edges length does not match num_genes*(num_genes-1)/2.");
    }

    if (num_genes == 0) return;

    float* d_theta = nullptr;
    float* d_mean = nullptr;
    float* d_median = nullptr;

    d_theta = RcppVectorToDevice<float>(theta_edges, expected_edges);

    PROPR_CUDA_CHECK(cudaMalloc(&d_mean, static_cast<size_t>(num_genes) * sizeof(float)));
    PROPR_CUDA_CHECK(cudaMalloc(&d_median, static_cast<size_t>(num_genes) * sizeof(float)));

    const int block = Config::BLK_X;
    const int grid = num_genes;
    const int cache_cap = Config::CACHE_CAP_VALUES;
    const size_t shared_bytes = static_cast<size_t>(cache_cap) * sizeof(float);

    {
        PROPR_PROFILE_CUDA("kernel", context.stream);
        propr::dispatch::cuda::genewise_theta_stats<float><<<grid, block, shared_bytes, context.stream>>>(
            d_theta,
            num_genes,
            d_mean,
            d_median,
            cache_cap);
        PROPR_CUDA_CHECK(cudaGetLastError());
        PROPR_STREAM_SYNCHRONIZE(context);
    }

    copyToNumericVector(d_mean, out_mean, num_genes);
    copyToNumericVector(d_median, out_median, num_genes);

    PROPR_CUDA_CHECK(cudaFree(d_theta));
    PROPR_CUDA_CHECK(cudaFree(d_mean));
    PROPR_CUDA_CHECK(cudaFree(d_median));
}
