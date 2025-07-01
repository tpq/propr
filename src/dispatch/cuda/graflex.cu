#include <array>

#include <cooperative_groups.h>

#include <cub/device/device_scan.cuh>

#include <thrust/device_ptr.h>
#include <thrust/execution_policy.h>
#include <thrust/transform_reduce.h>
#include <thrust/device_vector.h>

#include <propr/utils.hpp>
#include <propr/kernels/cuda/detail/graflex.cuh>
#include <propr/kernels/cuda/dispatch/graflex.cuh>



using namespace Rcpp;
using namespace propr;

void
dispatch::cuda::getOR(NumericVector& out, 
                      const IntegerMatrix& A, 
                      const IntegerMatrix& G, 
                      propr::propr_context context) {

    using scan_tile_state_t = cub::ScanTileState<int4>;
    
    const int n = A.ncol();

    int g_stride; int* d_G; d_G = RcppMatrixToDevice<int>(G, g_stride);
    int a_stride; int* d_A; d_A = RcppMatrixToDevice<int>(A, a_stride);

    const int numPairs      = (n * (n - 1)) / 2;
    const int blockSize     = 512;
    const int numBlocks     = cub::DivideAndRoundUp(numPairs, blockSize);
    const int numInitBlocks = cub::DivideAndRoundUp(numBlocks, blockSize);

    std::size_t device_partials_size{};
    scan_tile_state_t::AllocationSize(numBlocks, device_partials_size);

    std::size_t aligned_size = ((device_partials_size + 15) / 16) * 16;
    void *d_temp_storage;
    CUDA_CHECK(cudaMalloc(&d_temp_storage, aligned_size));
    
    scan_tile_state_t tile_state;
    tile_state.Init(
        numBlocks,
        d_temp_storage,
        device_partials_size
    );

    int4 *d_acc;
    CUDA_CHECK(cudaMalloc(&d_acc, sizeof(int4)));
    CUDA_CHECK(cudaMemset(d_acc, 0, sizeof(int4)));

    void *kernelArgs[] = {
        (void *)&tile_state,
        (void *)&d_A, 
        (void *)&a_stride,
        (void *)&d_G, 
        (void *)&g_stride,
        (void *)&n,
        (void *)&d_acc
    };
    cudaLaunchCooperativeKernel((void *)detail::cuda::compute_odd_ratio, numBlocks, blockSize, kernelArgs, 0, NULL);
    CUDA_CHECK(cudaDeviceSynchronize());
    CUDA_CHECK(cudaGetLastError());

    int4 h_acc;
    CUDA_CHECK(cudaMemcpy(&h_acc, d_acc, sizeof(int4), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaFree(d_acc));

    int h_a = h_acc.x;
    int h_b = h_acc.y;
    int h_c = h_acc.z;
    int h_d = h_acc.w;

    double odds_ratio = static_cast<double>(h_a * h_d) / (h_b * h_c);
    double log_odds_ratio = std::log(odds_ratio);
    out[0] = h_a;
    out[1] = h_b;
    out[2] = h_c;
    out[3] = h_d;
    out[4] = odds_ratio;
    out[5] = log_odds_ratio;
    out[6] = R_NaN;
    out[7] = R_NaN;
    
    CUDA_CHECK(cudaFree(d_G));
    CUDA_CHECK(cudaFree(d_A));
    CUDA_CHECK(cudaFree(d_temp_storage));
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
    const int n = permuted.size();
    double* d_permuted;
    CUDA_CHECK(cudaMalloc(&d_permuted, n * sizeof(double)));
    CUDA_CHECK(cudaMemcpy(
        d_permuted, permuted.begin(), n * sizeof(double),
        cudaMemcpyHostToDevice
    ));

    thrust::device_ptr<double> dev_permuted(d_permuted);
    auto policy = thrust::cuda::par.on(context.stream);
    constexpr int2 init{0,0};
    int2 result = thrust::transform_reduce(
        policy,
        dev_permuted,  dev_permuted + n,
        [=] __host__ __device__ (double x) {
            return int2{ x >= actual, x <= actual };
        },
        init,
        [] __host__ __device__ (const int2& a, const int2& b) {
            return int2{ a.x + b.x, a.y + b.y };
        }
    );

    CUDA_CHECK(cudaStreamSynchronize(context.stream));
    CUDA_CHECK(cudaFree(d_permuted));

    double fdr_over = static_cast<double>(result.x) / n;
    double fdr_under = static_cast<double>(result.y) / n;
    out["over"] = fdr_over;
    out["under"] = fdr_under;
}

void
dispatch::cuda::getG(IntegerMatrix& out, const IntegerVector& Gk, propr::propr_context context) {
    Rcpp::stop("Not implemented in CUDA for getG.");
}

void
dispatch::cuda::graflex(NumericVector& out, const IntegerMatrix& A, const IntegerVector& Gk, int p, propr::propr_context context) {
    Rcpp::stop("Not implemented in CUDA for graflex.");   
}