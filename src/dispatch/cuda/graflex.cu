#include <array>

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

    using scan_tile_state_t = cub::ScanTileState<uint4>;

    const int n = A.ncol();

    int a_stride; unsigned char* d_A; d_A = RcppMatrixToDevice<unsigned char,INTSXP, false>(A, a_stride,1);
    int g_stride; unsigned char* d_G; d_G = RcppMatrixToDevice<unsigned char,INTSXP, false>(G, g_stride,1);

    const int numPairs      = (n * (n - 1)) / 2;
    const int blockSize     = 1024;
    const int numBlocks     = cub::DivideAndRoundUp(numPairs, blockSize);
    const int numInitBlocks = cub::DivideAndRoundUp(numBlocks, blockSize);
        
    std::size_t device_partials_size{};
    scan_tile_state_t::AllocationSize(numBlocks, device_partials_size);

    thrust::device_vector<std::uint8_t> temp_storage(device_partials_size);
    std::uint8_t *d_temp_storage = thrust::raw_pointer_cast(temp_storage.data());
    
    scan_tile_state_t tile_state;
    tile_state.Init(
        numBlocks,
        d_temp_storage,
        device_partials_size
    );

    uint4 *d_acc;
    CUDA_CHECK(cudaMalloc(&d_acc, sizeof(uint4)));
    CUDA_CHECK(cudaMemset(d_acc, 0, sizeof(uint4)));

    detail::cuda::compute_odd_ratio_init<<<numInitBlocks, blockSize>>>(
        tile_state, numBlocks
    );
    CUDA_CHECK(cudaDeviceSynchronize());
    CUDA_CHECK(cudaGetLastError());
    
    detail::cuda::compute_odd_ratio<<<numBlocks, blockSize>>>(
        tile_state,
        d_A, a_stride,
        d_G, g_stride,
        n,
        d_acc
    );
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    uint4 h_acc;
    CUDA_CHECK(cudaMemcpy(&h_acc, d_acc, sizeof(uint4), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaFree(d_acc));

    int h_a = h_acc.x;
    int h_b = h_acc.y;
    int h_c = h_acc.z;
    int h_d = h_acc.w;

    double odds_ratio = static_cast<double>(double(h_a) / double(h_b)) *(double(h_d) / double(h_c));
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