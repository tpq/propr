#include <array>

#include <cub/device/device_scan.cuh>

#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/transform_reduce.h>

#include <propr/data/types.h>
#include <propr/utils/cuda_checks.h>
#include <propr/utils/rcpp_checks.h>
#include <propr/utils/rcpp_cuda.cuh>

#include <propr/kernels/cuda/detail/graflex.cuh>
#include <propr/kernels/cuda/dispatch/graflex.cuh>
#include <propr/kernels/cuda/traits/graflex.cuh>



using namespace Rcpp;


void
propr::dispatch::cuda::getOR(NumericVector& out, 
                      const IntegerMatrix& A, 
                      const IntegerMatrix& G, 
                      propr::propr_context context) 
{
    using Config = propr::cuda::traits::getOR_config;

    constexpr int MAX_PAIRS = std::numeric_limits<std::int32_t>::max();

    using scan_tile_state_t = cub::ScanTileState<uint4>;

    const int n = A.ncol();

    offset_t a_stride; unsigned char* d_A; d_A = RcppMatrixToDevice<unsigned char,INTSXP, false>(A, a_stride,1);
    offset_t g_stride; unsigned char* d_G; d_G = RcppMatrixToDevice<unsigned char,INTSXP, false>(G, g_stride,1);

    const std::int64_t Nmax = static_cast<std::int64_t>(std::floor(std::sqrt(2.0 * MAX_PAIRS)));
    const std::int32_t k = (n + Nmax - 1) / Nmax; // elements per thread

    const int numPairs      = int32_t(((int64_t(n) * (n - 1)) / 2) / k);
    const int numBlocks     = cub::DivideAndRoundUp(numPairs, Config::BLK_X);
    const int numInitBlocks = cub::DivideAndRoundUp(numBlocks, Config::BLK_X);
        
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
    PROPR_CUDA_CHECK(cudaMalloc(&d_acc, sizeof(uint4)));
    PROPR_CUDA_CHECK(cudaMemset(d_acc, 0, sizeof(uint4)));

    detail::cuda::compute_odd_ratio_init<<<numInitBlocks, Config::BLK_X>>>(
        tile_state, numBlocks
    );
    PROPR_CUDA_CHECK(cudaDeviceSynchronize());
    PROPR_CUDA_CHECK(cudaGetLastError());
    
    detail::cuda::compute_odd_ratio<Config::BLK_X><<<numBlocks, Config::BLK_X>>>(
        tile_state,
        d_A, a_stride,
        d_G, g_stride,
        n,
        d_acc
    );
    PROPR_CUDA_CHECK(cudaGetLastError());
    PROPR_CUDA_CHECK(cudaDeviceSynchronize());

    uint4 h_acc;
    PROPR_CUDA_CHECK(cudaMemcpy(&h_acc, d_acc, sizeof(uint4), cudaMemcpyDeviceToHost));
    PROPR_CUDA_CHECK(cudaFree(d_acc));

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
    
    PROPR_CUDA_CHECK(cudaFree(d_G));
    PROPR_CUDA_CHECK(cudaFree(d_A));
}

void
propr::dispatch::cuda::getORperm(NumericVector& out, const IntegerMatrix& A, const IntegerMatrix& G, const IntegerVector& perm, propr::propr_context context) {
    using Config = propr::cuda::traits::getORperm_config;
    constexpr int MAX_PAIRS = std::numeric_limits<std::int32_t>::max();

    using scan_tile_state_t = cub::ScanTileState<uint4>;
    //TODO: donot forget to consider the number of elements per thread
    const int n = A.ncol();

    offset_t a_stride; unsigned char* d_A = RcppMatrixPermToDevice<unsigned char, INTSXP, false>(A, perm,a_stride,1);
    offset_t g_stride; unsigned char* d_G = RcppMatrixToDevice<unsigned char, INTSXP, false>(G, g_stride,1);

    
    const std::int64_t Nmax = static_cast<std::int64_t>(std::floor(std::sqrt(2.0 * MAX_PAIRS)));
    const std::int32_t k = (n + Nmax - 1) / Nmax; // elements per thread

    const int numPairs      = int32_t(((int64_t(n) * (n - 1)) / 2) / k);
    const int numBlocks     = cub::DivideAndRoundUp(numPairs, Config::BLK_X);
    const int numInitBlocks = cub::DivideAndRoundUp(numBlocks, Config::BLK_X);
        
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
    PROPR_CUDA_CHECK(cudaMalloc(&d_acc, sizeof(uint4)));
    PROPR_CUDA_CHECK(cudaMemset(d_acc, 0, sizeof(uint4)));

    detail::cuda::compute_odd_ratio_init<<<numInitBlocks, Config::BLK_X>>>(
        tile_state, numBlocks
    );
    PROPR_CUDA_CHECK(cudaDeviceSynchronize());
    PROPR_CUDA_CHECK(cudaGetLastError());
    
    detail::cuda::compute_odd_ratio<Config::BLK_X><<<numBlocks, Config::BLK_X>>>(
        tile_state,
        d_A, a_stride,
        d_G, g_stride,
        n,
        d_acc
    );
    PROPR_CUDA_CHECK(cudaGetLastError());
    PROPR_CUDA_CHECK(cudaDeviceSynchronize());

    uint4 h_acc;
    PROPR_CUDA_CHECK(cudaMemcpy(&h_acc, d_acc, sizeof(uint4), cudaMemcpyDeviceToHost));
    PROPR_CUDA_CHECK(cudaFree(d_acc));

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
    
    PROPR_CUDA_CHECK(cudaFree(d_G));
    PROPR_CUDA_CHECK(cudaFree(d_A));
}

void
propr::dispatch::cuda::permuteOR(NumericMatrix& out, const IntegerMatrix& A, const IntegerMatrix& G, int p, propr::propr_context context) {
    PROPR_CHECK_MATRIX_DIMS(out, p, 8);
    int ncol = A.ncol();
    NumericVector or_tmp(8);
    for (int i = 0; i < p; ++i) {
        IntegerVector perm = sample(ncol, ncol, false) - 1;
        propr::dispatch::cuda::getORperm(or_tmp, A, G, perm);
        for (int j = 0; j < 8; ++j) {
            out(i, j) = or_tmp[j];
        }
    }
}

void
propr::dispatch::cuda::getFDR(List& out, double actual, const NumericVector& permuted, propr::propr_context context) {
    const int n = permuted.size();
    double* d_permuted;
    PROPR_CUDA_CHECK(cudaMalloc(&d_permuted, n * sizeof(double)));
    PROPR_CUDA_CHECK(cudaMemcpy(
        d_permuted, permuted.begin(), n * sizeof(double),
        cudaMemcpyHostToDevice
    ));

    auto policy = thrust::cuda::par.on(context.stream);
    constexpr int2 init{0,0};
    int2 result = thrust::transform_reduce(
        policy,
        d_permuted,  d_permuted + n,
        [=] __host__ __device__ (double x) {
            return int2{ x >= actual, x <= actual };
        },
        init,
        [] __host__ __device__ (const int2& a, const int2& b) {
            return int2{ a.x + b.x, a.y + b.y };
        }
    );

    PROPR_CUDA_CHECK(cudaStreamSynchronize(context.stream));
    PROPR_CUDA_CHECK(cudaFree(d_permuted));

    double fdr_over = static_cast<double>(result.x) / n;
    double fdr_under = static_cast<double>(result.y) / n;
    out["over"] = fdr_over;
    out["under"] = fdr_under;
}

void
propr::dispatch::cuda::getG(IntegerMatrix& out, const IntegerVector& Gk, propr::propr_context context) {
    using T  = char;
    size_t n = Gk.size();
    PROPR_CHECK_MATRIX_DIMS(out, n, n);

    T   *d_G = RcppVectorToDevice<T, INTSXP>(Gk, n);
    int *d_C  = nullptr;
    PROPR_CUDA_CHECK(cudaMalloc(&d_C, n * n * sizeof(*d_C)));

    cublasHandle_t handle;
    cublasCreate(&handle);
    float alpha = 1.0f, beta = 0.0f;

    // C = G * G^T
    PROPR_CUBLAS_CHECK(cublasGemmEx(handle,
                 CUBLAS_OP_N, CUBLAS_OP_T,
                 n, n, 1,
                 &alpha,
                 d_G, CUDA_R_8I, n,
                 d_G, CUDA_R_8I, n,
                 &beta,
                 d_C, CUDA_R_32F, n,
                 CUBLAS_COMPUTE_32F,
                 CUBLAS_GEMM_DEFAULT));

    
    std::vector<float> h_out;
    h_out.resize(static_cast<size_t>(n) * n);
    PROPR_CUDA_CHECK(cudaMemcpy(h_out.data(), d_C, n*n*sizeof(*d_C), cudaMemcpyDeviceToHost));
    for (int idx = 0; idx < n * n; ++idx) out(idx / n, idx % n) = int(h_out[idx]);
    cublasDestroy(handle);
    cudaFree(d_G);
    cudaFree(d_C);
}

void
propr::dispatch::cuda::graflex(NumericVector& out, const IntegerMatrix& A, const IntegerVector& Gk, int p, propr::propr_context context) {
    IntegerMatrix G_tmp(Gk.size(), Gk.size());
    propr::dispatch::cuda::getG(G_tmp, Gk);

    NumericVector actual_tmp(8);
    propr::dispatch::cuda::getOR(actual_tmp, A, G_tmp);

    PROPR_CHECK_VECTOR_SIZE(out, actual_tmp.length());
    for (int i = 0; i < actual_tmp.length(); ++i) {
        out[i] = actual_tmp[i];
    }

    if (!std::isnan(actual_tmp(4))) {
        NumericMatrix permuted_tmp(p, 8);
        propr::dispatch::cuda::permuteOR(permuted_tmp, A, G_tmp, p);

        List fdr_tmp; 
        propr::dispatch::cuda::getFDR(fdr_tmp, actual_tmp(4), permuted_tmp(_, 4));
        out(6) = as<double>(fdr_tmp["under"]);
        out(7) = as<double>(fdr_tmp["over"]);
    }   
}