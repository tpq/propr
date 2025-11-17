#include <Rcpp.h>

#include <propr/data/types.h>
#include <propr/utils/rcpp_checks.h>
#include <propr/utils/cuda_checks.h>
#include <propr/utils/rcpp_cuda.cuh>
#include <propr/utils/cuda_helpers.cuh>

#include <propr/kernels/cuda/dispatch/ctzRcpp.cuh>
#include <propr/kernels/cuda/detail/ctz.cuh>
#include <propr/kernels/cuda/traits/ctz.cuh>


using namespace Rcpp;
using namespace propr;

void
dispatch::cuda::ctzRcpp(NumericVector& out,
                        NumericMatrix& X,
                        propr_context context) 
{
    using Config = propr::cuda::traits::ctzRcpp_config;

    int nfeats = X.ncol();
    int nsubjs = X.nrow();
    size_t llt = size_t(nfeats) * (nfeats - 1) / 2;
    PROPR_CHECK_VECTOR_SIZE(out, llt);

    offset_t X_stride; float* d_X = RcppMatrixToDevice<float>(X, X_stride);

    int* d_zeroes;
    PROPR_CUDA_CHECK(cudaMalloc(&d_zeroes, nfeats * sizeof(int)));

    const int grid1 = nfeats;
    detail::cuda::count_per_feature<Config::PHASE_ONE_BLK_X><<<grid1, Config::PHASE_ONE_BLK_X, 0, context.stream>>>(
        d_X, X_stride, nsubjs, nfeats, d_zeroes
    );
    PROPR_CUDA_CHECK(cudaGetLastError());
    PROPR_STREAM_SYNCHRONIZE(context);

    int* d_result;
    PROPR_CUDA_CHECK(cudaMalloc(&d_result, llt * sizeof(int)));

    dim3 blockDim2(Config::PHASE_TWO_BLK_X, Config::PHASE_TWO_BLK_Y);
    dim3 gridDim2(propr::ceil_div(nfeats, Config::PHASE_TWO_BLK_X), propr::ceil_div(nfeats, Config::PHASE_TWO_BLK_Y));
    
    detail::cuda::count_joint_zeros<<<gridDim2, blockDim2, 0, context.stream>>>(
        d_zeroes, 1, nfeats, d_result
    );
    PROPR_CUDA_CHECK(cudaGetLastError());
    PROPR_STREAM_SYNCHRONIZE(context);

    copyToNumericVector(d_result, out, llt);

    PROPR_CUDA_CHECK(cudaFree(d_X));
    PROPR_CUDA_CHECK(cudaFree(d_zeroes));
    PROPR_CUDA_CHECK(cudaFree(d_result));
}
