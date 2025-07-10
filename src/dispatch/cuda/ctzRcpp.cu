#include <Rcpp.h>

#include <propr/data/types.h>
#include <propr/utils/rcpp_checks.h>
#include <propr/utils/cuda_checks.h>
#include <propr/utils/rcpp_cuda.cuh>

#include <propr/kernels/cuda/dispatch/ctzRcpp.cuh>
#include <propr/kernels/cuda/detail/ctz.cuh>

using namespace Rcpp;
using namespace propr;

void
dispatch::cuda::ctzRcpp(NumericVector& out,
                        NumericMatrix& X,
                        propr_context context) {
    int nfeats = X.ncol();
    int nsubjs = X.nrow();
    int llt    = nfeats * (nfeats - 1) / 2;
    CHECK_VECTOR_SIZE(out, llt);

    offset_t X_stride; float* d_X = RcppMatrixToDevice<float>(X, X_stride);

    int* d_zeroes;
    CUDA_CHECK(cudaMalloc(&d_zeroes, nfeats * sizeof(int)));

    const int BLK = 256;
    int grid1 = nfeats;
    detail::cuda::count_per_feature<BLK><<<grid1, BLK, 0, context.stream>>>(
        d_X, X_stride, nsubjs, nfeats, d_zeroes
    );
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaStreamSynchronize(context.stream));

    int* d_result;
    CUDA_CHECK(cudaMalloc(&d_result, llt * sizeof(int)));

    dim3 blockDim2(16, 16);
    dim3 gridDim2((nfeats + blockDim2.x - 1) / blockDim2.x,
                  (nfeats + blockDim2.y - 1) / blockDim2.y);
    detail::cuda::count_joint_zeros<<<gridDim2, blockDim2, 0, context.stream>>>(
        d_zeroes, 1, nfeats, d_result
    );
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaStreamSynchronize(context.stream));

    copyToNumericVector(d_result, out, llt);

    CUDA_CHECK(cudaFree(d_X));
    CUDA_CHECK(cudaFree(d_zeroes));
    CUDA_CHECK(cudaFree(d_result));
}
