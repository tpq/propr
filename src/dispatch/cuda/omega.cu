#include <Rcpp.h>
#include <propr/utils/rcpp_checks.h>
#include <propr/utils/cuda_checks.h>
#include <propr/utils/rcpp_cuda.cuh>


#include <propr/context.h>
#include <propr/kernels/cuda/dispatch/omega.cuh>
#include <propr/kernels/cuda/detail/omega.cuh>



using namespace Rcpp;
using namespace propr;


void printRMatrix(Rcpp::NumericMatrix mat) {
  int nrow = mat.nrow();
  int ncol = mat.ncol();
  for (int i = 0; i < nrow; ++i) {
    for (int j = 0; j < ncol; ++j) {
      std::cout << mat(i, j) << " ";
    } std::cout << std::endl;
  } std::cout << std::endl;
}

void 
propr::dispatch::cuda::dof_global(NumericVector& out, const NumericMatrix& W, propr_context context) {
    using Converter = cutlass::NumericConverter<cute::half_t, float>;
    using Config     = kernels::cutlass_impl::OmegaConfig;

    int nfeats  = W.ncol();
    int samples = W.nrow();
    
    std::cout << samples << std::endl;
    std::cout << nfeats  << std::endl;

    NumericMatrix Wl(64,4);

    for (int i = 0; i < Wl.nrow(); ++i) {
      for (int j = 0; j < Wl.ncol(); ++j) {
        Wl(i, j) = float( 1.0f);
      }
    }

    nfeats  = Wl.ncol();
    samples = Wl.nrow();

    const size_t llt = static_cast<size_t>(nfeats) * (nfeats - 1) / 2;
    // CHECK_VECTOR_SIZE(out, llt);
    int alignment = 16; //16
    offset_t W_stride;
    auto *d_W = RcppMatrixToDevice<cute::half_t, REALSXP>(Wl, W_stride, alignment);

    cute::half_t* d_out = nullptr;
    offset_t dout_stride = ((nfeats  + alignment - 1) / alignment) * alignment;
    CUDA_CHECK(cudaMalloc(&d_out, nfeats * dout_stride * sizeof(cute::half_t)));

    const int BX = (nfeats + Config::BLK_M - 1) / Config::BLK_M;
    const int BY = (nfeats + Config::BLK_M - 1) / Config::BLK_M;

    dim3 block(cute::size(Config::TiledMMA{}));
    dim3 grid(BX, BY);
    static constexpr int shm_size_AB = cute::cosize(Config::SmemLayoutA{}) + cute::cosize(Config::SmemLayoutB{});
    static constexpr int kShmSize = shm_size_AB * sizeof(__half);
    int shm_size = kShmSize;
    const auto fptr = propr::kernels::cutlass_impl::omega_kernel<Config,cute::half_t>;

    cudaFuncSetAttribute(fptr, cudaFuncAttributeMaxDynamicSharedMemorySize, shm_size);
    cudaFuncSetAttribute(fptr, cudaFuncAttributePreferredSharedMemoryCarveout, 100);
    fptr<<<grid, block, shm_size, context.stream>>>( nfeats, 
                                                     samples, 
                                                     d_W, W_stride, d_out, dout_stride);
    CUDA_CHECK(cudaStreamSynchronize(context.stream));

    auto h_full= new std::vector<cute::half_t> (nfeats * dout_stride);
    CUDA_CHECK(cudaMemcpy(
        h_full->data(),
        d_out,
        nfeats * dout_stride * sizeof(cute::half_t),
        cudaMemcpyDeviceToHost
    ));

    Converter cast_away;
    size_t counter = 0;
    double* out_ptr = REAL(out);
    std::cout << "[GPU]: ";
    for (int i = 1; i < nfeats; ++i) {
        for (int j = 0; j < i; ++j) {
            cute::half_t v = h_full->at(size_t(i) * dout_stride + j);
            // out_ptr[counter++] = static_cast<double>(cast_away(v));
            std::cout << cast_away(v) << " ";
        }
    }
    std::cout << std::endl;

    CUDA_CHECK(cudaFree(d_W));
    CUDA_CHECK(cudaFree(d_out));
}
 
void 
propr::dispatch::cuda::dof_population(NumericVector& out, const NumericMatrix& W, propr_context context) {
    // using Converter = cutlass::NumericConverter<cute::half_t, float>;
    // using Config     = kernels::cutlass_impl::OmegaConfig;

    // int nfeats   = W.ncol();
    // int samples = W.nrow();

    // const size_t llt = static_cast<size_t>(nfeats) * (nfeats - 1) / 2;
    // CHECK_VECTOR_SIZE(out, llt);

    // offset_t W_stride;
    // auto *d_W = RcppMatrixToDevice<cute::half_t, REALSXP, true>(W, W_stride);

    // cute::half_t* d_out = nullptr;
    // CUDA_CHECK(cudaMalloc(&d_out, nfeats * nfeats * sizeof(cute::half_t)));

    // const int BX = (nfeats + Config::BLK_M - 1) / Config::BLK_M;
    // const int BY = (nfeats + Config::BLK_M - 1) / Config::BLK_M;

    // dim3 block(cute::size(Config::TiledMMA{}));
    // dim3 grid(BX, BY);
    // static constexpr int shm_size_AB = cute::cosize(Config::SmemLayoutA{}) + cute::cosize(Config::SmemLayoutB{});
    // static constexpr int kShmSize = shm_size_AB * sizeof(__half);
    // int shm_size = kShmSize;
    // const auto fptr = propr::kernels::cutlass_impl::omega_kernel<Config,cute::half_t>;

    // cudaFuncSetAttribute(fptr, cudaFuncAttributeMaxDynamicSharedMemorySize, shm_size);
    // cudaFuncSetAttribute(fptr, cudaFuncAttributePreferredSharedMemoryCarveout, 100);
    // fptr<<<grid, block, shm_size, context.stream>>>(d_W, d_out, samples, W_stride);
    // CUDA_CHECK(cudaStreamSynchronize(context.stream));

    // std::vector<cute::half_t> h_full(nfeats * nfeats);
    // CUDA_CHECK(cudaMemcpy(
    //     h_full.data(),
    //     d_out,
    //     W_stride * W_stride * sizeof(cute::half_t),
    //     cudaMemcpyDeviceToHost
    // ));
    // Converter cast_away;
    // size_t counter = 0;
    // double* out_ptr = REAL(out);
    // for (int i = 1; i < nfeats; ++i) {
    //     for (int j = 0; j < i; ++j) {
    //         cute::half_t v = h_full[size_t(i) * W_stride + j];
    //         out_ptr[counter++] = static_cast<double>(cast_away(v));
    //     }
    // }

    // CUDA_CHECK(cudaFree(d_W));
    // CUDA_CHECK(cudaFree(d_out));
}