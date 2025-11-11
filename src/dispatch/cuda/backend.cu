#include <Rcpp.h>
#include <math.h>

#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/system/cuda/execution_policy.h>

#include <propr/data/types.h>

#include <propr/kernels/cuda/traits/backend.cuh>
#include <propr/kernels/cuda/dispatch/backend.cuh>
#include <propr/kernels/cuda/detail/backend.cuh>


#include <propr/utils/rcpp_cuda.cuh>
#include <propr/utils/rcpp_checks.h>
#include <propr/utils/rcpp_helpers.h>
#include <propr/utils/cuda_checks.h>
#include <propr/utils/cuda_helpers.cuh>



using namespace Rcpp;
using namespace propr;

void 
dispatch::cuda::wtmRcpp(double& out, const NumericVector& x, const NumericVector& w, propr_context context){
  const int BLK = 1024;
  PROPR_CHECK_VECTOR_SIZE(x, w.size());

  const int n = x.size();
  float* d_x = RcppVectorToDevice<float>(x, n);
  float* d_w = RcppVectorToDevice<float>(w, n);
  
  float h_mean  = 0;
  float *d_mean = nullptr;
  PROPR_CUDA_CHECK(cudaMalloc(&d_mean, sizeof(float)));
  detail::cuda::wtm<BLK><<<1, BLK, 0, context.stream>>>(d_mean, d_x,d_w, n);
  PROPR_CUDA_CHECK(cudaStreamSynchronize(context.stream));
  PROPR_CUDA_CHECK(cudaMemcpy(&h_mean, d_mean, sizeof(float), cudaMemcpyDeviceToHost));
  out = h_mean;
  PROPR_CUDA_CHECK(cudaFree(d_x));
  PROPR_CUDA_CHECK(cudaFree(d_w));
}

void 
dispatch::cuda::wtvRcpp(double& out, const NumericVector& x, const NumericVector& w, propr_context context) {
  const int BLK = 1024;
  PROPR_CHECK_VECTOR_SIZE(x, w.size());

  const int n = x.size();
  float* d_x = RcppVectorToDevice<float>(x, n);
  float* d_w = RcppVectorToDevice<float>(w, n);
  
  float h_var  = 0;
  float *d_var = nullptr;
  PROPR_CUDA_CHECK(cudaMalloc(&d_var, sizeof(float)));
  detail::cuda::wtv<BLK><<<1, BLK, 0, context.stream>>>(d_var, d_x,d_w, n);
  PROPR_CUDA_CHECK(cudaStreamSynchronize(context.stream));
  PROPR_CUDA_CHECK(cudaMemcpy(&h_var, d_var, sizeof(float), cudaMemcpyDeviceToHost));
  out = h_var;
  PROPR_CUDA_CHECK(cudaFree(d_x));
  PROPR_CUDA_CHECK(cudaFree(d_w));
}

void 
centerNumericMatrix(NumericMatrix& out, const NumericMatrix & X, propr_context context){
  using Config = propr::cuda::traits::centerNumericMatrix_config;

  PROPR_CHECK_MATRIX_DIMS(out, X.nrow(), X.ncol());
  
  offset_t d_out_stride; offset_t d_x_stride;
  auto *d_x   = RcppMatrixToDevice<float, REALSXP, true>(X  , d_x_stride  );
  auto *d_out = RcppMatrixToDevice<float, REALSXP, true>(out, d_out_stride);

  int block = Config::BLK_X;
  int grid= propr::ceil_div(X.ncol(), block);

  propr::detail::cuda::centerNumericMatrix<Config::BLK_X><<<grid,block,0,context.stream>>>(d_out, d_out_stride, d_x, d_x_stride, X.nrow(), X.ncol());
  PROPR_CUDA_CHECK(cudaStreamSynchronize(context.stream));

  int ncols = X.ncol();
  int nrows = X.nrow();
  float *centered_mat = new float[nrows * d_out_stride];
  PROPR_CUDA_CHECK(cudaMemcpy(
      centered_mat,
      d_out,
      nrows * d_out_stride * sizeof(float),
      cudaMemcpyDeviceToHost
  ));

  double *outptr = REAL(out);
  for (size_t j = 0; j < ncols; ++j) {
    for (size_t i = 0; i < nrows; ++i) {
      outptr[i + j * nrows] = static_cast<double>(centered_mat[i * d_out_stride + j]);
    }
  }

  delete[] centered_mat;
  centered_mat = nullptr;

  PROPR_CUDA_CHECK(cudaFree(d_x));
  PROPR_CUDA_CHECK(cudaFree(d_out));
}


void
dispatch::cuda::corRcpp(NumericMatrix& out, const NumericMatrix & X, propr_context context) {
  PROPR_CHECK_MATRIX_DIMS(out, X.ncol(), X.ncol());
  using Config = propr::cuda::traits::cor_config;
  int M = X.ncol();   // features
  int K = X.nrow();   // samples

  int M_pad = round_up(M, Config::BLK_M);
  int K_pad = round_up(K, Config::BLK_K);

  offset_t X_stride_src;
  float *d_X_src = RcppMatrixToDevice<float>(X, X_stride_src);

  float *d_X = nullptr;
  PROPR_CUDA_CHECK(cudaMalloc(&d_X, size_t(M_pad) * size_t(K_pad) * sizeof(float)));
  PROPR_CUDA_CHECK(cudaMemsetAsync(d_X, 0, size_t(M_pad) * size_t(K_pad) * sizeof(float), context.stream));

  size_t src_pitch = size_t(X_stride_src) * sizeof(float);
  size_t dst_pitch = size_t(K_pad) * sizeof(float);
  size_t width     = size_t(K) * sizeof(float);
  size_t height    = size_t(M);
  PROPR_CUDA_CHECK(cudaMemcpy2DAsync( d_X, dst_pitch, d_X_src, src_pitch, width, height, cudaMemcpyDeviceToDevice, context.stream));

  float* d_out = nullptr;
  offset_t dout_stride = M_pad;
  PROPR_CUDA_CHECK(cudaMalloc(&d_out, size_t(M_pad) * size_t(M_pad) * sizeof(*d_out)));
  PROPR_CUDA_CHECK(cudaMemsetAsync(d_out, 0, size_t(M_pad) * size_t(M_pad) * sizeof(*d_out), context.stream));

  dim3 block(Config::BLK_M / Config::TH_X, Config::BLK_M / Config::TH_Y);
  dim3 grid(M_pad / Config::BLK_M, M_pad / Config::BLK_M);

  propr::detail::cuda::corRcpp<Config><<<grid, block, 0, context.stream>>>( d_out, dout_stride, d_X, K_pad, /*rows=*/M, /*cols=*/K);
  PROPR_CUDA_CHECK(cudaStreamSynchronize(context.stream));

  auto h_full = new float[size_t(M) * size_t(M)];
  PROPR_CUDA_CHECK(cudaMemcpy2D(
    h_full, 
    size_t(M) * sizeof(float),
    d_out, 
    size_t(dout_stride) * sizeof(float),
    size_t(M) * sizeof(float), 
    size_t(M),
    cudaMemcpyDeviceToHost)
  );

  double *outptr = REAL(out);
  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < M; ++j) {
      outptr[i + j * M] = h_full[i * M + j];
    }
  }
  delete[] h_full;
  h_full = nullptr;
  PROPR_CUDA_CHECK(cudaFree(d_X_src));
  PROPR_CUDA_CHECK(cudaFree(d_X));
  PROPR_CUDA_CHECK(cudaFree(d_out));
}

void 
dispatch::cuda::covRcpp(NumericMatrix& out, const NumericMatrix & X, const int norm_type, propr_context context) {
  PROPR_CHECK_MATRIX_DIMS(out, X.ncol(), X.ncol());

  using Config = propr::cuda::traits::cov_config;
  int nfeats  = X.ncol();   // M
  int samples = X.nrow();   // K
  int M= X.ncol();
  int K= X.nrow();

  const int M_pad = propr::round_up(nfeats,  Config::BLK_M); // Pad rows to BLK_M so every block writes a full 128x128 tile with no if guards
  const int K_pad = propr::round_up(samples, Config::BLK_K); // Pad columns to BLK_K so the last tile can read a full BLK_K 'k-chunk' safely

  offset_t X_stride_pad = K_pad;
  float* d_Xpad = nullptr;
  PROPR_CUDA_CHECK(cudaMalloc(&d_Xpad, (size_t)M_pad * X_stride_pad * sizeof(float)));
  PROPR_CUDA_CHECK(cudaMemset(d_Xpad, 0, (size_t)M_pad * X_stride_pad * sizeof(float)));

  offset_t X_stride_orig;
  const float* d_X = RcppMatrixToDevice<float>(X, X_stride_orig); // returns row-major; X_stride_orig == samples
  PROPR_CUDA_CHECK(cudaMemcpy2D(
      d_Xpad,                         // dst base
      X_stride_pad * sizeof(float),   // dst pitch in bytes
      d_X,                            // src base
      X_stride_orig * sizeof(float),  // src pitch in bytes
      K * sizeof(float),              // width in bytes (valid columns)
      M,                              // number of rows
      cudaMemcpyDeviceToDevice));
  PROPR_CUDA_CHECK(cudaFree((void*)d_X));
  
  float* d_out = nullptr;
  offset_t dout_stride = M_pad;
  PROPR_CUDA_CHECK(cudaMalloc(&d_out, (size_t)M_pad * dout_stride * sizeof(float)));
  PROPR_CUDA_CHECK(cudaMemset(d_out, 0, (size_t)M_pad * dout_stride * sizeof(float)));

  dim3 block(Config::BLK_M / Config::TH_X, Config::BLK_M / Config::TH_Y);
  dim3 grid(M_pad / Config::BLK_M, M_pad / Config::BLK_M);

  propr::detail::cuda::covRcpp<Config><<<grid, block, 0, context.stream>>>(
    norm_type, d_out, dout_stride, d_Xpad, X_stride_pad, /*rows*/ M_pad, /*cols*/ K
  );

  PROPR_CUDA_CHECK(cudaStreamSynchronize(context.stream));

  auto h_full = std::vector<float>((size_t)M_pad * dout_stride);
  PROPR_CUDA_CHECK(cudaMemcpy(h_full.data(), d_out,
                              (size_t)M_pad * dout_stride * sizeof(float),
                              cudaMemcpyDeviceToHost));
  double* outptr = REAL(out);
  for (int i = 0; i < nfeats; ++i)
    for (int j = 0; j < nfeats; ++j)
      outptr[i + j * nfeats] = h_full[i * dout_stride + j];

  PROPR_CUDA_CHECK(cudaFree(d_Xpad));
  PROPR_CUDA_CHECK(cudaFree(d_out));
}

void 
dispatch::cuda::clrRcpp(NumericMatrix& out, const NumericMatrix & X, propr_context context){
    const int rows = X.nrow();
    const int cols = X.ncol();
    PROPR_CHECK_MATRIX_DIMS(out, rows, cols);

    offset_t d_out_stride; 
    auto *d_out = RcppMatrixToDevice<float, REALSXP>(out, d_out_stride);

    offset_t d_x_stride;
    auto *d_x = RcppMatrixToDevice<float, REALSXP>(X, d_x_stride);

    constexpr int BLK_X = 128;
    constexpr int BLK_Y = 4;
    int block = BLK_X * BLK_Y;
    int grid = propr::ceil_div(rows, BLK_Y);
    propr::detail::cuda::clrRcpp<BLK_X, BLK_Y, false><<<grid, block, 0, context.stream>>>(
        d_out, d_out_stride, d_x, d_x_stride, rows, cols
    );
    PROPR_CUDA_CHECK(cudaStreamSynchronize(context.stream));

    float *out_host = new float[cols * d_out_stride];
    PROPR_CUDA_CHECK(cudaMemcpy(
        out_host,
        d_out,
        cols * d_out_stride * sizeof(float),
        cudaMemcpyDeviceToHost
    ));

    double *outptr = REAL(out);
    for (int j = 0; j < cols; ++j) {
        for (int i = 0; i < rows; ++i) {
            outptr[i + j * rows] = static_cast<double>(out_host[i + j * d_out_stride]);
        }
    }

    delete [] out_host;
    PROPR_CUDA_CHECK(cudaFree(d_out));
    PROPR_CUDA_CHECK(cudaFree(d_x));
}


void 
dispatch::cuda::alrRcpp(NumericMatrix& out, const NumericMatrix & X, const int ivar, propr_context context){
    using Config = propr::cuda::traits::alrRcpp_config;
    if (ivar == 0) Rcpp::stop("Select non-zero ivar for alrRcpp.");
    const int nrows = X.nrow();
    const int ncols = X.ncol();
    if (ivar < 1 || ivar > ncols) {
        Rcpp::stop("ivar out of range: must be between 1 and number of columns (%d).", ncols);
    }
    PROPR_CHECK_MATRIX_DIMS(out, nrows, ncols);

    offset_t d_out_stride; auto *d_out = RcppMatrixToDevice<float, REALSXP>(out, d_out_stride);
    offset_t d_x_stride  ; auto *d_x   = RcppMatrixToDevice<float, REALSXP>(X, d_x_stride);

    int block = Config::BLK_X;
    int grid = propr::ceil_div(ncols, block);

    propr::detail::cuda::alrRcpp<Config::BLK_X><<<grid, block, 0, context.stream>>>(
        ivar,
        d_out,
        d_out_stride,
        d_x,
        d_x_stride,
        nrows,
        ncols
    );

    PROPR_CUDA_CHECK(cudaStreamSynchronize(context.stream));

    const size_t total_cols = ncols;
    const size_t total_elems_per_col = d_out_stride;
    const size_t host_elems = total_cols * total_elems_per_col;

    float *out_host = new float[host_elems];
    PROPR_CUDA_CHECK(cudaMemcpy(
        out_host,
        d_out,
        host_elems * sizeof(float),
        cudaMemcpyDeviceToHost
    ));

    double *outptr = REAL(out);
    for (int j = 0; j < total_cols; ++j) {
        for (int i = 0; i < nrows; ++i) {
            outptr[i + j * nrows] =
                static_cast<double>( out_host[i + j * total_elems_per_col] );
        }
    }

    delete[] out_host;
    PROPR_CUDA_CHECK(cudaFree(d_out));
    PROPR_CUDA_CHECK(cudaFree(d_x));
}

void 
dispatch::cuda::symRcpp(NumericMatrix& out, const NumericMatrix & X, propr_context context) {
  using Config = propr::cuda::traits::sym_config;
  PROPR_CHECK_MATRIX_DIMS(out, X.nrow(), X.ncol());
  int nrow = X.nrow();
  int ncol = X.ncol(); 
  
  offset_t X_stride; 
  auto *d_X = RcppMatrixToDevice<float>(X, X_stride, 1);

  offset_t dout_stride; 
  auto *d_out = RcppMatrixToDevice<float>(out, dout_stride, 1);

  dim3 block(Config::TILE, Config::BLK_N);
  dim3 grid(ceil_div(nrow, Config::TILE),ceil_div(nrow, Config::TILE));
  
  propr::detail::cuda::symRcpp<Config><<<grid, block, 0, context.stream>>>(d_out, dout_stride, d_X, X_stride, nrow, ncol);
  PROPR_CUDA_CHECK(cudaStreamSynchronize(context.stream));
  PROPR_CUDA_CHECK(cudaPeekAtLastError());
  auto h_full = new float[nrow * ncol ];
  PROPR_CUDA_CHECK(cudaMemcpy(h_full, d_out, nrow * ncol * sizeof(float), cudaMemcpyDeviceToHost));
  double *outptr = REAL(out);
  for (int j = 0; j < ncol; ++j) {
    for (int i = 0; i < nrow; ++i) {
        outptr[i + j * ncol] = h_full[i + j * ncol];
    }
  }

  delete[] h_full;
  h_full = nullptr;
  
  PROPR_CUDA_CHECK(cudaFree(d_X  ));
  PROPR_CUDA_CHECK(cudaFree(d_out));
}

void 
dispatch::cuda::vlrRcpp(Rcpp::NumericMatrix& out, const Rcpp::NumericMatrix & X, propr_context context){
  using Config = propr::cuda::traits::vlr_config;

  int M = X.ncol();   // features
  int K = X.nrow();   // samples

  int M_pad = ((M + Config::BLK_M - 1) / Config::BLK_M) * Config::BLK_M; // multiple of 128
  int K_pad = ((K + 4 - 1) / 4) * 4;

  int padTop = 0, padLeft = 0;
  int padBottom = M_pad - M;
  int padRight  = K_pad - K;
  Rcpp::NumericMatrix X_pad = rcpp::helpers::pad_matrix(X, padTop, padBottom, padLeft, padRight, 1.0);

  offset_t X_stride; 
  float *d_X = RcppMatrixToDevice<float>(X_pad, X_stride);

  float* d_out = nullptr;
  offset_t dout_stride = M_pad;
  PROPR_CUDA_CHECK(cudaMalloc(&d_out, (size_t)M_pad * dout_stride * sizeof(*d_out)));
  dim3 block(Config::BLK_M / Config::TH_X, Config::BLK_M / Config::TH_Y);
  dim3 grid(ceil_div(M_pad, Config::BLK_M), ceil_div(M_pad, Config::BLK_M));

  propr::detail::cuda::vlrRcpp<Config><<<grid, block, 0, context.stream>>>( d_out, dout_stride, d_X, X_stride, M, K);
  PROPR_CUDA_CHECK(cudaGetLastError());
  PROPR_CUDA_CHECK(cudaStreamSynchronize(context.stream));

  float * h_full =  new float[(size_t)M * M];
  PROPR_CUDA_CHECK(cudaMemcpy2D(
      h_full,                              // dst
      (size_t)M * sizeof(float),           // dst pitch (bytes)
      d_out,                               // src
      (size_t)dout_stride * sizeof(float), // src pitch (bytes)
      (size_t)M * sizeof(float),           // width in bytes
      (size_t)M,                           // height (rows)
      cudaMemcpyDeviceToHost));

  PROPR_CHECK_MATRIX_DIMS(out, M, M);
  double *outptr = REAL(out);
  for (int j = 0; j < M; ++j) {
      for (int i = 0; i < M; ++i) {
          outptr[i + j * M] = h_full[i * M + j];
      }
  }
  delete[] h_full;
  h_full = nullptr;
  PROPR_CUDA_CHECK(cudaFree(d_X));
  PROPR_CUDA_CHECK(cudaFree(d_out));
}

void 
dispatch::cuda::phiRcpp(NumericMatrix& out, NumericMatrix &X, const bool sym, propr_context context) {
    using Config = propr::cuda::traits::phi_config;
    PROPR_CHECK_MATRIX_DIMS(out, X.ncol(), X.ncol());
    int nfeats  = X.ncol(); 
    int samples = X.nrow();

    size_t N = static_cast<size_t>(nfeats);
    //size_t M = static_cast<size_t>(samples);

    offset_t X_stride; 
    auto *d_X = RcppMatrixToDevice<float>(X, X_stride);

    offset_t dout_stride = nfeats;
    float* d_out = nullptr;

    size_t d_out_elems = N * static_cast<size_t>(dout_stride);  // == N*N
    PROPR_CUDA_CHECK(cudaMalloc(&d_out, d_out_elems * sizeof(*d_out)));

    float* row_sums = nullptr;
    PROPR_CUDA_CHECK(cudaMalloc(&row_sums, N * sizeof(*row_sums)));

    float* mu_sum = nullptr;
    PROPR_CUDA_CHECK(cudaMalloc(&mu_sum, sizeof(*mu_sum)));

    int* gbar = nullptr;
    dim3 block(Config::BLK_M / Config::TH_X, Config::BLK_M / Config::TH_Y);
    dim3 grid (ceil_div(nfeats, Config::BLK_M), ceil_div(nfeats, Config::BLK_M));
    
    size_t gbar_len = static_cast<size_t>(grid.x) * grid.y;
    PROPR_CUDA_CHECK(cudaMalloc(&gbar, gbar_len * sizeof(*gbar)));

    PROPR_CUDA_CHECK(cudaMemset(row_sums, 0, N * sizeof(*row_sums)));
    PROPR_CUDA_CHECK(cudaMemset(mu_sum,   0, sizeof(*mu_sum)));
    PROPR_CUDA_CHECK(cudaMemset(gbar,     0, gbar_len * sizeof(*gbar)));

    void* args[] = {static_cast<void*>(const_cast<bool *>(&sym)),
                    static_cast<void*>(&d_out), static_cast<void*>(&dout_stride),
                    static_cast<void*>(&d_X)  , static_cast<void*>(&X_stride),
                    static_cast<void*>(&row_sums), static_cast<void*>(&mu_sum),
                    static_cast<void*>(&nfeats), static_cast<void*>(&samples),
                };

    // std::cout << "<<<(" << grid.x << "," << grid.y << "),(" << block.x << "," << block.y << ")>>>"<< std::endl;
    PROPR_CUDA_CHECK(cudaLaunchCooperativeKernel(propr::detail::cuda::phiRcpp<Config>, grid, block, args, 0, context.stream));
    PROPR_CUDA_CHECK(cudaStreamSynchronize(context.stream));
    auto h_full = new float[d_out_elems];
    PROPR_CUDA_CHECK(cudaMemcpy(h_full, d_out,
                        d_out_elems * sizeof(float),
                        cudaMemcpyDeviceToHost));


   double *outptr = REAL(out);
   for (int i = 0; i < nfeats; ++i) {
        for (int j = 0; j < nfeats; ++j) {
            outptr[i + j * nfeats] = h_full[i * dout_stride  + j];
        }
    }

  delete[] h_full;
  h_full = nullptr;
  
  PROPR_CUDA_CHECK(cudaFree(d_X  ));
  PROPR_CUDA_CHECK(cudaFree(d_out));
  PROPR_CUDA_CHECK(cudaFree(mu_sum));
}

void 
dispatch::cuda::rhoRcpp(NumericMatrix& out,
                        const NumericMatrix &X,
                        const NumericMatrix &lr,
                        const int ivar,
                        propr_context context)
{
  using Config = propr::cuda::traits::rho_config;
  const int nfeats  = lr.ncol();   // M
  const int samples = lr.nrow();   // K

  const int M_pad = propr::round_up(nfeats,  Config::BLK_M);
  const int K_pad = propr::round_up(samples, Config::BLK_K);

  offset_t x_stride_orig;
  const float* d_x_in = RcppMatrixToDevice<float>(X, x_stride_orig);

  float* d_x = nullptr;
  offset_t x_stride = K_pad;
  PROPR_CUDA_CHECK(cudaMalloc(&d_x, (size_t)M_pad * x_stride * sizeof(float)));
  PROPR_CUDA_CHECK(cudaMemset(d_x, 0, (size_t)M_pad * x_stride * sizeof(float)));
  PROPR_CUDA_CHECK(cudaMemcpy2D(
      d_x, x_stride * sizeof(float),
      d_x_in, x_stride_orig * sizeof(float),
      samples * sizeof(float), nfeats,
      cudaMemcpyDeviceToDevice));
  PROPR_CUDA_CHECK(cudaFree((void*)d_x_in));

  offset_t lr_stride_orig;
  const float* d_lr_in = RcppMatrixToDevice<float>(lr, lr_stride_orig);

  float* d_lr = nullptr;
  offset_t lr_stride = K_pad;
  PROPR_CUDA_CHECK(cudaMalloc(&d_lr, (size_t)M_pad * lr_stride * sizeof(float)));
  PROPR_CUDA_CHECK(cudaMemset(d_lr, 0, (size_t)M_pad * lr_stride * sizeof(float)));
  PROPR_CUDA_CHECK(cudaMemcpy2D(
      d_lr, lr_stride * sizeof(float),
      d_lr_in, lr_stride_orig * sizeof(float),
      samples * sizeof(float), nfeats,
      cudaMemcpyDeviceToDevice));
  PROPR_CUDA_CHECK(cudaFree((void*)d_lr_in));

  float* d_out = nullptr;
  offset_t dout_stride = M_pad;
  PROPR_CUDA_CHECK(cudaMalloc(&d_out, (size_t)M_pad * dout_stride * sizeof(float)));
  PROPR_CUDA_CHECK(cudaMemset(d_out, 0, (size_t)M_pad * dout_stride * sizeof(float)));

  dim3 block(Config::BLK_M / Config::TH_X, Config::BLK_M / Config::TH_Y);
  dim3 grid(M_pad / Config::BLK_M, M_pad / Config::BLK_M);

  propr::detail::cuda::rhoRcpp<Config><<<grid, block, 0, context.stream>>>(
      ivar, d_out, dout_stride, d_x, x_stride, d_lr, lr_stride,
      /*rows*/ M_pad, /*cols*/ samples /* true K */
  );
  PROPR_CUDA_CHECK(cudaStreamSynchronize(context.stream));

  auto h_full = std::vector<float>((size_t)M_pad * dout_stride);
  PROPR_CUDA_CHECK(cudaMemcpy(h_full.data(), d_out,
                              (size_t)M_pad * dout_stride * sizeof(float),
                              cudaMemcpyDeviceToHost));

  double* outptr = REAL(out);
  for (int i = 0; i < nfeats; ++i)
    for (int j = 0; j < nfeats; ++j)
      outptr[i + j * nfeats] = h_full[i * dout_stride + j];

  PROPR_CUDA_CHECK(cudaFree(d_out));
  PROPR_CUDA_CHECK(cudaFree(d_lr));
  PROPR_CUDA_CHECK(cudaFree(d_x));
}


void 
dispatch::cuda::indexPairs(std::vector<int>& out,
                           const NumericMatrix & X,
                           const String op,
                           const double ref,
                           propr_context context) {
  out.clear();
  int nfeats = X.nrow();
  for(int i = 1; i < nfeats; i++){
    for(int j = 0; j < i; j++){
      if(op == "==" || op == "="){
        if(X(i, j) == ref){
          out.push_back(j * nfeats + i + 1);
        }
      }else if(op == ">"){
        if(X(i, j) > ref){
          out.push_back(j * nfeats + i + 1);
        }
      }else if(op == ">="){
        if(X(i, j) >= ref){
          out.push_back(j * nfeats + i + 1);
        }
      }else if(op == "<"){
        if(X(i, j) < ref){
          out.push_back(j * nfeats + i + 1);
        }
      }else if(op == "<="){
        if(X(i, j) <= ref){
          out.push_back(j * nfeats + i + 1);
        }
      }else if(op == "!="){
        if(X(i, j) != ref){
          out.push_back(j * nfeats + i + 1);
        }
      }else if(op == "all"){
        out.push_back(j * nfeats + i + 1);
      }else{
        stop("Operator not found.");
      }
    }
  }
}

void 
dispatch::cuda::indexToCoord(List& out, IntegerVector V, int N, propr_context context) {
  using Config =  propr::cuda::traits::indexToCoord_config;

  const size_t len = V.length();

    IntegerVector rows(len);
    IntegerVector cols(len);

    if (len == 0) {
        out["feat1"] = rows;
        out["feat2"] = cols;
        return;
    }

    int* d_V = RcppVectorToDevice<int, INTSXP>(V, len);

    int *d_row = nullptr;
    int *d_col = nullptr;
    const size_t bytes = len * sizeof(int);
    PROPR_CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&d_row), bytes));
    PROPR_CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&d_col), bytes));
    PROPR_CUDA_CHECK(cudaMemset(d_row, 0, bytes));
    PROPR_CUDA_CHECK(cudaMemset(d_col, 0, bytes));

    const int block = Config::BLK_X;
    const int grid = propr::ceil_div(len,block);
    propr::detail::cuda::indexToCoord<<<grid, block, 0, context.stream>>>(
        N, d_V, d_row, d_col, len
    );
    PROPR_CUDA_CHECK(cudaStreamSynchronize(context.stream));

    copyToNumericVector<int, INTSXP>(d_row, rows, len);
    copyToNumericVector<int, INTSXP>(d_col, cols, len);

    out["feat1"] = rows;
    out["feat2"] = cols;

    PROPR_CUDA_CHECK(cudaFree(d_V));
    PROPR_CUDA_CHECK(cudaFree(d_row));
    PROPR_CUDA_CHECK(cudaFree(d_col));
}


void 
dispatch::cuda::coordToIndex(
    IntegerVector& out,
    IntegerVector row,
    IntegerVector col,
    int N,
    propr_context context
) {
    using Config = propr::cuda::traits::coordToIndex_config;

    PROPR_CHECK_VECTOR_SIZE(out, row.length());
    if (static_cast<size_t>(col.length()) != static_cast<size_t>(row.length())) {
        Rcpp::stop("coordToIndex: 'row' and 'col' must have the same length");
    }

    const size_t len = row.length();
    if (len == 0) return; 

    int *d_out = RcppVectorToDevice<int, INTSXP>(out, len);
    int *d_row = RcppVectorToDevice<int, INTSXP>(row, len);
    int *d_col = RcppVectorToDevice<int, INTSXP>(col, len);

    const int block = Config::BLK_X;
    const int grid  = propr::ceil_div(len, block);
    
    propr::detail::cuda::coordToIndex<<<grid, block, 0, context.stream>>>(
        N, d_out, d_row, d_col, len
    );
    PROPR_CUDA_CHECK(cudaStreamSynchronize(context.stream));
    copyToNumericVector(d_out, out, len);
    PROPR_CUDA_CHECK(cudaFree(d_out));
    PROPR_CUDA_CHECK(cudaFree(d_row));
    PROPR_CUDA_CHECK(cudaFree(d_col));
}


void 
dispatch::cuda::linRcpp(NumericMatrix& out, const NumericMatrix & rho, const NumericMatrix &lr, propr_context context){
    // PROPR_CHECK_MATRIX_DIMS(out, rho.ncol(), rho.ncol());

    using Config = propr::cuda::traits::lin_config;
    size_t nfeats  = lr.ncol();
    size_t samples = lr.nrow();

    offset_t lr_stride;
    auto *d_lr = RcppMatrixToDevice<float>(lr, lr_stride);

    offset_t rho_stride;
    auto *d_rho = RcppMatrixToDevice<float>(rho, rho_stride);

    float* d_out = nullptr;
    offset_t dout_stride = nfeats;
    PROPR_CUDA_CHECK(cudaMalloc(&d_out, nfeats * dout_stride * sizeof(*d_out)));
    
    dim3 block(Config::BLK_M / Config::TH_X, Config::BLK_M / Config::TH_Y);
    dim3 grid(ceil_div(nfeats, Config::BLK_M), ceil_div(nfeats, Config::BLK_M));
    
    propr::detail::cuda::linRcpp<Config><<<grid, block, 0, context.stream>>>(d_out, dout_stride, d_rho, rho_stride, d_lr, lr_stride, nfeats, samples);
    PROPR_CUDA_CHECK(cudaStreamSynchronize(context.stream));

    auto h_full = new float[nfeats * nfeats];
    PROPR_CUDA_CHECK(cudaMemcpy(
        h_full,
        d_out,
        nfeats * dout_stride * sizeof(float),
        cudaMemcpyDeviceToHost
    ));

    double *outptr = REAL(out);
    for (size_t i = 0; i < nfeats; ++i) {
            for (size_t j = 0; j < nfeats; ++j) {
                outptr[i + j * nfeats] = h_full[i * dout_stride  + j];
            }
        }

    delete[] h_full;
    h_full = nullptr;
    PROPR_CUDA_CHECK(cudaFree(d_out));
    PROPR_CUDA_CHECK(cudaFree(d_lr));
    PROPR_CUDA_CHECK(cudaFree(d_rho));

}

void 
dispatch::cuda::lltRcpp(NumericVector& out, const NumericMatrix & X, propr_context context){
    using Config = propr::cuda::traits::lltRcpp_config;
    int nfeats = X.nrow();
    int llt = nfeats * (nfeats - 1) / 2;
    // PROPR_CHECK_VECTOR_SIZE(out, llt);

    auto* d_out = RcppVectorToDevice<float>(out, llt);
    offset_t d_x_stride;
    auto *d_x = RcppMatrixToDevice<float, REALSXP>(X, d_x_stride);
    
    int block = Config::BLK_X;
    int grid= propr::ceil_div(llt,block);
    propr::detail::cuda::lltRcpp<<<grid, block,0,context.stream>>>(d_out, llt, d_x, d_x_stride);
    PROPR_CUDA_CHECK(cudaStreamSynchronize(context.stream));
    copyToNumericVector(d_out, out, llt);
    PROPR_CUDA_CHECK(cudaFree(d_x));
    PROPR_CUDA_CHECK(cudaFree(d_out));
}

void dispatch::cuda::urtRcpp(NumericVector& out, const NumericMatrix & X, propr_context context){
    using Config = propr::cuda::traits::urtRcpp_config;
    int nfeats = X.nrow();
    int llt = nfeats * (nfeats - 1) / 2;
    PROPR_CHECK_VECTOR_SIZE(out, llt);

    auto* d_out = RcppVectorToDevice<float>(out, llt);
    offset_t d_x_stride;
    auto *d_x = RcppMatrixToDevice<float, REALSXP>(X, d_x_stride);
    int block = Config::BLK_X;
    int grid= propr::ceil_div(llt, block);
    propr::detail::cuda::lltRcpp<<<grid, block,0,context.stream>>>(d_out, llt, d_x, d_x_stride);
    PROPR_CUDA_CHECK(cudaStreamSynchronize(context.stream));
    copyToNumericVector(d_out, out, llt);
    PROPR_CUDA_CHECK(cudaFree(d_x));
    PROPR_CUDA_CHECK(cudaFree(d_out));
}

void dispatch::cuda::labRcpp(List & out, int nfeats, propr_context context){
  using Config = propr::cuda::traits::labRcpp_config;
  int llt = nfeats * (nfeats - 1) / 2;

  int *d_partner; int *d_pair;
  PROPR_CUDA_CHECK(cudaMalloc(&d_partner, sizeof(*d_partner) * llt));
  PROPR_CUDA_CHECK(cudaMalloc(&d_pair, sizeof(*d_pair) * llt ));

  int block = Config::BLK_X;
  int grid= propr::ceil_div(llt, block);
  propr::detail::cuda::labRcpp<<<grid, block,0,context.stream>>>(d_partner, d_pair, nfeats);
  PROPR_CUDA_CHECK(cudaStreamSynchronize(context.stream));

  out["Partner"] = Rcpp::IntegerVector(llt);
  out["Pair"] = Rcpp::IntegerVector(llt);

  Rcpp::IntegerVector partner_ref = out["Partner"];
  Rcpp::IntegerVector pair_ref = out["Pair"];

  copyToNumericVector(d_partner,partner_ref, llt);
  copyToNumericVector(d_pair, pair_ref, llt);

  PROPR_CUDA_CHECK(cudaFree(d_partner));
  PROPR_CUDA_CHECK(cudaFree(d_pair));
}

void 
dispatch::cuda::half2mat(NumericMatrix& out, const NumericVector & X, propr_context context){
    using Config = propr::cuda::traits::half2mat_config;
    size_t nfeats = static_cast<int>(std::round(std::sqrt(2.0 * static_cast<double>(X.size()) + 0.25) + 0.5));
    PROPR_CHECK_MATRIX_DIMS(out, nfeats, nfeats);
    const size_t total_pairs = nfeats * static_cast<size_t>(nfeats - 1) / 2;

    if (static_cast<size_t>(X.size()) != total_pairs) {
        Rcpp::stop("half2mat: length(X) != nfeats*(nfeats-1)/2 (recomputed nfeats=%d, expected pairs=%zu, got=%zu)",
                   nfeats, total_pairs, static_cast<size_t>(X.size()));
    }

    offset_t d_out_stride;
    float *d_out = RcppMatrixToDevice<float, REALSXP>(out, d_out_stride);
    float *d_X   = RcppVectorToDevice<float, REALSXP>(X, total_pairs);

    const size_t block = Config::BLK_X;
    const int grid     = static_cast<int>(propr::ceil_div(total_pairs, block));
    propr::detail::cuda::half2mat<<<grid, block, 0, context.stream>>>(d_out, d_out_stride, d_X, nfeats);
    PROPR_CUDA_CHECK(cudaStreamSynchronize(context.stream));

    const size_t total_elems = d_out_stride * nfeats;
    const size_t total_bytes = total_elems * sizeof(float);
    float *out_host = new float[total_elems];
    PROPR_CUDA_CHECK(cudaMemcpy(out_host, d_out, total_bytes, cudaMemcpyDeviceToHost));

    double *outptr = REAL(out);
    for (size_t col = 0; col < nfeats; ++col) {
        for (size_t row = 0; row < nfeats; ++row) {
            const size_t host_idx = row + col * static_cast<size_t>(d_out_stride);
            const size_t r_idx    = row + col * nfeats;
            outptr[r_idx] = static_cast<double>(out_host[host_idx]);
        }
    }

    delete[] out_host;

    PROPR_CUDA_CHECK(cudaFree(d_out));
    PROPR_CUDA_CHECK(cudaFree(d_X));
}


void 
dispatch::cuda::vector2mat(
    NumericMatrix& out,
    const NumericVector & X,
    const IntegerVector & i,
    const IntegerVector & j,
    int nfeats,
    propr_context context
){
    using Config  = propr::cuda::traits::vector2mat_config;

    int nX = X.length();
    int ni = i.length();
    int nj = j.length();
    if (ni != nj) Rcpp::stop("i and j must be the same length.");
    if (ni != nX) Rcpp::stop("i, j, and X must be the same length.");
    PROPR_CHECK_MATRIX_DIMS(out, nfeats, nfeats);
    offset_t d_out_stride;
    auto *d_out = RcppMatrixToDevice<float, REALSXP>(out, d_out_stride);
    auto *d_X   = RcppVectorToDevice<float, REALSXP>(X, ni);
    auto *d_i   = RcppVectorToDevice<int, INTSXP>(i, ni);
    auto *d_j   = RcppVectorToDevice<int, INTSXP>(j, ni);

    const int block = Config::BLK_X;
    const int grid = propr::ceil_div(ni,block);

    propr::detail::cuda::vector2mat<<<grid, block, 0, context.stream>>>(
        d_out,
        d_out_stride,
        d_X,
        d_i,
        d_j,
        ni
    );
    PROPR_CUDA_CHECK(cudaStreamSynchronize(context.stream));

    const size_t total_elems = static_cast<size_t>(d_out_stride) * static_cast<size_t>(nfeats);
    const size_t total_bytes = total_elems * sizeof(float);
    float *out_host = new float[total_elems];
    PROPR_CUDA_CHECK(cudaMemcpy(
        out_host,
        d_out,
        total_bytes,
        cudaMemcpyDeviceToHost
    ));

    double *outptr = REAL(out);
    for (size_t col = 0; col < static_cast<size_t>(nfeats); ++col) {
        for (size_t row = 0; row < static_cast<size_t>(nfeats); ++row) {
            const size_t host_idx = row + col * static_cast<size_t>(d_out_stride);
            const size_t r_idx    = row + col * static_cast<size_t>(nfeats);
            outptr[r_idx] = static_cast<double>(out_host[host_idx]);
        }
    }

    delete[] out_host;
    PROPR_CUDA_CHECK(cudaFree(d_out));
    PROPR_CUDA_CHECK(cudaFree(d_X));
    PROPR_CUDA_CHECK(cudaFree(d_i));
    PROPR_CUDA_CHECK(cudaFree(d_j));
}


void 
dispatch::cuda::ratiosRcpp(NumericMatrix & out, const NumericMatrix & X, propr_context context){
    using Config = propr::cuda::traits::ratiosRcpp_config;

    int nfeats = X.ncol();
    int nsamps = X.nrow();
    int llt = nfeats * (nfeats - 1) / 2;
    PROPR_CHECK_MATRIX_DIMS(out, nsamps, llt);

    offset_t d_out_stride; auto *d_out = RcppMatrixToDevice<float>(out, d_out_stride);
    offset_t d_x_stride  ; auto *d_x   = RcppMatrixToDevice<float>(X, d_x_stride);
    
    int block = Config::BLK_X;
    int grid= propr::ceil_div(llt * nsamps, block);
    propr::detail::cuda::ratiosRcpp<<<grid,block,0,context.stream>>>(d_out, d_out_stride, d_x, d_x_stride, nfeats, nsamps);
    PROPR_CUDA_CHECK(cudaStreamSynchronize(context.stream));

    float *out_host = new float[llt * d_out_stride];
    PROPR_CUDA_CHECK(cudaMemcpy(
        out_host,
        d_out,
        llt * d_out_stride * sizeof(float),
        cudaMemcpyDeviceToHost
    ));

    double *outptr = REAL(out);
    for (size_t j = 0; j < llt; ++j) {
      for (size_t i = 0; i < nsamps; ++i) {
        outptr[i + j * nsamps] = static_cast<double>(out_host[i + j * d_out_stride]);
      }
    }
    delete[] out_host;
    PROPR_CUDA_CHECK(cudaFree(d_out));
    PROPR_CUDA_CHECK(cudaFree(d_x));
}

void dispatch::cuda::results2matRcpp(NumericMatrix & out, const DataFrame& results, int n, double diagonal, propr_context context){
    PROPR_CHECK_MATRIX_DIMS(out, n, n);
    Rcpp::stop("results2matRcpp is not implemented in CUDA. Falling back to CPU via dispatcher.");
}