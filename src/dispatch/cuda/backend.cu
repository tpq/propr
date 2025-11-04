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
  CHECK_VECTOR_SIZE(x, w.size());

  const int n = x.size();
  float* d_x = RcppVectorToDevice<float>(x, n);
  float* d_w = RcppVectorToDevice<float>(w, n);
  
  float h_mean  = 0;
  float *d_mean = nullptr;
  CUDA_CHECK(cudaMalloc(&d_mean, sizeof(float)));
  detail::cuda::wtm<BLK><<<1, BLK, 0, context.stream>>>(d_mean, d_x,d_w, n);
  CUDA_CHECK(cudaStreamSynchronize(context.stream));
  CUDA_CHECK(cudaMemcpy(&h_mean, d_mean, sizeof(float), cudaMemcpyDeviceToHost));
  out = h_mean;
  CUDA_CHECK(cudaFree(d_x));
  CUDA_CHECK(cudaFree(d_w));
}

void 
dispatch::cuda::wtvRcpp(double& out, const NumericVector& x, const NumericVector& w, propr_context context) {
  const int BLK = 1024;
  CHECK_VECTOR_SIZE(x, w.size());

  const int n = x.size();
  float* d_x = RcppVectorToDevice<float>(x, n);
  float* d_w = RcppVectorToDevice<float>(w, n);
  
  float h_var  = 0;
  float *d_var = nullptr;
  CUDA_CHECK(cudaMalloc(&d_var, sizeof(float)));
  detail::cuda::wtv<BLK><<<1, BLK, 0, context.stream>>>(d_var, d_x,d_w, n);
  CUDA_CHECK(cudaStreamSynchronize(context.stream));
  CUDA_CHECK(cudaMemcpy(&h_var, d_var, sizeof(float), cudaMemcpyDeviceToHost));
  out = h_var;
  CUDA_CHECK(cudaFree(d_x));
  CUDA_CHECK(cudaFree(d_w));
}

void 
centerNumericMatrix(NumericMatrix& out, const NumericMatrix & X, propr_context context){
  CHECK_MATRIX_DIMS(out, X.nrow(), X.ncol());
  
  offset_t d_out_stride; offset_t d_x_stride;
  auto *d_x   = RcppMatrixToDevice<float, REALSXP, true>(X  , d_x_stride  );
  auto *d_out = RcppMatrixToDevice<float, REALSXP, true>(out, d_out_stride);

  int block = 256;
  int grid = ceil_div(X.ncol(), block);

  propr::detail::cuda::centerNumericMatrix<256><<<grid,block,0,context.stream>>>(d_out, d_out_stride, d_x, d_x_stride, X.nrow(), X.ncol());
  CUDA_CHECK(cudaStreamSynchronize(context.stream));

  int ncols = X.ncol();
  int nrows = X.nrow();
  float *centered_mat = new float[nrows * d_out_stride];
  CUDA_CHECK(cudaMemcpy(
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

  CUDA_CHECK(cudaFree(d_x));
  CUDA_CHECK(cudaFree(d_out));
}

void 
dispatch::cuda::corRcpp(NumericMatrix& out, const NumericMatrix & X, propr_context context) {
  CHECK_MATRIX_DIMS(out, X.ncol(), X.ncol());
  using Config = propr::cuda::traits::cor_config;
  int nfeats  = X.ncol();
  int samples = X.nrow();

  offset_t X_stride;
  auto *d_X = RcppMatrixToDevice<float>(X, X_stride);

  float* d_out = nullptr;
  offset_t dout_stride = nfeats;
  CUDA_CHECK(cudaMalloc(&d_out, nfeats * dout_stride * sizeof(*d_out)));
  
  dim3 block(Config::BLK_M / Config::TH_X, Config::BLK_M / Config::TH_Y);
  dim3 grid(ceil_div(nfeats, Config::BLK_M), ceil_div(nfeats, Config::BLK_M));
  propr::detail::cuda::corRcpp<Config><<<grid, block, 0, context.stream>>>(d_out, dout_stride, d_X, X_stride, nfeats, samples);
  CUDA_CHECK(cudaStreamSynchronize(context.stream));

   auto h_full = new float[nfeats * dout_stride];
  CUDA_CHECK(cudaMemcpy(
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
  CUDA_CHECK(cudaFree(d_X));
  CUDA_CHECK(cudaFree(d_out));
}

void 
dispatch::cuda::covRcpp(NumericMatrix& out, const NumericMatrix & X, const int norm_type, propr_context context) {
  CHECK_MATRIX_DIMS(out, X.ncol(), X.ncol());
  using Config = propr::cuda::traits::cov_config;
  int nfeats  = X.ncol();
  int samples = X.nrow();

  offset_t X_stride;
  auto *d_X = RcppMatrixToDevice<float>(X, X_stride);

  float* d_out = nullptr;
  offset_t dout_stride = nfeats;
  CUDA_CHECK(cudaMalloc(&d_out, nfeats * dout_stride * sizeof(*d_out)));
  
  dim3 block(Config::BLK_M / Config::TH_X, Config::BLK_M / Config::TH_Y);
  dim3 grid(ceil_div(nfeats, Config::BLK_M), ceil_div(nfeats, Config::BLK_M));
  
  propr::detail::cuda::covRcpp<Config><<<grid, block, 0, context.stream>>>(norm_type, d_out, dout_stride, d_X, X_stride, nfeats, samples);
  CUDA_CHECK(cudaStreamSynchronize(context.stream));

   auto h_full = new float[nfeats * dout_stride];
  CUDA_CHECK(cudaMemcpy(
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
  CUDA_CHECK(cudaFree(d_X));
  CUDA_CHECK(cudaFree(d_out));
}

void 
dispatch::cuda::clrRcpp(NumericMatrix& out, const NumericMatrix & X, propr_context context){
    const int rows = X.nrow();
    const int cols = X.ncol();
    CHECK_MATRIX_DIMS(out, rows, cols);

    offset_t d_out_stride; 
    auto *d_out = RcppMatrixToDevice<float, REALSXP>(out, d_out_stride);

    offset_t d_x_stride;
    auto *d_x = RcppMatrixToDevice<float, REALSXP>(X, d_x_stride);

    constexpr int BLK_X = 128;
    constexpr int BLK_Y = 4;
    int block = BLK_X * BLK_Y;
    int grid  = ceil_div(rows, BLK_Y);
    propr::detail::cuda::clrRcpp<BLK_X, BLK_Y, false><<<grid, block, 0, context.stream>>>(
        d_out, d_out_stride, d_x, d_x_stride, rows, cols
    );
    CUDA_CHECK(cudaStreamSynchronize(context.stream));

    float *out_host = new float[cols * d_out_stride];
    CUDA_CHECK(cudaMemcpy(
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
    CUDA_CHECK(cudaFree(d_out));
    CUDA_CHECK(cudaFree(d_x));
}


void 
dispatch::cuda::alrRcpp(NumericMatrix& out, const NumericMatrix & X, const int ivar, propr_context context){
    if (ivar == 0) Rcpp::stop("Select non-zero ivar for alrRcpp.");
    const int nrows = X.nrow();
    const int ncols = X.ncol();
    if (ivar < 1 || ivar > ncols) {
        Rcpp::stop("ivar out of range: must be between 1 and number of columns (%d).", ncols);
    }
    CHECK_MATRIX_DIMS(out, nrows, ncols);

    offset_t d_out_stride; auto *d_out = RcppMatrixToDevice<float, REALSXP>(out, d_out_stride);
    offset_t d_x_stride  ; auto *d_x   = RcppMatrixToDevice<float, REALSXP>(X, d_x_stride);

    int block = 256;
    int grid  = ceil_div(ncols, block);

    propr::detail::cuda::alrRcpp<256><<<grid, block, 0, context.stream>>>(
        ivar,
        d_out,
        d_out_stride,
        d_x,
        d_x_stride,
        nrows,
        ncols
    );

    CUDA_CHECK(cudaStreamSynchronize(context.stream));

    const size_t total_cols = ncols;
    const size_t total_elems_per_col = d_out_stride;
    const size_t host_elems = total_cols * total_elems_per_col;

    float *out_host = new float[host_elems];
    CUDA_CHECK(cudaMemcpy(
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
    CUDA_CHECK(cudaFree(d_out));
    CUDA_CHECK(cudaFree(d_x));
}

void 
dispatch::cuda::symRcpp(NumericMatrix& out, const NumericMatrix & X, propr_context context) {
  using Config = propr::cuda::traits::sym_config;
  CHECK_MATRIX_DIMS(out, X.nrow(), X.ncol());
  int nrow = X.nrow();
  int ncol = X.ncol(); 
  
  offset_t X_stride; 
  auto *d_X = RcppMatrixToDevice<float>(X, X_stride, 1);

  offset_t dout_stride; 
  auto *d_out = RcppMatrixToDevice<float>(out, dout_stride, 1);

  dim3 block(Config::TILE, Config::BLK_N);
  dim3 grid(ceil_div(nrow, Config::TILE),ceil_div(nrow, Config::TILE));
  
  propr::detail::cuda::symRcpp<Config><<<grid, block, 0, context.stream>>>(d_out, dout_stride, d_X, X_stride, nrow, ncol);
  CUDA_CHECK(cudaStreamSynchronize(context.stream));
  CUDA_CHECK(cudaPeekAtLastError());
  auto h_full = new float[nrow * ncol ];
  CUDA_CHECK(cudaMemcpy(h_full, d_out, nrow * ncol * sizeof(float), cudaMemcpyDeviceToHost));
  double *outptr = REAL(out);
  for (int j = 0; j < ncol; ++j) {
    for (int i = 0; i < nrow; ++i) {
        outptr[i + j * ncol] = h_full[i + j * ncol];
    }
  }

  delete[] h_full;
  h_full = nullptr;
  
  CUDA_CHECK(cudaFree(d_X  ));
  CUDA_CHECK(cudaFree(d_out));
}

void 
dispatch::cuda::vlrRcpp(NumericMatrix& out, const NumericMatrix & X, propr_context context){
  CHECK_MATRIX_DIMS(out, X.ncol(), X.ncol());
  using Config = propr::cuda::traits::vlr_config;
  int nfeats  = X.ncol();
  int samples = X.nrow();

  offset_t X_stride;
  auto *d_X = RcppMatrixToDevice<float>(X, X_stride);

  float* d_out = nullptr;
  offset_t dout_stride = nfeats;
  CUDA_CHECK(cudaMalloc(&d_out, nfeats * dout_stride * sizeof(*d_out)));
  
  dim3 block(Config::BLK_M / Config::TH_X, Config::BLK_M / Config::TH_Y);
  dim3 grid(ceil_div(nfeats, Config::BLK_M), ceil_div(nfeats, Config::BLK_M));

  propr::detail::cuda::vlrRcpp<Config><<<grid, block, 0, context.stream>>>(d_out, dout_stride, d_X, X_stride, nfeats, samples);
  CUDA_CHECK(cudaStreamSynchronize(context.stream));

  auto h_full = new float[nfeats * nfeats];
  CUDA_CHECK(cudaMemcpy(
      h_full,
      d_out,
      nfeats * dout_stride * sizeof(float),
      cudaMemcpyDeviceToHost
  ));

   double *outptr = REAL(out);
   for (int i = 0; i < nfeats; ++i) {
        for (int j = 0; j < nfeats; ++j) {
            outptr[i + j * nfeats] = h_full[i * dout_stride  + j];
        }
    }

  delete[] h_full;
  h_full = nullptr;
  CUDA_CHECK(cudaFree(d_X));
  CUDA_CHECK(cudaFree(d_out));
}

void 
dispatch::cuda::phiRcpp(NumericMatrix& out, NumericMatrix &X, const bool sym, propr_context context) {
    using Config = propr::cuda::traits::phi_config;
    CHECK_MATRIX_DIMS(out, X.ncol(), X.ncol());
    int nfeats  = X.ncol(); 
    int samples = X.nrow();

    size_t N = static_cast<size_t>(nfeats);
    //size_t M = static_cast<size_t>(samples);

    offset_t X_stride; 
    auto *d_X = RcppMatrixToDevice<float>(X, X_stride);

    offset_t dout_stride = nfeats;
    float* d_out = nullptr;

    size_t d_out_elems = N * static_cast<size_t>(dout_stride);  // == N*N
    CUDA_CHECK(cudaMalloc(&d_out, d_out_elems * sizeof(*d_out)));

    float* row_sums = nullptr;
    CUDA_CHECK(cudaMalloc(&row_sums, N * sizeof(*row_sums)));

    float* mu_sum = nullptr;
    CUDA_CHECK(cudaMalloc(&mu_sum, sizeof(*mu_sum)));

    int* gbar = nullptr;
    dim3 block(Config::BLK_M / Config::TH_X, Config::BLK_M / Config::TH_Y);
    dim3 grid (ceil_div(nfeats, Config::BLK_M), ceil_div(nfeats, Config::BLK_M));
    
    size_t gbar_len = static_cast<size_t>(grid.x) * grid.y;
    CUDA_CHECK(cudaMalloc(&gbar, gbar_len * sizeof(*gbar)));

    CUDA_CHECK(cudaMemset(row_sums, 0, N * sizeof(*row_sums)));
    CUDA_CHECK(cudaMemset(mu_sum,   0, sizeof(*mu_sum)));
    CUDA_CHECK(cudaMemset(gbar,     0, gbar_len * sizeof(*gbar)));

    void* args[] = {static_cast<void*>(const_cast<bool *>(&sym)),
                    static_cast<void*>(&d_out), static_cast<void*>(&dout_stride),
                    static_cast<void*>(&d_X)  , static_cast<void*>(&X_stride),
                    static_cast<void*>(&row_sums), static_cast<void*>(&mu_sum),
                    static_cast<void*>(&nfeats), static_cast<void*>(&samples),
                };

    // std::cout << "<<<(" << grid.x << "," << grid.y << "),(" << block.x << "," << block.y << ")>>>"<< std::endl;
    CUDA_CHECK(cudaLaunchCooperativeKernel(propr::detail::cuda::phiRcpp<Config>, grid, block, args, 0, context.stream));
    CUDA_CHECK(cudaStreamSynchronize(context.stream));
    auto h_full = new float[d_out_elems];
    CUDA_CHECK(cudaMemcpy(h_full, d_out,
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
  
  CUDA_CHECK(cudaFree(d_X  ));
  CUDA_CHECK(cudaFree(d_out));
  CUDA_CHECK(cudaFree(mu_sum));
}

void 
dispatch::cuda::rhoRcpp(NumericMatrix& out, const NumericMatrix &X, const NumericMatrix &lr, const int ivar, propr_context context){
    using Config = propr::cuda::traits::rho_config;
    size_t nfeats  = lr.ncol();
    size_t samples = lr.nrow();

    offset_t lr_stride;
    auto *d_lr = RcppMatrixToDevice<float>(lr, lr_stride);

    offset_t x_stride;
    auto *d_x = RcppMatrixToDevice<float>(X, x_stride);

    float* d_out = nullptr;
    offset_t dout_stride = nfeats;
    CUDA_CHECK(cudaMalloc(&d_out, nfeats * dout_stride * sizeof(*d_out)));
    
    dim3 block(Config::BLK_M / Config::TH_X, Config::BLK_M / Config::TH_Y);
    dim3 grid(ceil_div(nfeats, Config::BLK_M), ceil_div(nfeats, Config::BLK_M));
    
    propr::detail::cuda::rhoRcpp<Config><<<grid, block, 0, context.stream>>>(ivar, d_out, dout_stride, d_x, x_stride, d_lr, lr_stride, nfeats, samples);
    CUDA_CHECK(cudaStreamSynchronize(context.stream));

    auto h_full = new float[nfeats * nfeats];
    CUDA_CHECK(cudaMemcpy(
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
    CUDA_CHECK(cudaFree(d_out));
    CUDA_CHECK(cudaFree(d_lr));
    CUDA_CHECK(cudaFree(d_x));
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
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&d_row), bytes));
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&d_col), bytes));
    CUDA_CHECK(cudaMemset(d_row, 0, bytes));
    CUDA_CHECK(cudaMemset(d_col, 0, bytes));

    const int block = 256;
    const int grid  = ceil_div(len,block);
    propr::detail::cuda::indexToCoord<<<grid, block, 0, context.stream>>>(
        N, d_V, d_row, d_col, len
    );
    CUDA_CHECK(cudaStreamSynchronize(context.stream));

    copyToNumericVector<int, INTSXP>(d_row, rows, len);
    copyToNumericVector<int, INTSXP>(d_col, cols, len);

    out["feat1"] = rows;
    out["feat2"] = cols;

    CUDA_CHECK(cudaFree(d_V));
    CUDA_CHECK(cudaFree(d_row));
    CUDA_CHECK(cudaFree(d_col));
}


void dispatch::cuda::coordToIndex(
    IntegerVector& out,
    IntegerVector row,
    IntegerVector col,
    int N,
    propr_context context
) {
    CHECK_VECTOR_SIZE(out, row.length());

    if (static_cast<size_t>(col.length()) != static_cast<size_t>(row.length())) {
        Rcpp::stop("coordToIndex: 'row' and 'col' must have the same length");
    }

    const size_t len = row.length();
    if (len == 0) return; 

    int *d_out = RcppVectorToDevice<int, INTSXP>(out, len);
    int *d_row = RcppVectorToDevice<int, INTSXP>(row, len);
    int *d_col = RcppVectorToDevice<int, INTSXP>(col, len);

    const int block = 256;
    const int grid  = ceil_div(len,block);
    
    propr::detail::cuda::coordToIndex<<<grid, block, 0, context.stream>>>(
        N, d_out, d_row, d_col, len
    );
    CUDA_CHECK(cudaStreamSynchronize(context.stream));
    copyToNumericVector(d_out, out, len);
    CUDA_CHECK(cudaFree(d_out));
    CUDA_CHECK(cudaFree(d_row));
    CUDA_CHECK(cudaFree(d_col));
}


void 
dispatch::cuda::linRcpp(NumericMatrix& out, const NumericMatrix & rho, const NumericMatrix &lr, propr_context context){
    // CHECK_MATRIX_DIMS(out, rho.ncol(), rho.ncol());

    using Config = propr::cuda::traits::lin_config;
    size_t nfeats  = lr.ncol();
    size_t samples = lr.nrow();

    offset_t lr_stride;
    auto *d_lr = RcppMatrixToDevice<float>(lr, lr_stride);

    offset_t rho_stride;
    auto *d_rho = RcppMatrixToDevice<float>(rho, rho_stride);

    float* d_out = nullptr;
    offset_t dout_stride = nfeats;
    CUDA_CHECK(cudaMalloc(&d_out, nfeats * dout_stride * sizeof(*d_out)));
    
    dim3 block(Config::BLK_M / Config::TH_X, Config::BLK_M / Config::TH_Y);
    dim3 grid(ceil_div(nfeats, Config::BLK_M), ceil_div(nfeats, Config::BLK_M));
    
    propr::detail::cuda::linRcpp<Config><<<grid, block, 0, context.stream>>>(d_out, dout_stride, d_rho, rho_stride, d_lr, lr_stride, nfeats, samples);
    CUDA_CHECK(cudaStreamSynchronize(context.stream));

    auto h_full = new float[nfeats * nfeats];
    CUDA_CHECK(cudaMemcpy(
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
    CUDA_CHECK(cudaFree(d_out));
    CUDA_CHECK(cudaFree(d_lr));
    CUDA_CHECK(cudaFree(d_rho));

}

void 
dispatch::cuda::lltRcpp(NumericVector& out, const NumericMatrix & X, propr_context context){
    int nfeats = X.nrow();
    int llt = nfeats * (nfeats - 1) / 2;
    // CHECK_VECTOR_SIZE(out, llt);

    auto* d_out = RcppVectorToDevice<float>(out, llt);
    offset_t d_x_stride;
    auto *d_x = RcppMatrixToDevice<float, REALSXP>(X, d_x_stride);
    
    int block = 256;
    int grid = ceil_div(llt,block);
    propr::detail::cuda::lltRcpp<<<grid, block,0,context.stream>>>(d_out, llt, d_x, d_x_stride);
    CUDA_CHECK(cudaStreamSynchronize(context.stream));
    copyToNumericVector(d_out, out, llt);
    CUDA_CHECK(cudaFree(d_x));
    CUDA_CHECK(cudaFree(d_out));
}

void dispatch::cuda::urtRcpp(NumericVector& out, const NumericMatrix & X, propr_context context){
    int nfeats = X.nrow();
    int llt = nfeats * (nfeats - 1) / 2;
    CHECK_VECTOR_SIZE(out, llt);

    auto* d_out = RcppVectorToDevice<float>(out, llt);
    offset_t d_x_stride;
    auto *d_x = RcppMatrixToDevice<float, REALSXP>(X, d_x_stride);
    int block = 256;
    int grid = ceil_div(llt, block);
    propr::detail::cuda::lltRcpp<<<grid, block,0,context.stream>>>(d_out, llt, d_x, d_x_stride);
    CUDA_CHECK(cudaStreamSynchronize(context.stream));
    copyToNumericVector(d_out, out, llt);
    CUDA_CHECK(cudaFree(d_x));
    CUDA_CHECK(cudaFree(d_out));
}

void dispatch::cuda::labRcpp(List & out, int nfeats, propr_context context){
  int llt = nfeats * (nfeats - 1) / 2;

  int *d_partner; int *d_pair;
  CUDA_CHECK(cudaMalloc(&d_partner, sizeof(*d_partner) * llt));
  CUDA_CHECK(cudaMalloc(&d_pair, sizeof(*d_pair) * llt ));

  int block = 256;
  int grid = ceil_div(llt, block);
  propr::detail::cuda::labRcpp<<<grid, block,0,context.stream>>>(d_partner, d_pair, nfeats);
  CUDA_CHECK(cudaStreamSynchronize(context.stream));

  out["Partner"] = Rcpp::IntegerVector(llt);
  out["Pair"] = Rcpp::IntegerVector(llt);

  Rcpp::IntegerVector partner_ref = out["Partner"];
  Rcpp::IntegerVector pair_ref = out["Pair"];

  copyToNumericVector(d_partner,partner_ref, llt);
  copyToNumericVector(d_pair, pair_ref, llt);

  CUDA_CHECK(cudaFree(d_partner));
  CUDA_CHECK(cudaFree(d_pair));
}

void 
dispatch::cuda::half2mat(NumericMatrix& out, const NumericVector & X, propr_context context){
    size_t nfeats = static_cast<int>(std::round(std::sqrt(2.0 * static_cast<double>(X.size()) + 0.25) + 0.5));
    CHECK_MATRIX_DIMS(out, nfeats, nfeats);
    const size_t total_pairs = nfeats * static_cast<size_t>(nfeats - 1) / 2;

    if (static_cast<size_t>(X.size()) != total_pairs) {
        Rcpp::stop("half2mat: length(X) != nfeats*(nfeats-1)/2 (recomputed nfeats=%d, expected pairs=%zu, got=%zu)",
                   nfeats, total_pairs, static_cast<size_t>(X.size()));
    }

    offset_t d_out_stride;
    float *d_out = RcppMatrixToDevice<float, REALSXP>(out, d_out_stride);
    float *d_X   = RcppVectorToDevice<float, REALSXP>(X, total_pairs);

    const size_t block = 256;
    const int grid     = static_cast<int>(ceil_div(total_pairs, block));
    propr::detail::cuda::half2mat<<<grid, block, 0, context.stream>>>(d_out, d_out_stride, d_X, nfeats);
    CUDA_CHECK(cudaStreamSynchronize(context.stream));

    const size_t total_elems = d_out_stride * nfeats;
    const size_t total_bytes = total_elems * sizeof(float);
    float *out_host = new float[total_elems];
    CUDA_CHECK(cudaMemcpy(out_host, d_out, total_bytes, cudaMemcpyDeviceToHost));

    double *outptr = REAL(out);
    for (size_t col = 0; col < nfeats; ++col) {
        for (size_t row = 0; row < nfeats; ++row) {
            const size_t host_idx = row + col * static_cast<size_t>(d_out_stride);
            const size_t r_idx    = row + col * nfeats;
            outptr[r_idx] = static_cast<double>(out_host[host_idx]);
        }
    }

    delete[] out_host;

    CUDA_CHECK(cudaFree(d_out));
    CUDA_CHECK(cudaFree(d_X));
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
    int nX = X.length();
    int ni = i.length();
    int nj = j.length();
    if (ni != nj) Rcpp::stop("i and j must be the same length.");
    if (ni != nX) Rcpp::stop("i, j, and X must be the same length.");
    CHECK_MATRIX_DIMS(out, nfeats, nfeats);
    offset_t d_out_stride;
    auto *d_out = RcppMatrixToDevice<float, REALSXP>(out, d_out_stride);
    auto *d_X   = RcppVectorToDevice<float, REALSXP>(X, ni);
    auto *d_i   = RcppVectorToDevice<int, INTSXP>(i, ni);
    auto *d_j   = RcppVectorToDevice<int, INTSXP>(j, ni);

    const int block = 256;
    const int grid  = ceil_div(ni,block);

    propr::detail::cuda::vector2mat<<<grid, block, 0, context.stream>>>(
        d_out,
        d_out_stride,
        d_X,
        d_i,
        d_j,
        ni
    );
    CUDA_CHECK(cudaStreamSynchronize(context.stream));

    const size_t total_elems = static_cast<size_t>(d_out_stride) * static_cast<size_t>(nfeats);
    const size_t total_bytes = total_elems * sizeof(float);
    float *out_host = new float[total_elems];
    CUDA_CHECK(cudaMemcpy(
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
    CUDA_CHECK(cudaFree(d_out));
    CUDA_CHECK(cudaFree(d_X));
    CUDA_CHECK(cudaFree(d_i));
    CUDA_CHECK(cudaFree(d_j));
}


void dispatch::cuda::ratiosRcpp(NumericMatrix & out, const NumericMatrix & X, propr_context context){
    int nfeats = X.ncol();
    int nsamps = X.nrow();
    int llt = nfeats * (nfeats - 1) / 2;
    CHECK_MATRIX_DIMS(out, nsamps, llt);

    offset_t d_out_stride; auto *d_out = RcppMatrixToDevice<float>(out, d_out_stride);
    offset_t d_x_stride  ; auto *d_x   = RcppMatrixToDevice<float>(X, d_x_stride);
    
    int block = 256;
    int grid = ceil_div(llt * nsamps, block);
    propr::detail::cuda::ratiosRcpp<<<grid,block,0,context.stream>>>(d_out, d_out_stride, d_x, d_x_stride, nfeats, nsamps);
    CUDA_CHECK(cudaStreamSynchronize(context.stream));

    float *out_host = new float[llt * d_out_stride];
    CUDA_CHECK(cudaMemcpy(
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
    CUDA_CHECK(cudaFree(d_out));
    CUDA_CHECK(cudaFree(d_x));
}

void dispatch::cuda::results2matRcpp(NumericMatrix & out, const DataFrame& results, int n, double diagonal, propr_context context){
    CHECK_MATRIX_DIMS(out, n, n);
    Rcpp::stop("results2matRcpp is not implemented in CUDA. Falling back to CPU via dispatcher.");
}