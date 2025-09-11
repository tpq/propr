#include <Rcpp.h>
#include <math.h>

#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/system/cuda/execution_policy.h>

#include <propr/data/types.h>

#include <propr/kernels/cuda/dispatch/backend.cuh>
#include <propr/kernels/cuda/detail/backend.cuh>


#include <propr/utils/rcpp_cuda.cuh>
#include <propr/utils/rcpp_checks.h>
#include <propr/utils/rcpp_helpers.h>
#include <propr/utils/cuda_checks.h>


using namespace Rcpp;
using namespace propr;

void dispatch::cuda::wtmRcpp(double& out, const NumericVector& x, const NumericVector& w, propr_context context){
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

void dispatch::cuda::wtvRcpp(double& out, const NumericVector& x, const NumericVector& w, propr_context context) {
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

void centerNumericMatrix(NumericMatrix& out, const NumericMatrix & X, propr_context context){
  CHECK_MATRIX_DIMS(out, X.nrow(), X.ncol());
  
  offset_t d_out_stride; offset_t d_x_stride;
  auto *d_x   = RcppMatrixToDevice<float, REALSXP, true>(X  , d_x_stride  );
  auto *d_out = RcppMatrixToDevice<float, REALSXP, true>(out, d_out_stride);

  int block = 512;
  int grid = (X.ncol() + block - 1) / block;
  propr::detail::cuda::centerNumericMatrix<<<grid,block,0,context.stream>>>(d_out, d_out_stride, d_x, d_x_stride, X.nrow(), X.ncol());
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

  delete centered_mat;
  centered_mat = nullptr;

  CUDA_CHECK(cudaFree(d_x));
  CUDA_CHECK(cudaFree(d_out));
}

void dispatch::cuda::corRcpp(NumericMatrix& out, const NumericMatrix & X, propr_context context) {
  CHECK_MATRIX_DIMS(out, X.ncol(), X.ncol());
  using Config = propr::detail::cuda::cor_config;
  int t = 128;
  auto X_padded = rcpp::helpers::pad_matrix(X, 0, ((X.nrow() + t - 1)/t)*t - X.nrow(), 0, ((X.ncol() + t - 1)/t)*t - X.ncol());
  int nfeats  = X_padded.ncol();
  int samples = X_padded.nrow();

  offset_t X_stride;
  auto *d_X = RcppMatrixToDevice<float>(X_padded, X_stride);

  float* d_out = nullptr;
  offset_t dout_stride = nfeats;
  CUDA_CHECK(cudaMalloc(&d_out, nfeats * dout_stride * sizeof(*d_out)));
  
  dim3 block(Config::BLK_M / Config::TH_X, Config::BLK_M / Config::TH_Y);
  dim3 grid(nfeats / Config::BLK_M, nfeats / Config::BLK_M);
  
  propr::detail::cuda::corRcpp<Config><<<grid, block, 0, context.stream>>>(d_X, d_out, nfeats, samples);
  CUDA_CHECK(cudaStreamSynchronize(context.stream));

  auto h_full = new float[nfeats * nfeats];
  CUDA_CHECK(cudaMemcpy(
      h_full,
      d_out,
      nfeats * nfeats * sizeof(float),
      cudaMemcpyDeviceToHost
  ));

  delete h_full;
  h_full = nullptr;
  CUDA_CHECK(cudaFree(d_W));
  CUDA_CHECK(cudaFree(d_out));
}

void dispatch::cuda::covRcpp(NumericMatrix& out, const NumericMatrix & X, const int norm_type, propr_context context) {
  CHECK_MATRIX_DIMS(out, X.ncol(), X.ncol());
  Rcpp::stop("covRcpp is not implemented in CUDA. Falling back to CPU via dispatcher.");
}

void dispatch::cuda::vlrRcpp(NumericMatrix& out, const NumericMatrix & X, propr_context context){
  CHECK_MATRIX_DIMS(out, X.ncol(), X.ncol());
  Rcpp::stop("vlrRcpp is not implemented in CUDA. Falling back to CPU via dispatcher.");
}

void dispatch::cuda::clrRcpp(NumericMatrix& out, const NumericMatrix & X, propr_context context){
  CHECK_MATRIX_DIMS(out, X.nrow(), X.ncol());
  Rcpp::stop("clrRcpp is not implemented in CUDA. Falling back to CPU via dispatcher.");
}

void dispatch::cuda::alrRcpp(NumericMatrix& out, const NumericMatrix & X, const int ivar, propr_context context){
  if(ivar == 0) Rcpp::stop("Select non-zero ivar for alrRcpp.");
  CHECK_MATRIX_DIMS(out, X.nrow(), X.ncol());
  Rcpp::stop("alrRcpp is not implemented in CUDA. Falling back to CPU via dispatcher.");
}

void dispatch::cuda::symRcpp(NumericMatrix& out, const NumericMatrix & X, propr_context context) {
  CHECK_MATRIX_DIMS(out, X.nrow(), X.ncol());
  Rcpp::stop("symRcpp is not implemented in CUDA. Falling back to CPU via dispatcher.");
}

void dispatch::cuda::phiRcpp(NumericMatrix& out, const NumericMatrix &X, const bool sym, propr_context context) {
  CHECK_MATRIX_DIMS(out, X.ncol(), X.ncol());
  Rcpp::stop("phiRcpp is not implemented in CUDA. Falling back to CPU via dispatcher.");
}

void dispatch::cuda::rhoRcpp(NumericMatrix& out, const NumericMatrix &X, const NumericMatrix &lr, const int ivar, propr_context context){
  CHECK_MATRIX_DIMS(out, X.ncol(), X.ncol());
  Rcpp::stop("rhoRcpp is not implemented in CUDA. Falling back to CPU via dispatcher.");
}


void dispatch::cuda::indexPairs(std::vector<int>& out,
                                const NumericMatrix & X,
                                const String op,
                                const double ref,
                                propr_context context) {
}

void dispatch::cuda::indexToCoord(List& out, IntegerVector V, int N, propr_context context){
    Rcpp::stop("indexToCoord is not implemented in CUDA. Falling back to CPU via dispatcher.");
}

void dispatch::cuda::coordToIndex(IntegerVector& out, IntegerVector row, IntegerVector col, int N, propr_context context){
    CHECK_VECTOR_SIZE(out, row.length());
    Rcpp::stop("coordToIndex is not implemented in CUDA. Falling back to CPU via dispatcher.");
}

void dispatch::cuda::linRcpp(NumericMatrix& out, const NumericMatrix & rho, const NumericMatrix &lr, propr_context context){
    CHECK_MATRIX_DIMS(out, rho.nrow(), rho.ncol());
    Rcpp::stop("linRcpp is not implemented in CUDA. Falling back to CPU via dispatcher.");
}

void dispatch::cuda::lltRcpp(NumericVector& out, const NumericMatrix & X, propr_context context){
    int nfeats = X.nrow();
    int llt = nfeats * (nfeats - 1) / 2;
    CHECK_VECTOR_SIZE(out, llt);

    auto* d_out = RcppVectorToDevice<float>(out, llt);
    offset_t d_x_stride;
    auto *d_x = RcppMatrixToDevice<float, REALSXP>(X, d_x_stride);
    
    int block = 512;
    int grid = (llt + block - 1) / block;
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
    int block = 512;
    int grid = (llt + block - 1) / block;
    propr::detail::cuda::lltRcpp<<<grid, block,0,context.stream>>>(d_out, llt, d_x, d_x_stride);
    CUDA_CHECK(cudaStreamSynchronize(context.stream));
    copyToNumericVector(d_out, out, llt);
    CUDA_CHECK(cudaFree(d_x));
    CUDA_CHECK(cudaFree(d_out));
}

void dispatch::cuda::labRcpp(List & out, int nfeats, propr_context context){
  int llt = nfeats * (nfeats - 1) / 2;

  int *d_partner; int *d_pair;
  CUDA_CHECK(cudaMalloc(&d_partner, sizeof(*d_partner)));
  CUDA_CHECK(cudaMalloc(&d_pair, sizeof(*d_pair)));

  int block = 512;
  int grid = (llt + block - 1) / block;
  propr::detail::cuda::labRcpp<<<grid, block,0,context.stream>>>(d_partner, d_pair, llt);
  CUDA_CHECK(cudaStreamSynchronize(context.stream));

  out["Partner"] = Rcpp::IntegerVector(llt);
  out["Pair"] = Rcpp::IntegerVector(llt);

  Rcpp::IntegerVector partner_ref = out["Partner"];
  Rcpp::IntegerVector pair_ref = out["Pair"];

  copyToNumericVector(d_partner,partner_ref, llt);
  copyToNumericVector(pair_ref, pair_ref, llt);

  CUDA_CHECK(cudaFree(d_partner));
  CUDA_CHECK(cudaFree(d_pair));
}

void dispatch::cuda::half2mat(NumericMatrix& out, const NumericVector & X, propr_context context){
    int nfeats = round(sqrt(2 * X.length() + 0.25) + 0.5); // Re-calculate nfeats from X.length()
    CHECK_MATRIX_DIMS(out, nfeats, nfeats);
    Rcpp::stop("half2mat is not implemented in CUDA. Falling back to CPU via dispatcher.");
}

void dispatch::cuda::vector2mat(NumericMatrix& out, const NumericVector & X, const IntegerVector & i, const IntegerVector & j, int nfeats, propr_context context){
    CHECK_MATRIX_DIMS(out, nfeats, nfeats);
    Rcpp::stop("vector2mat is not implemented in CUDA. Falling back to CPU via dispatcher.");
}

void dispatch::cuda::ratiosRcpp(NumericMatrix & out, const NumericMatrix & X, propr_context context){
    int nfeats = X.ncol();
    int nsamps = X.nrow();
    int llt = nfeats * (nfeats - 1) / 2;
    CHECK_MATRIX_DIMS(out, nsamps, llt);

    offset_t d_out_stride; auto *d_out = RcppMatrixToDevice<float, REALSXP>(out, d_out_stride);
    offset_t d_x_stride  ; auto *d_x   = RcppMatrixToDevice<float, REALSXP>(X, d_x_stride);
    
    int block = 512;
    int grid = (llt * nsamps + block - 1) / block;
    propr::detail::cuda::ratiosRcpp<<<grid,block,0,context.stream>>>(d_out, d_out_stride, dx, d_x_stride, nfeat, nsamps);
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
    delete out_host;
    CUDA_CHECK(cudaFree(d_out));
    CUDA_CHECK(cudaFree(d_x));
}

void dispatch::cuda::results2matRcpp(NumericMatrix & out, const DataFrame& results, int n, double diagonal, propr_context context){
    CHECK_MATRIX_DIMS(out, n, n);
    Rcpp::stop("results2matRcpp is not implemented in CUDA. Falling back to CPU via dispatcher.");
}