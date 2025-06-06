#pragma once

#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <vector>
#include <Rcpp.h>


// __constant__ const float FLT_MAX = 3.402823466e+38f;
// __constant__ const float FLT_MIN = 1.175494351e-38f;


#define CHECK_MATRIX_DIMS(mat, expected_rows, expected_cols) \
    if (mat.nrow() != expected_rows || mat.ncol() != expected_cols) { \
        Rcpp::stop("Error: Output matrix dimensions do not match expected dimensions."); \
    }

#define CHECK_VECTOR_SIZE(vec, expected_size) \
    if (vec.length() != expected_size) { \
        Rcpp::stop("Error: Output vector size does not match expected size."); \
    }


#define CUDA_CHECK(call) \
    do { \
        cudaError_t err = call; \
        if (err != cudaSuccess) { \
            fprintf(stderr, "CUDA error at %s:%d - %s\n", __FILE__, __LINE__, cudaGetErrorString(err)); \
            exit(EXIT_FAILURE); \
        } \
    } while (0)


static void checkCublas(cublasStatus_t st, const char* msg) {
    if (st != CUBLAS_STATUS_SUCCESS) {
        fprintf(stderr, "[cuBLAS ERROR] %s\n", msg);
        exit(-1);
    }
}


inline float* RcppNumericMatrixToDeviceFloat(Rcpp::NumericMatrix& mat, int& num_rows, int& num_cols) {
    num_rows = mat.nrow();
    num_cols = mat.ncol();
    float* d_ptr;
    CUDA_CHECK(cudaMalloc(&d_ptr, num_rows * num_cols * sizeof(float)));
    std::vector<float> h_vec(num_rows * num_cols);
    for (int j = 0; j < num_cols; ++j) {
        for (int i = 0; i < num_rows; ++i) {
            h_vec[i + j * num_rows] = static_cast<float>(mat(i, j));
        }
    }
    CUDA_CHECK(cudaMemcpy(d_ptr, h_vec.data(), num_rows * num_cols * sizeof(float), cudaMemcpyHostToDevice));
    return d_ptr;
}

inline void copyFloatToNumericVector(float* d_src, Rcpp::NumericVector& h_dest, int size) {
    std::vector<float> h_temp(size);
    CUDA_CHECK(cudaMemcpy(h_temp.data(), d_src, size * sizeof(float), cudaMemcpyDeviceToHost));
    for (int i = 0; i < size; ++i) {
        h_dest[i] = static_cast<double>(h_temp[i]);
    }
}