#pragma once

#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <vector>
#include <Rcpp.h>

#define IS_POWER_OF_2(x) (((x) != 0) && (((x) & ((x) - 1)) == 0))

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

inline void check_alignment(const void* ptr, int alignment) {
    if ((uintptr_t)ptr % (alignment * sizeof(float)) != 0) {
        printf("ERROR: Misaligned access at %p, required alignment: %d bytes\n", 
               ptr, alignment * (int)sizeof(float));
    }
}

inline float* RcppNumericMatrixToDeviceFloat(Rcpp::NumericMatrix& mat, int& memory_stride, int alignment = 16) {
    int original_rows = mat.nrow();
    int num_cols      = mat.ncol(); 
    
    int padded_rows = ((original_rows + alignment - 1) / alignment) * alignment;
    size_t total_size = padded_rows * num_cols * sizeof(float);
    float* d_ptr;
    CUDA_CHECK(cudaMalloc(&d_ptr, total_size));
    CUDA_CHECK(cudaMemset(d_ptr, 0, total_size));    
    std::vector h_vec(padded_rows * num_cols, 0.0f);
    for (int j = 0; j < num_cols; ++j) {
        for (int i = 0; i < original_rows; ++i) {
            h_vec[i + j * padded_rows] = static_cast<float>(mat(i, j));
        }
    }
    CUDA_CHECK(cudaMemcpy(d_ptr, h_vec.data(), total_size, cudaMemcpyHostToDevice));    
    memory_stride = padded_rows;
    check_alignment(d_ptr, alignment);
    return d_ptr;
}

template<typename T>
void copyToNumericVector(const T* d_src,
                                Rcpp::NumericVector& h_dest,
                                int size) {
    std::vector<T> h_temp(size);
    CUDA_CHECK(cudaMemcpy(
        h_temp.data(), d_src,
        size * sizeof(T),
        cudaMemcpyDeviceToHost
    ));
    for (int i = 0; i < size; ++i) {
        h_dest[i] = static_cast<double>(h_temp[i]);
    }
}
