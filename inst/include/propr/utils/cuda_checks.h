#pragma once

#include <Rcpp.h>

#include <cublas_v2.h> 
#include <cuda_runtime.h>


#define CUDA_CHECK(call) \
    do { \
        cudaError_t err = call; \
        if (err != cudaSuccess) { \
            fprintf(stderr, "CUDA error at %s:%d - %s\n", __FILE__, __LINE__, cudaGetErrorString(err)); \
            exit(EXIT_FAILURE); \
        } \
    } while (0)


#define CUBLAS_CHECK(call)                                                                        \
  do {                                                                                            \
    cublasStatus_t status = (call);                                                               \
    if (status != CUBLAS_STATUS_SUCCESS) {                                                        \
      std::string error_msg;                                                                      \
      switch (status) {                                                                           \
        case CUBLAS_STATUS_SUCCESS:          error_msg = "CUBLAS_STATUS_SUCCESS"; break;          \
        case CUBLAS_STATUS_NOT_INITIALIZED:  error_msg = "CUBLAS_STATUS_NOT_INITIALIZED"; break;  \
        case CUBLAS_STATUS_ALLOC_FAILED:     error_msg = "CUBLAS_STATUS_ALLOC_FAILED"; break;     \
        case CUBLAS_STATUS_INVALID_VALUE:    error_msg = "CUBLAS_STATUS_INVALID_VALUE"; break;    \
        case CUBLAS_STATUS_ARCH_MISMATCH:    error_msg = "CUBLAS_STATUS_ARCH_MISMATCH"; break;    \
        case CUBLAS_STATUS_MAPPING_ERROR:    error_msg = "CUBLAS_STATUS_MAPPING_ERROR"; break;    \
        case CUBLAS_STATUS_EXECUTION_FAILED: error_msg = "CUBLAS_STATUS_EXECUTION_FAILED"; break; \
        case CUBLAS_STATUS_INTERNAL_ERROR:   error_msg = "CUBLAS_STATUS_INTERNAL_ERROR"; break;   \
        case CUBLAS_STATUS_NOT_SUPPORTED:    error_msg = "CUBLAS_STATUS_NOT_SUPPORTED"; break;    \
        case CUBLAS_STATUS_LICENSE_ERROR:    error_msg = "CUBLAS_STATUS_LICENSE_ERROR"; break;    \
        default:                             error_msg = "Unknown cuBLAS error"; break;           \
      }                                                                                           \
      std::cerr << "cuBLAS error at " << __FILE__ << ":" << __LINE__ << ": "                      \
                << error_msg << " (" << status << ")" << std::endl;                               \
      std::exit(EXIT_FAILURE);                                                                    \
    }                                                                                             \
  } while (0)



static void checkCublas(cublasStatus_t st, const char* msg) {
    if (st != CUBLAS_STATUS_SUCCESS) {
        fprintf(stderr, "[cuBLAS ERROR] %s\n", msg);
        exit(-1);
    }
}