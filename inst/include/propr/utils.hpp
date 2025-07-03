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

inline void check_alignment(const void* ptr, int alignment) {
    if ((uintptr_t)ptr % (alignment * sizeof(float)) != 0) {
        printf("ERROR: Misaligned access at %p, required alignment: %d bytes\n", 
               ptr, alignment * (int)sizeof(float));
    }
}

template <typename OutT, int RTYPE, const bool RowMajor=false>
inline OutT* RcppMatrixToDevice(
    const Rcpp::Matrix<RTYPE>& mat,
    int& memory_stride,
    int alignment = 16
) {
    const int nrows = mat.nrow();
    const int ncols = mat.ncol();

    constexpr bool ColMajor = !RowMajor;
    const int slow_orig = ColMajor ? nrows : ncols;
    const int fast_orig = ColMajor ? ncols : nrows;

    const int slow_padded = ((slow_orig + alignment - 1) / alignment) * alignment;
    const size_t total_elems = size_t(slow_padded) * fast_orig;
    const size_t total_bytes = total_elems * sizeof(OutT);

    OutT* d_ptr = nullptr;
    CUDA_CHECK(cudaMalloc(&d_ptr, total_bytes));
    CUDA_CHECK(cudaMemset(d_ptr, 0, total_bytes));

    std::vector<OutT> h_buf(total_elems, OutT(0));
    if constexpr (ColMajor) {
        for (int j = 0; j < ncols; ++j)
            for (int i = 0; i < nrows; ++i)
                h_buf[i + j * slow_padded] = static_cast<OutT>(mat(i, j));
    } else {
        for (int i = 0; i < nrows; ++i)
            for (int j = 0; j < ncols; ++j)
                h_buf[j + i * slow_padded] = static_cast<OutT>(mat(i, j));
    }

    CUDA_CHECK(cudaMemcpy(d_ptr, h_buf.data(), total_bytes, cudaMemcpyHostToDevice));
    memory_stride = slow_padded;
    check_alignment(d_ptr, alignment);
    return d_ptr;
}



template <typename OutT, int RTYPE, const bool RowMajor=false>
inline OutT* RcppMatrixPermToDevice(
    const Rcpp::Matrix<RTYPE>& mat,
    const Rcpp::IntegerVector& perm,
    int& memory_stride,
    int alignment = 16
) {
    const int nrows = mat.nrow();
    const int ncols = mat.ncol();

    constexpr bool ColMajor = !RowMajor;
    const int slow_orig = ColMajor ? nrows : ncols;
    const int fast_orig = ColMajor ? ncols : nrows;

    const int slow_padded = ((slow_orig + alignment - 1) / alignment) * alignment;
    const size_t total_elems = size_t(slow_padded) * fast_orig;
    const size_t total_bytes = total_elems * sizeof(OutT);

    OutT* d_ptr = nullptr;
    CUDA_CHECK(cudaMalloc(&d_ptr, total_bytes));
    CUDA_CHECK(cudaMemset(d_ptr, 0, total_bytes));

    std::vector<OutT> h_buf(total_elems, OutT(0));
    if constexpr (ColMajor) {
        for (int j = 0; j < ncols; ++j)
            for (int i = 0; i < nrows; ++i)
                h_buf[i + j * slow_padded] = static_cast<OutT>(mat(perm[i], perm[j]));
    } else {
        for (int i = 0; i < nrows; ++i)
            for (int j = 0; j < ncols; ++j)
                h_buf[j + i * slow_padded] = static_cast<OutT>(mat(perm[i], perm[j]));
    }

    CUDA_CHECK(cudaMemcpy(d_ptr, h_buf.data(), total_bytes, cudaMemcpyHostToDevice));
    memory_stride = slow_padded;
    check_alignment(d_ptr, alignment);
    return d_ptr;
}



template<typename T, int RTYPE>
void copyToNumericVector(const T* d_src, Rcpp::Vector<RTYPE>& h_dest, int size) {
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

template<typename T, int RTYPE>
T* RcppVectorToDevice(const Rcpp::Vector<RTYPE>& h_src, int size) {
    using SrcType = typename Rcpp::traits::storage_type<RTYPE>::type;
    
    T* d_ptr     = nullptr;
    size_t bytes = static_cast<size_t>(size) * sizeof(T);
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&d_ptr), bytes));
    if constexpr (std::is_same<T, SrcType>::value) {
        CUDA_CHECK(cudaMemcpy(d_ptr,
                              static_cast<const void*>(h_src.begin()),
                              bytes,
                              cudaMemcpyHostToDevice));
    } else {
        std::vector<T> temp(size);
        for (int i = 0; i < size; ++i) {
            temp[i] = static_cast<T>(h_src[i]);
        }
        CUDA_CHECK(cudaMemcpy(d_ptr, temp.data(), bytes,
                              cudaMemcpyHostToDevice));
    }

    return d_ptr;
}