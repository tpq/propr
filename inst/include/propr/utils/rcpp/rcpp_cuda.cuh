#pragma once

#include <Rcpp.h>

#include <propr/data/numeric_conversion.cuh>
#include <propr/data/types.h>
#include <propr/utils/cuda/cuda_checks.h>
#include <propr/utils/cuda/alignment.h>

template <
  typename OutT,
  int RTYPE,
  bool RowMajor = false
>
inline OutT* RcppMatrixToDevice(
    const Rcpp::Matrix<RTYPE>& mat,
    offset_t&                  memory_stride,
    int                        alignment = 16) {

    propr::convert::NumericConverter<OutT, typename Rcpp::Matrix<RTYPE>::stored_type> CastOp;

    const offset_t nrows = mat.nrow();
    const offset_t ncols = mat.ncol();

    constexpr bool ColMajor = !RowMajor;
    const offset_t slow_orig = ColMajor ? nrows : ncols;
    const offset_t fast_orig = ColMajor ? ncols : nrows;

    const offset_t slow_padded = ((slow_orig + alignment - 1) / alignment) * alignment;
    const offset_t total_elems  = slow_padded * fast_orig;
    const size_t  total_bytes  = total_elems * sizeof(OutT);

    OutT* d_ptr = nullptr;
    PROPR_CUDA_CHECK(cudaMalloc(&d_ptr, total_bytes));
    PROPR_CUDA_CHECK(cudaMemset(d_ptr, 0, total_bytes));

    OutT* h_buf = static_cast<OutT*>(std::malloc(total_bytes));
    if (!h_buf) throw std::bad_alloc();
    std::memset(h_buf, 0, total_bytes);

    if constexpr (ColMajor) {
        for (offset_t j = 0; j < ncols; ++j) {
            for (offset_t i = 0; i < nrows; ++i) {
                const offset_t idx = i + j * slow_padded;
                h_buf[idx] = CastOp( mat(i, j) );
            }
        }
    } else {
        for (offset_t i = 0; i < nrows; ++i) {
            for (offset_t j = 0; j < ncols; ++j) {
                const offset_t idx = j + i * slow_padded;
                h_buf[idx] = CastOp( mat(i, j) );
            }
        }
    }

    PROPR_CUDA_CHECK(cudaMemcpy(d_ptr, h_buf, total_bytes, cudaMemcpyHostToDevice));
    std::free(h_buf);

    memory_stride = slow_padded;
    propr::cuda::check_pointer_alignment(d_ptr, alignment);
    return d_ptr;
}

template <typename OutT, int RTYPE, const bool RowMajor=false>
inline OutT* RcppMatrixPermToDevice(
    const Rcpp::Matrix<RTYPE>& mat,
    const Rcpp::IntegerVector& perm,
    offset_t& memory_stride,
    int alignment = 16
) {

    propr::convert::NumericConverter<OutT, typename Rcpp::Matrix<RTYPE>::stored_type> CastOp;

    const offset_t nrows = mat.nrow();
    const offset_t ncols = mat.ncol();

    constexpr bool ColMajor = !RowMajor;
    const offset_t slow_orig = ColMajor ? nrows : ncols;
    const offset_t fast_orig = ColMajor ? ncols : nrows;

    const offset_t slow_padded =
        ((slow_orig + alignment - 1) / alignment) * alignment;
    const offset_t total_elems = slow_padded * fast_orig;
    const size_t  total_bytes = total_elems * sizeof(OutT);

    OutT* d_ptr = nullptr;
    PROPR_CUDA_CHECK(cudaMalloc(&d_ptr, total_bytes));
    PROPR_CUDA_CHECK(cudaMemset(d_ptr, 0, total_bytes));

    OutT* h_buf = static_cast<OutT*>(std::malloc(total_bytes));
    if (!h_buf) throw std::bad_alloc();
    std::memset(h_buf, 0, total_bytes);

    if constexpr (ColMajor) {
        for (offset_t j = 0; j < ncols; ++j) {
            for (offset_t i = 0; i < nrows; ++i) {
                const offset_t idx = i + j * slow_padded;
                h_buf[idx] = CastOp(mat(perm[i], perm[j]));
            }
        }
    } else {
        for (offset_t i = 0; i < nrows; ++i) {
            for (offset_t j = 0; j < ncols; ++j) {
                const offset_t idx = j + i * slow_padded;
                h_buf[idx] = CastOp(mat(perm[i], perm[j]));
            }
        }
    }

    PROPR_CUDA_CHECK(cudaMemcpy(d_ptr, h_buf, total_bytes, cudaMemcpyHostToDevice));
    std::free(h_buf);

    memory_stride = slow_padded;
    propr::cuda::check_pointer_alignment(d_ptr, alignment);
    return d_ptr;
}


template<typename T, int RTYPE>
void copyToNumericVector(
    const T* d_src,
    Rcpp::Vector<RTYPE>& h_dest,
    const size_t size
) {
    using DstType = typename Rcpp::traits::storage_type<RTYPE>::type;
    
    const size_t bytes = size * sizeof(T);
    T* h_temp = static_cast<T*>(std::malloc(bytes));
    if (!h_temp) throw std::bad_alloc();

    PROPR_CUDA_CHECK(cudaMemcpy(h_temp, d_src, bytes, cudaMemcpyDeviceToHost));
    for (offset_t i = 0; i < size; ++i) {
        h_dest[i] = static_cast<DstType>(h_temp[i]);
    }
    std::free(h_temp);
}


template<typename T, int RTYPE>
T* RcppVectorToDevice(const Rcpp::Vector<RTYPE>& h_src, size_t size) {
    using SrcType = typename Rcpp::traits::storage_type<RTYPE>::type;
    propr::convert::NumericConverter<T, typename Rcpp::Vector<RTYPE>::stored_type> CastOp;
    T* d_ptr     = nullptr;
    const size_t bytes = size * sizeof(T);
    PROPR_CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&d_ptr), bytes));
    if constexpr (std::is_same_v<T, SrcType>) {
        PROPR_CUDA_CHECK(cudaMemcpy(d_ptr,
                              static_cast<const void*>(h_src.begin()),
                              bytes,
                              cudaMemcpyHostToDevice));
    } else {
        T* h_temp = static_cast<T*>(std::malloc(bytes));
        for (int i = 0; i < size; ++i) h_temp[i] = CastOp(h_src[i]);
        PROPR_CUDA_CHECK(cudaMemcpy(d_ptr, h_temp, bytes, cudaMemcpyHostToDevice));
        std::free(h_temp);
    }

    return d_ptr;
}
