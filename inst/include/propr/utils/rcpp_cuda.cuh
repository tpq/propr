#include <Rcpp.h>
#include <cuda_fp16.h>

#include <propr/data/types.h>
#include <propr/utils/cuda_checks.h>

template<typename OutT, typename InT>
struct CastOperator {
    OutT operator()(InT x) const { return static_cast<OutT>(x); }
};

template<> struct CastOperator<__half,float> {
    __half operator()(const float x) const { return __float2half(x); }
};
template<> struct CastOperator<float,__half> {
    float operator()(const __half h) const { return __half2float(h); }
};



template<> struct CastOperator<__half,double> {
    __half operator()(const double x) const { return __float2half(static_cast<float>(x)); }
};
template<> struct CastOperator<double,__half> {
    double operator()(const __half h) const { return static_cast<double>(__half2float(h)); }
};



inline void check_alignment(const void* ptr, int alignment) {
    if (reinterpret_cast<uintptr_t>(ptr) % (alignment * sizeof(float)) != 0) {
        printf("ERROR: Misaligned access at %p, required alignment: %d bytes\n",  ptr, alignment * (int)sizeof(float));
    }
}

template <
  typename OutT,
  int RTYPE,
  bool RowMajor = false
>
inline OutT* RcppMatrixToDevice(
    const Rcpp::Matrix<RTYPE>& mat,
    offset_t&                  memory_stride,
    int                        alignment = 16) {

    CastOperator<OutT, typename Rcpp::Matrix<RTYPE>::stored_type> CastOp;

    const offset_t nrows = mat.nrow();
    const offset_t ncols = mat.ncol();

    constexpr bool ColMajor = !RowMajor;
    const offset_t slow_orig = ColMajor ? nrows : ncols;
    const offset_t fast_orig = ColMajor ? ncols : nrows;

    const offset_t slow_padded = ((slow_orig + alignment - 1) / alignment) * alignment;
    const offset_t total_elems  = slow_padded * fast_orig;
    const size_t  total_bytes  = total_elems * sizeof(OutT);

    OutT* d_ptr = nullptr;
    CUDA_CHECK(cudaMalloc(&d_ptr, total_bytes));
    CUDA_CHECK(cudaMemset(d_ptr, 0, total_bytes));

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

    CUDA_CHECK(cudaMemcpy(d_ptr, h_buf, total_bytes, cudaMemcpyHostToDevice));
    std::free(h_buf);

    memory_stride = slow_padded;
    check_alignment(d_ptr, alignment);
    return d_ptr;
}

template <typename OutT, int RTYPE, const bool RowMajor=false>
inline OutT* RcppMatrixPermToDevice(
    const Rcpp::Matrix<RTYPE>& mat,
    const Rcpp::IntegerVector& perm,
    offset_t& memory_stride,
    int alignment = 16
) {

    CastOperator<OutT, typename Rcpp::Matrix<RTYPE>::stored_type> CastOp;

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
    CUDA_CHECK(cudaMalloc(&d_ptr, total_bytes));
    CUDA_CHECK(cudaMemset(d_ptr, 0, total_bytes));

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

    CUDA_CHECK(cudaMemcpy(d_ptr, h_buf, total_bytes, cudaMemcpyHostToDevice));
    std::free(h_buf);

    memory_stride = slow_padded;
    check_alignment(d_ptr, alignment);
    return d_ptr;
}


template<typename T, int RTYPE>
void copyToNumericVector(
    const T* d_src,
    Rcpp::Vector<RTYPE>& h_dest,
    const size_t size
) {
    const size_t bytes = size * sizeof(T);
    T* h_temp = static_cast<T*>(std::malloc(bytes));
    if (!h_temp) throw std::bad_alloc();

    CUDA_CHECK(cudaMemcpy(h_temp, d_src, bytes, cudaMemcpyDeviceToHost));
    for (offset_t i = 0; i < size; ++i) {
        h_dest[i] = static_cast<double>(h_temp[i]);
    }
    std::free(h_temp);
}


template<typename T, int RTYPE>
T* RcppVectorToDevice(const Rcpp::Vector<RTYPE>& h_src, size_t size) {
    using SrcType = typename Rcpp::traits::storage_type<RTYPE>::type;
    CastOperator<T, typename Rcpp::Vector<RTYPE>::stored_type> CastOp;

    T* d_ptr     = nullptr;
    const size_t bytes = size * sizeof(T);
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&d_ptr), bytes));
    if constexpr (std::is_same_v<T, SrcType>) {
        CUDA_CHECK(cudaMemcpy(d_ptr,
                              static_cast<const void*>(h_src.begin()),
                              bytes,
                              cudaMemcpyHostToDevice));
    } else {
        std::vector<T> temp(size);
        for (int i = 0; i < size; ++i) temp[i] = CastOp(h_src[i]);
        CUDA_CHECK(cudaMemcpy(d_ptr, temp.data(), bytes,
                              cudaMemcpyHostToDevice));
    }

    return d_ptr;
}