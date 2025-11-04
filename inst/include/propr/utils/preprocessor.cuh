
#pragma once

#ifdef __CUDA_ARCH__
    #define PROPR_HOST_DEVICE  __host__ __device__
    #define PROPR_HOST         __host__
    #define PROPR_DEVICE       __device__
#else
    #define PROPR_HOST_DEVICE
    #define PROPR_HOST
    #define PROPR_DEVICE
#endif


#define PROPR_INLINE       inline
#define PROPR_FORCE_INLINE __forceinline__

#if !defined(__CUDACC_RTC__) && (defined(__CUDA_ARCH__) || defined(_NVHPC_CUDA))
    #define PROPR_UNROLL    #pragma unroll
    #define PROPR_NO_UNROLL #pragma unroll 1
#elif defined(__CUDACC_RTC__)
    #define PROPR_UNROLL    _Pragma("unroll")
    #define PROPR_NO_UNROLL _Pragma("unroll 1")
#else
    #define PROPR_UNROLL
    #define PROPR_NO_UNROLL
#endif  


// #define FETCH_FLOAT4(pointer) (reinterpret_cast<float4*>(&(pointer))[0])
#define OFFSET(row, col, ld) ((row) * (ld) + (col))

//cub::ThreadLoad<cub::LOAD_DEFAULT>((float4 *) 