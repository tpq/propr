
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