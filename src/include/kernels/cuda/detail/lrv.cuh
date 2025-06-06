#pragma once
#include <cuda_runtime.h>

namespace propr{
    namespace detail {
        namespace cuda {

            __global__
            void lrv_basic(float* __restrict__ d_Y, 
                           float* __restrict__ d_variances, 
                           int nb_samples, int nb_genes);


            __global__
            void lrv_weighted(float* __restrict__ d_Y, 
                              float* __restrict__ d_W,
                              float* __restrict__ d_variances, 
                              int nb_samples, int nb_genes);


            __global__
            void lrv_alpha(float* __restrict__ d_Y,
                           float* __restrict__ d_Yfull,
                           float a,
                           float* __restrict__ d_variances,
                           int nb_samples, int nb_genes);


            __global__
            void lrv_alpha_weighted(float* __restrict__ d_Y,
                                    float* __restrict__ d_Yfull,
                                    float* __restrict__ d_W,
                                    float* __restrict__ d_Wfull,
                                    float a,
                                    float* __restrict__ d_variances,
                                    int nb_samples, int nb_genes);

        }
    }
}