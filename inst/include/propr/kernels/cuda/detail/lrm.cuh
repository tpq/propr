#pragma once

#include <cuda_runtime.h>
#include <propr/utils.hpp>

namespace propr {
    namespace detail {
        namespace cuda {

            extern "C" __global__
            void lrm_basic(float* __restrict__ d_Y, 
                                float* __restrict__ d_mean, 
                                int nb_samples, 
                                int nb_genes);


            extern "C" __global__
            void lrm_weighted(float* __restrict__ d_Y, 
                              float* __restrict__ d_W,
                              float* __restrict__ d_mean, 
                              int nb_samples, 
                              int nb_genes);


            extern "C" __global__
            void lrm_alpha(float* __restrict__ d_Y,
                            float* __restrict__ d_Yfull,
                            int N1, int NT,
                            float a,
                            float* __restrict__ d_means,
                            int nb_samples, 
                            int nb_genes);



            extern "C" __global__
            void lrm_alpha_weighted(float* __restrict__ d_Y,
                                    float* __restrict__ d_Yfull,
                                    float* __restrict__ d_W,
                                    float* __restrict__ d_Wfull,
                                    int N1, int NT,
                                    float a,
                                    float* __restrict__ d_means,
                                    int nb_genes);
            
        }
    }
}