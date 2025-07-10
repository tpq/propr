#pragma once

#include <cuda_runtime.h>
#include <propr/data/types.h>

namespace propr {
    namespace detail {
        namespace cuda {

            __global__
            void
            lrm_basic(float* __restrict__ d_Y, offset_t d_Y_stride,
                      float* __restrict__ d_mean,
                      int nb_samples,
                      int nb_genes) {
                int i = blockIdx.x * blockDim.x + threadIdx.x;
                int j = blockIdx.y * blockDim.y + threadIdx.y;
                if (i >= nb_genes || j >= i) return;

                float4 accum = {0.0f, 0.0f, 0.0f, 0.0f};
                int k = 0;
                #pragma unroll
                for (; k < (nb_samples/4)*4; k += 4) {
                    float4 y_i = *reinterpret_cast<float4*>(&d_Y[k + i * d_Y_stride]);
                    float4 y_j = *reinterpret_cast<float4*>(&d_Y[k + j * d_Y_stride]);
                    accum.x = __fmaf_rn(1.0f, __logf(__fdividef(y_i.x, y_j.x)), accum.x);
                    accum.y = __fmaf_rn(1.0f, __logf(__fdividef(y_i.y, y_j.y)), accum.y);
                    accum.z = __fmaf_rn(1.0f, __logf(__fdividef(y_i.z, y_j.z)), accum.z);
                    accum.w = __fmaf_rn(1.0f, __logf(__fdividef(y_i.w, y_j.w)), accum.w);
                }

                accum.x = accum.x + accum.y + accum.z + accum.w;
                for (; k < nb_samples; ++k) {
                    float yi = d_Y[k + i * d_Y_stride];
                    float yj = d_Y[k + j * d_Y_stride];
                    accum.x  = __fmaf_rn(1.0f, __logf(__fdividef(yi, yj)), accum.x);
                }

                float inv_n = __frcp_rn(static_cast<float>(nb_samples));
                float mean  = accum.x * inv_n;
                int pair_index = (i * (i - 1)) / 2 + j;
                d_mean[pair_index] = mean;
            }


            __global__
            void
            lrm_weighted(float* __restrict__ d_Y, offset_t d_Y_stride,
                         float* __restrict__ d_W, offset_t d_W_stride,
                         float* __restrict__ d_mean,
                         int nb_samples,
                         int nb_genes) {
                int i = blockIdx.x * blockDim.x + threadIdx.x;
                int j = blockIdx.y * blockDim.y + threadIdx.y;
                if (i >= nb_genes || j >= i) return;

                // accum.x = w_sum, accum.y = mean
                float2 accum = make_float2(0.0f, 0.0f);
                int k = 0;
                #pragma unroll
                for (; k < (nb_samples/4)*4; k += 4) {
                    float4 y_i = *reinterpret_cast<float4*>(&d_Y[k + i * d_Y_stride]);
                    float4 y_j = *reinterpret_cast<float4*>(&d_Y[k + j * d_Y_stride]);

                    float4 w_i = *reinterpret_cast<float4*>(&d_W[k + i * d_W_stride]);
                    float4 w_j = *reinterpret_cast<float4*>(&d_W[k + j * d_W_stride]);

                    #pragma unroll
                    for (int m = 0; m < 4; ++m) {
                        float mean_old = accum.y;
                        float w_k      = (&w_i.x)[m] * (&w_j.x)[m];
                        float w        = w_k;
                        accum.x        = __fmaf_rn(1.0f, w, accum.x);

                        float log_val = __logf(__fdividef((&y_i.x)[m], (&y_j.x)[m]));
                        float delta   = log_val - mean_old;
                        float w_ratio = __fdividef(w, accum.x);
                        accum.y       = __fmaf_rn(w_ratio, delta, mean_old);
                    }
                }

                for (; k < nb_samples; ++k) {
                    float y_ik = d_Y[k + i * d_Y_stride];
                    float y_jk = d_Y[k + j * d_Y_stride];

                    float w_ik = d_W[k + i * d_W_stride];
                    float w_jk = d_W[k + j * d_W_stride];

                    float w_k = w_ik * w_jk;

                    float ratio    = __fdividef(y_ik, y_jk);
                    float log_val  = __logf(ratio);
                    float mean_old = accum.y;
                    accum.x       += w_k;

                    float delta   = log_val - mean_old;
                    float w_ratio = w_k / accum.x;
                    accum.y      += w_ratio * delta;
                }

                int pair_index = (i * (i - 1)) / 2 + j;
                d_mean[pair_index] = accum.y;
                
            }


            __global__
            void
            lrm_alpha(float* __restrict__     d_Y, offset_t d_Y_stride,
                      float* __restrict__ d_Yfull, offset_t d_Yfull_stride,
                      int N1, int NT,
                      float a,
                      float* __restrict__ d_means,
                      int nb_samples,
                      int nb_genes) {
                int i = blockIdx.x * blockDim.x + threadIdx.x;
                int j = blockIdx.y * blockDim.y + threadIdx.y;
                if (i >= nb_genes || j >= i) return;

                float mu_full_i = 0.0f, mu_full_j = 0.0f;
                float S_i = 0.0f, S_j = 0.0f;

                float T = 0.0f, D = 0.0f;
                int n, N, k;

                k = 0; N = 0;
                #pragma unroll
                for (; k < (NT/4) * 4; k += 4) {
                    float4 yfull_i = *reinterpret_cast<float4*>(&d_Yfull[k + i * d_Yfull_stride]);
                    float4 yfull_j = *reinterpret_cast<float4*>(&d_Yfull[k + j * d_Yfull_stride]);

                    #pragma unroll
                    for (int m = 0; m < 4; ++m) {
                        float inv_N    = __frcp_rn(static_cast<float>(++N));
                        float X_full_i = __powf(reinterpret_cast<float*>(&yfull_i)[m], a);
                        float X_full_j = __powf(reinterpret_cast<float*>(&yfull_j)[m], a);
                        float delta_mu_i = X_full_i - mu_full_i;
                        float delta_mu_j = X_full_j - mu_full_j;
                        mu_full_i = __fmaf_rn(delta_mu_i, inv_N,mu_full_i);
                        mu_full_j = __fmaf_rn(delta_mu_j, inv_N,mu_full_j);
                        T  = T + X_full_i - X_full_j;
                    }
                }

                for (; k < NT; ++k) {
                    float yfull_i = d_Yfull[k + i * d_Yfull_stride];
                    float yfull_j = d_Yfull[k + j * d_Yfull_stride];
                    float inv_N    = __frcp_rn(static_cast<float>(++N));
                    float X_full_i = __powf(yfull_i, a);
                    float X_full_j = __powf(yfull_j, a);
                    float delta_mu_i = X_full_i - mu_full_i;
                    float delta_mu_j = X_full_j - mu_full_j;
                    mu_full_i =__fmaf_rn(delta_mu_i, inv_N,mu_full_i);
                    mu_full_j =__fmaf_rn(delta_mu_j, inv_N,mu_full_j);
                    T  = T + (X_full_i - X_full_j);
                } T *= NT;
                
                k = 0; n = 0;
                #pragma unroll
                for (; k < (N1/4) * 4; k += 4) {
                    float4 y_i = *reinterpret_cast<float4*>(&d_Y[k + i * d_Y_stride]);
                    float4 y_j = *reinterpret_cast<float4*>(&d_Y[k + j * d_Y_stride]);

                    #pragma unroll
                    for (int m = 0; m < 4; ++m) {
                        float X_i   = __powf(reinterpret_cast<float*>(&y_i)[m], a);
                        float X_j   = __powf(reinterpret_cast<float*>(&y_j)[m], a);
                        S_i = S_i + X_i; 
                        S_j = S_j + X_j;
                        D += (X_i - X_j);
                    }
                }

                for (; k < N1; ++k) {
                    float y_i   = d_Y[k + i * d_Y_stride];
                    float y_j   = d_Y[k + j * d_Y_stride];
                    float X_i   = __powf(y_i, a);
                    float X_j   = __powf(y_j, a);
                    S_i = S_i + X_i; S_j = S_j + X_j;
                    D += (X_i - X_j);
                }


                float complement_term = float((N1 < NT)) * (T - D) / (NT - n);
                float C_z = D / n + complement_term;
                float M_z = (S_i / (N1 * mu_full_i)) - (S_j / (N1 * mu_full_j));

                int pair_index = (i * (i - 1)) / 2 + j;
                d_means[pair_index] =  ( (C_z / 2) + M_z ) / a;

            }



            __global__
            void
            lrm_alpha_weighted( float* __restrict__ d_Y    , offset_t Y_stride,
                                float* __restrict__ d_Yfull, offset_t Yfull_stride,
                                float* __restrict__ d_W    , offset_t W_stride,
                                float* __restrict__ d_Wfull, offset_t Wfull_stride,
                                int N1, int NT,
                                float a,
                                float* __restrict__ d_means,
                                int nb_genes) {
                int i = blockIdx.x * blockDim.x + threadIdx.x;
                int j = blockIdx.y * blockDim.y + threadIdx.y;
                if (i >= nb_genes || j >= i) return;

                // Phase 1:
                float sum_w_full    = 0.0f;
                float sum_wx_full_i = 0.0f;
                float sum_wx_full_j = 0.0f;
                int k = 0;
                for (; k < (NT/4)*4; k += 4) {
                    float4 yfull_i4 = *reinterpret_cast<float4*>(&d_Yfull[k + i * Yfull_stride]);
                    float4 yfull_j4 = *reinterpret_cast<float4*>(&d_Yfull[k + j * Yfull_stride]);
                    float4 wfull_i4 = *reinterpret_cast<float4*>(&d_Wfull[k + i * Wfull_stride]);
                    float4 wfull_j4 = *reinterpret_cast<float4*>(&d_Wfull[k + j * Wfull_stride]);

                    for (int m = 0; m < 4; ++m) {
                        float y_i = reinterpret_cast<float*>(&yfull_i4)[m];
                        float y_j = reinterpret_cast<float*>(&yfull_j4)[m];
                        float w_i = reinterpret_cast<float*>(&wfull_i4)[m];
                        float w_j = reinterpret_cast<float*>(&wfull_j4)[m];

                        float w_ij = w_i * w_j;
                        float X_i = __powf(y_i, a);
                        float X_j = __powf(y_j, a);

                        sum_w_full    += w_ij;
                        sum_wx_full_i += w_ij * X_i;
                        sum_wx_full_j += w_ij * X_j;
                    }
                }
                for (; k < NT; ++k) {
                    float y_i = d_Yfull[k + i * Yfull_stride];
                    float y_j = d_Yfull[k + j * Yfull_stride];
                    float w_i = d_Wfull[k + i * Wfull_stride];
                    float w_j = d_Wfull[k + j * Wfull_stride];

                    float w_ij = w_i * w_j;
                    float X_i = __powf(y_i, a);
                    float X_j = __powf(y_j, a);

                    sum_w_full += w_ij;
                    sum_wx_full_i += w_ij * X_i;
                    sum_wx_full_j += w_ij * X_j;
                }

                float mu_i_full = 0.0f, mu_j_full = 0.0f;
                float T_full = 0.0f;
                if (sum_w_full > 1e-10) {
                    mu_i_full = sum_wx_full_i / sum_w_full;
                    mu_j_full = sum_wx_full_j / sum_w_full;
                    T_full = sum_wx_full_i - sum_wx_full_j;
                }

                // Phase 2:
                float sum_w_current    = 0.0f;
                float sum_wx_current_i = 0.0f;   
                float sum_wx_current_j = 0.0f;   

                k = 0;
                for (; k < (N1/4)*4; k += 4) {
                    float4 y_i4 = *reinterpret_cast<float4*>(&d_Y[k + i * Y_stride]);
                    float4 y_j4 = *reinterpret_cast<float4*>(&d_Y[k + j * Y_stride]);
                    float4 w_i4 = *reinterpret_cast<float4*>(&d_W[k + i * W_stride]);
                    float4 w_j4 = *reinterpret_cast<float4*>(&d_W[k + j * W_stride]);

                    for (int m = 0; m < 4; ++m) {
                        float y_i = reinterpret_cast<float*>(&y_i4)[m];
                        float y_j = reinterpret_cast<float*>(&y_j4)[m];
                        float w_i = reinterpret_cast<float*>(&w_i4)[m];
                        float w_j = reinterpret_cast<float*>(&w_j4)[m];

                        float w_ij = w_i * w_j;
                        float X_i = __powf(y_i, a);
                        float X_j = __powf(y_j, a);

                        sum_w_current += w_ij;
                        sum_wx_current_i += w_ij * X_i;
                        sum_wx_current_j += w_ij * X_j;
                    }
                }

                for (; k < N1; ++k) {
                    float y_i = d_Y[k + i * Y_stride];
                    float y_j = d_Y[k + j * Y_stride];
                    float w_i = d_W[k + i * W_stride];
                    float w_j = d_W[k + j * W_stride];

                    float w_ij = w_i * w_j;
                    float X_i = __powf(y_i, a);
                    float X_j = __powf(y_j, a);

                    sum_w_current += w_ij;
                    sum_wx_current_i += w_ij * X_i;
                    sum_wx_current_j += w_ij * X_j;
                }

                float T_current = sum_wx_current_i - sum_wx_current_j;

                float complement_term = 0.0f;
                float denom_complement = sum_w_full - sum_w_current;
                if (denom_complement > 1e-10) {
                    complement_term = (T_full - T_current) / denom_complement;
                }

                float C_z = 0.0f;
                if (sum_w_current > 1e-10) {
                    C_z = (T_current / sum_w_current) + complement_term;
                } else if (denom_complement > 1e-10) {
                    C_z = T_full / sum_w_full;
                }

                float M_z = 0.0f;
                if (sum_w_current > 1e-10 && mu_i_full > 1e-10 && mu_j_full > 1e-10) {
                    M_z = (sum_wx_current_i / mu_i_full - sum_wx_current_j / mu_j_full) / sum_w_current;
                }

                float result = ((C_z / 2.0f) + M_z) / a;

                int pair_index = (i * (i - 1)) / 2 + j;
                d_means[pair_index] = result;
            }

        }
    }
}