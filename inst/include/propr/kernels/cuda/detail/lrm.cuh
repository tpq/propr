#pragma once

#include <cuda_runtime.h>
#include <propr/data/types.h>
#include <propr/utils/preprocessor.cuh>
#include <propr/internal/device/cuda/thread/mem_ops.cuh>


using namespace propr::cuda::internal;


namespace propr {
    namespace detail {
        namespace cuda {

            template <class Config>
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
                PROPR_UNROLL
                for (; k < (nb_samples/4)*4; k += 4) {
                    float4 y_i = thread::load<Config::LoadModifer,float4>(&d_Y[k + i * d_Y_stride]);
                    float4 y_j = thread::load<Config::LoadModifer,float4>(&d_Y[k + j * d_Y_stride]);                    
                    
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

            template<class Config>
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

                PROPR_UNROLL
                for (; k < (nb_samples / 4) * 4; k += 4) {
                    float4 y_i = thread::load<Config::LoadModifer,float4>(&d_Y[k + i * d_Y_stride]);
                    float4 y_j = thread::load<Config::LoadModifer,float4>(&d_Y[k + j * d_Y_stride]);

                    float4 w_i = thread::load<Config::LoadModifer,float4>(&d_W[k + i * d_W_stride]);
                    float4 w_j = thread::load<Config::LoadModifer,float4>(&d_W[k + j * d_W_stride]);

                    PROPR_UNROLL
                    for (int m = 0; m < 4; ++m) {
                        float mean_old = accum.y;

                        float w_im  = (&w_i.x)[m];
                        float w_jm  = (&w_j.x)[m];
                        float denom = w_im + w_jm;
                        float w     = (denom > 0.0f) ? (2.0f * w_im * w_jm / denom) : 0.0f;

                        accum.x = __fmaf_rn(1.0f, w, accum.x);

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

                    float denom = w_ik + w_jk;
                    float w_k   = (denom > 0.0f) ? (2.0f * w_ik * w_jk / denom) : 0.0f;

                    float ratio    = __fdividef(y_ik, y_jk);
                    float log_val  = __logf(ratio);
                    float mean_old = accum.y;

                    accum.x += w_k;

                    float delta   = log_val - mean_old;
                    float w_ratio = w_k / accum.x;
                    accum.y      += w_ratio * delta;
                }

                int pair_index = (i * (i - 1)) / 2 + j;
                d_mean[pair_index] = accum.y;
            }

            template<class Config>
            __global__
            void
            lrm_alpha(float* __restrict__     d_Y,         size_t d_Y_stride,
                    float* __restrict__     d_Yfull,     size_t d_Yfull_stride,
                    int                      N1,
                    int                      NT,
                    float                    a,
                    float* __restrict__     d_means,
                    int                      nb_samples,
                    int                      nb_genes) {
                int i = blockIdx.x * blockDim.x + threadIdx.x;
                int j = blockIdx.y * blockDim.y + threadIdx.y;
                if (i >= nb_genes || j >= i) return;

                float mu_full_i = 0.0f, mu_full_j = 0.0f;
                float S_i = 0.0f, S_j = 0.0f;

                float T = 0.0f, D = 0.0f;
                int N = 0, k = 0;

                // --- full-group running means + T sum ---
                PROPR_UNROLL
                for (; k + 3 < NT; k += 4) {
                    float4 yfull_i = thread::load<Config::LoadModifer,float4>(&d_Yfull[k + i * d_Yfull_stride]);
                    float4 yfull_j = thread::load<Config::LoadModifer,float4>(&d_Yfull[k + j * d_Yfull_stride]);
                    PROPR_UNROLL
                    for (int m = 0; m < 4; ++m) {
                        float inv_N    = __frcp_rn(static_cast<float>(++N));
                        float X_full_i = __powf(reinterpret_cast<float*>(&yfull_i)[m], a);
                        float X_full_j = __powf(reinterpret_cast<float*>(&yfull_j)[m], a);
                        mu_full_i = __fmaf_rn(X_full_i - mu_full_i, inv_N, mu_full_i);
                        mu_full_j = __fmaf_rn(X_full_j - mu_full_j, inv_N, mu_full_j);
                        T += (X_full_i - X_full_j);
                    }
                }
                for (; k < NT; ++k) {
                    float inv_N    = __frcp_rn(static_cast<float>(++N));
                    float X_full_i = __powf(d_Yfull[k + i * d_Yfull_stride], a);
                    float X_full_j = __powf(d_Yfull[k + j * d_Yfull_stride], a);
                    mu_full_i = __fmaf_rn(X_full_i - mu_full_i, inv_N, mu_full_i);
                    mu_full_j = __fmaf_rn(X_full_j - mu_full_j, inv_N, mu_full_j);
                    T += (X_full_i - X_full_j);
                }

                k = 0;
                PROPR_UNROLL
                for (; k + 3 < N1; k += 4) {
                    float4 y_i = thread::load<Config::LoadModifer,float4>(&d_Y[k + i * d_Y_stride]);
                    float4 y_j = thread::load<Config::LoadModifer,float4>(&d_Y[k + j * d_Y_stride]);
                    PROPR_UNROLL
                    for (int m = 0; m < 4; ++m) {
                        float X_i = __powf(reinterpret_cast<float*>(&y_i)[m], a);
                        float X_j = __powf(reinterpret_cast<float*>(&y_j)[m], a);
                        S_i += X_i; 
                        S_j += X_j;
                        D   += (X_i - X_j);
                    }
                }
                for (; k < N1; ++k) {
                    float X_i = __powf(d_Y[k + i * d_Y_stride], a);
                    float X_j = __powf(d_Y[k + j * d_Y_stride], a);
                    S_i += X_i; 
                    S_j += X_j;
                    D   += (X_i - X_j);
                }

                // --- final combination using n == N1 ---
                float complement_term = float((N1 < NT)) * (T - D) / (NT - N1);
                float C_z = D / float(N1) + complement_term;
                float M_z = (S_i/mu_full_i - S_j/mu_full_j) / float(N1);

                int pair_index = (i * (i - 1)) / 2 + j;
                d_means[pair_index] = ((C_z / 2) + M_z) / a;
            }


            template<class Config>
            __global__
            void
            lrm_alpha_weighted( float* __restrict__ d_Y    , offset_t Y_stride,
                                float* __restrict__ d_Yfull, offset_t Yfull_stride,
                                float* __restrict__ d_W    , offset_t W_stride,
                                float* __restrict__ d_Wfull, offset_t Wfull_stride,
                                int N1, int NT,
                                float a,
                                float* __restrict__ d_means,
                                int nb_genes) 
            {
                int i = blockIdx.x * blockDim.x + threadIdx.x;
                int j = blockIdx.y * blockDim.y + threadIdx.y;
                if (i >= nb_genes || j >= i) return;

                // =====================
                // Phase 1: FULL (Wfullij)
                // =====================
                float sum_w_full    = 0.0f;
                float sum_wx_full_i = 0.0f;
                float sum_wx_full_j = 0.0f;
                int k = 0;

                PROPR_UNROLL
                for (; k < (NT/4)*4; k += 4) {
                    float4 yfull_i4 = thread::load<Config::LoadModifer,float4>(&d_Yfull[k + i * Yfull_stride]);
                    float4 yfull_j4 = thread::load<Config::LoadModifer,float4>(&d_Yfull[k + j * Yfull_stride]);
                    float4 wfull_i4 = thread::load<Config::LoadModifer,float4>(&d_Wfull[k + i * Wfull_stride]);
                    float4 wfull_j4 = thread::load<Config::LoadModifer,float4>(&d_Wfull[k + j * Wfull_stride]);

                    PROPR_UNROLL
                    for (int m = 0; m < 4; ++m) {
                        float y_i = reinterpret_cast<float*>(&yfull_i4)[m];
                        float y_j = reinterpret_cast<float*>(&yfull_j4)[m];
                        float w_i = reinterpret_cast<float*>(&wfull_i4)[m];
                        float w_j = reinterpret_cast<float*>(&wfull_j4)[m];

                        // Wfullij_k = 2 * w_i * w_j / (w_i + w_j)
                        float denom = w_i + w_j;
                        float w_ij  = (denom > 0.0f) ? (2.0f * w_i * w_j / denom) : 0.0f;

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

                    float denom = w_i + w_j;
                    float w_ij  = (denom > 0.0f) ? (2.0f * w_i * w_j / denom) : 0.0f;

                    float X_i = __powf(y_i, a);
                    float X_j = __powf(y_j, a);

                    sum_w_full    += w_ij;
                    sum_wx_full_i += w_ij * X_i;
                    sum_wx_full_j += w_ij * X_j;
                }

                float mu_i_full = 0.0f, mu_j_full = 0.0f;
                float T_full = 0.0f;  // sum(Wfullij * (Xfull_i - Xfull_j))
                if (sum_w_full > 1e-10f) {
                    mu_i_full = sum_wx_full_i / sum_w_full;  // mean_Xfull_i
                    mu_j_full = sum_wx_full_j / sum_w_full;  // mean_Xfull_j
                    T_full    = sum_wx_full_i - sum_wx_full_j;
                }

                // =====================
                // Phase 2: CURRENT (Wij)
                // =====================
                float sum_w_current    = 0.0f;
                float sum_wx_current_i = 0.0f;
                float sum_wx_current_j = 0.0f;

                k = 0;
                PROPR_UNROLL
                for (; k < (N1/4)*4; k += 4) {
                    float4 y_i4 = thread::load<Config::LoadModifer,float4>(&d_Y[k + i * Y_stride]);
                    float4 y_j4 = thread::load<Config::LoadModifer,float4>(&d_Y[k + j * Y_stride]);
                    float4 w_i4 = thread::load<Config::LoadModifer,float4>(&d_W[k + i * W_stride]);
                    float4 w_j4 = thread::load<Config::LoadModifer,float4>(&d_W[k + j * W_stride]);

                    PROPR_UNROLL
                    for (int m = 0; m < 4; ++m) {
                        float y_i = reinterpret_cast<float*>(&y_i4)[m];
                        float y_j = reinterpret_cast<float*>(&y_j4)[m];
                        float w_i = reinterpret_cast<float*>(&w_i4)[m];
                        float w_j = reinterpret_cast<float*>(&w_j4)[m];

                        // Wij_k = 2 * w_i * w_j / (w_i + w_j)
                        float denom = w_i + w_j;
                        float w_ij  = (denom > 0.0f) ? (2.0f * w_i * w_j / denom) : 0.0f;

                        float X_i = __powf(y_i, a);
                        float X_j = __powf(y_j, a);

                        sum_w_current    += w_ij;
                        sum_wx_current_i += w_ij * X_i;
                        sum_wx_current_j += w_ij * X_j;
                    }
                }

                for (; k < N1; ++k) {
                    float y_i = d_Y[k + i * Y_stride];
                    float y_j = d_Y[k + j * Y_stride];
                    float w_i = d_W[k + i * W_stride];
                    float w_j = d_W[k + j * W_stride];

                    float denom = w_i + w_j;
                    float w_ij  = (denom > 0.0f) ? (2.0f * w_i * w_j / denom) : 0.0f;

                    float X_i = __powf(y_i, a);
                    float X_j = __powf(y_j, a);

                    sum_w_current    += w_ij;
                    sum_wx_current_i += w_ij * X_i;
                    sum_wx_current_j += w_ij * X_j;
                }

                // T_current = sum(Wij * (X_i - X_j))
                float T_current = sum_wx_current_i - sum_wx_current_j;

                // -------- C_z term --------
                float complement_term   = 0.0f;
                float denom_complement  = sum_w_full - sum_w_current;
                if (denom_complement > 1e-10f) {
                    complement_term = (T_full - T_current) / denom_complement;
                }

                float C_z = 0.0f;
                if (sum_w_current > 1e-10f) {
                    C_z = (T_current / sum_w_current) + complement_term;
                } else if (denom_complement > 1e-10f) {
                    // no "current" part, fall back to full
                    C_z = T_full / sum_w_full;
                }

                // -------- M_z term --------
                float M_z = 0.0f;
                if (sum_w_current > 1e-10f && mu_i_full > 1e-10f && mu_j_full > 1e-10f) {
                    // (sum(Wij * X_i)/mu_i_full - sum(Wij * X_j)/mu_j_full) / sum(Wij)
                    M_z = (sum_wx_current_i / mu_i_full - sum_wx_current_j / mu_j_full) / sum_w_current;
                }

                float result = ((C_z / 2.0f) + M_z) / a;

                int pair_index = (i * (i - 1)) / 2 + j;
                d_means[pair_index] = result;
            }

        }
    }
}