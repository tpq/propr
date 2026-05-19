#pragma once

#include <cuda_runtime.h>
#include <limits>
#include <propr/data/types.h>
#include <propr/utils/common/preprocessor.cuh>
#include <propr/internal/device/cuda/thread/mem_ops.cuh>


using namespace propr::cuda::internal;


namespace propr {
    namespace detail {
        namespace cuda {

            template <class Config>
            __global__
            void
            lrm_basic_phase_1(float* __restrict__ d_Y,
                               offset_t d_Y_stride,
                               float* __restrict__ d_mean_log,
                               int nb_samples,
                               int nb_genes) {
                const auto EPS = std::numeric_limits<float>::epsilon();
                const int g = blockIdx.x * blockDim.x + threadIdx.x;
                if (g >= nb_genes) return;

                const offset_t g_offset = static_cast<offset_t>(g) * d_Y_stride;

                float s0 = 0.0;
                float s1 = 0.0;
                float s2 = 0.0;
                float s3 = 0.0;
                int k = 0;

                PROPR_UNROLL
                for (; k < (nb_samples / 4) * 4; k += 4) {
                    const float4 y = thread::load<Config::LoadModifer, float4>(&d_Y[g_offset + k]);
                    s0 += logf(fmaxf(y.x, EPS));
                    s1 += logf(fmaxf(y.y, EPS));
                    s2 += logf(fmaxf(y.z, EPS));
                    s3 += logf(fmaxf(y.w, EPS));
                }

                float sum = (s0 + s1) + (s2 + s3);
                for (; k < nb_samples; ++k) {
                    const float y = thread::load<Config::LoadModifer, float>(&d_Y[g_offset + k]);
                    sum += logf(fmaxf(y, EPS));
                }

                const float mean_log = sum / static_cast<float>(nb_samples);
                thread::store<Config::StoreModifer, float>(&d_mean_log[g], mean_log);
            }

            template <class Config>
            __global__
            void
            lrm_basic_phase_2(float* __restrict__ d_mean_log,
                              float* __restrict__ d_mean,
                              int nb_genes) {
                using P2_Layout = typename Config::P2_Layout;
                static_assert(P2_Layout::BLK_X == P2_Layout::BLK_Y, "Tile size must be square");
                constexpr int TILE_G = P2_Layout::BLK_X;

                const int li = threadIdx.x;
                const int lj = threadIdx.y;

                const int gi = blockIdx.x * TILE_G + li;
                const int gj = blockIdx.y * TILE_G + lj;

                if (blockIdx.y > blockIdx.x) return;

                __shared__ float sh_i[TILE_G], sh_j[TILE_G];

                if (lj == 0) {
                    sh_i[li] = (gi < nb_genes)
                        ? thread::load<Config::LoadModifer, float>(&d_mean_log[gi])
                        : 0.0f;
                }

                if (li == 0) {
                    sh_j[lj] = (gj < nb_genes)
                        ? thread::load<Config::LoadModifer, float>(&d_mean_log[gj])
                        : 0.0f;
                }

                __syncthreads();

                if (gi < nb_genes && gj < nb_genes && gj < gi) {
                    const offset_t pair_index =
                        (static_cast<offset_t>(gi) * static_cast<offset_t>(gi - 1)) / 2 +
                        static_cast<offset_t>(gj);
                    thread::store<Config::StoreModifer, float>(&d_mean[pair_index], sh_i[li] - sh_j[lj]);
                }
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

                        accum.x += w;

                        float log_val = logf((&y_i.x)[m] / (&y_j.x)[m]);
                        float delta   = log_val - mean_old;
                        float w_ratio = w / accum.x;
                        accum.y       = fmaf(w_ratio, delta, mean_old);
                    }
                }

                for (; k < nb_samples; ++k) {
                    float y_ik = d_Y[k + i * d_Y_stride];
                    float y_jk = d_Y[k + j * d_Y_stride];

                    float w_ik = d_W[k + i * d_W_stride];
                    float w_jk = d_W[k + j * d_W_stride];

                    float denom = w_ik + w_jk;
                    float w_k   = (denom > 0.0f) ? (2.0f * w_ik * w_jk / denom) : 0.0f;

                    float log_val  = logf(y_ik / y_jk);
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
            lrm_alpha_phase_1(float* __restrict__ d_Y,
                              offset_t d_Y_stride,
                              float* __restrict__ d_Yfull,
                              offset_t d_Yfull_stride,
                              int N1,
                              int NT,
                              float a,
                              float* __restrict__ d_h,
                              int nb_genes) {
                const auto EPS = std::numeric_limits<float>::epsilon();
                const int g = blockIdx.x * blockDim.x + threadIdx.x;
                if (g >= nb_genes) return;

                const offset_t y_offset = static_cast<offset_t>(g) * d_Y_stride;
                const offset_t yfull_offset = static_cast<offset_t>(g) * d_Yfull_stride;

                float U = 0.0;
                float S = 0.0;
                int k = 0;

                PROPR_UNROLL
                for (; k + 3 < NT; k += 4) {
                    const float4 y = thread::load<Config::LoadModifer, float4>(&d_Yfull[yfull_offset + k]);

                    const float y0 = fmaxf(y.x, EPS);
                    const float y1 = fmaxf(y.y, EPS);
                    const float y2 = fmaxf(y.z, EPS);
                    const float y3 = fmaxf(y.w, EPS);

                    U += powf(y0, a);
                    U += powf(y1, a);
                    U += powf(y2, a);
                    U += powf(y3, a);
                }

                for (; k < NT; ++k) {
                    const float y = fmaxf(thread::load<Config::LoadModifer, float>(&d_Yfull[yfull_offset + k]), EPS);
                    U += powf(y, a);
                }

                k = 0;
                PROPR_UNROLL
                for (; k + 3 < N1; k += 4) {
                    const float4 y = thread::load<Config::LoadModifer, float4>(&d_Y[y_offset + k]);

                    const float y0 = fmaxf(y.x, EPS);
                    const float y1 = fmaxf(y.y, EPS);
                    const float y2 = fmaxf(y.z, EPS);
                    const float y3 = fmaxf(y.w, EPS);

                    S += powf(y0, a);
                    S += powf(y1, a);
                    S += powf(y2, a);
                    S += powf(y3, a);
                }

                for (; k < N1; ++k) {
                    const float y = fmaxf(thread::load<Config::LoadModifer, float>(&d_Y[y_offset + k]), EPS);
                    S += powf(y, a);
                }

                const float inv_N1 = 1.0f / N1;

                float A = S * inv_N1;
                if (N1 < NT) {
                    A += (U - S) / (NT - N1);
                }

                const float U_safe = (EPS > 0.0f) ? fmax(U, EPS) : U;
                const float B = (NT * S) / (N1 * U_safe);

                const float h = (0.5f * A + B) / a;
                thread::store<Config::StoreModifer, float>(&d_h[g], h);
            }

            template<class Config>
            __global__
            void
            lrm_alpha_phase_2(float* __restrict__ d_h,
                              float* __restrict__ d_means,
                              int nb_genes) {
                using P2_Layout = typename Config::P2_Layout;
                static_assert(P2_Layout::BLK_X == P2_Layout::BLK_Y, "Tile size must be square");
                constexpr int TILE_G = P2_Layout::BLK_X;

                const int li = threadIdx.x;
                const int lj = threadIdx.y;

                const int gi = blockIdx.x * TILE_G + li;
                const int gj = blockIdx.y * TILE_G + lj;

                if (blockIdx.y > blockIdx.x) return;

                __shared__ float sh_i[TILE_G];
                __shared__ float sh_j[TILE_G];

                if (lj == 0) {
                    sh_i[li] = (gi < nb_genes) ? thread::load<Config::LoadModifer, float>(&d_h[gi]) : 0.0f;
                }

                if (li == 0) {
                    sh_j[lj] = (gj < nb_genes)? thread::load<Config::LoadModifer, float>(&d_h[gj]) : 0.0f;
                }

                __syncthreads();

                if (gi < nb_genes && gj < nb_genes && gj < gi) {
                    const offset_t pair_index =
                        (static_cast<offset_t>(gi) * static_cast<offset_t>(gi - 1)) / 2 +
                        static_cast<offset_t>(gj);
                    thread::store<Config::StoreModifer, float>(&d_means[pair_index], sh_i[li] - sh_j[lj]);
                }
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

                        float X_i = powf(y_i, a);
                        float X_j = powf(y_j, a);

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

                    float X_i = powf(y_i, a);
                    float X_j = powf(y_j, a);

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

                        float X_i = powf(y_i, a);
                        float X_j = powf(y_j, a);

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

                    float X_i = powf(y_i, a);
                    float X_j = powf(y_j, a);

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
