#pragma once

#include <cuda_runtime.h>
#include <propr/data/types.h>
#include <propr/utils/common/preprocessor.cuh>
#include <propr/internal/device/cuda/thread/mem_ops.cuh>


using namespace propr::cuda::internal;

namespace propr{
    namespace detail {
        namespace cuda {

            template <class Config>
            __global__
            void lrv_basic(float* __restrict__ d_Y, offset_t stride,
                           float* __restrict__ d_variances,
                           int nb_samples, int nb_genes) {
                int i = blockIdx.x * blockDim.x + threadIdx.x;
                int j = blockIdx.y * blockDim.y + threadIdx.y;
                if (i >= nb_genes || j >= i) return;

                float sum_log_ratios  = 0.0f;
                float sum_log_ratios2 = 0.0f;

                int k = 0;
                PROPR_UNROLL
                for (; k < (nb_samples / 4) * 4; k += 4) {
                    float4 y_i = thread::load<Config::LoadModifer,float4>(&d_Y[k + i * stride]);
                    float4 y_j = thread::load<Config::LoadModifer,float4>(&d_Y[k + j * stride]);

                    PROPR_UNROLL
                    for (int m = 0; m < 4; m++) {
                        float log_val   = logf((&y_i.x)[m] / (&y_j.x)[m]);
                        sum_log_ratios  += log_val;
                        sum_log_ratios2 += log_val * log_val;
                    }
                }

                PROPR_UNROLL
                for (; k < nb_samples; ++k) {
                    const float yi = d_Y[k + i * stride];
                    const float yj = d_Y[k + j * stride];

                    const float log_val = logf(yi / yj);
                    sum_log_ratios  += log_val;
                    sum_log_ratios2 += log_val * log_val;
                }

                float inv_n    = 1.0f / static_cast<float>(nb_samples);
                float mean     = sum_log_ratios * inv_n;
                float variance = (sum_log_ratios2 - nb_samples * mean * mean) / static_cast<float>(nb_samples - 1);

                int pair_index = (i * (i - 1)) / 2 + j;
                d_variances[pair_index] = variance;
            }


            template <class Config>
            __global__
            void lrv_weighted(float* __restrict__ d_Y, offset_t Y_stride,
                            float* __restrict__ d_W, offset_t W_stride,
                            float* __restrict__ d_variances,
                            int nb_samples, int nb_genes) {
                int i = blockIdx.x * blockDim.x + threadIdx.x;
                int j = blockIdx.y * blockDim.y + threadIdx.y;
                if (i >= nb_genes || j >= i) return;

                // accum.x = S_w      (sum of weights)
                // accum.y = S_w2     (sum of squared weights)
                // accum.z = mean     (weighted mean of log-ratio)
                // accum.w = M2       (sum w * (z - mean)^2)
                float4 accum = make_float4(0.0f, 0.0f, 0.0f, 0.0f);
                int k = 0;

                PROPR_UNROLL
                for (; k < (nb_samples / 4) * 4; k += 4) {
                    float4 y_i = thread::load<Config::LoadModifer,float4>(&d_Y[k + i * Y_stride]);
                    float4 y_j = thread::load<Config::LoadModifer,float4>(&d_Y[k + j * Y_stride]);
                    float4 w_i = thread::load<Config::LoadModifer,float4>(&d_W[k + i * W_stride]);
                    float4 w_j = thread::load<Config::LoadModifer,float4>(&d_W[k + j * W_stride]);

                    PROPR_UNROLL
                    for (int m = 0; m < 4; ++m) {
                        float yi = (&y_i.x)[m];
                        float yj = (&y_j.x)[m];
                        float wi = (&w_i.x)[m];
                        float wj = (&w_j.x)[m];

                        float w_sum = wi + wj;
                        if (w_sum <= 0.0f) continue;  // skip if degenerate

                        // Harmonic weight: 2 wi wj / (wi + wj)
                        float w = 2.0f * wi * wj / w_sum;

                        float log_val = logf(yi / yj);

                        float mean_old = accum.z;

                        // Update S_w and S_w2
                        accum.x += w;              // S_w += w
                        accum.y += w * w;          // S_w2 += w^2

                        // Update mean
                        float delta   = log_val - mean_old;
                        float w_ratio = w / accum.x;          // w / S_w
                        accum.z       = fmaf(w_ratio, delta, mean_old);

                        // Update M2
                        float delta_new = log_val - accum.z;
                        accum.w = fmaf(w, delta * delta_new, accum.w);
                    }
                }

                // Tail loop
                PROPR_UNROLL
                for (; k < nb_samples; ++k) {
                    float y_ik = d_Y[k + i * Y_stride];
                    float y_jk = d_Y[k + j * Y_stride];
                    float w_ik = d_W[k + i * W_stride];
                    float w_jk = d_W[k + j * W_stride];

                    float w_sum = w_ik + w_jk;
                    if (w_sum <= 0.0f) continue;

                    // Harmonic combined weight
                    float w_k = 2.0f * w_ik * w_jk / w_sum;
                    float log_val = logf(y_ik / y_jk);
                    float mean_old = accum.z;

                    // Update S_w and S_w2
                    accum.x += w_k;           // S_w
                    accum.y += w_k * w_k;     // S_w2

                    // Update mean
                    float delta   = log_val - mean_old;
                    float w_ratio = w_k / accum.x;   // w_k / S_w
                    accum.z      += w_ratio * delta;

                    // Update M2
                    float delta_new = log_val - accum.z;
                    accum.w += w_k * delta * delta_new;
                }

                float S_total   = accum.x;   // S_w
                float S_total2  = accum.y;   // S_w2
                float M2        = accum.w;

                float lrv = 0.0f;
                float denom   = S_total * S_total - S_total2;
                int   valid_i = (S_total > 0.0f) & (denom > 0.0f);
                float valid   = (float)valid_i;
                float safe_denom = valid * denom + (1.0f - valid) * (1.0f + fabsf(denom));

                float factor = S_total / safe_denom;
                lrv = M2 * factor * valid;

                int pair_index = (i * (i - 1)) / 2 + j;
                d_variances[pair_index] = lrv;
            }


            template <class Config>
            __global__
            void lrv_alpha(float* __restrict__ d_Y    , offset_t Y_stride,
                           float* __restrict__ d_Yfull, offset_t Yfull_stride,
                           float a,
                           float* __restrict__ d_variances,
                           int nb_samples,
                           int nb_samples_full,
                           int nb_genes) {

                    int i = blockIdx.x * blockDim.x + threadIdx.x;
                    int j = blockIdx.y * blockDim.y + threadIdx.y;
                    if (i >= nb_genes || j >= i) return;

                    float2 sum_full = make_float2(0.0f,0.0f); //compute mean of Yfull^a for columns i and j ---
                    int k = 0;

                    PROPR_UNROLL
                    for (; k + 4 <= nb_samples_full; k += 4) {
                        float4 yfull_i = thread::load<Config::LoadModifer,float4>(&d_Yfull[k + i * Yfull_stride]);
                        float4 yfull_j = thread::load<Config::LoadModifer,float4>(&d_Yfull[k + j * Yfull_stride]);
                        PROPR_UNROLL
                        for (int m=0;m<4;++m) {
                            float Xf_i = powf(reinterpret_cast<float*>(&yfull_i)[m], a);
                            float Xf_j = powf(reinterpret_cast<float*>(&yfull_j)[m], a);
                            sum_full.x += Xf_i;
                            sum_full.y += Xf_j;
                        }
                    }

                    for (; k < nb_samples_full; ++k) {
                        float Xf_i = powf(d_Yfull[k + i * Yfull_stride], a);
                        float Xf_j = powf(d_Yfull[k + j * Yfull_stride], a);
                        sum_full.x += Xf_i;
                        sum_full.y += Xf_j;
                    }

                    float mu_full_i = (nb_samples_full > 0) ? (sum_full.x / static_cast<float>(nb_samples_full)) : 0.0f;
                    float mu_full_j = (nb_samples_full > 0) ? (sum_full.y / static_cast<float>(nb_samples_full)) : 0.0f;

                    float a_i = (mu_full_i != 0.0f) ? (1.0f / mu_full_i) : 0.0f;
                    float a_j = (mu_full_j != 0.0f) ? (1.0f / mu_full_j) : 0.0f;
                    float ai_sq = a_i * a_i;
                    float aj_sq = a_j * a_j;
                    float aij   = a_i * a_j;

                    float4 mu_m = make_float4(0.0f, 0.0f, 0.0f, 0.0f); // mu_m.x = mean Xi, mu_m.y = mean Xj
                    float C = 0.0f;
                    float acc_x = 0.0f; // sum of Xi^2
                    float acc_y = 0.0f; // sum of Xj^2
                    int n = 0;
                    k = 0;
                    PROPR_UNROLL
                    for (; k + 4 <= nb_samples; k += 4) {
                        float4 y_i = *reinterpret_cast<float4*>(&d_Y[k + i * Y_stride]);
                        float4 y_j = *reinterpret_cast<float4*>(&d_Y[k + j * Y_stride]);
                        for (int m=0;m<4;++m) {
                            n++;
                            float inv_n = 1.0f / static_cast<float>(n);
                            float X_i = powf(reinterpret_cast<float*>(&y_i)[m], a);
                            float X_j = powf(reinterpret_cast<float*>(&y_j)[m], a);

                            float prev_mu_i = mu_m.x;
                            float dx_i = X_i - prev_mu_i;
                            mu_m.x = fmaf(dx_i, inv_n, prev_mu_i);

                            float prev_mu_j = mu_m.y;
                            float dx_j = X_j - prev_mu_j;
                            mu_m.y = fmaf(dx_j, inv_n, prev_mu_j);

                            float dxj_muj = X_j - mu_m.y;
                            C = fmaf(dx_i, dxj_muj, C);

                            acc_x = fmaf(X_i, X_i, acc_x);
                            acc_y = fmaf(X_j, X_j, acc_y);
                        }
                    }

                    for (; k < nb_samples; ++k) {
                        n++;
                        float inv_n = 1.0f / static_cast<float>(n);
                        float X_i = powf(d_Y[k + i * Y_stride], a);
                        float X_j = powf(d_Y[k + j * Y_stride], a);

                        float prev_mu_i = mu_m.x;
                        float dx_i = X_i - prev_mu_i;
                        mu_m.x = fmaf(dx_i, inv_n, prev_mu_i);

                        float prev_mu_j = mu_m.y;
                        float dx_j = X_j - prev_mu_j;
                        mu_m.y = fmaf(dx_j, inv_n, prev_mu_j);

                        float dxj_muj = X_j - mu_m.y;
                        C = fmaf(dx_i, dxj_muj, C);

                        acc_x = fmaf(X_i, X_i, acc_x);
                        acc_y = fmaf(X_j, X_j, acc_y);
                    }

                    float n_mui_sq = n * mu_m.x * mu_m.x;
                    float n_muj_sq = n * mu_m.y * mu_m.y;
                    float sum_sq_i = acc_x - n_mui_sq;
                    float sum_sq_j = acc_y - n_muj_sq;

                    // combined numerator S = sum_sq_i/ mu_full_i^2 + sum_sq_j/mu_full_j^2 - 2*C/(mu_full_i*mu_full_j)
                    float term1 = sum_sq_i * ai_sq;
                    float term2 = sum_sq_j * aj_sq;
                    float term3 = 2.0f * aij * C;
                    float S         = term1 + term2 - term3;
                    float a_sq      = a * a;
                    float denom     = (n > 1) ? (a_sq * static_cast<float>(n - 1)) : 1.0f;
                    float lrv_value = S / denom;

                    int pair_index = (i * (i - 1)) / 2 + j;
                    d_variances[pair_index] = lrv_value;
            }


            template <class Config>
            __global__
            void
            lrv_alpha_weighted(
                float* __restrict__ d_Y    , offset_t Y_stride,
                float* __restrict__ d_Yfull, offset_t Yfull_stride,
                float* __restrict__ d_W    , offset_t W_stride,
                float* __restrict__ d_Wfull, offset_t Wfull_stride,
                float a,
                float* __restrict__ d_variances,
                int nb_samples,
                int nb_samples_full,
                int nb_genes)
            {
                int i = blockIdx.x * blockDim.x + threadIdx.x;
                int j = blockIdx.y * blockDim.y + threadIdx.y;
                if (i >= nb_genes || j >= i) return;

                // -------------------------
                // 1) compute weighted sums on Yfull^a with harmonic weights Wfull_ij
                // -------------------------
                float sum_w_full = 0.0f;
                float sum_w_full_X_full_i = 0.0f;
                float sum_w_full_X_full_j = 0.0f;

                int k = 0;
                PROPR_UNROLL
                for (; k + 4 <= nb_samples_full; k += 4) {
                    float4 yfull_i = thread::load<Config::LoadModifer,float4>(&d_Yfull[k + i * Yfull_stride]);
                    float4 yfull_j = thread::load<Config::LoadModifer,float4>(&d_Yfull[k + j * Yfull_stride]);
                    float4 wfull_i = thread::load<Config::LoadModifer,float4>(&d_Wfull[k + i * Wfull_stride]);
                    float4 wfull_j = thread::load<Config::LoadModifer,float4>(&d_Wfull[k + j * Wfull_stride]);

                    PROPR_UNROLL
                    for (int m = 0; m < 4; ++m) {
                        float X_full_i = powf(reinterpret_cast<float*>(&yfull_i)[m], a);
                        float X_full_j = powf(reinterpret_cast<float*>(&yfull_j)[m], a);

                        float wi_full = reinterpret_cast<float*>(&wfull_i)[m];
                        float wj_full = reinterpret_cast<float*>(&wfull_j)[m];
                        float wsum_full = wi_full + wj_full;

                        // harmonic combined weight: 2 * wi * wj / (wi + wj)
                        float w_full = (2.0f * wi_full * wj_full) /
                            (wsum_full + FLT_MIN);  // avoid /0

                        sum_w_full += w_full;
                        sum_w_full_X_full_i = fmaf(w_full, X_full_i, sum_w_full_X_full_i);
                        sum_w_full_X_full_j = fmaf(w_full, X_full_j, sum_w_full_X_full_j);
                    }
                }

                for (; k < nb_samples_full; ++k) {
                    float X_full_i = powf(d_Yfull[k + i * Yfull_stride], a);
                    float X_full_j = powf(d_Yfull[k + j * Yfull_stride], a);

                    float wi_full = d_Wfull[k + i * Wfull_stride];
                    float wj_full = d_Wfull[k + j * Wfull_stride];
                    float wsum_full = wi_full + wj_full;

                    float w_full = (2.0f * wi_full * wj_full) /
                        (wsum_full + FLT_MIN);

                    sum_w_full += w_full;
                    sum_w_full_X_full_i = fmaf(w_full, X_full_i, sum_w_full_X_full_i);
                    sum_w_full_X_full_j = fmaf(w_full, X_full_j, sum_w_full_X_full_j);
                }

                // inv_sum_w_full = 1/(sum_w_full + FLT_MIN)
                float inv_sum_w_full = 1.0f / (sum_w_full + FLT_MIN);
                float mu_full_i = sum_w_full_X_full_i * inv_sum_w_full;
                float mu_full_j = sum_w_full_X_full_j * inv_sum_w_full;

                // mu_mask ~ mu/(mu + FLT_MIN) => ~1 if mu >> FLT_MIN, ~0 if mu==0
                float mu_mask_i = mu_full_i / (mu_full_i + FLT_MIN);
                float mu_mask_j = mu_full_j / (mu_full_j + FLT_MIN);

                // -------------------------
                // 2) pass over Y/W to compute weighted sums for X = Y^a (within-group)
                // -------------------------
                float sum_w = 0.0f;
                float sum_w_sq = 0.0f;

                float sum_wX_i = 0.0f;
                float sum_wX_j = 0.0f;

                float sum_wX_i_sq = 0.0f;
                float sum_wX_j_sq = 0.0f;

                float sum_wX_iX_j = 0.0f;

                k = 0;
                for (; k + 4 <= nb_samples; k += 4) {
                    float4 y_i = thread::load<Config::LoadModifer,float4>(&d_Y[k + i * Y_stride]);
                    float4 y_j = thread::load<Config::LoadModifer,float4>(&d_Y[k + j * Y_stride]);
                    float4 w_i = thread::load<Config::LoadModifer,float4>(&d_W[k + i * W_stride]);
                    float4 w_j = thread::load<Config::LoadModifer,float4>(&d_W[k + j * W_stride]);

                    for (int m = 0; m < 4; ++m) {
                        float X_i = powf(reinterpret_cast<float*>(&y_i)[m], a);
                        float X_j = powf(reinterpret_cast<float*>(&y_j)[m], a);

                        float wi = reinterpret_cast<float*>(&w_i)[m];
                        float wj = reinterpret_cast<float*>(&w_j)[m];
                        float wsum = wi + wj;

                        // harmonic combined weight: 2 * wi * wj / (wi + wj)
                        float w = (2.0f * wi * wj) /
                            (wsum + FLT_MIN);

                        float X_i_sq = X_i * X_i;
                        float X_j_sq = X_j * X_j;
                        float X_iX_j = X_i * X_j;

                        sum_w    += w;
                        sum_w_sq = fmaf(w, w, sum_w_sq);

                        sum_wX_i = fmaf(w, X_i, sum_wX_i);
                        sum_wX_j = fmaf(w, X_j, sum_wX_j);

                        sum_wX_i_sq  = fmaf(w, X_i_sq, sum_wX_i_sq);
                        sum_wX_j_sq  = fmaf(w, X_j_sq, sum_wX_j_sq);

                        sum_wX_iX_j = fmaf(w, X_iX_j, sum_wX_iX_j);
                    }
                }
                // tail
                for (; k < nb_samples; ++k) {
                    float X_i = powf(d_Y[k + i * Y_stride], a);
                    float X_j = powf(d_Y[k + j * Y_stride], a);

                    float wi = d_W[k + i * W_stride];
                    float wj = d_W[k + j * W_stride];
                    float wsum = wi + wj;

                    float w = (2.0f * wi * wj) /
                        (wsum + FLT_MIN);

                    float X_i_sq = X_i * X_i;
                    float X_j_sq = X_j * X_j;
                    float X_iX_j = X_i * X_j;

                    sum_w    += w;
                    sum_w_sq = fmaf(w, w, sum_w_sq);

                    sum_wX_i = fmaf(w, X_i, sum_wX_i);
                    sum_wX_j = fmaf(w, X_j, sum_wX_j);

                    sum_wX_i_sq  = fmaf(w, X_i_sq, sum_wX_i_sq);
                    sum_wX_j_sq  = fmaf(w, X_j_sq, sum_wX_j_sq);

                    sum_wX_iX_j = fmaf(w, X_iX_j, sum_wX_iX_j);
                }

                // -------------------------
                // weighted central sums and denominator
                // -------------------------
                float inv_sum_w = 1.0f / (sum_w + FLT_MIN);

                // central sums
                float sum_sq_i = sum_wX_i_sq - sum_wX_i * sum_wX_i * inv_sum_w;
                float sum_sq_j = sum_wX_j_sq - sum_wX_j * sum_wX_j * inv_sum_w;
                float C = sum_wX_iX_j - sum_wX_i * sum_wX_j * inv_sum_w;

                // effective degrees-of-freedom term (weighted)
                float denom_term = sum_w - sum_w_sq * inv_sum_w;

                // positive-part
                float denom_pos = fmaxf(denom_term, 0.0f);
                float denom_mask = denom_pos / (denom_pos + FLT_MIN);

                // mu masks: mu_mask_i, mu_mask_j from above
                float valid_mask = mu_mask_i * mu_mask_j * denom_mask;

                // inverses for scaled numerator
                float inv_mu_full_i = 1.0f / (mu_full_i + FLT_MIN);
                float inv_mu_full_j = 1.0f / (mu_full_j + FLT_MIN);

                float inv_mu_full_i_sq = inv_mu_full_i * inv_mu_full_i;
                float inv_mu_full_j_sq = inv_mu_full_j * inv_mu_full_j;
                float inv_mu_full_ij   = inv_mu_full_i * inv_mu_full_j;

                float term1 = sum_sq_i * inv_mu_full_i_sq;
                float term2 = sum_sq_j * inv_mu_full_j_sq;
                float term3 = 2.0f * inv_mu_full_ij * C;
                float numerator = term1 + term2 - term3;

                float a_sq = a * a;

                // denom = a^2 * (denom_pos + FLT_MIN)
                float denom = a_sq * (denom_pos + FLT_MIN);

                float raw = numerator / denom;
                float lrv_value = raw * valid_mask;

                int pair_index = (i * (i - 1)) / 2 + j;
                d_variances[pair_index] = lrv_value;
            }

        }
    }
}
