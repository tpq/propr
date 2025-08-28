#pragma once

#include <cuda_runtime.h>
#include <propr/data/types.h>
#include <propr/utils/preprocessor.cuh>

namespace propr{
    namespace detail {
        namespace cuda {

            template<typename InT=float, typename OutT=float, typename AccT=float>
            __global__ void
            lrv_basic(InT* __restrict__ d_Y, offset_t stride,
                      OutT* __restrict__ d_variances,
                      int nb_samples, int nb_genes) {
                int i = blockIdx.x * blockDim.x + threadIdx.x;
                int j = blockIdx.y * blockDim.y + threadIdx.y;
                if (i >= nb_genes || j >= i) return;

                AccT sum_log_ratios  = 0.0f;
                AccT sum_log_ratios2 = 0.0f;

                int k = 0;
                PROPR_UNROLL
                for (; k < (nb_samples / 4) * 4; k += 4) {
                    float4 y_i = *reinterpret_cast<float4*>(&d_Y[k + i * stride]);
                    float4 y_j = *reinterpret_cast<float4*>(&d_Y[k + j * stride]);

                    PROPR_UNROLL
                    for (int m = 0; m < 4; m++) {
                        float ratio     = __fdividef((&y_i.x)[m], (&y_j.x)[m]);
                        float log_val   = __logf(ratio);
                        sum_log_ratios  = __fmaf_rn(1.0f, log_val,    sum_log_ratios);
                        sum_log_ratios2 = __fmaf_rn(log_val, log_val, sum_log_ratios2); 
                    }
                }
                
                PROPR_UNROLL
                for (; k < nb_samples; ++k) {
                    const float yi = d_Y[k + i * stride];
                    const float yj = d_Y[k + j * stride];

                    const float ratio   = __fdividef(yi, yj);
                    const float log_val = __logf(ratio);

                    sum_log_ratios  = __fmaf_rn(1.0f, log_val,    sum_log_ratios);
                    sum_log_ratios2 = __fmaf_rn(log_val, log_val, sum_log_ratios2);
                }

                float inv_n    = __frcp_rn(static_cast<float>(nb_samples));
                float mean     = sum_log_ratios * inv_n;
                float variance = (sum_log_ratios2 - __fmul_rn(nb_samples, __fmul_rn(mean, mean))) * __frcp_rn(static_cast<float>(nb_samples - 1));

                int pair_index = (i * (i - 1)) / 2 + j;
                d_variances[pair_index] = variance;
            }


            __global__
            void
            lrv_weighted(float* __restrict__ d_Y, offset_t Y_stride,
                         float* __restrict__ d_W, offset_t W_stride,
                         float* __restrict__ d_variances, 
                         int nb_samples, int nb_genes) {
                int i = blockIdx.x * blockDim.x + threadIdx.x;
                int j = blockIdx.y * blockDim.y + threadIdx.y;
                if (i >= nb_genes || j >= i) return;

                // accum.x = w_sum, accum.y = w_sum2, accum.z = mean, accum.w = numerator (S)
                float4 accum = make_float4(0.0f, 0.0f, 0.0f, 0.0f);
                int k = 0;

                PROPR_UNROLL
                for (; k < (nb_samples/4)*4 ; k += 4) {
                    float4 y_i = *reinterpret_cast<float4*>(&d_Y[k + i * Y_stride]);
                    float4 y_j = *reinterpret_cast<float4*>(&d_Y[k + j * Y_stride]);
                    float4 w_i = *reinterpret_cast<float4*>(&d_W[k + i * W_stride]);
                    float4 w_j = *reinterpret_cast<float4*>(&d_W[k + j * W_stride]);
                
                    PROPR_UNROLL
                    for (int m = 0; m < 4; ++m) {
                        float ratio = __fdividef((&y_i.x)[m], (&y_j.x)[m]);
                        float log_val = __logf(ratio);
                        float mean_old = accum.z;
                        float w = (&w_i.x)[m] * (&w_j.x)[m];
                
                        accum.x = __fmaf_rn(1.0f, w, accum.x);
                        accum.y = __fmaf_rn(w, w, accum.y);
                
                        float delta = log_val - mean_old;
                        float w_ratio = __fdividef(w, accum.x);
                        accum.z = __fmaf_rn(w_ratio, delta, mean_old);
                
                        float delta_new = log_val - accum.z;
                        accum.w = __fmaf_rn(w, delta * delta_new, accum.w);
                    }
                }

                PROPR_UNROLL
                for (; k < nb_samples; ++k) {
                    float y_ik = d_Y[k + i * Y_stride];
                    float y_jk = d_Y[k + j * Y_stride];
                    float w_ik = d_W[k + i * W_stride];
                    float w_jk = d_W[k + j * W_stride];
                    float w_k = w_ik * w_jk;

                    float ratio = __fdividef(y_ik, y_jk);
                    float log_val = __logf(ratio);
                    float mean_old = accum.z;

                    // Update w_sum and w_sum2
                    accum.x += w_k;
                    accum.y += w_k * w_k;

                    // Compute new mean
                    float delta = log_val - mean_old;
                    float w_ratio = w_k / accum.x;
                    accum.z += w_ratio * delta;

                    // Update numerator (S)
                    float delta_new = log_val - accum.z;
                    accum.w += w_k * delta * delta_new;
                }

                float S_total = accum.x;
                float Q_total = accum.y;
                float numerator = accum.w;

                float denominator = S_total - (Q_total / S_total);
                float lrv = numerator / denominator;

                int pair_index = (i * (i - 1)) / 2 + j;
                d_variances[pair_index] = lrv;
            }


            __global__
            void
            lrv_alpha(float* __restrict__ d_Y    , offset_t Y_stride,
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
                        float4 yfull_i = *reinterpret_cast<float4*>(&d_Yfull[k + i * Yfull_stride]);
                        float4 yfull_j = *reinterpret_cast<float4*>(&d_Yfull[k + j * Yfull_stride]);
                        for (int m=0;m<4;++m) {
                            float Xf_i = __powf(reinterpret_cast<float*>(&yfull_i)[m], a);
                            float Xf_j = __powf(reinterpret_cast<float*>(&yfull_j)[m], a);
                            sum_full.x += Xf_i;
                            sum_full.y += Xf_j;
                        }
                    }
                    
                    for (; k < nb_samples_full; ++k) {
                        float Xf_i = __powf(d_Yfull[k + i * Yfull_stride], a);
                        float Xf_j = __powf(d_Yfull[k + j * Yfull_stride], a);
                        sum_full.x += Xf_i;
                        sum_full.y += Xf_j;
                    }

                    float mu_full_i = (nb_samples_full > 0) ? (sum_full.x / static_cast<float>(nb_samples_full)) : 0.0f;
                    float mu_full_j = (nb_samples_full > 0) ? (sum_full.y / static_cast<float>(nb_samples_full)) : 0.0f;

                    float a_i = (mu_full_i != 0.0f) ? __frcp_rn(mu_full_i) : 0.0f;
                    float a_j = (mu_full_j != 0.0f) ? __frcp_rn(mu_full_j) : 0.0f;
                    float ai_sq = __fmul_rn(a_i, a_i);
                    float aj_sq = __fmul_rn(a_j, a_j);
                    float aij   = __fmul_rn(a_i, a_j);

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
                            float inv_n = __frcp_rn(static_cast<float>(n));
                            float X_i = __powf(reinterpret_cast<float*>(&y_i)[m], a);
                            float X_j = __powf(reinterpret_cast<float*>(&y_j)[m], a);

                            float prev_mu_i = mu_m.x;
                            float dx_i = __fsub_rn(X_i, prev_mu_i);
                            mu_m.x = __fmaf_rn(dx_i, inv_n, prev_mu_i);

                            float prev_mu_j = mu_m.y;
                            float dx_j = __fsub_rn(X_j, prev_mu_j);
                            mu_m.y = __fmaf_rn(dx_j, inv_n, prev_mu_j);

                            float dxj_muj = __fsub_rn(X_j, mu_m.y);
                            C = __fmaf_rn(dx_i, dxj_muj, C);

                            acc_x = __fmaf_rn(X_i, X_i, acc_x);
                            acc_y = __fmaf_rn(X_j, X_j, acc_y);
                        }
                    }

                    for (; k < nb_samples; ++k) {
                        n++;
                        float inv_n = __frcp_rn(static_cast<float>(n));
                        float X_i = __powf(d_Y[k + i * Y_stride], a);
                        float X_j = __powf(d_Y[k + j * Y_stride], a);

                        float prev_mu_i = mu_m.x;
                        float dx_i = __fsub_rn(X_i, prev_mu_i);
                        mu_m.x = __fmaf_rn(dx_i, inv_n, prev_mu_i);

                        float prev_mu_j = mu_m.y;
                        float dx_j = __fsub_rn(X_j, prev_mu_j);
                        mu_m.y = __fmaf_rn(dx_j, inv_n, prev_mu_j);

                        float dxj_muj = __fsub_rn(X_j, mu_m.y);
                        C = __fmaf_rn(dx_i, dxj_muj, C);

                        acc_x = __fmaf_rn(X_i, X_i, acc_x);
                        acc_y = __fmaf_rn(X_j, X_j, acc_y);
                    }

                    float n_mui_sq = __fmul_rn(n, __fmul_rn(mu_m.x, mu_m.x));
                    float n_muj_sq = __fmul_rn(n, __fmul_rn(mu_m.y, mu_m.y));
                    float sum_sq_i = __fsub_rn(acc_x, n_mui_sq);
                    float sum_sq_j = __fsub_rn(acc_y, n_muj_sq);

                    // combined numerator S = sum_sq_i/ mu_full_i^2 + sum_sq_j/mu_full_j^2 - 2*C/(mu_full_i*mu_full_j)
                    float term1 = __fmul_rn(sum_sq_i, ai_sq);
                    float term2 = __fmul_rn(sum_sq_j, aj_sq);
                    float term3 = __fmul_rn(2.0f, __fmul_rn(aij, C));
                    float S         = __fadd_rn(term1, __fsub_rn(term2, term3));
                    float a_sq      = __fmul_rn(a, a);
                    float denom     = (n > 1) ? __fmul_rn(a_sq, static_cast<float>(n - 1)) : 1.0f;
                    float lrv_value = __fdiv_rn(S, denom);

                    int pair_index = (i * (i - 1)) / 2 + j;
                    d_variances[pair_index] = lrv_value;
            }

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
                // 1) compute weighted sums on Yfull^a with weights Wfull (separate pass)
                // -------------------------
                float sum_w_full = 0.0f;
                float sum_w_full_X_full_i = 0.0f;
                float sum_w_full_X_full_j = 0.0f;

                int k = 0;
                PROPR_UNROLL
                for (; k + 4 <= nb_samples_full; k += 4) {
                    float4 yfull_i = *reinterpret_cast<float4*>(&d_Yfull[k + i * Yfull_stride]);
                    float4 yfull_j = *reinterpret_cast<float4*>(&d_Yfull[k + j * Yfull_stride]);
                    float4 wfull_i = *reinterpret_cast<float4*>(&d_Wfull[k + i * Wfull_stride]);
                    float4 wfull_j = *reinterpret_cast<float4*>(&d_Wfull[k + j * Wfull_stride]);

                    for (int m = 0; m < 4; ++m) {
                        float X_full_i = __powf(reinterpret_cast<float*>(&yfull_i)[m], a);
                        float X_full_j = __powf(reinterpret_cast<float*>(&yfull_j)[m], a);
                        float w_full = __fmul_rn(reinterpret_cast<float*>(&wfull_i)[m],
                                                reinterpret_cast<float*>(&wfull_j)[m]);

                        sum_w_full = __fadd_rn(sum_w_full, w_full);
                        sum_w_full_X_full_i = __fmaf_rn(w_full, X_full_i, sum_w_full_X_full_i);
                        sum_w_full_X_full_j = __fmaf_rn(w_full, X_full_j, sum_w_full_X_full_j);
                    }
                }
                
                for (; k < nb_samples_full; ++k) {
                    float X_full_i = __powf(d_Yfull[k + i * Yfull_stride], a);
                    float X_full_j = __powf(d_Yfull[k + j * Yfull_stride], a);
                    float w_full   = __fmul_rn(d_Wfull[k + i * Wfull_stride],
                                            d_Wfull[k + j * Wfull_stride]);

                    sum_w_full = __fadd_rn(sum_w_full, w_full);
                    sum_w_full_X_full_i = __fmaf_rn(w_full, X_full_i, sum_w_full_X_full_i);
                    sum_w_full_X_full_j = __fmaf_rn(w_full, X_full_j, sum_w_full_X_full_j);
                }

                // branchless compute mu_full_i, mu_full_j:
                // inv_sum_w_full = 1/(sum_w_full + FLT_MIN)
                float inv_sum_w_full = __frcp_rn(__fadd_rn(sum_w_full, FLT_MIN));
                float mu_full_i = __fmul_rn(sum_w_full_X_full_i, inv_sum_w_full);
                float mu_full_j = __fmul_rn(sum_w_full_X_full_j, inv_sum_w_full);

                // create multiplicative masks in branchless way:
                // mu_mask ~ mu/(mu + FLT_MIN) => ~1 if mu >> FLT_MIN, ~0 if mu==0
                float mu_mask_i = __fdiv_rn(mu_full_i, __fadd_rn(mu_full_i, FLT_MIN));
                float mu_mask_j = __fdiv_rn(mu_full_j, __fadd_rn(mu_full_j, FLT_MIN));

                // -------------------------
                // 2) pass over Y/W to compute weighted sums for X = Y^a
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
                    float4 y_i = *reinterpret_cast<float4*>(&d_Y[k + i * Y_stride]);
                    float4 y_j = *reinterpret_cast<float4*>(&d_Y[k + j * Y_stride]);
                    float4 w_i = *reinterpret_cast<float4*>(&d_W[k + i * W_stride]);
                    float4 w_j = *reinterpret_cast<float4*>(&d_W[k + j * W_stride]);

                    for (int m = 0; m < 4; ++m) {
                        float X_i = __powf(reinterpret_cast<float*>(&y_i)[m], a);
                        float X_j = __powf(reinterpret_cast<float*>(&y_j)[m], a);
                        float w   = __fmul_rn(reinterpret_cast<float*>(&w_i)[m],
                                            reinterpret_cast<float*>(&w_j)[m]);

                        float X_i_sq = __fmul_rn(X_i, X_i);
                        float X_j_sq = __fmul_rn(X_j, X_j);
                        float X_iX_j  = __fmul_rn(X_i, X_j);

                        sum_w    = __fadd_rn(sum_w, w);
                        sum_w_sq = __fmaf_rn(w, w, sum_w_sq);

                        sum_wX_i = __fmaf_rn(w, X_i, sum_wX_i);
                        sum_wX_j = __fmaf_rn(w, X_j, sum_wX_j);

                        sum_wX_i_sq  = __fmaf_rn(w, X_i_sq, sum_wX_i_sq);
                        sum_wX_j_sq  = __fmaf_rn(w, X_j_sq, sum_wX_j_sq);

                        sum_wX_iX_j = __fmaf_rn(w, X_iX_j, sum_wX_iX_j);
                    }
                }
                // tail
                for (; k < nb_samples; ++k) {
                    float X_i = __powf(d_Y[k + i * Y_stride], a);
                    float X_j = __powf(d_Y[k + j * Y_stride], a);
                    float w   = __fmul_rn(d_W[k + i * W_stride],
                                        d_W[k + j * W_stride]);

                    float X_i_sq = __fmul_rn(X_i, X_i);
                    float X_j_sq = __fmul_rn(X_j, X_j);
                    float X_iX_j  = __fmul_rn(X_i, X_j);

                    sum_w    = __fadd_rn(sum_w, w);
                    sum_w_sq = __fmaf_rn(w, w, sum_w_sq);

                    sum_wX_i = __fmaf_rn(w, X_i, sum_wX_i);
                    sum_wX_j = __fmaf_rn(w, X_j, sum_wX_j);

                    sum_wX_i_sq  = __fmaf_rn(w, X_i_sq, sum_wX_i_sq);
                    sum_wX_j_sq  = __fmaf_rn(w, X_j_sq, sum_wX_j_sq);

                    sum_wX_iX_j = __fmaf_rn(w, X_iX_j, sum_wX_iX_j);
                }

                // -------------------------
                // weighted central sums and denominator
                // -------------------------
                float inv_sum_w = __frcp_rn(__fadd_rn(sum_w, FLT_MIN));

                // central sums (numerators)
                float sum_sq_i = __fsub_rn(sum_wX_i_sq,
                                        __fmul_rn(__fmul_rn(sum_wX_i, sum_wX_i), inv_sum_w));
                float sum_sq_j = __fsub_rn(sum_wX_j_sq,
                                        __fmul_rn(__fmul_rn(sum_wX_j, sum_wX_j), inv_sum_w));
                float C = __fsub_rn(sum_wX_iX_j,
                                    __fmul_rn(__fmul_rn(sum_wX_i, sum_wX_j), inv_sum_w));

                // effective degrees-of-freedom term (weighted)
                float denom_term = __fsub_rn(sum_w, __fmul_rn(sum_w_sq, inv_sum_w));

                // positive-part: denom_pos = max(goo, 0)
                float denom_pos = fmaxf(denom_term, 0.0f);

                // denom_mask = denom_pos/(denom_pos + FLT_MIN) => ~1 if denom_pos >> FLT_MIN, ~0 if denom_pos == 0
                float denom_mask = __fdiv_rn(denom_pos, __fadd_rn(denom_pos, FLT_MIN));

                // mu masks computed earlier: mu_mask_i, mu_mask_j
                float valid_mask = __fmul_rn(__fmul_rn(mu_mask_i, mu_mask_j), denom_mask);

                // inverses for scaled numerator (guarded by FLT_MIN)
                float inv_mu_full_i = __frcp_rn(__fadd_rn(mu_full_i, FLT_MIN));
                float inv_mu_full_j = __frcp_rn(__fadd_rn(mu_full_j, FLT_MIN));

                float inv_mu_full_i_sq = __fmul_rn(inv_mu_full_i, inv_mu_full_i);
                float inv_mu_full_j_sq = __fmul_rn(inv_mu_full_j, inv_mu_full_j);
                float inv_mu_full_ij   = __fmul_rn(inv_mu_full_i, inv_mu_full_j);

                float term1 = __fmul_rn(sum_sq_i, inv_mu_full_i_sq);
                float term2 = __fmul_rn(sum_sq_j, inv_mu_full_j_sq);
                float term3 = __fmul_rn(2.0f, __fmul_rn(inv_mu_full_ij, C));
                float numerator = __fadd_rn(term1, __fsub_rn(term2, term3));

                float a_sq = __fmul_rn(a, a);

                // denom = a^2 * (denom_pos + FLT_MIN) => avoid /0; denom_pos + FLT_MIN >= FLT_MIN
                float denom = __fmul_rn(a_sq, __fadd_rn(denom_pos, FLT_MIN));

                // raw value (may be huge if denom FLT_MIN), but will be scaled down by valid_mask -> branchless zeroing
                float raw = __fdiv_rn(numerator, denom);

                // final branchless output
                float lrv_value = __fmul_rn(raw, valid_mask);

                // store
                int pair_index = (i * (i - 1)) / 2 + j;
                d_variances[pair_index] = lrv_value;
            }

        }
    }
}