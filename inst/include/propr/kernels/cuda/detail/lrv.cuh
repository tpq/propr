#pragma once

#include <cuda_runtime.h>
#include <propr/data/math.cuh>
#include <propr/data/traits.cuh>
#include <propr/data/types.h>
#include <propr/utils/common/cuda_helpers.cuh>
#include <propr/utils/common/preprocessor.cuh>
#include <propr/internal/device/cuda/thread/indexing.cuh>
#include <propr/internal/device/cuda/thread/mem_ops.cuh>


using namespace propr::cuda::internal;

namespace propr{
    namespace detail {
        namespace cuda {

            template <typename Real, class Config>
            __global__
            void lrv_basic(Real* __restrict__ d_Y, offset_t stride,
                           Real* __restrict__ d_variances,
                           int nb_samples, int nb_genes) {
                using Wide = cuda_wide_vector_t<Real>;
                constexpr int Lanes = cuda_wide_lanes_v<Real>;
                int i = blockIdx.x * blockDim.x + threadIdx.x;
                int j = blockIdx.y * blockDim.y + threadIdx.y;
                if (i >= nb_genes || j >= i) return;

                Real sum_log_ratios  = Real(0);
                Real sum_log_ratios2 = Real(0);

                int k = 0;
                PROPR_UNROLL
                for (; k < (nb_samples / Lanes) * Lanes; k += Lanes) {
                    Wide y_i = thread::load<Config::LoadModifer, Wide>(&d_Y[k + i * stride]);
                    Wide y_j = thread::load<Config::LoadModifer, Wide>(&d_Y[k + j * stride]);

                    PROPR_UNROLL
                    for (int m = 0; m < Lanes; m++) {
                        Real log_val = propr::math::log_t(lane_at(y_i, m) / lane_at(y_j, m));
                        sum_log_ratios  += log_val;
                        sum_log_ratios2 += log_val * log_val;
                    }
                }

                PROPR_UNROLL
                for (; k < nb_samples; ++k) {
                    const Real yi = d_Y[k + i * stride];
                    const Real yj = d_Y[k + j * stride];

                    const Real log_val = propr::math::log_t(yi / yj);
                    sum_log_ratios  += log_val;
                    sum_log_ratios2 += log_val * log_val;
                }

                Real inv_n    = Real(1) / static_cast<Real>(nb_samples);
                Real mean     = sum_log_ratios * inv_n;
                Real variance = (sum_log_ratios2 - static_cast<Real>(nb_samples) * mean * mean) / static_cast<Real>(nb_samples - 1);

                int pair_index = (i * (i - 1)) / 2 + j;
                d_variances[pair_index] = variance;
            }


            template <typename Real, class Config>
            __global__
            void lrv_weighted(Real* __restrict__ d_Y, offset_t Y_stride,
                            Real* __restrict__ d_W, offset_t W_stride,
                            Real* __restrict__ d_variances,
                            int nb_samples, int nb_genes) {
                using Wide = cuda_wide_vector_t<Real>;
                constexpr int Lanes = cuda_wide_lanes_v<Real>;
                int i = blockIdx.x * blockDim.x + threadIdx.x;
                int j = blockIdx.y * blockDim.y + threadIdx.y;
                if (i >= nb_genes || j >= i) return;

                Real sum_w = Real(0);
                Real sum_w2 = Real(0);
                Real mean = Real(0);
                Real M2 = Real(0);
                int k = 0;

                PROPR_UNROLL
                for (; k < (nb_samples / Lanes) * Lanes; k += Lanes) {
                    Wide y_i = thread::load<Config::LoadModifer, Wide>(&d_Y[k + i * Y_stride]);
                    Wide y_j = thread::load<Config::LoadModifer, Wide>(&d_Y[k + j * Y_stride]);
                    Wide w_i = thread::load<Config::LoadModifer, Wide>(&d_W[k + i * W_stride]);
                    Wide w_j = thread::load<Config::LoadModifer, Wide>(&d_W[k + j * W_stride]);

                    PROPR_UNROLL
                    for (int m = 0; m < Lanes; ++m) {
                        Real yi = lane_at(y_i, m);
                        Real yj = lane_at(y_j, m);
                        Real wi = lane_at(w_i, m);
                        Real wj = lane_at(w_j, m);

                        Real w_sum = wi + wj;
                        if (w_sum <= Real(0)) continue;

                        Real w = Real(2) * wi * wj / w_sum;

                        Real log_val = propr::math::log_t(yi / yj);

                        Real mean_old = mean;

                        sum_w += w;
                        sum_w2 += w * w;

                        Real delta = log_val - mean_old;
                        mean += (w / sum_w) * delta;

                        Real delta_new = log_val - mean;
                        M2 += w * delta * delta_new;
                    }
                }

                // Tail loop
                PROPR_UNROLL
                for (; k < nb_samples; ++k) {
                    Real y_ik = d_Y[k + i * Y_stride];
                    Real y_jk = d_Y[k + j * Y_stride];
                    Real w_ik = d_W[k + i * W_stride];
                    Real w_jk = d_W[k + j * W_stride];

                    Real w_sum = w_ik + w_jk;
                    if (w_sum <= Real(0)) continue;

                    Real w_k = Real(2) * w_ik * w_jk / w_sum;
                    Real log_val = propr::math::log_t(y_ik / y_jk);
                    Real mean_old = mean;

                    sum_w += w_k;
                    sum_w2 += w_k * w_k;

                    Real delta = log_val - mean_old;
                    mean += (w_k / sum_w) * delta;

                    Real delta_new = log_val - mean;
                    M2 += w_k * delta * delta_new;
                }

                Real lrv = Real(0);
                Real denom = sum_w * sum_w - sum_w2;
                if (sum_w > Real(0) && denom > Real(0)) {
                    lrv = M2 * sum_w / denom;
                }

                int pair_index = (i * (i - 1)) / 2 + j;
                d_variances[pair_index] = lrv;
            }


            template <typename Real, class Config>
            __global__
            void lrv_alpha(Real* __restrict__ d_Y    , offset_t Y_stride,
                           Real* __restrict__ d_Yfull, offset_t Yfull_stride,
                           Real a,
                           Real* __restrict__ d_variances,
                           int nb_samples,
                           int nb_samples_full,
                           int nb_genes) {

                    int i = blockIdx.x * blockDim.x + threadIdx.x;
                    int j = blockIdx.y * blockDim.y + threadIdx.y;
                    if (i >= nb_genes || j >= i) return;

                    Real sum_full_i = Real(0);
                    Real sum_full_j = Real(0);
                    int k = 0;

                    for (; k < nb_samples_full; ++k) {
                        Real Xf_i = propr::math::pow_t(d_Yfull[k + i * Yfull_stride], a);
                        Real Xf_j = propr::math::pow_t(d_Yfull[k + j * Yfull_stride], a);
                        sum_full_i += Xf_i;
                        sum_full_j += Xf_j;
                    }

                    Real mu_full_i = (nb_samples_full > 0) ? (sum_full_i / static_cast<Real>(nb_samples_full)) : Real(0);
                    Real mu_full_j = (nb_samples_full > 0) ? (sum_full_j / static_cast<Real>(nb_samples_full)) : Real(0);

                    Real a_i = (mu_full_i != Real(0)) ? (Real(1) / mu_full_i) : Real(0);
                    Real a_j = (mu_full_j != Real(0)) ? (Real(1) / mu_full_j) : Real(0);
                    Real ai_sq = a_i * a_i;
                    Real aj_sq = a_j * a_j;
                    Real aij   = a_i * a_j;

                    Real mu_i = Real(0);
                    Real mu_j = Real(0);
                    Real C = Real(0);
                    Real acc_x = Real(0);
                    Real acc_y = Real(0);
                    int n = 0;
                    k = 0;

                    for (; k < nb_samples; ++k) {
                        n++;
                        Real inv_n = Real(1) / static_cast<Real>(n);
                        Real X_i = propr::math::pow_t(d_Y[k + i * Y_stride], a);
                        Real X_j = propr::math::pow_t(d_Y[k + j * Y_stride], a);

                        Real prev_mu_i = mu_i;
                        Real dx_i = X_i - prev_mu_i;
                        mu_i += dx_i * inv_n;

                        Real prev_mu_j = mu_j;
                        Real dx_j = X_j - prev_mu_j;
                        mu_j += dx_j * inv_n;

                        Real dxj_muj = X_j - mu_j;
                        C += dx_i * dxj_muj;

                        acc_x += X_i * X_i;
                        acc_y += X_j * X_j;
                    }

                    Real n_mui_sq = static_cast<Real>(n) * mu_i * mu_i;
                    Real n_muj_sq = static_cast<Real>(n) * mu_j * mu_j;
                    Real sum_sq_i = acc_x - n_mui_sq;
                    Real sum_sq_j = acc_y - n_muj_sq;

                    Real term1 = sum_sq_i * ai_sq;
                    Real term2 = sum_sq_j * aj_sq;
                    Real term3 = Real(2) * aij * C;
                    Real S         = term1 + term2 - term3;
                    Real a_sq      = a * a;
                    Real denom     = (n > 1) ? (a_sq * static_cast<Real>(n - 1)) : Real(1);
                    Real lrv_value = S / denom;

                    int pair_index = (i * (i - 1)) / 2 + j;
                    d_variances[pair_index] = lrv_value;
            }


            template <typename Real, class Config>
            __global__
            void
            lrv_alpha_weighted(
                Real* __restrict__ d_Y    , offset_t Y_stride,
                Real* __restrict__ d_Yfull, offset_t Yfull_stride,
                Real* __restrict__ d_W    , offset_t W_stride,
                Real* __restrict__ d_Wfull, offset_t Wfull_stride,
                Real a,
                Real* __restrict__ d_variances,
                int nb_samples,
                int nb_samples_full,
                int nb_genes)
            {
                int i = blockIdx.x * blockDim.x + threadIdx.x;
                int j = blockIdx.y * blockDim.y + threadIdx.y;
                if (i >= nb_genes || j >= i) return;

                Real sum_w_full = Real(0);
                Real sum_w_full_X_full_i = Real(0);
                Real sum_w_full_X_full_j = Real(0);

                int k = 0;
                for (; k < nb_samples_full; ++k) {
                    Real X_full_i = propr::math::pow_t(d_Yfull[k + i * Yfull_stride], a);
                    Real X_full_j = propr::math::pow_t(d_Yfull[k + j * Yfull_stride], a);

                    Real wi_full = d_Wfull[k + i * Wfull_stride];
                    Real wj_full = d_Wfull[k + j * Wfull_stride];
                    Real wsum_full = wi_full + wj_full;

                    Real w_full = (wsum_full > Real(0)) ? (Real(2) * wi_full * wj_full / wsum_full) : Real(0);

                    sum_w_full += w_full;
                    sum_w_full_X_full_i += w_full * X_full_i;
                    sum_w_full_X_full_j += w_full * X_full_j;
                }

                const Real eps = propr::math::eps<Real>();
                Real inv_sum_w_full = (sum_w_full > eps) ? Real(1) / sum_w_full : Real(0);
                Real mu_full_i = sum_w_full_X_full_i * inv_sum_w_full;
                Real mu_full_j = sum_w_full_X_full_j * inv_sum_w_full;

                Real sum_w = Real(0);
                Real sum_w_sq = Real(0);

                Real sum_wX_i = Real(0);
                Real sum_wX_j = Real(0);

                Real sum_wX_i_sq = Real(0);
                Real sum_wX_j_sq = Real(0);

                Real sum_wX_iX_j = Real(0);

                k = 0;
                for (; k < nb_samples; ++k) {
                    Real X_i = propr::math::pow_t(d_Y[k + i * Y_stride], a);
                    Real X_j = propr::math::pow_t(d_Y[k + j * Y_stride], a);

                    Real wi = d_W[k + i * W_stride];
                    Real wj = d_W[k + j * W_stride];
                    Real wsum = wi + wj;

                    Real w = (wsum > Real(0)) ? (Real(2) * wi * wj / wsum) : Real(0);

                    Real X_i_sq = X_i * X_i;
                    Real X_j_sq = X_j * X_j;
                    Real X_iX_j = X_i * X_j;

                    sum_w    += w;
                    sum_w_sq += w * w;

                    sum_wX_i += w * X_i;
                    sum_wX_j += w * X_j;

                    sum_wX_i_sq += w * X_i_sq;
                    sum_wX_j_sq += w * X_j_sq;

                    sum_wX_iX_j += w * X_iX_j;
                }

                Real inv_sum_w = (sum_w > eps) ? Real(1) / sum_w : Real(0);

                Real sum_sq_i = sum_wX_i_sq - sum_wX_i * sum_wX_i * inv_sum_w;
                Real sum_sq_j = sum_wX_j_sq - sum_wX_j * sum_wX_j * inv_sum_w;
                Real C = sum_wX_iX_j - sum_wX_i * sum_wX_j * inv_sum_w;

                Real denom_term = sum_w - sum_w_sq * inv_sum_w;

                Real denom_pos = propr::math::max_t(denom_term, Real(0));

                if (mu_full_i <= eps || mu_full_j <= eps || denom_pos <= eps) {
                    d_variances[(i * (i - 1)) / 2 + j] = Real(0);
                    return;
                }

                Real inv_mu_full_i = Real(1) / mu_full_i;
                Real inv_mu_full_j = Real(1) / mu_full_j;

                Real inv_mu_full_i_sq = inv_mu_full_i * inv_mu_full_i;
                Real inv_mu_full_j_sq = inv_mu_full_j * inv_mu_full_j;
                Real inv_mu_full_ij   = inv_mu_full_i * inv_mu_full_j;

                Real term1 = sum_sq_i * inv_mu_full_i_sq;
                Real term2 = sum_sq_j * inv_mu_full_j_sq;
                Real term3 = Real(2) * inv_mu_full_ij * C;
                Real numerator = term1 + term2 - term3;

                Real a_sq = a * a;

                Real denom = a_sq * denom_pos;

                Real lrv_value = numerator / denom;

                int pair_index = (i * (i - 1)) / 2 + j;
                d_variances[pair_index] = lrv_value;
            }

        }
    }
}
