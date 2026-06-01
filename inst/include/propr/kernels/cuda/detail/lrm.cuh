#pragma once

#include <cuda_runtime.h>
#include <limits>
#include <propr/data/math.cuh>
#include <propr/data/traits.cuh>
#include <propr/data/types.h>
#include <propr/utils/common/cuda_helpers.cuh>
#include <propr/utils/common/preprocessor.cuh>
#include <propr/internal/device/cuda/thread/mem_ops.cuh>


using namespace propr::cuda::internal;


namespace propr {
    namespace detail {
        namespace cuda {

            template <typename Real, class Config>
            __global__
            void
            lrm_basic_phase_1(Real* __restrict__ d_Y,
                               offset_t d_Y_stride,
                               Real* __restrict__ d_mean_log,
                               int nb_samples,
                               int nb_genes) {
                using Wide = cuda_wide_vector_t<Real>;
                constexpr int Lanes = cuda_wide_lanes_v<Real>;
                const Real EPS = propr::math::eps<Real>();
                const int g = blockIdx.x * blockDim.x + threadIdx.x;
                if (g >= nb_genes) return;

                const offset_t g_offset = static_cast<offset_t>(g) * d_Y_stride;

                Real sum = Real(0);
                int k = 0;

                PROPR_UNROLL
                for (; k < (nb_samples / Lanes) * Lanes; k += Lanes) {
                    const Wide y = thread::load<Config::LoadModifer, Wide>(&d_Y[g_offset + k]);
                    PROPR_UNROLL
                    for (int lane = 0; lane < Lanes; ++lane) {
                        sum += propr::math::log_t(propr::math::max_t(lane_at(y, lane), EPS));
                    }
                }

                for (; k < nb_samples; ++k) {
                    const Real y = thread::load<Config::LoadModifer, Real>(&d_Y[g_offset + k]);
                    sum += propr::math::log_t(propr::math::max_t(y, EPS));
                }

                const Real mean_log = sum / static_cast<Real>(nb_samples);
                thread::store<Config::StoreModifer, Real>(&d_mean_log[g], mean_log);
            }

            template <typename Real, class Config>
            __global__
            void
            lrm_basic_phase_2(Real* __restrict__ d_mean_log,
                              Real* __restrict__ d_mean,
                              int nb_genes) {
                using P2_Layout = typename Config::P2_Layout;
                static_assert(P2_Layout::BLK_X == P2_Layout::BLK_Y, "Tile size must be square");
                constexpr int TILE_G = P2_Layout::BLK_X;

                const int li = threadIdx.x;
                const int lj = threadIdx.y;

                const int gi = blockIdx.x * TILE_G + li;
                const int gj = blockIdx.y * TILE_G + lj;

                if (blockIdx.y > blockIdx.x) return;

                __shared__ Real sh_i[TILE_G], sh_j[TILE_G];

                if (lj == 0) {
                    sh_i[li] = (gi < nb_genes)
                        ? thread::load<Config::LoadModifer, Real>(&d_mean_log[gi])
                        : Real(0);
                }

                if (li == 0) {
                    sh_j[lj] = (gj < nb_genes)
                        ? thread::load<Config::LoadModifer, Real>(&d_mean_log[gj])
                        : Real(0);
                }

                __syncthreads();

                if (gi < nb_genes && gj < nb_genes && gj < gi) {
                    const offset_t pair_index =
                        (static_cast<offset_t>(gi) * static_cast<offset_t>(gi - 1)) / 2 +
                        static_cast<offset_t>(gj);
                    thread::store<Config::StoreModifer, Real>(&d_mean[pair_index], sh_i[li] - sh_j[lj]);
                }
            }


            template<typename Real, class Config>
            __global__
            void
            lrm_weighted(Real* __restrict__ d_Y, offset_t d_Y_stride,
                         Real* __restrict__ d_W, offset_t d_W_stride,
                         Real* __restrict__ d_mean,
                         int nb_samples,
                         int nb_genes) {
                using Wide = cuda_wide_vector_t<Real>;
                constexpr int Lanes = cuda_wide_lanes_v<Real>;
                int i = blockIdx.x * blockDim.x + threadIdx.x;
                int j = blockIdx.y * blockDim.y + threadIdx.y;
                if (i >= nb_genes || j >= i) return;

                Real sum_w = Real(0);
                Real mean = Real(0);
                int k = 0;

                PROPR_UNROLL
                for (; k < (nb_samples / Lanes) * Lanes; k += Lanes) {
                    Wide y_i = thread::load<Config::LoadModifer, Wide>(&d_Y[k + i * d_Y_stride]);
                    Wide y_j = thread::load<Config::LoadModifer, Wide>(&d_Y[k + j * d_Y_stride]);

                    Wide w_i = thread::load<Config::LoadModifer, Wide>(&d_W[k + i * d_W_stride]);
                    Wide w_j = thread::load<Config::LoadModifer, Wide>(&d_W[k + j * d_W_stride]);

                    PROPR_UNROLL
                    for (int m = 0; m < Lanes; ++m) {
                        Real mean_old = mean;

                        Real w_im  = lane_at(w_i, m);
                        Real w_jm  = lane_at(w_j, m);
                        Real denom = w_im + w_jm;
                        Real w     = (denom > Real(0)) ? (Real(2) * w_im * w_jm / denom) : Real(0);

                        sum_w += w;

                        Real log_val = propr::math::log_t(lane_at(y_i, m) / lane_at(y_j, m));
                        Real delta   = log_val - mean_old;
                        Real w_ratio = (sum_w > Real(0)) ? w / sum_w : Real(0);
                        mean += w_ratio * delta;
                    }
                }

                for (; k < nb_samples; ++k) {
                    Real y_ik = d_Y[k + i * d_Y_stride];
                    Real y_jk = d_Y[k + j * d_Y_stride];

                    Real w_ik = d_W[k + i * d_W_stride];
                    Real w_jk = d_W[k + j * d_W_stride];

                    Real denom = w_ik + w_jk;
                    Real w_k   = (denom > Real(0)) ? (Real(2) * w_ik * w_jk / denom) : Real(0);

                    Real log_val  = propr::math::log_t(y_ik / y_jk);
                    Real mean_old = mean;

                    sum_w += w_k;

                    Real delta   = log_val - mean_old;
                    Real w_ratio = (sum_w > Real(0)) ? w_k / sum_w : Real(0);
                    mean += w_ratio * delta;
                }

                int pair_index = (i * (i - 1)) / 2 + j;
                d_mean[pair_index] = mean;
            }

            template<typename Real, class Config>
            __global__
            void
            lrm_alpha_phase_1(Real* __restrict__ d_Y,
                              offset_t d_Y_stride,
                              Real* __restrict__ d_Yfull,
                              offset_t d_Yfull_stride,
                              int N1,
                              int NT,
                              Real a,
                              Real* __restrict__ d_h,
                              int nb_genes) {
                using Wide = cuda_wide_vector_t<Real>;
                constexpr int Lanes = cuda_wide_lanes_v<Real>;
                const Real EPS = propr::math::eps<Real>();
                const int g = blockIdx.x * blockDim.x + threadIdx.x;
                if (g >= nb_genes) return;

                const offset_t y_offset = static_cast<offset_t>(g) * d_Y_stride;
                const offset_t yfull_offset = static_cast<offset_t>(g) * d_Yfull_stride;

                Real U = Real(0);
                Real S = Real(0);
                int k = 0;

                PROPR_UNROLL
                for (; k < (NT / Lanes) * Lanes; k += Lanes) {
                    const Wide y = thread::load<Config::LoadModifer, Wide>(&d_Yfull[yfull_offset + k]);
                    PROPR_UNROLL
                    for (int lane = 0; lane < Lanes; ++lane) {
                        U += propr::math::pow_t(propr::math::max_t(lane_at(y, lane), EPS), a);
                    }
                }

                for (; k < NT; ++k) {
                    const Real y = propr::math::max_t(thread::load<Config::LoadModifer, Real>(&d_Yfull[yfull_offset + k]), EPS);
                    U += propr::math::pow_t(y, a);
                }

                k = 0;
                PROPR_UNROLL
                for (; k < (N1 / Lanes) * Lanes; k += Lanes) {
                    const Wide y = thread::load<Config::LoadModifer, Wide>(&d_Y[y_offset + k]);
                    PROPR_UNROLL
                    for (int lane = 0; lane < Lanes; ++lane) {
                        S += propr::math::pow_t(propr::math::max_t(lane_at(y, lane), EPS), a);
                    }
                }

                for (; k < N1; ++k) {
                    const Real y = propr::math::max_t(thread::load<Config::LoadModifer, Real>(&d_Y[y_offset + k]), EPS);
                    S += propr::math::pow_t(y, a);
                }

                const Real inv_N1 = Real(1) / static_cast<Real>(N1);

                Real A = S * inv_N1;
                if (N1 < NT) {
                    A += (U - S) / (NT - N1);
                }

                const Real U_safe = propr::math::max_t(U, EPS);
                const Real B = (static_cast<Real>(NT) * S) / (static_cast<Real>(N1) * U_safe);

                const Real h = (Real(0.5) * A + B) / a;
                thread::store<Config::StoreModifer, Real>(&d_h[g], h);
            }

            template<typename Real, class Config>
            __global__
            void
            lrm_alpha_phase_2(Real* __restrict__ d_h,
                              Real* __restrict__ d_means,
                              int nb_genes) {
                using P2_Layout = typename Config::P2_Layout;
                static_assert(P2_Layout::BLK_X == P2_Layout::BLK_Y, "Tile size must be square");
                constexpr int TILE_G = P2_Layout::BLK_X;

                const int li = threadIdx.x;
                const int lj = threadIdx.y;

                const int gi = blockIdx.x * TILE_G + li;
                const int gj = blockIdx.y * TILE_G + lj;

                if (blockIdx.y > blockIdx.x) return;

                __shared__ Real sh_i[TILE_G];
                __shared__ Real sh_j[TILE_G];

                if (lj == 0) {
                    sh_i[li] = (gi < nb_genes) ? thread::load<Config::LoadModifer, Real>(&d_h[gi]) : Real(0);
                }

                if (li == 0) {
                    sh_j[lj] = (gj < nb_genes)? thread::load<Config::LoadModifer, Real>(&d_h[gj]) : Real(0);
                }

                __syncthreads();

                if (gi < nb_genes && gj < nb_genes && gj < gi) {
                    const offset_t pair_index =
                        (static_cast<offset_t>(gi) * static_cast<offset_t>(gi - 1)) / 2 +
                        static_cast<offset_t>(gj);
                    thread::store<Config::StoreModifer, Real>(&d_means[pair_index], sh_i[li] - sh_j[lj]);
                }
            }

            template<typename Real, class Config>
            __global__
            void
            lrm_alpha_weighted( Real* __restrict__ d_Y    , offset_t Y_stride,
                                Real* __restrict__ d_Yfull, offset_t Yfull_stride,
                                Real* __restrict__ d_W    , offset_t W_stride,
                                Real* __restrict__ d_Wfull, offset_t Wfull_stride,
                                int N1, int NT,
                                Real a,
                                Real* __restrict__ d_means,
                                int nb_genes)
            {
                using Wide = cuda_wide_vector_t<Real>;
                constexpr int Lanes = cuda_wide_lanes_v<Real>;
                int i = blockIdx.x * blockDim.x + threadIdx.x;
                int j = blockIdx.y * blockDim.y + threadIdx.y;
                if (i >= nb_genes || j >= i) return;

                // =====================
                // Phase 1: FULL (Wfullij)
                // =====================
                Real sum_w_full    = Real(0);
                Real sum_wx_full_i = Real(0);
                Real sum_wx_full_j = Real(0);
                int k = 0;

                PROPR_UNROLL
                for (; k < (NT / Lanes) * Lanes; k += Lanes) {
                    Wide yfull_i4 = thread::load<Config::LoadModifer, Wide>(&d_Yfull[k + i * Yfull_stride]);
                    Wide yfull_j4 = thread::load<Config::LoadModifer, Wide>(&d_Yfull[k + j * Yfull_stride]);
                    Wide wfull_i4 = thread::load<Config::LoadModifer, Wide>(&d_Wfull[k + i * Wfull_stride]);
                    Wide wfull_j4 = thread::load<Config::LoadModifer, Wide>(&d_Wfull[k + j * Wfull_stride]);

                    PROPR_UNROLL
                    for (int m = 0; m < Lanes; ++m) {
                        Real y_i = lane_at(yfull_i4, m);
                        Real y_j = lane_at(yfull_j4, m);
                        Real w_i = lane_at(wfull_i4, m);
                        Real w_j = lane_at(wfull_j4, m);

                        Real denom = w_i + w_j;
                        Real w_ij  = (denom > Real(0)) ? (Real(2) * w_i * w_j / denom) : Real(0);

                        Real X_i = propr::math::pow_t(y_i, a);
                        Real X_j = propr::math::pow_t(y_j, a);

                        sum_w_full    += w_ij;
                        sum_wx_full_i += w_ij * X_i;
                        sum_wx_full_j += w_ij * X_j;
                    }
                }

                for (; k < NT; ++k) {
                    Real y_i = d_Yfull[k + i * Yfull_stride];
                    Real y_j = d_Yfull[k + j * Yfull_stride];
                    Real w_i = d_Wfull[k + i * Wfull_stride];
                    Real w_j = d_Wfull[k + j * Wfull_stride];

                    Real denom = w_i + w_j;
                    Real w_ij  = (denom > Real(0)) ? (Real(2) * w_i * w_j / denom) : Real(0);

                    Real X_i = propr::math::pow_t(y_i, a);
                    Real X_j = propr::math::pow_t(y_j, a);

                    sum_w_full    += w_ij;
                    sum_wx_full_i += w_ij * X_i;
                    sum_wx_full_j += w_ij * X_j;
                }

                Real mu_i_full = Real(0), mu_j_full = Real(0);
                Real T_full = Real(0);
                const Real eps = propr::math::eps<Real>();
                if (sum_w_full > eps) {
                    mu_i_full = sum_wx_full_i / sum_w_full;  // mean_Xfull_i
                    mu_j_full = sum_wx_full_j / sum_w_full;  // mean_Xfull_j
                    T_full    = sum_wx_full_i - sum_wx_full_j;
                }

                // =====================
                // Phase 2: CURRENT (Wij)
                // =====================
                Real sum_w_current    = Real(0);
                Real sum_wx_current_i = Real(0);
                Real sum_wx_current_j = Real(0);

                k = 0;
                PROPR_UNROLL
                for (; k < (N1 / Lanes) * Lanes; k += Lanes) {
                    Wide y_i4 = thread::load<Config::LoadModifer, Wide>(&d_Y[k + i * Y_stride]);
                    Wide y_j4 = thread::load<Config::LoadModifer, Wide>(&d_Y[k + j * Y_stride]);
                    Wide w_i4 = thread::load<Config::LoadModifer, Wide>(&d_W[k + i * W_stride]);
                    Wide w_j4 = thread::load<Config::LoadModifer, Wide>(&d_W[k + j * W_stride]);

                    PROPR_UNROLL
                    for (int m = 0; m < Lanes; ++m) {
                        Real y_i = lane_at(y_i4, m);
                        Real y_j = lane_at(y_j4, m);
                        Real w_i = lane_at(w_i4, m);
                        Real w_j = lane_at(w_j4, m);

                        Real denom = w_i + w_j;
                        Real w_ij  = (denom > Real(0)) ? (Real(2) * w_i * w_j / denom) : Real(0);

                        Real X_i = propr::math::pow_t(y_i, a);
                        Real X_j = propr::math::pow_t(y_j, a);

                        sum_w_current    += w_ij;
                        sum_wx_current_i += w_ij * X_i;
                        sum_wx_current_j += w_ij * X_j;
                    }
                }

                for (; k < N1; ++k) {
                    Real y_i = d_Y[k + i * Y_stride];
                    Real y_j = d_Y[k + j * Y_stride];
                    Real w_i = d_W[k + i * W_stride];
                    Real w_j = d_W[k + j * W_stride];

                    Real denom = w_i + w_j;
                    Real w_ij  = (denom > Real(0)) ? (Real(2) * w_i * w_j / denom) : Real(0);

                    Real X_i = propr::math::pow_t(y_i, a);
                    Real X_j = propr::math::pow_t(y_j, a);

                    sum_w_current    += w_ij;
                    sum_wx_current_i += w_ij * X_i;
                    sum_wx_current_j += w_ij * X_j;
                }

                // T_current = sum(Wij * (X_i - X_j))
                Real T_current = sum_wx_current_i - sum_wx_current_j;

                // -------- C_z term --------
                Real complement_term   = Real(0);
                Real denom_complement  = sum_w_full - sum_w_current;
                if (denom_complement > eps) {
                    complement_term = (T_full - T_current) / denom_complement;
                }

                Real C_z = Real(0);
                if (sum_w_current > eps) {
                    C_z = (T_current / sum_w_current) + complement_term;
                } else if (denom_complement > eps) {
                    C_z = T_full / sum_w_full;
                }

                // -------- M_z term --------
                Real M_z = Real(0);
                if (sum_w_current > eps && mu_i_full > eps && mu_j_full > eps) {
                    // (sum(Wij * X_i)/mu_i_full - sum(Wij * X_j)/mu_j_full) / sum(Wij)
                    M_z = (sum_wx_current_i / mu_i_full - sum_wx_current_j / mu_j_full) / sum_w_current;
                }

                Real result = ((C_z / Real(2)) + M_z) / a;

                int pair_index = (i * (i - 1)) / 2 + j;
                d_means[pair_index] = result;
            }

        }
    }
}
