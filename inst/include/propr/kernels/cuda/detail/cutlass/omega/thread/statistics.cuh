#pragma once

#include <cute/tensor.hpp>

#include <propr/utils/preprocessor.cuh>

namespace propr {
    namespace kernels {
        namespace cutlass {
            namespace thread{

                template <typename OmegaConfig>
                struct OmegaHarmonic {
                    typename OmegaConfig::TiledMMA tiled_mma;
                    using DLayout = decltype(cute::make_layout(cute::Int<OmegaConfig::BLK_M>{}, cute::Int<OmegaConfig::BLK_M>{}));
                    using TD      =  cute::ViewEngine<cute::gmem_ptr<cute::half_t *>>;

                    cute::Tensor<TD, DLayout> &D;

                    decltype(tiled_mma.get_thread_slice(0u)) thread_mma;
                    decltype(thread_mma.partition_fragment_C(D)) tDrD_n;
                    decltype(thread_mma.partition_fragment_C(D)) tDrD_s;

                    PROPR_DEVICE
                    OmegaHarmonic(cute::Tensor<TD, DLayout>& D_)
                       : D(D_),
                         thread_mma(tiled_mma.get_thread_slice(threadIdx.x)),
                         tDrD_n(thread_mma.partition_fragment_C(D)),
                         tDrD_s(thread_mma.partition_fragment_C(D)) {
                        cute::clear(tDrD_n); cute::clear(tDrD_s);
                    }

                    PROPR_DEVICE
                    template <class TA, class ALayout, class TB, class BLayout>
                    void operator()(cute::Tensor<TA, ALayout>& a_view, cute::Tensor<TB, BLayout>& b_view) {
                        CUTE_UNROLL
                        for (int i = 0; i < size(a_view); ++i) {
                            auto a = a_view(i); auto b = b_view(i);
                            a_view(i) = (2.0f / (a + b)) * a;
                        }
                        cute::gemm(tiled_mma, tDrD_n, a_view, b_view, tDrD_n);
                        CUTE_UNROLL
                        for (int i = 0; i < size(a_view); ++i) a_view(i) = a_view(i) * a_view(i);
                        CUTE_UNROLL
                        for (int i = 0; i < size(b_view); ++i) b_view(i) = b_view(i) * b_view(i);
                        cute::gemm(tiled_mma, tDrD_s, a_view, b_view, tDrD_s);
                    }

                    PROPR_DEVICE
                    float get(int i) {
                        return  (1.0f - (tDrD_n(i) == 0.0f)) * (tDrD_n(i) - tDrD_s(i) / (tDrD_n(i) + (tDrD_n(i) == 0.0f)));
                    }

                    PROPR_DEVICE
                    auto size() { return cute::size(tDrD_n); }
                };



                template <typename OmegaConfig>
                struct OmegaPopulation {
                    using DLayout = decltype(cute::make_layout(cute::Int<OmegaConfig::BLK_M>{}, cute::Int<OmegaConfig::BLK_M>{}));
                    using TD      =  cute::ViewEngine<cute::gmem_ptr<cute::half_t *>>;
                    typename OmegaConfig::TiledMMA tiled_mma;

                    cute::Tensor<TD, DLayout> &D;

                    decltype(tiled_mma.get_thread_slice(0u)) thread_mma;
                    decltype(thread_mma.partition_fragment_C(D)) tDrD_n;

                    PROPR_DEVICE
                    OmegaPopulation(cute::Tensor<TD, DLayout>& D_) : D(D_),
                         thread_mma(tiled_mma.get_thread_slice(threadIdx.x)),
                         tDrD_n(thread_mma.partition_fragment_C(D)) {
                        cute::clear(tDrD_n);
                    }

                    PROPR_DEVICE
                    template <class TA, class ALayout, class TB, class BLayout>
                    void operator()(cute::Tensor<TA, ALayout>& a_view, cute::Tensor<TB, BLayout>& b_view) {
                        CUTE_UNROLL
                        for (int i = 0; i < size(a_view); ++i) {
                            auto a = a_view(i); auto b = b_view(i);
                            a_view(i) = (2.0f / (a + b)) * a;
                        }
                        cute::gemm(tiled_mma, tDrD_n, a_view, b_view, tDrD_n);
                    }

                    PROPR_DEVICE
                    float get(int i) { return tDrD_n(i); }

                    PROPR_DEVICE
                    auto size() { return cute::size(tDrD_n); }
                };

            }
        }
    }
}
