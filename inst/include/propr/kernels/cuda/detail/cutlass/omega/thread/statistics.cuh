#pragma once

#include <cute/tensor.hpp>

namespace propr {
    namespace kernels {
        namespace cutlass {
            namespace thread{

                template <class TA, class ALayout,
                          class TB, class BLayout,
                          class FrgTensorAccm,
                          class TiledMMA>
                CUTE_HOST_DEVICE
                void harmonic_mean(TiledMMA const& tiled_mma,
                                   cute::Tensor<TA, ALayout>& a_view,
                                   cute::Tensor<TB, BLayout>& b_view,
                                   FrgTensorAccm const& num_acc,
                                   FrgTensorAccm const& den_acc) {
                    CUTE_UNROLL
                    for (int i = 0; i < size(a_view); ++i) {
                        auto a = a_view(i);
                        auto b = b_view(i);
                        a_view(i) = (2.0f / (a + b)) * a;
                    }
                    cute::gemm(tiled_mma, num_acc, a_view, b_view, num_acc);
                    CUTE_UNROLL
                    for (int i = 0; i < size(a_view); ++i) a_view(i) = a_view(i) * a_view(i);
                    CUTE_UNROLL
                    for (int i = 0; i < size(b_view); ++i) b_view(i) = b_view(i) * b_view(i);
                    cute::gemm(tiled_mma, den_acc, a_view, b_view, den_acc);
                }

                void
                CUTE_HOST_DEVICE
            }
        }
    }
}
