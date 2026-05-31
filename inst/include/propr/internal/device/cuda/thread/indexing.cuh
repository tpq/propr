#pragma once

#include <cmath>
#include <cstddef>

#include <propr/utils/common/preprocessor.cuh>

namespace propr {
    namespace cuda {
        namespace internal {
            namespace thread {

                PROPR_HOST_DEVICE
                PROPR_FORCE_INLINE
                void lower_triangle_pair_from_index(std::size_t k, int& i, int& j) {
                    const double t = sqrt(1.0 + 8.0 * static_cast<double>(k));
                    i = static_cast<int>((1.0 + t) / 2.0);
                    j = static_cast<int>( k - (static_cast<std::size_t>(i) * static_cast<std::size_t>(i - 1)) / 2);
                }

                // TODO: add the upper tri-indexing

                // TODO: add warp_id using uniform registers

            } // namespace thread
        } // namespace internal
    } // namespace cuda
} // namespace propr
