#pragma once

#include <propr/utils/common/preprocessor.cuh>
#include <propr/utils/common/constants.h>

namespace propr {
    template<typename T>
    struct ReduceSum {
        PROPR_NO_DISCARD
        PROPR_DEVICE 
        PROPR_FORCE_INLINE T operator()(T a, T b) const { return a + b; }
        static constexpr T identity = T(0);
    };
}

namespace propr {
    namespace cuda {
        namespace internal {
            namespace warp {
                template<typename ReduceOp, typename T>
                PROPR_NO_DISCARD
                PROPR_DEVICE 
                PROPR_FORCE_INLINE T warp_reduce(T v, unsigned mask, ReduceOp op = ReduceOp()) {
                    const int lane = static_cast<int>(threadIdx.x) & (PROPR_WARP_SIZE - 1);
                    PROPR_UNROLL
                    for (int offset = PROPR_WARP_SIZE / 2; offset > 0; offset >>= 1) {
                        const int src_lane = lane + offset;
                        const bool src_active = (src_lane < PROPR_WARP_SIZE) && ((mask & (1u << src_lane)) != 0u);
                        const T other = __shfl_down_sync(mask, v, offset);
                        if (src_active) v = op(v, other);
                    }
                    return v;
                }

                template<typename ReduceOp, typename T>
                PROPR_NO_DISCARD
                PROPR_DEVICE
                PROPR_FORCE_INLINE T warp_reduce(T v, ReduceOp op = ReduceOp()) {
                    return warp_reduce(v, __activemask(), op);
                }

            } // namespace warp
        } // namespace internal
    } // namespace cuda
} // namespace propr
