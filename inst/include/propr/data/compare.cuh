#pragma once

#include <propr/utils/common/preprocessor.cuh>

namespace propr {
    namespace compare {

        template <typename T>
        struct less_than {
            PROPR_FORCE_INLINE
            PROPR_HOST_DEVICE
            bool operator()(T lhs, T rhs) const {
                return lhs < rhs;
            }
        };

        template <typename T>
        struct less_equal_than {
            PROPR_FORCE_INLINE
            PROPR_HOST_DEVICE
            bool operator()(T lhs, T rhs) const {
                return lhs <= rhs;
            }
        };

        template <typename T>
        struct greater_than {
            PROPR_FORCE_INLINE
            PROPR_HOST_DEVICE
            bool operator()(T lhs, T rhs) const {
                return lhs > rhs;
            }
        };

        template <typename T>
        struct greater_equal_than {
            PROPR_FORCE_INLINE
            PROPR_HOST_DEVICE
            bool operator()(T lhs, T rhs) const {
                return lhs >= rhs;
            }
        };

        template <typename T>
        struct equal {
            PROPR_FORCE_INLINE
            PROPR_HOST_DEVICE
            bool operator()(T lhs, T rhs) const {
                return lhs == rhs;
            }
        };

        template <typename T>
        struct not_equal {
            PROPR_FORCE_INLINE
            PROPR_HOST_DEVICE
            bool operator()(T lhs, T rhs) const {
                return lhs != rhs;
            }
        };

    } // namespace compare
} // namespace propr
