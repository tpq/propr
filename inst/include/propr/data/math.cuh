#pragma once

#include <algorithm>
#include <cmath>
#include <math.h>

#include <propr/utils/common/preprocessor.cuh>

namespace propr {
    namespace math {

        static PROPR_DEVICE PROPR_FORCE_INLINE
        bool is_nan(float value) {
            return isnan(value);
        }

        static PROPR_DEVICE PROPR_FORCE_INLINE
        bool is_nan(double value) {
            return isnan(value);
        }

        inline bool
        nearly_equal(double a, double b, double relative_tolerance = 1e-10) {
            const double scale = std::max(std::max(std::fabs(a), std::fabs(b)), 1.0);
            return std::fabs(a - b) <= relative_tolerance * scale;
        }

    } // namespace math
} // namespace propr
