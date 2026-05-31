#pragma once

#include <algorithm>
#include <cfloat>
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

        template <typename T>
        PROPR_HOST_DEVICE PROPR_FORCE_INLINE
        T zero() {
            return T(0);
        }

        template <typename T>
        PROPR_HOST_DEVICE PROPR_FORCE_INLINE
        T one() {
            return T(1);
        }

        template <typename T>
        PROPR_HOST_DEVICE PROPR_FORCE_INLINE
        T eps();

        template <>
        PROPR_HOST_DEVICE PROPR_FORCE_INLINE
        float eps<float>() {
            return FLT_EPSILON;
        }

        template <>
        PROPR_HOST_DEVICE PROPR_FORCE_INLINE
        double eps<double>() {
            return DBL_EPSILON;
        }

        PROPR_HOST_DEVICE PROPR_FORCE_INLINE
        float log_t(float x) {
            return logf(x);
        }

        PROPR_HOST_DEVICE PROPR_FORCE_INLINE
        double log_t(double x) {
            return log(x);
        }

        PROPR_HOST_DEVICE PROPR_FORCE_INLINE
        float pow_t(float x, float a) {
            return powf(x, a);
        }

        PROPR_HOST_DEVICE PROPR_FORCE_INLINE
        double pow_t(double x, double a) {
            return pow(x, a);
        }

        PROPR_HOST_DEVICE PROPR_FORCE_INLINE
        float sqrt_t(float x) {
            return sqrtf(x);
        }

        PROPR_HOST_DEVICE PROPR_FORCE_INLINE
        double sqrt_t(double x) {
            return sqrt(x);
        }

        PROPR_DEVICE PROPR_FORCE_INLINE
        float rsqrt_t(float x) {
            return rsqrtf(x);
        }

        PROPR_DEVICE PROPR_FORCE_INLINE
        double rsqrt_t(double x) {
            return rsqrt(x);
        }

        PROPR_HOST_DEVICE PROPR_FORCE_INLINE
        float min_t(float a, float b) {
            return fminf(a, b);
        }

        PROPR_HOST_DEVICE PROPR_FORCE_INLINE
        double min_t(double a, double b) {
            return fmin(a, b);
        }

        PROPR_HOST_DEVICE PROPR_FORCE_INLINE
        float max_t(float a, float b) {
            return fmaxf(a, b);
        }

        PROPR_HOST_DEVICE PROPR_FORCE_INLINE
        double max_t(double a, double b) {
            return fmax(a, b);
        }

        template <typename T>
        PROPR_HOST_DEVICE PROPR_FORCE_INLINE
        T clamp_t(T x, T lo, T hi) {
            return max_t(lo, min_t(hi, x));
        }

        inline bool
        nearly_equal(double a, double b, double relative_tolerance = 1e-10) {
            const double scale = std::max(std::max(std::fabs(a), std::fabs(b)), 1.0);
            return std::fabs(a - b) <= relative_tolerance * scale;
        }

    } // namespace math
} // namespace propr
