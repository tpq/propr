#pragma once

#include <cuda_fp16.h>
#include <propr/utils/common/preprocessor.cuh>

namespace propr {
    namespace convert {

        template<typename OutT, typename InT>
        struct NumericConverter {
            PROPR_HOST_DEVICE
            PROPR_FORCE_INLINE
            OutT operator()(InT x) const { return static_cast<OutT>(x); }
        };

        template<>
        struct NumericConverter<__half, float> {
            PROPR_HOST_DEVICE
            PROPR_FORCE_INLINE
            __half operator()(const float x) const { return __float2half(x); }
        };

        template<>
        struct NumericConverter<float, __half> {
            PROPR_HOST_DEVICE
            PROPR_FORCE_INLINE
            float operator()(const __half h) const { return __half2float(h); }
        };

        template<>
        struct NumericConverter<__half, double> {
            PROPR_HOST_DEVICE
            PROPR_FORCE_INLINE
            __half operator()(const double x) const { return __float2half(static_cast<float>(x)); }
        };

        template<>
        struct NumericConverter<double, __half> {
            PROPR_HOST_DEVICE
            PROPR_FORCE_INLINE
            double operator()(const __half h) const { return static_cast<double>(__half2float(h)); }
        };

    } // namespace convert
} // namespace propr
