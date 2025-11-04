#pragma once

#include <propr/utils/preprocessor.cuh>


namespace propr {
        template <class A, class B,
                class = std::enable_if_t<std::is_integral_v<A> && std::is_integral_v<B>>>
        PROPR_HOST_DEVICE
        PROPR_FORCE_INLINE
        constexpr std::make_unsigned_t<std::common_type_t<A,B>>
        ceil_div(A a, B b) {
            using R = std::make_unsigned_t<std::common_type_t<A,B>>;
            R ra = static_cast<R>(a);
            R rb = static_cast<R>(b);
            if (rb == 0) return 0;              // or: assert(false)
            if (ra == 0) return 0;              // avoids underflow in (ra-1)
            return 1 + ((ra - 1) / rb);         // overflow-safe vs (ra + rb - 1)
        }
}

