#pragma once

#include <cstdint>

namespace propr {

    template <typename To, typename From>
    static __device__ __forceinline__ To bit_cast(From v) {
        static_assert(sizeof(To) == sizeof(From)); 
        union { From f; To t; } u;  u.f = v; 
        return u.t;
    }
    template <> __device__ __forceinline__ uint32_t bit_cast(float v)    { return __float_as_uint(v);  }
    template <> __device__ __forceinline__ float    bit_cast(uint32_t v) { return __uint_as_float(v);  }
    template <> __device__ __forceinline__ uint64_t bit_cast(double v)   { return __double_as_longlong(v); }
    template <> __device__ __forceinline__ double   bit_cast(uint64_t v) { return __longlong_as_double(v); }

    template <typename T>
    static __device__ __forceinline__ T  get_bitfield(T v, int pos, int len) {
        // retrieves a contiguous sequence of bits from value (v), starting at bit position (pos) and spanning len bits (R2L). 
        // the extracted bits are returned right-aligned (shifted to the least-significant position)
        T mask = (len == (int)(sizeof(T) * 8)) ? ~T(0) : ((T(1) << len) - 1);
        return (v >> pos) & mask;
    }
            
    template <typename T>
    static __device__ __forceinline__ T set_bitfield(T val, T to_insert, int pos, int len) {
        T fieldMask = (len == (int)(sizeof(T) * 8)) ? ~T(0) : ((T(1) << len) - 1);
        T mask = fieldMask << pos;
        return (val & ~mask) | ((to_insert << pos) & mask);
    }
}