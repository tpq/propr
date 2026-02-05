#pragma once

#include <cstdint>
#include <type_traits>


#include <cuda_runtime.h>


#include <propr/data/bitwise.cuh>
#include <propr/data/traits.cuh>

namespace propr {
    namespace radix {
        
        template <typename T>
        static __device__ __forceinline__ 
        unsigned_type_t<T> radix_convert(T v) {
            using U = unsigned_type_t<T>;
            constexpr U S = U(1) << (sizeof(U) * 8 - 1);
            if constexpr (is_floating_v<T>)        { U x = bit_cast<U>(v); return x ^ ((x & S) ? ~U(0) : S); }
            else if constexpr (is_signed_int_v<T>)  return bit_cast<U>(v) ^ S; 
            else return v;
        }

        template <typename T>
        static __device__ __forceinline__ 
        T radix_deconvert(unsigned_type_t<T> v) {
            using U = unsigned_type_t<T>;
            constexpr U S = U(1) << (sizeof(U) * 8 - 1);
            if constexpr (is_floating_v<T>)        return bit_cast<T>(v ^ ((v & S) ? S : ~U(0))); 
            else if constexpr (is_signed_int_v<T>) return bit_cast<T>(v ^ S); 
            else return v;
        }

    } // namespace radix
}// namespace propr