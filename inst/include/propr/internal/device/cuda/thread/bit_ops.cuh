#pragma once

#include <propr/utils/common/preprocessor.cuh>

namespace propr {
    namespace internal {
        namespace cuda {
            namespace thread {                
                template <typename T>
                PROPR_DEVICE PROPR_FORCE_INLINE 
                T get_bitfield(T v, int pos, int len) {
                    // retrieves a contiguous sequence of bits from value (v), starting at bit position (pos) and spanning len bits (R2L). 
                    // the extracted bits are returned right-aligned (shifted to the least-significant position)
                    T mask = (len == (int)(sizeof(T) * 8)) ? ~T(0) : ((T(1) << len) - 1);
                    return (v >> pos) & mask;
                }
                
                template <typename T>
                PROPR_DEVICE PROPR_FORCE_INLINE 
                T set_bitfield(T val, T to_insert, int pos, int len) {
                    T fieldMask = (len == (int)(sizeof(T) * 8)) ? ~T(0) : ((T(1) << len) - 1);
                    T mask = fieldMask << pos;
                    return (val & ~mask) | ((to_insert << pos) & mask);
                }
            }
        }
    }
}