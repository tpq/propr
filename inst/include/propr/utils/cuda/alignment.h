#pragma once

#include <cstdint>
#include <cstdio>

namespace propr {
    namespace cuda {

        inline void check_pointer_alignment(const void* ptr, int alignment) {
            if (reinterpret_cast<uintptr_t>(ptr) % (alignment * sizeof(float)) != 0) {
                printf("ERROR: Misaligned access at %p, required alignment: %d bytes\n", 
                    ptr, alignment * static_cast<int>(sizeof(float)));
            }
        }

    } // namespace cuda
} // namespace propr
