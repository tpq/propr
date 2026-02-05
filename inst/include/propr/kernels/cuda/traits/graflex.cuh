#pragma once

#include <propr/kernels/cuda/traits/common.cuh>

namespace propr {
    namespace cuda {
        namespace traits {
            struct getOR_config     : thread_layout_1d<1024> {};
            struct getORperm_config : thread_layout_1d<1024> {};
        }
    }
}