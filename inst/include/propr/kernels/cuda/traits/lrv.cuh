#pragma once

#include <propr/kernels/cuda/traits/common.cuh>

namespace propr {
    namespace cuda {
        namespace traits {
            struct lrv_basic          : thread_layout_2d<>{};
            struct lrv_weighted       : thread_layout_2d<>{};
            struct lrv_alpha          : thread_layout_2d<>{};
            struct lrv_alpha_weighted : thread_layout_2d<>{};
        }
    }
}