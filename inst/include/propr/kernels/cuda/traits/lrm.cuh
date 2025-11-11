#pragma once

#include <propr/kernels/cuda/traits/common.cuh>

namespace propr {
    namespace cuda {
        namespace traits {
            struct lrm_basic          : thread_layout_2d<>{};
            struct lrm_weighted       : thread_layout_2d<>{};
            struct lrm_alpha          : thread_layout_2d<>{};
            struct lrm_alpha_weighted : thread_layout_2d<>{};
        }
    }
}