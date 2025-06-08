#pragma once
#include <cuda_runtime.h>

namespace propr {
    typedef struct context {
        cudaStream_t stream;
    } propr_context;

    static const propr_context DEFAULT_GLOBAL_CONTEXT = { cudaStreamDefault };
}