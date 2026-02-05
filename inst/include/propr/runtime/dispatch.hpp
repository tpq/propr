#pragma once

#include <Rcpp.h>
#include <string>

namespace propr {
    namespace runtime {

        enum class Backend { CPU, CUDA };

        // Resolves "auto"/"cpu"/"cuda" to the actual backend to use.
        // - "auto": reads R option `propr.backend`, defaults to "cpu"
        // - "cuda": checks cuda_is_available(), warns and returns CPU if not
        // - "cpu": always returns CPU
        Backend resolve_backend(const Rcpp::String& requested);

        // Whether CUDA runtime is present and has devices.
        // Defined in cuda_available.cu (CUDA builds) or resolve_backend.cpp (CPU-only).
        bool cuda_is_available();

    }  // namespace runtime
}  // namespace propr
