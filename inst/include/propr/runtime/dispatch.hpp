#pragma once

#include <Rcpp.h>
#include <stdexcept>
#include <string>
#include <utility>

namespace propr {
    namespace runtime {

        enum class Backend { CPU, CUDA };
        enum class Precision { Float32, Float64 };

        template <typename T>
        struct precision_tag {
            using type = T;
        };

        // Resolves "auto"/"cpu"/"cuda" to the actual backend to use.
        // - "auto": reads R option `propr.backend`, defaults to "cpu"
        // - "cuda": checks cuda_is_available(), warns and returns CPU if not
        // - "cpu": always returns CPU
        Backend resolve_backend(const Rcpp::String& requested);

        // Reads R option `propr.precision`, defaults to "double".
        Precision resolve_precision();

        template <typename F>
        decltype(auto) with_precision(F&& f) {
            switch (resolve_precision()) {
                case Precision::Float32: return std::forward<F>(f)(precision_tag<float>{});
                case Precision::Float64: return std::forward<F>(f)(precision_tag<double>{});
            }
            throw std::logic_error("Internal error: unknown propr.precision.");
        }

        // Whether CUDA runtime is present and has devices.
        // Defined in cuda_available.cu (CUDA builds) or resolve_backend.cpp (CPU-only).
        bool cuda_is_available();

    }  // namespace runtime
}  // namespace propr
