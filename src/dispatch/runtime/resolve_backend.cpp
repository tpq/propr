#include <propr/runtime/dispatch.hpp>

#include <algorithm>
#include <string>

namespace propr::runtime {

    #if !defined(PROPR_HAS_CUDA) || !PROPR_HAS_CUDA
        bool cuda_is_available() { return false; }
    #endif

    Backend resolve_backend(const Rcpp::String& requested) {
        std::string req(requested.get_cstring());
        std::transform(req.begin(), req.end(), req.begin(), ::tolower);

        if (req == "auto") {
            Rcpp::Function getOption("getOption");
            SEXP val = getOption("propr.backend");
            if (!Rf_isNull(val)) {
                req = Rcpp::as<std::string>(val);
                std::transform(req.begin(), req.end(), req.begin(), ::tolower);
            } else {
                req = "cpu";
            }
        }

        if (req == "cuda" || req == "gpu" ) {
            if (cuda_is_available()) return Backend::CUDA;
            static bool warned = false;
            if (!warned) {
                warned = true;
                Rcpp::warning("CUDA backend requested but not available; falling back to CPU.");
            }
            return Backend::CPU;
        }
        return Backend::CPU;
    }

    Precision resolve_precision() {
        Rcpp::Function getOption("getOption");
        SEXP val = getOption("propr.precision", "double");

        if (Rf_isNull(val)) return Precision::Float64;
        if (!Rf_isString(val) || Rf_length(val) != 1) {
            Rcpp::stop("Option 'propr.precision' must be 'float' or 'double'.");
        }

        std::string req = Rcpp::as<std::string>(val);
        std::transform(req.begin(), req.end(), req.begin(), ::tolower);
        if (req == "float") return Precision::Float32;
        if (req == "double") return Precision::Float64;
        Rcpp::stop("Option 'propr.precision' must be 'float' or 'double'.");
    }

}  // namespace propr::runtime
