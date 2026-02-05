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

    if (req == "cuda") {
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

}  // namespace propr::runtime
