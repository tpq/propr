#pragma once

#include <Rcpp.h>

namespace propr {
    namespace dispatch {
        namespace cpu {
            __forceinline__
            void ctzRcpp(Rcpp::NumericVector& out, Rcpp::NumericMatrix& X);
        }
    }
}
