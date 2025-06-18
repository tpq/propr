#pragma once

#include <Rcpp.h>

namespace propr {
    namespace dispatch {
        namespace cpu {
            void ctzRcpp(Rcpp::NumericVector& out, Rcpp::NumericMatrix& X);
        }
    }
}
