#pragma once

#include <Rcpp.h>

namespace propr {
    namespace dispatch {
        namespace cpu    {
            
            __forceinline__
            void dof_global    (Rcpp::NumericVector& out, Rcpp::NumericMatrix& W);
            
            __forceinline__
            void dof_population(Rcpp::NumericVector& out, const Rcpp::NumericMatrix& W);
        }
    }
}
