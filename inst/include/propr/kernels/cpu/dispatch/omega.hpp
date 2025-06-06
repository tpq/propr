#pragma once

#include <Rcpp.h>

namespace propr {
    namespace dispatch {
        namespace cpu    {
            void dof_global    (Rcpp::NumericMatrix& W, Rcpp::NumericVector& out);
            void dof_population(const Rcpp::NumericMatrix& W, Rcpp::NumericVector& out);
        }
    }
}
