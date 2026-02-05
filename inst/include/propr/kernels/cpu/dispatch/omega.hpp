#pragma once

#include <Rcpp.h>

namespace propr {
    namespace dispatch {
        namespace cpu    {
            void dof_global    (Rcpp::NumericVector& out, Rcpp::NumericMatrix& W);
            void dof_population(Rcpp::NumericVector& out, const Rcpp::NumericMatrix& W);
        }
    }
}
