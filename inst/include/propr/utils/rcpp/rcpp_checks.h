#pragma once

#include <Rcpp.h>
#include <propr/data/types.h>

#define PROPR_CHECK_MATRIX_DIMS(mat, expected_rows, expected_cols) \
    if (mat.nrow() != expected_rows || mat.ncol() != expected_cols) { \
        Rcpp::stop("Error: Output matrix dimensions do not match expected dimensions."); \
    }

#define PROPR_CHECK_VECTOR_SIZE(vec, expected_size) \
    if (vec.length() != expected_size) { \
        Rcpp::stop("Error: Output vector size does not match expected size."); \
    }
