#pragma once

#include <random>
#include <Rcpp.h>

Rcpp::NumericMatrix generate_counts(int rows, int cols, int min_val = 0, int max_val = 1000, double zero_prob = 0.05) {
    Rcpp::NumericMatrix mat(rows, cols);
    std::default_random_engine generator;
    generator.seed(42);

    std::uniform_int_distribution<int> int_distribution(min_val, max_val);
    std::uniform_real_distribution<double> real_distribution(0.0, 1.0);
    for (int j = 0; j < cols; ++j) {
        for (int i = 0; i < rows; ++i) {
            if (min_val == 0 && real_distribution(generator) < zero_prob) {
                mat(i, j) = 0;
            } else {
                mat(i, j) = int_distribution(generator);
            }
        }
    }
    return mat;
}

Rcpp::NumericMatrix apply_zero_replacement(const Rcpp::NumericMatrix& counts_in) {
    Rcpp::NumericMatrix counts_out = Rcpp::clone(counts_in);
    bool has_zeros = false;
    for (int j = 0; j < counts_out.ncol(); ++j) {
        for (int i = 0; i < counts_out.nrow(); ++i) {
            if (counts_out(i, j) == 0) {
                has_zeros = true;
                break;
            }
        }
        if (has_zeros) break;
    }

    if (has_zeros) {
        double min_non_zero = std::numeric_limits<double>::infinity();
        bool found_non_zero = false;
        for (int j = 0; j < counts_out.ncol(); ++j) {
            for (int i = 0; i < counts_out.nrow(); ++i) {
                if (counts_out(i, j) > 0 && counts_out(i, j) < min_non_zero) {
                    min_non_zero = counts_out(i, j);
                    found_non_zero = true;
                }
            }
        }

        if (!found_non_zero)  min_non_zero = 1.0;
        for (int j = 0; j < counts_out.ncol(); ++j) {
            for (int i = 0; i < counts_out.nrow(); ++i) {
                if (counts_out(i, j) == 0) {
                    counts_out(i, j) = min_non_zero;
                }
            }
        }
    }
    return counts_out;
}

Rcpp::NumericMatrix generate_weights_matrix(int rows, int cols, double min_val = 0.1, double max_val = 1.0) {
    Rcpp::NumericMatrix mat(rows, cols);
    std::default_random_engine generator;
    generator.seed(42);
    std::uniform_real_distribution<double> distribution(min_val, max_val);
    for (int j = 0; j < cols; ++j) {
        for (int i = 0; i < rows; ++i) {
            mat(i, j) = distribution(generator);
        }
    }
    return mat;
}