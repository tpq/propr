#pragma once

#include <Rcpp.h>

namespace propr {
    namespace dispatch {
        namespace cpu {
              void genewise_connectivity(
                Rcpp::IntegerVector& per_gene_count,
                Rcpp::IntegerVector& per_gene_conn,
                Rcpp::NumericVector& per_gene_wconn,
                Rcpp::NumericVector& per_gene_fdr_sum,
                const Rcpp::IntegerVector& partner,
                const Rcpp::IntegerVector& pair,
                const Rcpp::NumericVector& theta,
                const Rcpp::NumericVector& fdr,
                int num_genes,
                double fdr_thresh);
        }
    }
}