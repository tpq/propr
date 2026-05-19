#pragma once

#include <Rcpp.h>
#include <propr/context.h>

namespace propr {
    namespace dispatch {
        namespace cuda {

            // partner/pair are 1-based feature indices as stored in propd@results.
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
                double fdr_thresh,
                propr_context context = DEFAULT_GLOBAL_CONTEXT);

            // theta_edges are expected in partner-major packed order:
            // idx = J*(J-1)/2 + I for I < J, 0-based.
            void genewise_theta_stats(
                Rcpp::NumericVector& out_mean,
                Rcpp::NumericVector& out_median,
                const Rcpp::NumericVector& theta_edges,
                int num_genes,
                propr_context context = DEFAULT_GLOBAL_CONTEXT);

        }
    }
}
