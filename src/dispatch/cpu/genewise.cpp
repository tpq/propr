#include<propr/kernels/cpu/dispatch/genewise.hpp>

using namespace Rcpp;
using namespace propr;

void
dispatch::cpu::genewise_connectivity(
                Rcpp::IntegerVector& per_gene_count,
                Rcpp::IntegerVector& per_gene_conn,
                Rcpp::NumericVector& per_gene_wconn,
                Rcpp::NumericVector& per_gene_fdr_sum,
                const Rcpp::IntegerVector& partner,
                const Rcpp::IntegerVector& pair,
                const Rcpp::NumericVector& theta,
                const Rcpp::NumericVector& fdr,
                int num_genes,
                double fdr_thresh) {
    
    for (int i = 0; i < num_genes; ++i) {
        const int p = partner[i];
        const int q = pair[i];

        if (p == NA_INTEGER || q == NA_INTEGER)               Rcpp::stop("partner/pair cannot contain NA values.");
        if (p < 1 || p > num_genes || q < 1 || q > num_genes) Rcpp::stop("partner/pair indices must be in [1, num_genes].");

        const int p0 = p - 1;
        const int q0 = q - 1;
        const double f = fdr[i];
        const bool has_f = !R_IsNA(f) && !R_IsNaN(f);

        if (has_f) {
            per_gene_count[p0]   += 1;
            per_gene_count[q0]   += 1;
            per_gene_fdr_sum[p0] += f;
            per_gene_fdr_sum[q0] += f;
        }

        if (has_f && f > 0.0 && f < fdr_thresh) {
            per_gene_conn[p0] += 1;
            per_gene_conn[q0] += 1;
            const double t = theta[i];
            if (!R_IsNA(t) && !R_IsNaN(t)) {
                const double w = 1.0 - t;
                per_gene_wconn[p0] += w;
                per_gene_wconn[q0] += w;
            }
        }
    }
}