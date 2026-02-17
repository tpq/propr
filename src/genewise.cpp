#include <Rcpp.h>

#include <propr/interface/genewise.hpp>
#include <propr/runtime/cuda_executor.hpp>
#include <propr/runtime/dispatch.hpp>

using namespace propr;

// [[Rcpp::export]]
Rcpp::List genewiseConnectivityRcpp(
    const Rcpp::IntegerVector& partner,
    const Rcpp::IntegerVector& pair,
    const Rcpp::NumericVector& theta,
    const Rcpp::NumericVector& fdr,
    int num_genes,
    double pairwise_fdr,
    Rcpp::String backend) {
    if (num_genes < 0) {
        Rcpp::stop("num_genes must be non-negative.");
    }

    const int num_edges = partner.size();
    if (pair.size() != num_edges || theta.size() != num_edges || fdr.size() != num_edges) {
        Rcpp::stop("partner, pair, theta, and fdr must have identical lengths.");
    }

    Rcpp::IntegerVector per_gene_count(num_genes);
    Rcpp::IntegerVector per_gene_conn(num_genes);
    Rcpp::NumericVector per_gene_wconn(num_genes);
    Rcpp::NumericVector per_gene_fdr_sum(num_genes);

    if (runtime::resolve_backend(backend) == runtime::Backend::CUDA) {
        runtime::cuda_executor::genewise_connectivity(
            per_gene_count,
            per_gene_conn,
            per_gene_wconn,
            per_gene_fdr_sum,
            partner,
            pair,
            theta,
            fdr,
            num_genes,
            pairwise_fdr);
    } else {
        for (int i = 0; i < num_edges; ++i) {
            const int p = partner[i];
            const int q = pair[i];

            if (p == NA_INTEGER || q == NA_INTEGER) {
                Rcpp::stop("partner/pair cannot contain NA values.");
            }
            if (p < 1 || p > num_genes || q < 1 || q > num_genes) {
                Rcpp::stop("partner/pair indices must be in [1, num_genes].");
            }

            const int p0 = p - 1;
            const int q0 = q - 1;
            const double f = fdr[i];
            const bool has_f = !R_IsNA(f) && !R_IsNaN(f);

            if (has_f) {
                per_gene_count[p0] += 1;
                per_gene_count[q0] += 1;
                per_gene_fdr_sum[p0] += f;
                per_gene_fdr_sum[q0] += f;
            }

            if (has_f && f > 0.0 && f < pairwise_fdr) {
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

    Rcpp::NumericVector fdr_mean(num_genes);
    for (int g = 0; g < num_genes; ++g) {
        const int denom = per_gene_count[g] + 1;  // +1 for diagonal zero in rowMeans(results_to_matrix(...), na.rm=TRUE)
        fdr_mean[g] = (denom > 0) ? (per_gene_fdr_sum[g] / static_cast<double>(denom)) : NA_REAL;
    }

    return Rcpp::List::create(
        Rcpp::Named("connectivity") = per_gene_conn,
        Rcpp::Named("wconnectivity") = per_gene_wconn,
        Rcpp::Named("FDR_mean") = fdr_mean);
}
