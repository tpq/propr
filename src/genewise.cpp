#include <Rcpp.h>

#include <propr/interface/genewise.hpp>
#include <propr/runtime/cuda_executor.hpp>
#include <propr/runtime/dispatch.hpp>
#include <propr/kernels/cpu/dispatch/genewise.hpp>

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
            partner, pair, theta, fdr, num_genes, pairwise_fdr);
    } else {
        dispatch::cpu::genewise_connectivity(
            per_gene_count,
            per_gene_conn,
            per_gene_wconn,
            per_gene_fdr_sum, 
            partner, pair, theta, fdr, num_genes, pairwise_fdr);
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
