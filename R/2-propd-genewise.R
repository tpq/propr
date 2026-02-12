#' Convert pairwise propd results into genewise
#' 
#' This function converts pairwise propd results into genewise results. The resulting
#' genewise results are direct indirectors of genes being differentially expressed.
#' 
#' @param propd A \code{\link{propd}} object, with FDR values from updateF. Note:
#' for the moment, only theta results with F-stats are supported. Later we will look 
#' on the option to get FDR values based on permutations.
#' @param metric A character string. The metric to use for genewise results. Options 
#' are "connectivity" and "wconnectivity".
#' @param pairwise_fdr 
#' @return A data frame with genewise results. 
#' 
#' @details The "connectivity" metric refers to the number of significant pairwise 
#' relationships a gene has, while the "wconnectivity" metric weights these 
#' relationships by their strength (i.e., theta values). The resulting genewise results 
#' can be used to identify genes that are central to the observed proportionality changes 
#' across groups (these genes tend to be differentially expressed).
#' 
#' @rdname propdGenewise
#' @export
propdGenewise <- function(propd,
                          metric = c("connectivity", "wconnectivity"),
                          pairwise_fdr = 0.05) {

  metric <- match.arg(metric)

  if (!"FDR" %in% colnames(propd@results)) { 
    stop("Please run updateF on the propd object to get FDR values before running propdGenewise.") 
  }

  if (metric == "connectivity") {
    get_connectivity(propd, pairwise_fdr)
  } else {
    get_weighted_connectivity(propd, pairwise_fdr)
  }
}

#' Get connectivity for each gene based on pairwise propd results
#' 
#' It computes the per-gene connectivity by counting the number of
#' significant pairwise relationships each gene has. It also provides
#' a genewise significance stat (for the moment we provide a fake one
#' by averaging all the pairwise FDR values for each gene, but later
#' we will look into more robust methods to get genewise p-values).
#' @inheritParams propdGenewise
#' @return A data frame with genewise connectivity results.
get_connectivity <- function(propd, pairwise_fdr = 0.05) {

  features <- colnames(propd@counts)

  # Build FDR matrix and adjacency matrix of significant pairs
  fdr_mat <- results_to_matrix(propd@results, what = "FDR", features = features)
  adj <- (fdr_mat > 0) & (fdr_mat < pairwise_fdr)

  # Connectivity: number of significant pairs per gene (row sums of adjacency)
  connectivity <- rowSums(adj)

  # Mean FDR across all pairs per gene
  fdr_mean <- rowMeans(fdr_mat)

  data.frame(
    Gene = features,
    connectivity = connectivity,
    FDR = fdr_mean,
    stringsAsFactors = FALSE
  )
}

#' Get weighted connectivity for each gene based on pairwise propd results 
#' 
#' It computes the per-gene weighted connectivity by summing the strength 
#' of significant pairwise relationships each gene has. The strength of a 
#' relationship can be defined #' as the inverse of the theta value (i.e., 1/theta) 
#' (i.e., 1/theta) for significant pairs, which gives more weight to stronger 
#' relationships. It also provides a genewise significance stat by averaging 
#' all the pairwise FDR values for each gene.
#' @inheritParams propdGenewise
#' @return A data frame with genewise weighted connectivity results.
get_weighted_connectivity <- function(propd, pairwise_fdr = 0.05) {

  features <- colnames(propd@counts)

  # Build FDR matrix and theta matrix
  fdr_mat <- results_to_matrix(propd@results, what = "FDR", features = features)
  theta_mat <- results_to_matrix(propd@results, what = "theta", features = features)

  # Adjacency matrix of significant pairs
  adj <- (fdr_mat > 0) & (fdr_mat < pairwise_fdr)

  # Weighted connectivity: sum of 1/theta for significant pairs per gene
  weight_mat <- ifelse(adj, 1 / theta_mat, 0)
  wconnectivity <- rowSums(weight_mat)

  # Mean FDR across all pairs per gene
  fdr_mean <- rowMeans(fdr_mat)

  data.frame(
    Gene = features,
    wconnectivity = wconnectivity,
    FDR = fdr_mean,
    stringsAsFactors = FALSE
  )
}