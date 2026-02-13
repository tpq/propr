#' Convert pairwise propd results into genewise
#' 
#' This function converts pairwise propd results into genewise results. The resulting
#' genewise results are direct indirectors of genes being differentially expressed.
#' 
#' @param propd A \code{\link{propd}} object, with FDR values from updateF. Note:
#' for the moment, only theta results with F-stats are supported. Later we will look 
#' on the option to get FDR values based on permutations.
#' @param pairwise_fdr FDR threashold to consider a pairwise relationship as significant.
#' Default is 0.05.
#' @return A data frame with genewise results that can be used to identify differentially
#' expressed genes. It contains the following columns:
#'  - "id": gene identifier
#'  - "lfc": Log Fold Change of the gene, using the geometric mean of all genes as reference.
#'  - "lrmD": Log Ratio Mean Difference of the gene. Equivalent to the LFC, but using a subset 
#'     of genes as reference (only the ones that are significantly connected to the gene).
#'  - "connectivity": number of significant pairwise relationships the gene has.
#'  - "wconnectivity": weighted connectivity, which sums the strength of significant pairwise 
#'     relationships (defined as 1/theta) for the gene.
#' - "FDR_mean": average FDR value across all pairwise comparisons for the gene, which can be 
#'     used as a genewise significance statistic. Note that this is a simple average and may 
#'     not be the most robust method for determining genewise significance, but it provides a 
#'     starting point for identifying genes of interest based on their pairwise relationships. 
#'     Future updates may include more sophisticated methods for calculating genewise p-values 
#'     or FDR values based on the pairwise results.
#' 
#' @rdname propdGenewise
#' @export
propdGenewise <- function(propd, pairwise_fdr = 0.05) {

  # for the moment it only works with theta results with F-stats,
  # but later we will look into the option to get FDR values based on permutations.
  if (!"FDR" %in% colnames(propd@results)) {
    stop("Please run updateF on the propd object to get FDR values before running propdGenewise.")
  }

  # not working for alpha != NA too
  if (!is.na(propd@alpha)) {
    stop("propdGenewise currently only works for alpha = NA. Future updates may include support for other alpha values.")
  }

  # not working for more than 2 groups too
  if (length(unique(propd@group)) > 2) {
    stop("propdGenewise currently only works for 2 groups. Future updates may include support for more than 2 groups.")
  }

  # get features and number of features
  features <- colnames(propd@counts)
  nfeatures <- length(features)

  ## ---- Build matrices needed for connectivity metrics ----
  fdr_mat <- results_to_matrix(propd@results, what = "FDR", features = features)
  theta_mat <- results_to_matrix(propd@results, what = "theta", features = features)
  adj <- (fdr_mat > 0) & (fdr_mat < pairwise_fdr)

  ## ---- Connectivity ----
  connectivity <- rowSums(adj)

  ## ---- Weighted connectivity ----
  weight_mat <- ifelse(adj, 1 / theta_mat, 0)
  wconnectivity <- rowSums(weight_mat)

  ## ---- FDR mean ----
  fdr_mean <- rowMeans(fdr_mat)

  ## ---- Build lrm matrices ----
  lrm1_all <- results_to_matrix(propd@results, what = "lrm1", features = features)
  lrm2_all <- results_to_matrix(propd@results, what = "lrm2", features = features)
  # results_to_matrix returns symmetric matrices, but lrm values are directed:
  # lrm(Partner, Pair) = mean(log(x_Partner / x_Pair)), with Partner > Pair.
  # Negate the upper triangle so that mat[g, j] = mean(log(x_g / x_j)) for all g, j.
  lrm1_all[upper.tri(lrm1_all)] <- -lrm1_all[upper.tri(lrm1_all)]
  lrm2_all[upper.tri(lrm2_all)] <- -lrm2_all[upper.tri(lrm2_all)]

  ## ---- LFC (CLR-based log fold change) ----
  # lrm differences represent log fold changes; averaging across all genes as
  # reference is equivalent to using the geometric mean (CLR transformation)
  lrm_diff <- lrm1_all - lrm2_all
  lfc <- rowMeans(lrm_diff) / log(2)

  ## ---- lrmD (LFC using only significant partners as reference) ----
  lrm_diff <- ifelse(adj, lrm_diff, NA) # keep only significant pairwise relationships
  lrmD <- apply(lrm_diff, 1, median, na.rm = TRUE) / log(2)

  ## ---- Compile results into a data frame ----
  data.frame(
    id = features,
    lfc = lfc,
    lrmD = lrmD,
    connectivity = connectivity,
    wconnectivity = wconnectivity,
    FDR_mean = fdr_mean,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}
