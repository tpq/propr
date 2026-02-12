#' Convert pairwise propd results into genewise
#' 
#' This function converts pairwise propd results into genewise results. The resulting
#' genewise results are direct indirectors of genes being differentially expressed.
#' 
#' @param propd A \code{\link{propd}} object. 
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

}

#' Get connectivity for each gene based on pairwise propd results
#' @inheritParams propdGenewise
#' @return A data frame with genewise connectivity results.
get_connectivity <- function(propd, pairwise_fdr = 0.05) {

}

#' Get weighted connectivity for each gene based on pairwise propd results 
#' @inheritParams propdGenewise
#' @return A data frame with genewise weighted connectivity results.
get_weighted_connectivity <- function(propd, pairwise_fdr = 0.05) { 

}