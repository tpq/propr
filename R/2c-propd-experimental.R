#' Get Per-Feature Theta
#'
#' This function calculates the differential proportionality
#'  between each feature and a set of normalization factors. When the
#'  normalization factors correctly remove the compositional bias, the
#'  resultant thetas indicate differential expression (DE). However, unlike
#'  other DE tests, the p-value for differential proportionality is
#'  not linked to the normalization factors. Here, normalization factors
#'  only affect the interpretation, not the statistics.
#'
#' @param object A \code{\link{propd}} object.
#' @param norm.factors A numeric vector. The effective library size
#'  normalization factors (e.g., from edgeR or DESeq2).
#' @return A numeric vector. A theta for each feature.
#' @export
runNormalization <- function(object, norm.factors) {
  if (!inherits(object, "propd")) {
    stop("Please provide a propd object.")
  }
  if (!identical(length(norm.factors), nrow(object@counts))) {
    stop("The norm factors should have one value for each subject.")
  }

  # compute thetas
  newCounts <- cbind(norm.factors, object@counts)
  newPD <-
    propd(
      newCounts,
      group = object@group,
      alpha = object@alpha,
      p = 0,
      weighted = object@weighted
    )
  if (object@active == "theta_mod") {
    newPD <- updateF(newPD, moderated = TRUE)
  }
  newPD <- setActive(newPD, object@active)

  # parse thetas for each gene
  rawRes <- newPD@results
  perFeature <- rawRes[rawRes$Pair == 1,]
  if (!identical(perFeature$Partner, 2:(ncol(newCounts))))
    stop("DEBUG ERROR #GET001.")
  thetas <- perFeature$theta
  names(thetas) <- colnames(object@counts)
  
  return(thetas)
}

#' Perform Post-Hoc Testing
#'
#' After running an ANOVA of more than 2 groups, it is useful
#'  to know which of the groups differ from the others. This
#'  question is often answered with post-hoc testing. This
#'  function implements post-hoc pairwise differential
#'  proportionality analyses for more than 2 groups.
#'
#' The ANOVA p-values are adjusted once (column-wise) during
#'  the original multi-group analysis. The post-hoc p-values
#'  are adjusted once (row-wise) for the number
#'  of post-hoc tests. The post-hoc adjustment
#'  is p times the number of post-hoc tests.
#'
#' Please note that a significant post-hoc test without
#'  a significant ANOVA test is not significant!
#'
#' @param object A \code{\link{propd}} object.
#' @return A \code{\link{propd}} object.
#' @export
runPostHoc <- function(object) {
  groups <- unique(object@group)
  if (!length(groups) > 2) {
    stop("This function requires more than 2 groups.")
  }

  if (!"Pval" %in% colnames(object@results)) {
    message("Alert: Calculating ANOVA p-values without moderation.")
    object <- updateF(object)
  }

  for (i in 1:length(groups)) {
    for (j in 1:length(groups)) {
      if (j < i) {
        group1 <- groups[i]
        group2 <- groups[j]

        index <- object@group == group1 | object@group == group2
        x.ij <- object@counts[index, ]
        y.ij <- object@group[index]
        object.ij <-
          suppressMessages(propd(
            x.ij,
            y.ij,
            alpha = object@alpha,
            weighted = object@weighted
          ))

        if (is.na(object@Fivar) | is.null(object@Fivar)) {
          mod <- FALSE
        } else{
          mod <- TRUE
        }

        object.ij <-
          suppressMessages(updateF(object.ij, moderated = mod))
        new_result <- data.frame(object.ij@results$Pval)
        colnames(new_result) <-
          paste0(group1, ".vs.", group2, ".adj")
        ntests <- length(groups) * (length(groups) - 1) / 2
        object@results <- cbind(object@results, new_result * ntests)
      }
    }
  }

  message("Alert: Use 'getResults' function to obtain post-hoc tests.")
  message("Alert: Use 'Pval' column for ANOVA significance.")
  return(object)
}