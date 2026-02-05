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
  NVTX_PUSH("runNormalization", 0)

  NVTX_PUSH("validate_input", 0)
  if (!inherits(object, "propd")) {
    NVTX_POP() # validate_input
    NVTX_POP() # runNormalization
    stop("Please provide a propd object.")
  }
  if (!identical(length(norm.factors), nrow(object@counts))) {
    NVTX_POP() # validate_input
    NVTX_POP() # runNormalization
    stop("The norm factors should have one value for each subject.")
  }
  NVTX_POP() # validate_input

  # compute thetas
  NVTX_PUSH("compute_thetas", 0)
  newCounts <- cbind(norm.factors, object@counts)

  NVTX_PUSH("propd", 0)
  newPD <-
    propd(
      newCounts,
      group = object@group,
      alpha = object@alpha,
      p = 0,
      weighted = object@weighted
    )
  NVTX_POP() # propd

  if (object@active == "theta_mod") {
    NVTX_PUSH("updateF_moderated", 0)
    newPD <- updateF(newPD, moderated = TRUE)
    NVTX_POP() # updateF_moderated
  }

  NVTX_PUSH("setActive", 0)
  newPD <- setActive(newPD, object@active)
  NVTX_POP() # setActive
  NVTX_POP() # compute_thetas

  # parse thetas for each gene
  NVTX_PUSH("parse_thetas", 0)
  rawRes <- newPD@results
  perFeature <- rawRes[rawRes$Pair == 1, ]
  if (!identical(perFeature$Partner, 2:(ncol(newCounts)))) {
    NVTX_POP() # parse_thetas
    NVTX_POP() # runNormalization
    stop("DEBUG ERROR #GET001.")
  }
  thetas <- perFeature$theta
  names(thetas) <- colnames(object@counts)
  NVTX_POP() # parse_thetas

  NVTX_POP() # runNormalization
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
  NVTX_PUSH("runPostHoc", 0)

  NVTX_PUSH("check_groups", 0)
  groups <- unique(object@group)
  if (!length(groups) > 2) {
    NVTX_POP() # check_groups
    NVTX_POP() # runPostHoc
    stop("This function requires more than 2 groups.")
  }
  NVTX_POP() # check_groups

  NVTX_PUSH("check_Pval_column", 0)
  if (!"Pval" %in% colnames(object@results)) {
    message("Alert: Calculating ANOVA p-values without moderation.")
    NVTX_PUSH("updateF_ANOVA", 0)
    object <- updateF(object)
    NVTX_POP() # updateF_ANOVA
  }
  NVTX_POP() # check_Pval_column

  NVTX_PUSH("loop_over_pairs", 0)
  for (i in 1:length(groups)) {
    for (j in 1:length(groups)) {
      if (j < i) {
        NVTX_PUSH("pair_computation", 0)

        group1 <- groups[i]
        group2 <- groups[j]

        NVTX_PUSH("build_pair_subset", 0)
        index <- object@group == group1 | object@group == group2
        x.ij <- object@counts[index, ]
        y.ij <- object@group[index]
        NVTX_POP() # build_pair_subset

        NVTX_PUSH("propd_pair", 0)
        object.ij <-
          suppressMessages(propd(
            x.ij,
            y.ij,
            alpha = object@alpha,
            weighted = object@weighted
          ))
        NVTX_POP() # propd_pair

        NVTX_PUSH("determine_moderation", 0)
        if (is.na(object@Fivar) | is.null(object@Fivar)) {
          mod <- FALSE
        } else {
          mod <- TRUE
        }
        NVTX_POP() # determine_moderation

        NVTX_PUSH("updateF_pair", 0)
        object.ij <-
          suppressMessages(updateF(object.ij, moderated = mod))
        NVTX_POP() # updateF_pair

        NVTX_PUSH("update_results", 0)
        new_result <- data.frame(object.ij@results$Pval)
        colnames(new_result) <-
          paste0(group1, ".vs.", group2, ".adj")
        ntests <- length(groups) * (length(groups) - 1) / 2
        object@results <- cbind(object@results, new_result * ntests)
        NVTX_POP() # update_results

        NVTX_POP() # pair_computation
      }
    }
  }
  NVTX_POP() # loop_over_pairs

  message("Alert: Use 'getResults' function to obtain post-hoc tests.")
  message("Alert: Use 'Pval' column for ANOVA significance.")

  NVTX_POP() # runPostHoc
  return(object)
}
