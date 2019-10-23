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
#'  a significant ANOVA test is never significant!
#'
#' @param object A \code{\link{propd}} object.
#' @return A \code{\link{propd}} object.
#' @export
posthoc <- function(object){

  groups <- unique(object@group)
  if(!length(groups) > 2){
    stop("This function requires more than 2 groups.")
  }

  if(!"Pval" %in% colnames(object@results)){
    message("Alert: Calculating ANOVA p-values without moderation.")
    object <- updateF(object)
  }

  for(i in 1:length(groups)){

    for(j in 1:length(groups)){

      if(j < i){

        group1 <- groups[i]
        group2 <- groups[j]

        index <- object@group == group1 | object@group == group2
        x.ij <- object@counts[index, ]
        y.ij <- object@group[index]
        object.ij <- suppressMessages(propd(x.ij, y.ij, alpha = object@alpha, weighted = object@weighted))

        if(is.na(object@Fivar) | is.null(object@Fivar)){
          mod <- FALSE
        }else{
          mod <- TRUE
        }

        object.ij <- suppressMessages(updateF(object.ij, moderated = mod))
        new_result <- data.frame(object.ij@results$Pval)
        colnames(new_result) <- paste0(group1, ".vs.", group2, ".adj")
        ntests <- length(groups)*(length(groups)-1)/2
        object@results <- cbind(object@results, new_result*ntests)
      }
    }
  }

  message("Alert: Use 'getResults' function to obtain post-hoc tests.")
  message("Alert: Use 'Pval' column for ANOVA significance.")
  return(object)
}
