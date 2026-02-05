#' Select Optimal Reference Component
#'
#' This function selects the optimal reference component from the log-ratio
#'  transformed data matrix based on the provided \code{ivar} (index variable)
#'  and \code{alpha} values.
#'
#' The function transforms the input \code{counts} matrix into log space using
#'  the \code{logratio} function. Then, it calculates the variance of each
#'  component and identifies the component with the minimum variance,
#'  which is considered the optimal reference.
#'
#' @inheritParams propr
#' @return The column name or index of the optimal reference component.
#' @examples
#' # Sample counts matrix
#' counts_matrix <- matrix(c(10, 20, 30, 40, 0, 50, 60, 70, 0), nrow = 3, byrow = TRUE)
#' colnames(counts_matrix) <- c("A", "B", "C")
#' rownames(counts_matrix) <- c("Sample1", "Sample2", "Sample3")
#'
#' # Select optimal reference component
#' selectReference(counts_matrix, ivar = "A", alpha = 0.5)
#'
#' @export
selectReference <- function(counts, ivar, alpha) {
  NVTX_PUSH("selectReference", 0)

  # replace zeros
  NVTX_PUSH("simple_zero_replacement", 0)
  counts <- simple_zero_replacement(counts)
  NVTX_POP()

  # Transform data into log space
  NVTX_PUSH("logratio", 0)
  lr <- logratio(counts, ivar, alpha)
  NVTX_POP()

  # Calculate var of each component
  NVTX_PUSH("variance_per_component", 0)
  vars <- apply(lr, 2, stats::var)
  NVTX_POP()

  NVTX_PUSH("select_min_var_reference", 0)
  res <- if (!is.null(colnames(counts))) {
    colnames(counts)[which.min(vars)]
  } else {
    which.min(vars)
  }
  NVTX_POP()

  NVTX_POP()  # selectReference
  return(res)
}
