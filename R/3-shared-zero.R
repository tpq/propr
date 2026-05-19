#' Simple Zero Replacement in a Count Matrix
#'
#' This function replaces zeros with the next smallest non-zero value in the
#'  input count matrix. If the matrix contains no zeros, it produces an
#'  informational message indicating that no replacements were made.
#'
#' @param ct A data matrix containing numerical values.
#' @return A matrix with zero values replaced by the next smallest non-zero value.
#' If no zeros are found, the function returns the original matrix.
#' @examples
#' # Sample input count data with zeros
#' data <- matrix(c(0, 2, 3, 4, 5, 0), nrow = 2, byrow = TRUE)
#' @export
simple_zero_replacement <- function(ct) {
  NVTX_PUSH("simple_zero_replacement", 0)
  
  if (any(ct == 0)) {
    NVTX_PUSH("zero_replacement_branch", 0)
    message("Alert: replacing zeros with minimun value.")
    zeros <- ct == 0
    NVTX_PUSH("compute_min_nonzero", 0)
    ct[zeros] <- min(ct[!zeros])
    NVTX_POP()  # compute_min_nonzero
    NVTX_POP()  # zero_replacement_branch
  } else {
    NVTX_PUSH("no_zero_branch", 0)
    message("Alert: No 0s found that need replacement.")
    NVTX_POP()  # no_zero_branch
  }
  
  NVTX_POP()  # simple_zero_replacement
  return(ct)
}
