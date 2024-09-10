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
  if (any(ct == 0)) {
    message("Alert: replacing zeros with minimun value.")
    zeros <- ct == 0
    ct[zeros] <- min(ct[!zeros])
  } else{
    message("Alert: No 0s found that need replacement.")
  }
  return(ct)
}
