#' Create permuted data
#'
#' This function creates p permuted data matrices
#'
#' This function wraps \code{updatePermutes.propr} and
#'  \code{updatePermutes.propd}.
#'
#' @param object A \code{propr} or \code{propd} object.
#' @param p The number of permutations to perform. Default is 100.
#' @param permutation_option A character string indicating if permute the data
#'  sample-wise or feature-wise. Default is "feature-wise". Note that this flag
#'  is only relevant for \code{propr} objects.
#' @return A \code{propr} or \code{propd} object with the permutes slot updated.
#' @export
updatePermutes <- function(object, p = 100, permutation_option = c("feature-wise", "sample-wise")) {
  NVTX_PUSH("updatePermutes", 0)
  
  if (inherits(object, "propr")) {
    res <- updatePermutes.propr(object, p, permutation_option)
    NVTX_POP()
    return(res)

  } else if (inherits(object, "propd")) {
    res <- updatePermutes.propd(object, p)
    NVTX_POP()
    return(res)

  } else {
    NVTX_POP()
    stop("Provided 'object' not recognized.")
  }
}

updatePermutes.propr <- function(object, p, permutation_option = c("feature-wise", "sample-wise")) {
  NVTX_PUSH("updatePermutes.propr", 0)
  
  message("Alert: Fixing permutations to active random seed.")
  ct <- object@counts
  
  NVTX_PUSH("validate_permutation_option", 0)
  option <- permutation_option[1]
  if (!option %in% c("feature-wise", "sample-wise")) {
    NVTX_POP()  # validate_permutation_option
    NVTX_POP()  # updatePermutes.propr
    stop("Invalid permutation option. Choose either 'feature-wise' or 'sample-wise'.")
  }
  NVTX_POP()  # validate_permutation_option
  
  NVTX_PUSH("generate_permutes", 0)
  permutes <- vector("list", p)
  for (ins in seq_len(p)) {
    if (option == "feature-wise") {
      # Permute features
      permutes[[ins]] <- apply(ct, 2, sample)
    } else { # option == "sample-wise"
      # Permute samples
      permutes[[ins]] <- t(apply(ct, 1, sample))
    }
  }
  NVTX_POP()  # generate_permutes
  
  NVTX_PUSH("update_object", 0)
  object@permutes <- permutes
  NVTX_POP()  # update_object
  
  NVTX_POP()  # updatePermutes.propr
  return(object)
}

updatePermutes.propd <- function(object, p) {
  NVTX_PUSH("updatePermutes.propd", 0)
  
  message("Alert: Fixing permutations to active random seed.")
  ct <- object@counts
  
  NVTX_PUSH("generate_permutes_matrix", 0)
  permutes <- as.data.frame(matrix(0, nrow = nrow(ct), ncol = p))
  for (col in seq_len(p)) {
    permutes[, col] <- sample(seq_len(nrow(ct)))
  }
  NVTX_POP()  # generate_permutes_matrix
  
  NVTX_PUSH("update_object", 0)
  object@permutes <- permutes
  NVTX_POP()  # update_object
  
  NVTX_POP()  # updatePermutes.propd
  return(object)
}
