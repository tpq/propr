#' @param counts A data matrix representing counts.
#'   It is assumed that the matrix contains numerical values only.
#' @param metric A character vector indicating the metric used for computing
#'  the association matrix. It can take the following values:
#'   - "rho": Propr matrix based on the rho coefficient.
#'   - "phi": Propr matrix based on the phi coefficient.
#'   - "phs": Propr matrix based on the symmetric phi coefficient.
#'   - "cor": Propr matrix based on the simple Pearson correlation coefficient.
#'   - "vlr": Propr matrix based on the variance of log-ratio (VLR).
#'   - "pcor": Propr matrix based on the partial correlation coefficient
#'    (using ppcor package).
#'   - "pcor.shrink": Propr matrix based on the shrinkage-estimated partial
#'    correlation coefficient (using corpcor package).
#'   - "pcor.bshrink": Propr matrix based on the partial correlation
#'    coefficient with basis shrinkage (ivar argument must be 'clr' or 'alr').
#' @param ivar An indicator specifying the method for log-ratio transformation.
#'  It can take the following values:
#'   - "clr" (default): Centered log-ratio transformation.
#'   - "alr": Additive log-ratio transformation ("pcor.bshrink" metric only).
#'   - "iqlr": Inter-quartile log-ratio transformation from ALDEx2.
#'   - The explicit name(s) or index(es) of variable(s) to use as a reference.
#'   - Use NA to skip log-ratio transformation and any other pre-processing, like
#'  zero replacement. This is useful when the input data is already pre-processed.
#' @param select A numeric vector representing the indices of features to be
#'  used for computing the Propr matrix. This argument is optional. If
#'  provided, it reduces the data size by using only the selected features.
#' @param symmetrize A logical value indicating whether to force symmetry in
#'  the output Propr matrix when the metric is "phi". If `TRUE`, the function
#'  will symmetrize the matrix; otherwise, it will return the original matrix.
#' @param alpha The alpha parameter used in the alpha log-ratio transformation.
#' @param p The number of permutations to perform for calculating the false
#'  discovery rate (FDR). The default is 0.
#' @param permutation_option A character string indicating if permute the data
#' sample-wise or feature-wise. Default is "feature-wise"
#'
#' @return A propr object containing the Propr matrix, associated log-ratio
#'  transformation, and other calculated statistics.
#'
#' @details The function performs log-ratio transformation and computes a
#'  Propr matrix using different measures of association.
#'
#' @examples
#' # Sample input count data
#' data <- matrix(c(10, 5, 15, 20, 30, 25), nrow = 2, byrow = TRUE)
#'
#' # Calculate Propr matrix using correlation coefficient
#' result_cor <- propr(data, metric = "cor", ivar = "clr")
#'
#' # Calculate Propr matrix using variance of log-ratio (VLR)
#' result_vlr <- propr(data, metric = "vlr", ivar = "clr")
#'
#' # Calculate Propr matrix using partial correlation coefficient
#' result_pcor <- propr(data, metric = "pcor", ivar = "clr")
#'
#' @rdname propr
#' @export
propr <- function(counts,
                  metric = c("rho",
                             "phi",
                             "phs",
                             "cor",
                             "vlr",
                             "ppcor",
                             "pcor",
                             "pcor.shrink",
                             "pcor.bshrink"),
                  ivar = "clr",
                  select = NA,
                  symmetrize = FALSE,
                  alpha = NA,
                  p = 0,
                  permutation_option = c("feature-wise", "sample-wise")) {
  ##############################################################################
  ### CLEAN UP ARGS
  ##############################################################################

  # Special handling for equivalent args
  if (identical(alpha, 0))
    alpha <- NA

  # Special handling for 'metric'
  metric <- metric[1]
  ivar_pcor <- NA
  if (metric == "pcor.bshrink") {
    if (length(ivar) != 1) {
      stop("For 'pcor.bshrink', the 'ivar' argument must be a single value.")
    }
    if (!ivar %in% c("clr", "alr")) {
      stop("For 'pcor.bshrink', the 'ivar' argument must be 'clr' or 'alr'.")
    }
    message("Alert: Log-ratio transform will be handled by 'bShrink'.")
    ivar_pcor <- ivar
    ivar <- NA # skips log-ratio transform
  } else {
    if (length(ivar) == 1 && !is.na(ivar) && ivar == "alr") {
      stop("Please give the index or name of the reference gene instead of 'alr'.")
    } else if (length(ivar) > 1 && any(c("alr","clr","iqlr") %in% ivar)) {
      stop("Please check the ivar argument is correct.")
    }
  }

  # handle permutation option
  permutation_option <- permutation_option[1]
  if (!permutation_option %in% c("feature-wise", "sample-wise"))
    stop("permutation_option must be 'feature-wise' or 'sample-wise'")

  ##############################################################################
  ### PERFORM ZERO REPLACEMENT AND LOG-RATIO TRANSFORM
  ##############################################################################

  # NOTE: counts are the original counts, while ct may have zeros replaced
  counts <- as_safe_matrix(counts)
  if (length(ivar) == 1 && is.na(ivar) && is.na(ivar_pcor)) {
    ct <- counts
  } else {
    ct <- simple_zero_replacement(counts)
  }
  lr <- logratio(ct, ivar, alpha)

  ##############################################################################
  ### OPTIONALLY REDUCE DATA SIZE BEFORE COMPUTING MATRIX
  ##############################################################################

  if (!is.na(select)) {
    message("Alert: Using 'select' may make permutation testing unreliable.")
    counts <- counts[, select]
    ct <- ct[, select]
    lr <- lr[, select]
  }

  ##############################################################################
  ### COMPUTE THE ASSOCIATION MATRIX TO RETURN
  ##############################################################################

  lrv <- lr2vlr(lr)
  lambda <- NULL

  if (metric == "rho") {
    mat <- lr2rho(lr)

  } else if (metric == "phi") {
    mat <- lr2phi(lr)
    if (symmetrize)
      symRcpp(mat) # optionally force symmetry

  } else if (metric == "phs") {
    mat <- lr2phs(lr)

  } else if (metric == "cor") {
    mat <- stats::cor(lr)

  } else if (metric == "vlr") {
    mat <- lrv

  } else if (metric == "ppcor") {
    packageCheck("ppcor")
    mat <- ppcor::pcor(lr)$estimate

  } else if (metric == "pcor") {
    packageCheck("corpcor")
    cov <- cov(lr)
    mat <- corpcor::cor2pcor(cov)
    class(mat) <- "matrix"

  } else if (metric == "pcor.shrink") {
    packageCheck("corpcor")
    cov <- corpcor::cov.shrink(lr)
    mat <- corpcor::cor2pcor(cov)
    mat <- matrix(mat, ncol=ncol(lr), nrow=ncol(lr))
    class(mat) <- "matrix"
    lambda <- attr(cov, "lambda")

  } else if (metric == "pcor.bshrink") {
    with(pcor.bshrink(ct, outtype = ivar_pcor), {
      mat <<- matrix
      lambda <<- lambda
    })

  } else {
    stop("Provided 'metric' not recognized.")
  }

  ##############################################################################
  ### BUILD propr OBJECT TO RETURN TO USER
  ##############################################################################

  # Build propr object
  result <- new("propr")
  result@counts <- as.data.frame(ct)
  result@alpha <- as.numeric(alpha)
  result@metric <- metric[1]
  result@ivar <- ivar
  result@logratio <- as.data.frame(lr)
  result@pairs <- vector("numeric")
  result@permutes <- list(NULL)
  result@lambda <- lambda
  result@direct <- ifelse(metric[1] %in% c("rho", "cor", "pcor", "pcor.shrink", "pcor.bshrink"), TRUE, FALSE)
  result@has_meaningful_negative_values <- ifelse(metric[1] %in% c("cor", "pcor", "pcor.shrink", "pcor.bshrink"), TRUE, FALSE)  # metrics like proportionality has negative values that are difficult to interpret, whereas correlation metrics have a clear interpretation
  result@permutation_option <- permutation_option

  # ivar should not be NA for pcor.bshrink, otherwise updateCutoffs does not work
  if (metric == 'pcor.bshrink') result@ivar <- ivar_pcor

  # Clean row and column names
  result@matrix <- mat
  colnames(result@matrix) <- colnames(result@logratio)
  rownames(result@matrix) <- colnames(result@logratio)

  # Set up @results
  labels <- labRcpp(ncol(lr))
  result@results <-
    data.frame(
      "Partner" = labels[[1]],
      "Pair" = labels[[2]],
      "lrv" = lltRcpp(lrv),
      "metric" = metric,
      "alpha" = alpha,
      "propr" = lltRcpp(mat),
      "Zeros" = ctzRcpp(counts)
    )

  # permute data
  if (p > 0) result <- updatePermutes(result, p, permutation_option)

  ##############################################################################
  ### GIVE HELPFUL MESSAGES TO USER
  ##############################################################################

  message("Alert: Use 'updateCutoffs' to calculate FDR.")

  return(result)
}