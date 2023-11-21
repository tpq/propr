#' Example Absolute mRNA
#'
#' Data generated with supplemental script provided by
#'  <DOI:10.1371/journal.pcbi.1004075>. Data originally
#'  sourced from <DOI:10.1016/j.cell.2012.09.019>.
#'  A time series of yeast mRNA abundance after removal
#'  of a key nutrient. Absolute abundance estimated
#'  by multiplying microarray signal (relative to first
#'  time point) by the initial nCounter-calibrated and
#'  copy-per-cell-adjusted RNA-seq abundance (averaged
#'  across two replicates). Divide absolute abundances
#'  by total sample abundance to make data relative.
#'
#' @usage data(marg.abs)
"marg.abs"

#' Ensure Matrix Has Dim Names
#'
#' Makes sure input data has correct format. For back-end use only.
#'
#' @param counts A data matrix representing counts.
#' @return A matrix with dim names.
as_safe_matrix <-
  function(counts) {
    if ("data.frame" %in% class(counts))
      counts <- as.matrix(counts)
    if (is.null(colnames(counts)))
      colnames(counts) <- as.character(1:ncol(counts))
    if (is.null(rownames(counts)))
      rownames(counts) <- as.character(1:nrow(counts))
    if (any(is.na(counts)))
      stop("Remove NAs from 'counts' before proceeding.")
    return(counts)
  }

#' Make Progress Bar
#'
#' Makes a progress bar. For back-end use only.
#'
#' @param i The current iteration.
#' @param k Total iterations.
#' @param numTicks The result of \code{progress}.
#' @return The next \code{numTicks} argument.
progress <- function(i, k, numTicks) {
  if (i == 1)
    numTicks <- 0

  if (numTicks == 0)
    cat("|-")

  while (i > numTicks * (k / 40)) {
    cat("-")
    if (numTicks == 10)
      cat("(25%)")
    if (numTicks == 20)
      cat("(50%)")
    if (numTicks == 30)
      cat("(75%)")
    numTicks <- numTicks + 1
  }

  if (i == k)
    cat("-|\n")

  return(numTicks)
}

#' Package Check
#'
#' Checks whether the user has the required package installed.
#'  For back-end use only.
#'
#' @param package A character string. An R package.
#' @return Returns TRUE if no error.
packageCheck <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(
      "Uh oh! This method depends on ",
      package,
      "! ",
      "\nTry running: install.packages('",
      package,
      "')",
      "\nor: BiocManager::install('",
      package,
      "')"
    )
  }
  TRUE
}
