#' @rdname propr
#' @export
setClass(
  "propr",
  slots = c(
    counts = "data.frame",
    alpha = "numeric",
    metric = "character",
    direct = "logical",
    has_meaningful_negative_values = "logical",
    permutation_option = "character",
    ivar = "ANY",
    lambda = "ANY",
    logratio = "data.frame",
    matrix = "matrix",
    pairs = "numeric",
    results = "data.frame",
    permutes = "list",
    fdr = "data.frame",
    tails = "character"
  )
)

#' @rdname propr
#' @section Methods (by generic):
#' \code{show:} Method to show \code{propr} object.
#' @export
setMethod("show", "propr",
          function(object) {
            cat(
              "Not weighted",
              "and",
              ifelse(
                is.na(object@alpha),
                "not alpha-transformed",
                "alpha-transformed"
              ),
              "\n"
            )

            cat(
              "@counts summary:",
              nrow(object@counts),
              "subjects by",
              ncol(object@counts),
              "features\n"
            )

            cat(
              "@logratio summary:",
              nrow(object@logratio),
              "subjects by",
              ncol(object@logratio),
              "features\n"
            )

            cat(
              "@matrix summary:",
              nrow(object@matrix),
              "features by",
              ncol(object@matrix),
              "features\n"
            )

            if (length(object@pairs) > 0 |
                nrow(object@matrix) == 0) {
              cat("@pairs summary:", length(object@pairs), "feature pairs\n")

            } else{
              cat("@pairs summary: index with `[` method\n")
            }

            cat("@fdr summary:",
                ncol(object@permutes), "iterations\n")

            if (nrow(object@fdr) > 0) {
              print(object@fdr)
            }

            cat("See ?propr for object methods\n")
          })
