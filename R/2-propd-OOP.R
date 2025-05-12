#' @rdname propd
#' @export
setClass(
  "propd",
  slots = c(
    counts = "data.frame",
    alpha = "numeric",
    group = "character",
    weighted = "logical",
    weights = "data.frame",
    active = "character",
    Fivar = "ANY",
    dfz = "numeric",
    results = "data.frame",
    permutes = "data.frame",
    fdr = "data.frame"
  )
)

#' @rdname propd
#' @section Methods (by generic):
#' \code{show:} Method to show \code{propd} object.
#' @export
setMethod("show", "propd",
          function(object) {
            cat(
              ifelse(object@weighted, "Weighted", "Not weighted"),
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
              "@group summary:",
              length(unique(object@group)),
              "unique groups (",
              paste(table(object@group), collapse = " x "),
              ")\n"
            )

            cat("@results summary:",
                nrow(object@results),
                "feature pairs (",
                object@active,
                ")\n")

            cat("@fdr summary:",
                ncol(object@permutes), "iterations\n")

            if (nrow(object@fdr) > 0) {
              print(object@fdr)
            }

            cat("See ?propd for object methods\n")
          })
