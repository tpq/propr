#' Simplify Indexed Matrix
#'
#' This convenience function takes an indexed \code{\link{propr}} object
#'  and subsets the object based on that index. Then, it populates the
#'  \code{@@pairs} slot of the new object with an updated version
#'  of the original index. You can call \code{simplify} from within the
#'  \code{[} method using the argument \code{tiny}.
#'
#' @inheritParams propr
#'
#' @return Returns a \code{propr} object.
#'
#' @export
simplify <- function(object){

  if(!class(object) == "propr" | length(object@pairs) == 0){

    stop("Uh oh! You can only simplify an indexed 'propr' object.")
  }

  # Call indexToCoord on indexed propr object
  coords <- indexToCoord(object@pairs, nrow(object@matrix))
  selection <- sort(union(coords[[1]], coords[[2]]))
  object@pairs <- vector("numeric")

  # Subset propr object based on index
  new <- subset(object, select = selection)

  # Repopulate the pairs slot
  for(i in 1:length(coords[[1]])){

    coords[[1]][i] <- which(selection == coords[[1]][i])
    coords[[2]][i] <- which(selection == coords[[2]][i])
  }

  new@pairs <- (coords[[2]] - 1) * nrow(new@matrix) + (coords[[1]] - 1) + 1

  return(new)
}

#' Make Adjacency Object
#'
#' This function uses pairs indexed in the \code{@@pairs}
#'  slot to build a symmetric adjacency matrix.
#'
#' @inheritParams propr
#'
#' @return Returns a \code{propr} object with the adjacency
#'  matrix saved to the \code{@@matrix} slot.
#'
#' @export
adjacent <- function(object){

  if(!class(object) == "propr" | length(object@pairs) == 0){

    stop("Uh oh! This function requires an indexed 'propr' object.")
  }

  N <- nrow(object@matrix)
  mat <- matrix(0, N, N)
  mat[object@pairs] <- 1
  diag(mat) <- 1
  symRcpp(mat)

  adj <- object
  adj@matrix <- mat

  return(adj)
}

#' Export Indexed Pairs
#'
#' This function exports a \code{data.frame} of indexed pairs
#'  along with their corresponding values from \code{object@matrix}.
#'  In doing so, this function provides a preview of the
#'  interaction network, built using \code{igraph}.
#'  We recommend using this result as input to a standalone
#'  network visualization tool like Cytoscape.
#'
#' @inheritParams propr
#' @param minPairs An integer scalar. Subsets the interaction
#'  network to exclude any pair without a node that participates
#'  in at least this many total pairs.
#'
#' @return Returns a \code{data.frame} of indexed pairs.
#'
#' @export
cytescape <- function(object, minPairs = 2){

  packageCheck("igraph")

  if(!class(object) == "propr" | length(object@pairs) == 0){

    stop("Uh oh! This function requires an indexed 'propr' object.")
  }

  # Prepare data
  rho <- object@matrix[object@pairs]
  coords <- indexToCoord(object@pairs, nrow(object@matrix))
  df <- data.frame("Partner" = coords[[1]], "Pair" = coords[[2]], rho)

  # Remove extraneous pairs
  keep <- which(table(c(df$Partner, df$Pair)) >= minPairs)
  sub <- df[df$Partner %in% keep | df$Pair %in% keep, ]

  # Build and color igraph
  g <- igraph::graph_from_data_frame(sub, directed = FALSE)
  igraph::E(g)$color <- ifelse(sub$rho > 0, "red", "blue")
  plot(g, vertex.size = 2, vertex.label = NA)

  # Retrieve node names
  names <- colnames(object@logratio)
  if(!is.null(names)){
    sub$Partner <- names[sub$Partner]
    sub$Pair <- names[sub$Pair]
  }

  return(sub)
}
