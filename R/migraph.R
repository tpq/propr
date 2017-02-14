#' \code{igraph} Helper Functions
#'
#' @param g An \code{igraph} object.
#' @param names1,names2 A character vector. The \code{names1}
#'  argument defines a first set of vertices. The \code{names2}
#'  argument defines a second set of vertices to which the
#'  first set connects (i.e., element-wise), thereby defining
#'  a set of edges.
#' @param col A character string. The color applied to all
#'  vertices (or edges) specified by the \code{names1} (or
#'  \code{names2}) argument.
#' @param force A boolean. If true, the function adds any
#'  missing vertcies before adding edges. If false, the
#'  function only adds edges that have both vertices
#'  already present.
#'
#' @return An \code{igraph} object.
#' @name migraph

#' @rdname migraph
#' @export
migraph.add <- function(g, names1, names2, force = TRUE){

  if(!requireNamespace("igraph", quietly = TRUE)){
    stop("Uh oh! This method depends on igraph. ",
         "Try running: install.packages('igraph')")
  }

  if(missing(names2)) names2 <- NULL

  if(force | is.null(names2)){

    # Add any missing vertices before adding edges
    all <- union(names1, names2)
    new <- all[!all %in% igraph::V(g)$name]
    if(length(new) > 0){
      g <- igraph::add.vertices(g, length(new), "name" = new)
    }

  }else{

    # Only add edges in which both vertices appear on graph
    keep <- names1 %in% igraph::V(g)$name & names2 %in% igraph::V(g)$name
    names1 <- names1[keep]
    names2 <- names2[keep]
    if(length(names1) == 0){
      stop("No new edges to add.")
    }
  }

  if(!is.null(names2)){

    if(length(names1) != length(names2)){

      stop("Argument 'names1' and 'names2' should have same length.")
    }

    edges <- unlist(
      lapply(1:length(names1),
             function(i) c(names1[i], names2[i])))

    g <- igraph::add.edges(g, edges, color = "black")
  }

  g <- igraph::as.undirected(g)
  g <- igraph::simplify(g)

  return(g)
}

#' @rdname migraph
#' @export
migraph.color <- function(g, names1, names2, col){

  if(!requireNamespace("igraph", quietly = TRUE)){
    stop("Uh oh! This method depends on igraph. ",
         "Try running: install.packages('igraph')")
  }

  if(is.null(igraph::V(g)$color)) igraph::V(g)$color <- "white"
  if(is.null(igraph::E(g)$color)) igraph::E(g)$color <- "black"

  if(missing(names2)){

    igraph::V(g)$color <- ifelse(igraph::V(g)$name %in% names1,
                                 col, igraph::V(g)$color)

  }else{

    if(length(names1) != length(names2)){

      stop("Argument 'names1' and 'names2' should have same length.")
    }

    names <- unlist(
      lapply(1:length(names1),
             function(i) paste(names1[i], names2[i], sep = "-")))
    also <- unlist(
      lapply(1:length(names1),
             function(i) paste(names2[i], names1[i], sep = "-")))

    edges <- apply(igraph::get.edgelist(g), 1, paste, collapse = "-")
    igraph::E(g)$color <- ifelse(edges %in% names | edges %in% also,
                                 col, igraph::E(g)$color)
  }

  return(g)
}

#' @rdname migraph
#' @export
migraph.clean <- function(g){

  if(!requireNamespace("igraph", quietly = TRUE)){
    stop("Uh oh! This method depends on igraph. ",
         "Try running: install.packages('igraph')")
  }

  igraph::V(g)$size <- 2
  igraph::V(g)$label <- NA

  return(g)
}
