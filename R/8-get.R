#' Get Results from Object
#'
#' This function provides a unified wrapper to retrieve results
#'  from a \code{propr} or \code{propd} object.
#'
#' @param object A \code{propr} or \code{propd} object.
#' @param cutoff This argument indicates the value at which to
#'  cutoff the results. For "rho" and "cor", the function
#'  returns pairs with a value greater than the cutoff.
#'  For "theta", "phi", and "phs", the function returns pairs
#'  with a value less than the cutoff. Leave the argument as
#'  \code{NA} to return all results.
#' @param include This argument indicates which features by
#'  name should belong to a pair for that pair to get included
#'  in the results. Subset performed by
#'  \code{Partner \%in\% subset | Pair \%in\% subset}.
#' @param or A boolean. If \code{FALSE}, \code{include} subsets
#'  by \code{Partner \%in\% subset & Pair \%in\% subset}.
#'
#' @return A \code{data.frame} of results.
#'
#' @export
getResults <- function(object, cutoff = NA, include = NA, or = TRUE){

  # Unify @results slot subset procedure
  if(class(object) == "propr"){

    if(object@metric == "rho" | object@metric == "cor"){

      outcome <- "propr"
      keep <- object@results[,outcome] >= cutoff

    }else if(object@metric == "phi" | object@metric == "phs"){

      outcome <- "propr"
      keep <- object@results[,outcome] <= cutoff

    }else{

      stop("Provided 'propr' metric not recognized.")
    }

  }else if(class(object) == "propd"){

    # Handle cutoff > 1 for backwards compatibility
    if(!is.na(cutoff)){
      if(cutoff > nrow(object@results)) cutoff <- nrow(object@results)
      if(cutoff > 1) cutoff <- sort(object@results$theta)[cutoff]
    }

    outcome <- "theta"
    keep <- object@results[,outcome] <= cutoff

  }else{

    stop("Provided 'object' not recognized.")
  }

  if(!is.na(cutoff)){

    # Apply numeric cutoff
    df <- object@results[keep,]

  }else{

    # Apply NA cutoff
    df <- object@results
  }

  # Name features of the new data.frame
  if(nrow(df) == 0) stop("No results remain after cutoff.")
  names <- colnames(object@counts)
  df$Partner <- names[df$Partner]
  df$Pair <- names[df$Pair]

  # Order results by outcome
  df <- df[order(df[,outcome]),]

  # Subset by 'include'
  if(!any(is.na(include))){

    if(class(include) != "character") stop("Provide 'include' as character.")

    if(or){
      index <- df$Partner %in% include | df$Pair %in% include
    }else{
      index <- df$Partner %in% include & df$Pair %in% include
    }

    df <- df[index,]
  }

  return(df)
}

#' Get Network from Object
#'
#' This function provides a unified wrapper to build networks
#'  from \code{propr} and \code{propd} objects.
#'
#' @param object Any \code{propr} or \code{propd} object.
#' @param cutoff A cutoff argument for \code{object}, passed
#'  to \code{\link{getResults}}.
#' @param propr.object,thetad.object,thetae.object A \code{propr}
#'  object or an appropriate \code{propd} object.
#' @param propr.cutoff,thetad.cutoff,thetae.cutoff A cutoff
#'  argument passed to \code{\link{getResults}}.
#' @inheritParams all
#'
#' @return A network object.
#'
#' @export
getNetwork <- function(object, cutoff = NA, propr.object, propr.cutoff = NA,
                       thetad.object, thetad.cutoff = NA,
                       thetae.object, thetae.cutoff = NA,
                       col1, col2, include = NA, or = TRUE, d3 = FALSE){

  if(!missing(object)){
    if(class(object) == "propr"){
      message("Alert: Treating 'object' as the proportionality network.")
      propr.object <- object
      propr.cutoff <- cutoff
    }else if(class(object) == "propd" & object@active == "theta_d"){
      message("Alert: Treating 'object' as the disjointed proportionality network.")
      thetad.object <- object
      thetad.cutoff <- cutoff
    }else if(class(object) == "propd" & object@active == "theta_e"){
      message("Alert: Treating 'object' as the emergent proportionality network.")
      thetae.object <- object
      thetae.cutoff <- cutoff
    }else{
      stop("Provide a valid object to the 'object' argument.")
    }
  }

  g <- igraph::make_empty_graph(directed = FALSE)

  # Add propr nodes to network
  if(!missing(propr.object)){

    if(class(propr.object) != "propr") stop("Provide a valid object to the 'propr.object' argument.")
    propr.df <- getResults(propr.object, propr.cutoff, include = include, or = or)
    g <- migraph.add(g, propr.df$Partner, propr.df$Pair)
  }

  # Add propd nodes to network
  if(!missing(thetad.object)){

    if(class(thetad.object) != "propd") stop("Provide a valid object to the 'thetad.object' argument.")
    if(thetad.object@active != "theta_d") stop("Provide a valid object to the 'thetad.object' argument.")
    thetad.group <- unique(thetad.object@group)
    thetad.df <- getResults(thetad.object, thetad.cutoff, include = include, or = or)
    g <- migraph.add(g, thetad.df$Partner, thetad.df$Pair)
  }

  # Add propd nodes to network
  if(!missing(thetae.object)){

    if(class(thetae.object) != "propd") stop("Provide a valid object to the 'thetae.object' argument.")
    if(thetae.object@active != "theta_e") stop("Provide a valid object to the 'thetae.object' argument.")
    thetae.group <- unique(thetae.object@group)
    thetae.df <- getResults(thetae.object, thetae.cutoff, include = include, or = or)
    g <- migraph.add(g, thetae.df$Partner, thetae.df$Pair)
  }

  # Add propr edges to network
  if(!missing(propr.object)){

    g <- migraph.color(g, propr.df$Partner, propr.df$Pair, "forestgreen")
    message("Green: Pair positively proportional across all samples.")

    invProp <- propr.df$propr < 0
    if(any(invProp)){

      g <- migraph.color(g, propr.df$Partner[invProp], propr.df$Pair[invProp], "burlywood4")
      message("Brown: Pair inversely proportional across all samples.")
    }
  }

  # Add propd edges to network
  if(!missing(thetad.object)){

    g <- migraph.color(g, thetad.df[thetad.df$lrm1 > thetad.df$lrm2, "Partner"],
                       thetad.df[thetad.df$lrm1 > thetad.df$lrm2, "Pair"], "coral1") # red
    g <- migraph.color(g, thetad.df[thetad.df$lrm1 < thetad.df$lrm2, "Partner"],
                       thetad.df[thetad.df$lrm1 < thetad.df$lrm2, "Pair"], "coral1") # blue
    message("Red: Pair has a different LRM in group ", thetad.group[1],
            " than in group ", thetad.group[2])
    # message("Blue: Pair has higher LRM in group ", thetad.group[2],
    #         " than in group ", thetad.group[1])
  }

  # Add propd edges to network
  if(!missing(thetae.object)){

    g <- migraph.color(g, thetae.df[thetae.df$lrv1 < thetae.df$lrv2, "Partner"],
                       thetae.df[thetae.df$lrv1 < thetae.df$lrv2, "Pair"], "gold2") # gold
    g <- migraph.color(g, thetae.df[thetae.df$lrv1 > thetae.df$lrv2, "Partner"],
                       thetae.df[thetae.df$lrv1 > thetae.df$lrv2, "Pair"], "blueviolet") # purple
    message("Gold: Nearly all of total LRV explained by ", thetae.group[2])
    message("Purple: Nearly all of total LRV explained by ", thetae.group[1])
  }

  # Finalize network
  if(!missing(col1)) g <- migraph.color(g, col1, col = "darkred")
  if(!missing(col2)) g <- migraph.color(g, col2, col = "darkslateblue")
  g <- migraph.clean(g)

  # Plot network
  if(d3){
    packageCheck("rgl")
    coords <- igraph::layout_with_fr(g, dim = 3)
    suppressWarnings(igraph::rglplot(g, layout = coords))
  }else{
    plot(g)
  }

  return(g)
}

#' Get (Log-)ratios from Object
#'
#' This function provides a unified wrapper to retrieve (log-)ratios
#'  from \code{propr} and \code{propd} objects.
#'
#' When the \code{alpha} argument is provided, this function returns
#'  the (log-)ratios as \code{(partner^alpha - pair^alpha) / alpha}.
#'
#' @inheritParams getResults
#' @param melt A boolean. Toggles whether to melt the results for
#'  visualization with \code{ggplot2}.
#'
#' @return A \code{data.frame} of (log-)ratios.
#'
#' @export
getRatios <- function(object, cutoff = NA, include = NA, or = TRUE, melt = FALSE){

  # Get results based on cutoff
  df <- getResults(object, cutoff, include = include, or = or)

  index <- colnames(object@counts) %in% union(df$Partner, df$Pair)
  ct <- object@counts[, index]
  alpha <- object@alpha

  # Get (log-)ratios [based on alpha]
  lr <- ratios(ct, alpha)

  # Subset after calling ratios() [to get only pairs in results]
  lr <- lr[, paste0(df$Partner, "/", df$Pair)]

  # For propd objects, define ratio so Group 1 is at top
  if(class(object) == "propd"){

    switchRatio <- function(x){
      text <- unlist(strsplit(x, "/"))
      paste0(text[2], "/", text[1])
    }

    for(i in 1:ncol(lr)){
      if(df$lrm2[i] > df$lrm1[i]){
        lr[,i] <- -1 * lr[,i]
        colnames(lr)[i] <- switchRatio(colnames(lr)[i])
      }
    }
  }

  # Melt data if appropriate
  if(melt){
    return(wide2long(lr))
  }else{
    return(lr)
  }
}

#' Get Reference to Approximate CLR
#'
#' This function finds the reference that is most
#'  proportional to the geometric mean of the samples.
#'  Using this reference for an alr-transformation
#'  will permit an analysis that resembles that of a
#'  clr-transformation, but arguably with greater
#'  interpretability.
#'
#' @inheritParams all
#'
#' @return A reference as a character or integer.
#'
#' @export
getReference <- function(counts, alpha = NA){

  counts <- as.matrix(counts)

  # Transform data into log space
  if(is.na(alpha)){

    if(any(counts == 0)){

      message("Alert: Replacing 0s with next smallest value.")
      zeros <- counts == 0
      counts[zeros] <- min(counts[!zeros])
    }

    logX <- log(counts)

  }else{

    logX <- (counts^alpha - 1)/alpha
  }

  # Sweep out g(x) center to get log(x/ref)
  ref <- rowMeans(logX)
  lr <- sweep(logX, 1, ref, "-")

  # Calculate var of each component
  vars <- apply(lr, 2, var)
  if(!is.null(colnames(counts))){
    colnames(counts)[which.min(vars)]
  }else{
    which.min(vars)
  }
}

#' Get Object as Adjacency Matrix
#'
#' This function uses \code{getNetwork} to
#'  return a \code{propr} or \code{propd}
#'  object as an adjacency matrix.
#'  The diagonal is set to 1.
#'
#' @inheritParams getResults
#'
#' @return An adjacency matrix.
#'
#' @export
getAdjacency <- function(object, cutoff = NA, include = NA, or = TRUE){

  g <- getNetwork(object, cutoff, include = include, or = or)
  a <- as.matrix(igraph::as_adjacency_matrix(g))
  diag(a) <- 1
  return(a)
}

#' Get Object as Symmetric Matrix
#'
#' This function returns a symmetric matrix
#'  of \code{propr} or \code{propd} values.
#'
#' @inheritParams getResults
#'
#' @return A symmetric matrix.
#'
#' @export
getMatrix <- function(object){

  if(class(object) == "propr"){

    outcome <- "propr"

  }else if(class(object) == "propd"){

    outcome <- "theta"

  }else{

    stop("Provided 'object' not recognized.")
  }

  if(object@results$Partner[1] != 2 |
     object@results$Pair[1] != 1){

    stop("Unexpected sorting of results slot.")
  }

  mat <- half2mat(object@results[,outcome])
  rownames(mat) <- colnames(object@counts)
  colnames(mat) <- colnames(object@counts)
  return(mat)
}
