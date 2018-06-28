#' Import \code{ALDEx2} Object
#'
#' This method constructs a \code{propr} object from an
#'  \code{aldex.clr} object. See Details.
#'
#' The \code{ALDEx2} package has two exceptional features useful
#'  in proportionality analysis too. First, \code{ALDEx2} offers
#'  a number of extra log-ratio transformations, toggled
#'  by the \code{denom} argument in \code{aldex.clr}. Second,
#'  \code{ALDEx2} estimates per-feature technical variation
#'  within each sample using Monte-Carlo instances drawn
#'  from the Dirichlet distribution.
#'
#' The \code{aldex2propr} function takes advantage of both
#'  of these features by constructing a \code{propr} object
#'  directly from an \code{aldex.clr} object. When interpreting
#'  the resultant \code{propr} object, keep in mind that
#'  \code{ALDEx2} adds 0.5 to all \code{@@counts} regardless
#'  of whether the counts contain any zeros. Otherwise,
#'  the \code{@@logratio} slot contains the log-ratio
#'  transformed counts as averaged across all Monte Carlo
#'  instances. Likewise, the \code{@@matrix} slot gets
#'  filled with the proportionality matrix as averaged
#'  across all Monte Carlo instances.
#'
#' The \code{select} argument subsets the feature matrix
#'  after log-ratio transformation but before calculating
#'  proportionality. This reduces the run-time and RAM
#'  overhead without impacting the final result. Removing
#'  lowly abundant features prior to log-ratio transformation
#'  could otherwise change the proportionality measure.
#'
#' @param aldex.clr An \code{aldex.clr} object.
#' @param how A character string. The proportionality metric
#'  used to build the \code{propr} object. Choose from
#'  "rho", "phi", or "phs".
#' @param select Optional. Use this to subset the final
#'  proportionality matrix without altering the result.
#'
#' @return Returns a \code{propr} object.
#'
#' @export
aldex2propr <- function(aldex.clr, how = "rho", select){

  packageCheck("ALDEx2")

  if(class(aldex.clr) != "aldex.clr"){

    stop("This method expects an 'aldex.clr' object.")
  }

  if(how %in% c("perb", "rho", "lr2rho")){

    how <- "lr2rho"

  }else if(how %in% c("phit", "phi", "lr2phi")){

    how <- "lr2phi"

  }else if(how %in% c("phis", "phis", "phs", "lr2phs")){

    how <- "lr2phs"

  }else{

    stop("Provided 'how' not recognized.")
  }

  # Keep a running sum of propr instances
  counts <- t(as.matrix(aldex.clr@reads))
  mc <- ALDEx2::getMonteCarloInstances(aldex.clr)
  k <- ALDEx2::numMCInstances(aldex.clr)
  logratio <- 0
  prop <- 0
  for(i in 1:k){

    numTicks <- progress(i, k, numTicks)

    # Extract i-th Monte Carlo instance
    mci_lr <- t(sapply(mc, function(x) x[, i]))

    # Subset log-ratio transformed data
    if(!missing(select)){

      if(i == 1){

        # Make select boolean (it's OK if it's integer)
        if(is.character(select)) select <- match(select, colnames(mci_lr))
        if(any(is.na(select))) stop("Uh oh! Provided select reference not found in data.")
        counts <- counts[, select]
      }

      mci_lr <- mci_lr[, select]
    }

    # Add i-th log-ratio transformation to cumulative sum
    logratio <- logratio + mci_lr

    # Add i-th proportionality matrix to cumulative sum
    prop.i <- do.call(how, list("lr" = mci_lr))
    prop <- prop + prop.i
  }

  propr <- new("propr")
  propr@counts <- as.data.frame(counts)
  propr@logratio <- as.data.frame(logratio) / k
  propr@matrix <- prop / k

  message("Alert: Using 'aldex2propr' is not compatible the @results table.")
  propr@results <- data.frame()

  message("Alert: Using 'aldex2propr' disables permutation testing.")
  propr@permutes <- list(NULL)

  return(propr)
}

#' Correlate CLR Data with a Continuous Measurement
#'
#' This function correlates each log-ratio transformed
#'  feature vector with a continuous numeric variable.
#'
#' @param lr A data.frame. A log-ratio transformed counts matrix.
#' @param conditions A numeric vector of a continuous variable.
#' @param ... Arguments passed to \code{cor.test}.
#'
#' @return Returns a data.frame of the estimated
#'  correlation statistic (e.g., \code{r}) and p-value
#'  (\code{p}) for each log-ratio transformed feature,
#'  with FDR appended as the \code{BH} column.
#'
#' @importFrom stats cor.test p.adjust
#' @export
lr2cor <- function(lr, conditions, ...){

  if(!is.numeric(conditions)){

    stop("Please provide a numeric 'conditions' argument.")
  }

  if(nrow(lr) != length(conditions)){

    stop("Incorrect length for 'conditions' argument.")
  }

  cors <- apply(lr, 2, function(x){
    cor.test(x, conditions, ...)
  })

  r <- sapply(cors, getElement, "statistic")
  p <- sapply(cors, getElement, "p.value")
  BH <- p.adjust(p, method = "BH")

  data.frame(r, p, BH,
             row.names = colnames(lr))
}

#' Correlate CLR Data with a Continuous Measurement
#'
#' This function uses the Monte Carlo instances from an
#'  \code{aldex.clr} object to correlate each log-ratio
#'  transformed feature vector with a continuous
#'  numeric variable. See \code{\link{lr2cor}}.
#'
#' @param clr An \code{aldex.clr} object.
#' @param conditions A numeric vector of a continuous variable.
#' @param ... Arguments passed to \code{cor.test}.
#'
#' @return Returns a data.frame of the average
#'  correlation statistic (e.g., \code{r}) and p-value
#'  (\code{p}) for each log-ratio transformed feature,
#'  with FDR appended as the \code{BH} column.
#'
#' @export
aldex.cor <- function(clr, conditions, ...){

  packageCheck("ALDEx2")

  # Keep a running sum of lr2cor instances
  mc <- ALDEx2::getMonteCarloInstances(clr)
  k <- ALDEx2::numMCInstances(clr)
  r <- 0
  for(i in 1:k){

    numTicks <- progress(i, k, numTicks)
    mci_lr <- t(sapply(mc, function(x) x[, i]))
    r <- r + lr2cor(mci_lr, conditions, ...)
  }

  r / k # return expected
}

#' Build a Generalized Linear Model from CLR Data
#'
#' This function builds a generalized linear model from
#'  log-ratio transformed data.
#'
#' @param lr A data.frame. A log-ratio transformed counts matrix.
#' @param model.matrix A \code{model.matrix}.
#' @param ... Arguments passed to \code{glm}.
#'
#' @return Returns a data.frame of the estimated
#'  coefficients and their p-values for each feature,
#'  with FDR appended as a \code{BH} column.
#'
#' @importFrom stats glm p.adjust coef
#' @export
lr2glm <- function(lr, model.matrix, ...){

  if(class(model.matrix) != "matrix" &
     !("assign" %in% names(attributes(model.matrix)))){

    stop("Please use the 'model.matrix' function to prepare 'model.matrix'.")
  }

  if(nrow(lr) != nrow(model.matrix)){

    stop("Input data and 'model.matrix' should have same number of rows.")
  }

  # Build the glm models
  model. <- model.matrix
  glms <- apply(lr, 2, function(x){
    glm(x ~ model., ...)
  })

  # Extract coefficients and p-values
  extract <- function(model){
    x <- coef(summary(model))
    coefs <- lapply(1:nrow(x), function(i){
      y <- x[i,,drop=FALSE]
      colnames(y) <- paste(rownames(y), colnames(y))
      y})
    do.call("cbind", coefs)
  }

  # Combine to make data.frame
  extracts <- lapply(glms, extract)
  df <- do.call("rbind", extracts)
  rownames(df) <- colnames(lr)
  df <- as.data.frame(df)

  # Create new data.frame for FDR
  pvals <- colnames(df)[grepl("Pr\\(>", colnames(df))]
  df.bh <- df[,pvals]
  colnames(df.bh) <- paste0(colnames(df.bh), ".BH")
  for(j in 1:ncol(df.bh)){
    df.bh[,j] <- p.adjust(df.bh[,j])
  }

  # Merge results with FDR
  cbind(df, df.bh)
}

#' Build a Generalized Linear Model from CLR Data
#'
#' This function uses the Monte Carlo instances from an
#'  \code{aldex.clr} object to build a generalized linear
#'  model from log-ratio transformed data.
#'  See \code{\link{lr2glm}}.
#'
#' @param clr An \code{aldex.clr} object.
#' @param model.matrix A \code{model.matrix}.
#' @param ... Arguments passed to \code{glm}.
#'
#' @return Returns a data.frame of the average
#'  coefficients and their p-values for each feature,
#'  with FDR appended as a \code{BH} column.
#'
#' @export
aldex.glm <- function(clr, model.matrix, ...){

  packageCheck("ALDEx2")

  # Keep a running sum of lr2glm instances
  mc <- ALDEx2::getMonteCarloInstances(clr)
  k <- ALDEx2::numMCInstances(clr)
  r <- 0
  for(i in 1:k){

    numTicks <- progress(i, k, numTicks)
    mci_lr <- t(sapply(mc, function(x) x[, i]))
    r <- r + lr2glm(mci_lr, model.matrix, ...)
  }

  r / k # return expected
}
