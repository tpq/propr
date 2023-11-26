updatePermutes <-
  function(object, p) {
    if (inherits(object, "propr")) {
      updatePermutes.propr(object, p)

    } else if (inherits(object, "propd")) {
      updatePermutes.propd(object, p)

    } else{
      stop("Provided 'object' not recognized.")
    }
  }

updatePermutes.propr <-
  function(object, p = 100) {
    # Shuffle row-wise so that per-sample CLR does not change
    message("Alert: Fixing permutations to active random seed.")
    ct <- object@counts
    permutes <- vector("list", p)
    for (ins in 1:length(permutes))
      permutes[[ins]] <- t(apply(ct, 1, sample))
    object@permutes <- permutes
    return(object)
  }

updatePermutes.propd <-
  function(object, p = 100) {
    message("Alert: Fixing permutations to active random seed.")
    ct <- object@counts
    permutes <- as.data.frame(matrix(0, nrow = nrow(ct), ncol = p))
    for (col in 1:ncol(permutes))
      permutes[, col] <- sample(1:nrow(ct))
    object@permutes <- permutes
    return(object)
  }
