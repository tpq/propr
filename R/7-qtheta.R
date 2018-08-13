#' Calculate a theta Cutoff
#'
#' This function uses the F distribution to calculate a cutoff of
#'  theta for a p-value given by the \code{pval} argument.
#'
#' @inheritParams all
#' @param pval A p-value at which to calculate a theta cutoff.
#'
#' @return A cutoff of theta from [0, 1].
#'
#' @export
qtheta <- function(propd, moderated = FALSE, pval = 0.05){

  if(pval < 0 | pval > 1) stop("Provide a p-value cutoff from [0, 1].")

  K <- length(unique(propd@group))
  N <- length(propd@group)

  if(moderated){

    propd <- suppressMessages(updateF(propd, moderated = TRUE))
    z.df <- propd@dfz

    Q <- qf(pval, K - 1, N + z.df - K, lower.tail = FALSE)
    # # Fstat <- (n1 + n2 + z.df - 2) * Fprime
    # # theta_mod <- 1 / (1 + Fprime)
    # # Q = Fstat
    # # Q = (n1 + n2 + z.df - 2) * Fprime
    # # Fprime = 1/theta_mod - 1
    R <- N - 2 + z.df
    # # Q = R * (1/theta_mod - 1)
    # # Q = R/theta_mod - R
    theta_a05 <- R/(Q+R)

  }else{

    Q <- qf(pval, K - 1, N - K, lower.tail = FALSE)
    # # Fstat <- (N - 2) * (1 - propd@theta$theta) / propd@theta$theta
    # # Q = Fstat
    # # Q = (N-2) * (1-theta) / theta
    # # Q / (N-2) = (1/theta) - 1
    # # 1/theta = Q / (N-2) + 1 = Q(N-2)/(N-2)
    # # theta = (N-2)/(Q+(N-2))
    theta_a05 <- (N-2)/(Q+(N-2))
  }

  return(theta_a05)
}
