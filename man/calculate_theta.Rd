% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2a-propd-backend.R
\name{calculate_theta}
\alias{calculate_theta}
\title{Calculate Theta and Related Statistics}
\usage{
calculate_theta(
  counts,
  group,
  alpha = NA,
  lrv = NA,
  only = "all",
  weighted = FALSE,
  weights = as.matrix(NA),
  shrink = FALSE
)
}
\arguments{
\item{counts}{A data matrix representing counts.
It is assumed that the matrix contains numerical values only.}

\item{group}{A character vector representing group labels indicating the
assignment of each count to different groups.}

\item{alpha}{The alpha parameter used in the alpha log-ratio transformation.}

\item{lrv}{If LRV is provided, it is not computed within the function.}

\item{only}{A character vector specifying the type of theta to calculate.}

\item{weighted}{A logical value indicating whether weighted calculations
should be performed.}

\item{weights}{A weight matrix.}

\item{shrink}{A logical value indicating whether to apply shrinkage}
}
\value{
A data frame containing the computed theta values and
 related statistics, depending on the `only` parameter.
}
\description{
This function calculates theta and related statistics based on the input
 count matrix and other parameters. The function provides various options
 for calculating theta (theta_d, theta_e, theta_f, theta_g).
}
\examples{
# Sample input count data and group assignments
data <- iris[1:100, 1:4]
group <- iris[1:100, 5]

# Calculate all theta types
result_all <- calculate_theta(data, group, alpha = 0.5)

# Calculate only theta_d
result_theta_d <- calculate_theta(data, group, alpha = 0.5, only = "theta_d")

}
