% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2a-propd-backend.R
\name{lrv_with_shrinkage}
\alias{lrv_with_shrinkage}
\title{Calculate Logratio Variance with shrinkage}
\usage{
lrv_with_shrinkage(ct, shrink = TRUE)
}
\arguments{
\item{ct}{A count matrix.}

\item{shrink}{A logical value indicating whether to apply shrinkage.}
}
\value{
A shrunk logratio variance matrix.
}
\description{
This function computes the logratio variance (LRV) with the option
 to apply shrinkage. It uses the `corpcor` package to compute a shrunk
 covariance matrix and then converts it to a logratio variance matrix.
}
