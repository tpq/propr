% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/3-shared-updatePermutes.R
\name{updatePermutes}
\alias{updatePermutes}
\title{Create permuted data}
\usage{
updatePermutes(object, p = 100)
}
\arguments{
\item{object}{A \code{propr} or \code{propd} object.}

\item{p}{The number of permutations to perform. Default is 100.}
}
\value{
A \code{propr} or \code{propd} object with the permutes slot updated.
}
\description{
This function creates p permuted data matrices
}
\details{
This function wraps \code{updatePermutes.propr} and
 \code{updatePermutes.propd}.
}
