% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/5-selectRatios.R
\name{search_tree}
\alias{search_tree}
\title{Search Tree Function}
\usage{
search_tree(data, Z, nclust = ncol(data)/10, nsearch = 1, lrm = NULL)
}
\arguments{
\item{data}{The input data matrix for clustering.}

\item{Z}{The matrix used to fit vegan model.}

\item{nclust}{The number of clusters to create during hierarchical clustering.
Default is calculated as ncol(data) / 10.}

\item{nsearch}{The number of best clusters to search for. Default is 1.}

\item{lrm}{The Log Ratio Matrix. Default is NULL.}
}
\value{
A numeric vector containing the percentage of variance explained
 by CCA for each cluster identified.
}
\description{
This function performs a hierarchical clustering on the given data and
 identifies the best clusters based on variance explained by
 Canonical Correspondence Analysis (CCA).
}
