% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/norm_clust_strict_weighted.R
\name{norm_clust_strict_weighted}
\alias{norm_clust_strict_weighted}
\title{Normalize latent factors and clustering}
\usage{
norm_clust_strict_weighted(U, Knn = 20)
}
\arguments{
\item{U}{Latent factors U before normalization.}

\item{Knn}{The maximum number of nearest neighbors to search (default 20).}
}
\value{
Normalized U and clustering results.
}
\description{
Quantile normalization and clustering for the estimated latent factors U. Use weights in Louvain.
}
