% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/covMat.R
\name{covMat}
\alias{covMat}
\title{Compute Sample Covariance Matrix}
\usage{
covMat(data, center = TRUE, scale = TRUE, scaled_matrix = FALSE)
}
\arguments{
\item{data}{The data for which to compute the sample covariance matrix.}

\item{center}{A \code{logical} indicating whether the target and background
data sets should be centered to mean zero.}

\item{scale}{A \code{logical} indicating whether the target and background
data sets should be scaled to unit variance.}

\item{scaled_matrix}{A \code{logical} indicating whether to output a
\code{\link[ScaledMatrix]{ScaledMatrix}} object. The centering and scaling
procedure is delayed until later, permitting more efficient matrix
multiplication and row or column sums downstream. However, this comes at the
at the cost of numerical precision. Defaults to \code{FALSE}.}
}
\value{
the covariance matrix of the data.
}
\description{
\code{covMat} computes the sample covariance matrix of a data
  set. If a variable in the dataset has zero variance, then its
  corresponding row and column in the covariance matrix are zero vectors.
}
\keyword{internal}
