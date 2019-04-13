#' Compute Sample Covariance Matrix
#'
#' @description \code{covMat} computes the sample covariance matrix of a data
#'   set. If a variable in the dataset has zero variance, then its corresponding
#'   row and column in the covariance matrix are zero vectors.
#'
#' @param data The data with which to compute the sample covariance matrix.
#' @param center Indicates if the data should be centered to have mean 0.
#'   Defaults to \code{TRUE}.
#' @param scale Indicates if the data should be scaled to have variance 1.
#'   Defaults to \code{TRUE}.
#'
#' @return the covariance matrix of the data.
#'
#' @importFrom stats cov
#'
#' @author Philippe Boileau, \email{philippe_Boileau@@berkeley.edu}
covMat <- function(data, center = TRUE, scale = TRUE){

  # convert data to a dataframe
  data <- as.data.frame(data)

  # center the data matrix if required
  if(center)
    data <- as.data.frame(scale(data, center = TRUE, scale = FALSE))

  # scale the matrix if required
  # if there are constant columns, replace NAs by 0s
  if(scale){

    # identify columns with zero variance
    no_var_idx <- which(sapply(data, var) == 0)

    # scale the data and replace the columns with no variation by zero vectors
    data <- scale(data, center = FALSE, scale = TRUE)
    data[, no_var_idx] <- rep(0, nrow(data))

  }

  # compute the covariance matrix of the data
  cov_mat <- stats::cov(data)

  return(cov_mat)
}
