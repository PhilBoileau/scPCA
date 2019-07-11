#' Compute Sample Covariance Matrix
#'
#' @description \code{covMat} computes the sample covariance matrix of a data
#'   set. If a variable in the dataset has zero variance, then its corresponding
#'   row and column in the covariance matrix are zero vectors.
#'
#' @param data The data with which to compute the sample covariance matrix.
#' @param center A \code{logical} indicating whether the data sets' columns
#'  should be centered so as to have mean zero.
#' @param scale A \code{logical} indicating whether the data sets' columns
#'  should be re-scaled to have unit variance.
#'
#' @return the covariance matrix of the data.
#'
#' @importFrom stats cov var
#'
#' @author Philippe Boileau, \email{philippe_boileau@@berkeley.edu}
#'
covMat <- function(data, center = TRUE, scale = TRUE) {

  # convert data to a data.frame
  # NOTE: consider replacing with DataFrame
  data <- as.data.frame(data)

  # center the data matrix if required
  if (center) {
    data <- as.data.frame(scale(data, center = TRUE, scale = FALSE))
  }

  # scale the matrix if required
  # NOTE: if there are constant columns, replace NAs by 0s
  if (scale) {
    # identify columns with zero variance
    no_var_idx <- which(sapply(data, stats::var) == 0)

    # scale the data and replace the columns with no variation by zero vectors
    data <- scale(data, center = FALSE, scale = TRUE)
    data[, no_var_idx] <- rep(0, nrow(data))
  }

  # compute the covariance matrix of the data
  cov_mat <- stats::cov(data)
  return(cov_mat)
}
