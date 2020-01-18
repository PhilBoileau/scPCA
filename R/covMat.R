#' Compute Sample Covariance Matrix
#'
#' @description \code{covMat} computes the sample covariance matrix of a data
#'   set. If a variable in the dataset has zero variance, then its
#'   corresponding row and column in the covariance matrix are zero vectors.
#'
#' @param data The data for which to compute the sample covariance matrix.
#' @param center A \code{logical} indicating whether the target and background
#'  data sets should be centered to mean zero.
#' @param scale A \code{logical} indicating whether the target and background
#'  data sets should be scaled to unit variance.
#'
#' @return the covariance matrix of the data.
#'
#' @importFrom stats cov var
#'
#' @keywords internal
covMat <- function(data, center = TRUE, scale = TRUE) {
  # convert data to a data.frame
  # NOTE: consider replacing with DataFrame
  data <- as.data.frame(data)

  # center the data matrix if required
  if (center) {
    data <- as.data.frame(safeColScale(data, center = TRUE, scale = FALSE))
  }

  # scale the matrix if required
  if (scale) {
    # scale the data and avoid columns with no variation by small perturbation
    data <- safeColScale(data, center = FALSE, scale = TRUE)
  }

  # compute the covariance matrix of the data
  cov_mat <- stats::cov(data)
  return(cov_mat)
}
