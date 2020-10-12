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
#' @importFrom coop covar
#'
#' @keywords internal
covMat <- function(data, center = TRUE, scale = TRUE) {

  # center and scale the data matrix if required
  if (center || scale) {
    data <- safeColScale(data, center = center, scale = scale)
  }
  
  # compute the covariance matrix of the data
  if (class(data)[1] == "DelayedMatrix" || class(data)[1] == "dgCMatrix" ||
      class(data)[1] == "dgeMatrix") {
    cov_mat <- 1/(nrow(data)-1) * t(data) %*% data
    cov_mat <- as.matrix(cov_mat)
  } else {
    cov_mat <- coop::covar(data)
  }
  return(cov_mat)
}
