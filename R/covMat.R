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
#' @importFrom Matrix crossprod
#' @importFrom tibble is_tibble
#'
#' @keywords internal
covMat <- function(data, center = TRUE, scale = TRUE) {

  # center and scale the data matrix if required
  if (scale) {
    data <- safeColScale(data, center = TRUE, scale = scale)
  } else {
    data <- safeColScale(data, center = TRUE, scale = FALSE)
  }
  
  # compute the covariance matrix of the data
  if (is.matrix(data) || is.data.frame(data) || is_tibble(data)) {
    cov_mat <- coop::covar(data)
  } else {
    cov_mat <- 1 / (nrow(data) - 1) * Matrix::crossprod(data)
    cov_mat <- as.matrix(cov_mat)
  }
  return(cov_mat)
}
