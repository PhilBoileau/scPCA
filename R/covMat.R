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
#' @param scaled_matrix A \code{logical} indicating whether to output a
#'  \code{\link[ScaledMatrix]{ScaledMatrix}} object. The centering and scaling
#'  procedure is delayed until later, permitting more efficient matrix
#'  multiplication and row or column sums downstream. However, this comes at the
#'  at the cost of numerical precision. Defaults to \code{FALSE}.
#'
#' @return the covariance matrix of the data.
#'
#' @importFrom coop covar
#' @importFrom Matrix crossprod
#' @importFrom tibble is_tibble
#'
#' @keywords internal
covMat <- function(data, center = TRUE, scale = TRUE, scaled_matrix = FALSE) {

  # center and scale the data matrix if required
  if (scale) {
    data <- safeColScale(data, center = TRUE, scale = TRUE,
                         scaled_matrix = scaled_matrix)
  } else {
    data <- safeColScale(data, center = TRUE, scale = FALSE,
                         scaled_matrix = scaled_matrix)
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
