#' Contrastive Covariance Matrices
#'
#' @description Compute the list of contrastive covariance matrices.
#'
#' @param target The target (experimental) data set, in a standard format such
#'  as a \code{data.frame} or \code{matrix}.
#' @param background The background data set, in a standard format such as a
#'  \code{data.frame} or \code{matrix}.
#' @param contrasts A \code{numeric} vector of the contrastive parameters.
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
#' @return A list of contrastive covariance matrices. Each element has an
#'   associated contrastive parameter in the \code{contrasts} vector.
#'
#' @keywords internal
contrastiveCov <- function(
  target, background, contrasts, center, scale, scaled_matrix = FALSE
) {
  # get the covariance matrices of the target and background
  c_target <- covMat(
    target, center = center, scale = scale, scaled_matrix = scaled_matrix
  )
  c_background <- covMat(
    background, center = center, scale = scale, scaled_matrix = scaled_matrix
  )

  # get the list of contrastive covariance matrices
  c_contrasts <- lapply(contrasts, function(x) {
    c_target - x * c_background
  })

  # output
  return(c_contrasts)
}

################################################################################

#' Parallelized Contrastive Covariance Matrices
#'
#' @description Compute the list of contrastive covariance matrices in parallel
#'   using \code{\link[BiocParallel]{bplapply}}.
#'
#' @param target The target (experimental) data set, in a standard format such
#'  as a \code{data.frame} or \code{matrix}.
#' @param background The background data set, in a standard format such as a
#'  \code{data.frame} or \code{matrix}.
#' @param contrasts A \code{numeric} vector of the contrastive parameters.
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
#' @return A list of contrastive covariance matrices. Each element has an
#'   associated contrastive parameter in the \code{contrasts} vector.
#'
#' @importFrom BiocParallel bplapply
#'
#' @keywords internal
bpContrastiveCov <- function(
  target, background, contrasts, center, scale, scaled_matrix = FALSE
) {

  # get the covariance matrices of the target and background
  c_target <- covMat(
    target, center = center, scale = scale, scaled_matrix = scaled_matrix
  )
  c_background <- covMat(
    background, center = center, scale = scale, scaled_matrix = scaled_matrix
  )

  # get the list of contrastive covariance matrices
  c_contrasts <- BiocParallel::bplapply(contrasts, function(x) {
    c_target - x * c_background
  })

  # output
  return(c_contrasts)
}
