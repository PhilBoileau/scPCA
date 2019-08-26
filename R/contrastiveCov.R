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
#'
#' @return A list of contrastive covariance matrices. Each element has an
#'   associated contrastive parameter in the \code{contrasts} vector.
#'
#' @keywords internal
#'
contrastiveCov <- function(target, background, contrasts, center, scale) {
  # get the covariance matrices of the target and background
  c_target <- covMat(target, center = center, scale = scale)
  c_background <- covMat(background, center = center, scale = scale)

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
#'
#' @return A list of contrastive covariance matrices. Each element has an
#'   associated contrastive parameter in the \code{contrasts} vector.
#'
#' @importFrom BiocParallel bplapply
#'
#' @keywords internal
#'
bpContrastiveCov <- function(target, background, contrasts, center, scale) {
  # get the covariance matrices of the target and background
  c_target <- covMat(target, center = center, scale = scale)
  c_background <- covMat(background, center = center, scale = scale)

  # get the list of contrastive covariance matrices
  c_contrasts <- BiocParallel::bplapply(contrasts, function(x) {
    c_target - x * c_background
  })

  # output
  return(c_contrasts)
}
