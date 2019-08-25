#' Contrastive Covariance Matrices
#'
#' @description Compute the list of contrastive covariance matrices.
#'
#' @param target The target data set.
#' @param background The background data set.
#' @param contrasts The vector of contrastive parameters.
#' @param center A \code{logical} indicating whether the data sets' columns
#'  should be centered so as to have mean zero.
#' @param scale A \code{logical} indicating whether the data sets' columns
#'  should be re-scaled to have unit variance.
#'
#' @author Philippe Boileau, \email{philippe_boileau@berkeley.edu}
#'
#' @return A list of contrastive covariance matrices. Each element has an
#'   associated contrastive parameter in the \code{contrasts} vector.
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
