#' Contrastive Covariance Matrices
#'
#' @description Compute the list of contrastive covariance matrices.
#'
#' @param target The target data set
#' @param background The bacground data set
#' @param contrasts The vector of contrastive parameters
#' @param center Whether the data sets' columns should be centered
#' @param scale Whetehr the data sets' columns should be scaled to have variance
#'   1
#'
#' @author Philippe Boileau, \email{philippe_Boileau@@berkeley.edu}
#'
#' @return A list of contrastive covariance matrices. Each element has an
#'   associated contrastive parameter in the \code{contrasts} vector.
contrastiveCov <- function(target, background, contrasts, center, scale) {

  # get the covariance matrices of the target and background
  c_target <- covMat(target, center = center, scale = scale)
  c_background <- covMat(background, center = center, scale = scale)

  # get the list of contrastive covariance matrices
  c_contrasts <- lapply(contrasts, function(x) {
    c_target - x * c_background
  })

  return(c_contrasts)
}
