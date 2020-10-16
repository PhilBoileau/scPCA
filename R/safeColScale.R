#' Safe Centering and Scaling of Columns
#'
#' @description \code{safeColScale} is a safe utility for centering and scaling
#'  an input matrix \code{X}. It is intended to avoid the drawback of using
#'  \code{\link[base]{scale}} on data with constant variance by inducing adding
#'  a small perturbation to truncate the values in such columns. It also takes
#'  the opportunity to be faster than \code{\link[base]{scale}} through relying
#'  on \pkg{matrixStats} or \pkg{DelayedMatrixStats}, depending on the type of
#'  matrix being processed, for a key internal computation.
#'
#' @param X An input \code{matrix} to be centered and/or scaled. If \code{X} is
#'  not of class \code{matrix} or \code{DelayedMatrix}, then it must be
#'  coercible to a \code{matrix}.
#' @param center A \code{logical} indicating whether to re-center the columns
#'  of the input \code{X}.
#' @param scale A \code{logical} indicating whether to re-scale the columns of
#'  the input \code{X}.
#' @param tol A tolerance level for the lowest column variance (or standard
#'  deviation) value to be tolerated when scaling is desired. The default is
#'  set to \code{double.eps} of machine precision \code{.Machine}.
#' @param eps The desired lower bound of the estimated variance for a given
#'  column. When the lowest estimate falls below \code{tol}, it is truncated
#'  to the value specified in this argument. The default is 0.01.
#'
#' @return A centered and/or scaled version of the input data.
#'
#' @importFrom matrixStats colSds
#' @importFrom Matrix t colMeans
#' @importFrom assertthat assert_that
#'
#' @keywords internal
safeColScale <- function(X,
                         center = TRUE,
                         scale = TRUE,
                         tol = .Machine$double.eps,
                         eps = 0.01) {
  
  # check argument types
  assertthat::assert_that(is.logical(center))
  assertthat::assert_that(is.logical(scale))
  assertthat::assert_that(is.numeric(tol))
  assertthat::assert_that(is.numeric(eps))

  # input X must be a matrix for matrixStats
  if (!is.matrix(X) && class(X)[1] != "dgCMatrix" &&
      class(X)[1] != "DelayedMatrix") {
    X <- as.matrix(X)
  }

  # compute column means
  if (center) {
    colMeansX <- Matrix::colMeans(X, na.rm = TRUE)
  } else {
    # just subtract off zero if not centering
    colMeansX <- rep(0, ncol(X))
  }

  # compute scaling; replace by one if not scaling
  if (scale) {
    if (is.matrix(X)) {
      colSdsX <- matrixStats::colSds(X, na.rm = TRUE)
    } else if (class(X)[1] == "dgCMatrix") {
      colSdsX <- sparseMatrixStats::colSds(X, na.rm = TRUE)
    } else if (class(X)[1] == "DelayedMatrix") {
      colSdsX <- DelayedMatrixStats::colSds(X, na.rm = TRUE)
    }
    colSdsX[colSdsX < tol] <- eps
  } else {
    colSdsX <- rep(1, length(colMeansX))
  }

  # compute re-centered and re-scaled output
  # NOTE: there might be a  _faster_ way to do this?
  stdX <- t((t(X)-colMeansX)/colSdsX)

  # return output
  return(stdX)
}
