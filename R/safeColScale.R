#' Safe Centering and Scaling of Columns
#'
#' This is a safer utility for centering and scaling an input matrix \code{X}.
#' It is intended to avoid some of the drawbacks of \code{\link[base]{scale}}
#' while simultaneously being faster through relying on \pkg{matrixStats} for
#' key internal computations.
#'
#' @param X An input \code{matrix} to be centered and/or scaled.
#' @param center A \code{logical} indicating whether to re-center the columns
#'  of the input matrix \code{X}.
#' @param scale A \code{logical} indicating whether to re-scale the columns of
#'  the input matrix \code{X}.
#' @param tolVar A tolerance level for the lowest column variance (or standard
#'  deviation) value to be tolerated when scaling is desired. The default is
#'  set to \code{double.eps} of machine precision \code{\link[base]{.Machine}}.
#' @param epsVar The desired lower bound of the estimated variance for a given
#'  column. When the lowest estimate falls below \code{tolVar}, it is truncated
#'  to the value specified in this argument. The default is 0.01.
#'
#' @importFrom matrixStats colSds
#' @importFrom assertthat assert_that
#'
#' @keywords internal
safeColScale <- function(X,
                         center = TRUE,
                         scale = TRUE,
                         tol = .Machine$double.eps,
                         eps = 0.01) {
  # check argument types
  assertthat::assert_that(is.matrix(X))
  assertthat::assert_that(is.logical(center))
  assertthat::assert_that(is.logical(scale))
  assertthat::assert_that(is.numeric(tol))
  assertthat::assert_that(is.numeric(eps))

  # compute column means
  colMeansX <- colMeans(X, na.rm = TRUE)

  # just subtract off zero if not centering
  if (!center) {
    colMeansX <- rep(0, length(colMeansX))
  }
  # compute scaling; replace by one if not scaling
  if (scale) {
    colSdsX <- matrixStats::colSds(X)
    colSdsX[colSdsX < tol] <- eps
  } else {
    colSdsX <- rep(1, length(colMeansX))
  }

  # compute re-centered and re-scaled output
  # NOTE: there might be a  _faster_ way to do this?
  stdX <- t((t(X) - colMeansX) / colSdsX)

  # tweak attributes to exactly match output of base::scale
  if (center) {
    attr(stdX, "scaled:center") <- colMeansX
  }
  if (scale) {
    names(colSdsX) <- names(colMeansX)
    attr(stdX, "scaled:scale") <- colSdsX
  }

  # return output
  return(stdX)
}
