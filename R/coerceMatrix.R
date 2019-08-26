#' Coerce Matrix Objects to Base Matrix Objects
#'
#' @description Coerces an object from classes in the \code{Matrix} package to
#'  the base matrix class.
#'
#' @param data The data, usually expected to be a \code{data.frame} or
#'  \code{matrix}, to be coerced to a \code{matrix} object if formatted as
#'  \code{dgCMatrix} or \code{dgeMatrix}.
#'
#' @importFrom methods is
#'
#' @keywords internal
#'
#' @return A coerced matrix object.
#'
coerceMatrix <- function(data) {
  if (is(data, "dgCMatrix") || is(data, "dgeMatrix")) {
    data <- as.matrix(data)
  }
  return(data)
}
