#' Coerce Matrix Objects to Base Matrix Objects
#'
#' @description Coerces an object from classes in the \code{Matrix} package to
#'  the base matrix class.
#'
#' @param data The data to be coerced to a matrix object.
#'
#' @keywords internal
#'
#' @return A coerced matrix object.
#'
coerceMatrix <- function(data) {
  if (is(data, "dgCMatrix") || is(data, "dgeMatrix")) {
    data <- as.matrix(data)
  } else {
    data
  }
  return(data)
}
