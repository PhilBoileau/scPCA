#' Coerce Matrix Objects to Base Matrix Objects
#'
#' @description Coerces an object from the Matrix class to a base matrix object.
#'
#' @param data The data to be coerced to a matrix object.
#'
#' @return returns a coerced matrix object.
coerceMatrix <- function(data) {

  if(class(data) %in% c("dgCMatrix", "dgeMatrix"))
    data <- as.matrix(data)
  else
    data
}
