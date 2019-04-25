#' Constructor for class scePCA
#'
#' @return class \code{scePCA} object, sub-classed from SingleCellExperiment.
#'
#' @import BiocGenerics
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom methods setClass
#'
#' @export .scePCA
#' @exportClass scePCA
#'
#' @rdname scePCA-class
#'
#' @examples
#' library(SingleCellExperiment)
#' library(scPCA)
#'
#' example_scePCA_class <- function(sce) {
#'   call <- match.call(expand.dots = TRUE)
#'   scepca <- .scePCA(
#'     SingleCellExperiment(
#'       assays = assay(sce),
#'       rowData = rowData(sce),
#'       colData = colData(sce)
#'     ),
#'     call = call,
#'     tmleOut = as.data.frame(matrix(NA, 10, 10)),
#'     topTable = as.data.frame(matrix(NA, 10, 10))
#'   )
#'   return(scepca)
#' }
#'
#' example_class <- example_scepca_class(sce = EXAMPLEDATA)
.scePCA <- methods::setClass(
  Class = "scePCA",
  slots = list(
    call = "call",
    newslot1 = "DataFrame",
    newslot2 = "DataFrame"
  ),
  contains = "SingleCellExperiment"
)
