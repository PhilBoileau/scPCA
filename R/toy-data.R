#' Simulated Target Data for cPCA and scPCA
#'
#' The toy data consisting of 400 observations and 31 variables was simulated as
#' follows:
#' \itemize{
#'   \item Each of the first 10 variables was drawn from $N(0, 10)$
#'   \item For group 1 and 2, variables 11 through 20 were drawn from $N(0, 1)$
#'   \item For group 3 and 4, variables 11 through 20 were drawn from $N(3, 1)$
#'   \item For group 1 and 3, variables 21 though 30 were drawn from $N(-3, 1)$
#'   \item For group 2 and 4, variables 21 though 30 were drawn from $N(0, 1)$
#'   \item The last column provides each observations group number
#'}
#'
#' @docType data
#'
#' @usage data(toy_df)
#'
#' @format A simple \code{data.frame}.
#'
#' @keywords datasets
#'
#'
#' @examples
#' data(toy_df)
"toy_df"
