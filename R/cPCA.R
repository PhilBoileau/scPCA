#' Constrastive Principal Component Analysis
#'
#' @description Given target and background dataframes or matrices, \code{cPCA}
#'   will perform the contrastive principal component analysis of the target
#'   data for a given number of eigenvectors and a predetermined value
#'   for the contrast parameter. For more information on this method, see
#'   \href{https://www.nature.com/articles/s41467-018-04608-8#ref-CR29}
#'   {Abid et al.}.
#'
#' @param target The target data. Either a dataframe or a matrix with
#'   observations as rows and features as columns.
#' @param background The background data. Either a dataframe or a matrix with
#'   observations as rows and features as columns. The number of features
#'   must match the number of features in the target data.
#' @param num_eigen The number of contrastive principal components to compute.
#' @param contrast The value of the contrastive parameter.
#'
#' @return A list containing the lower-dimensional representation of the target
#'   data as a data frame and the column matrix of contrastive principal
#'   components.
#'
#' @importFrom RSpectra eigs_sym
#'
#' @author Philippe Boileau, \email{philippe_Boileau@@berkeley.edu}
#'
#' @examples
cPCA <- function(target, background, num_eigen, contrast){

  # make sure that target and background are matrics or dataframes,
  # num_eigen is an integer value between 1 and the number of columns in
  # target and that contrast is a non-negative number.
  if(!(class(target) %in% c("tbl_df", "tbl", "data.frame", "matrix"))){
    stop("Is your target data in the proper format? Check the documentation.")
  } else if(!(class(background) %in%
              c("tbl_df", "tbl", "data.frame", "matrix")) ||
            ncol(target) != ncol(background)){
    stop(paste("Is your background data in the proper format?",
               "Check the documentation."))
  } else if(class(num_eigen) != "numeric" || length(num_eigen) != 1 ||
            num_eigen < 1 || num_eigen > ncol(target)){
    stop(paste("The num_eigen parameter must be a non-negative integer",
               "between 1 and the number of columns in the target data."))
  } else if(class(contrast) != "numeric" || length(contrast) != 1 ||
            contrast >= 0){
    stop("The contrast parameter must be a non-negative numeric valule.")
  }

  # center both dataframes
  target <- scale(target, center = TRUE, scale = FALSE)
  background <- scale(background, center = TRUE, scale = FALSE)

  # calculate the empirical covariance matrices, correct scalling factor
  len_target <- nrow(target)
  c_target <- (len_target-1)/len_target*var(target)
  len_background <- nrow(background)
  c_background <- (len_background-1)/len_background*var(background)

  # perform cPCA on the contrasted covariance matrices
  c_contrast <- c_target - contrast*c_background
  loadings_mat <- RSpectra::eigs_sym(A = c_contrast, k = num_eigen)$vectors

  # compute the contrastive principal components
  cpc_target <- target %*% loadings_mat

  return(cpc_target)

}
