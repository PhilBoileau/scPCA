#' Constrastive Principal Component Analysis
#'
#' @description Given target and background dataframes or matrices, \code{cPCA}
#'   will perform the contrastive principal component analysis of the target
#'   data for a given number of eigenvectors and a vector of real valued
#'   contrast parameters. For more information on this method, see
#'   \href{https://www.nature.com/articles/s41467-018-04608-8#ref-CR29}
#'   {Abid et al.}.
#'
#' @param target The target data. Either a numeric dataframe or a matrix with
#'   observations as rows and features as columns.
#' @param background The background data. Either a numeric dataframe or a matrix
#'   with observations as rows and features as columns. The number of features
#'   must match the number of features in the target data.
#' @param num_eigen The number of contrastive principal components to compute.
#'   Must be a non-negative integer between 1 and the number of columns in the
#'   target data. Default is 2.
#' @param contrasts The numeric vector of the contrastive parameters. Each
#'   element must be a unique non-negative real number. Defaults to 40
#'   logarithmically spaced values between 0.1 and 1000 when \code{start},
#'   \code{end} and \code{num_contrasts} parameters are \code{NULL}.
#'   Alternatively, see \code{start}, \code{end} and \code{num_contrasts}
#'   parameters.
#' @param start The starting value in the sequence of contrast parameters. Use
#'   \code{start}, \code{end} and \code{num_contrasts} parameters to create the
#'   sequence of contrast parameters instead of an explicit vector passed to
#'   \code{contrasts}. Defaults to \code{NULL}. Otherwise, must be an element of
#'   the positive reals.
#' @param end The ending value in the log sequence of contrast parameters. Use
#'   \code{start}, \code{end} and \code{num_contrasts} parameters to create the
#'   sequence of contrast parameters instead of an explicit vector passed to
#'   \code{contrasts}. Defaults to \code{NULL}. Otherwise, must be an element of
#'   the positive reals.
#' @param num_contrasts The number of elements in the log sequence of contrast
#'   parameters. Use \code{start}, \code{end} and \code{num_contrasts} parameters
#'   to create the log sequence of contrast parameters instead of an explicit
#'   vector passed to \code{contrasts}. Defaults to \code{NULL}. Otherswise, must
#'   be an integer larger than 1.
#' @param num_medoids The number of medoids to select during the spectral
#'   clustering of the constrastive subspaces.
#'
#' @return A list containing the vector of contrast parameters that were
#'   selected as medoids by spectral clustering, the list of eigenvector
#'   matrices associated with each of these contrast parameters and the list
#'   of reduced-dimension target data associated with these contrast parameters.
#'
#' @importFrom RSpectra eigs_sym
#' @importFrom kernlab specc as.kernelMatrix
#'
#' @author Philippe Boileau, \email{philippe_Boileau@@berkeley.edu}
#'
#' @examples
cPCA <- function(target, background, num_eigen = 2,
                 contrasts = exp(seq(log(0.1), log(1000), length.out = 40)),
                 start = NULL, end = NULL, num_contrasts = NULL){

  # make sure that all parameters are input properly
  if(!(class(target) %in% c("tbl_df", "tbl", "data.frame", "matrix"))){
    stop("Is your target data in the proper format? Check the documentation.")
  } else if(!(class(background) %in%
              c("tbl_df", "tbl", "data.frame", "matrix")) ||
            ncol(target) != ncol(background)){
    stop(paste("Is your background data in the proper format?",
               "Check the documentation."))
  } else if(class(num_eigen) != "numeric" || length(num_eigen) != 1 ||
            num_eigen < 1 || num_eigen > ncol(target) || num_eigen%%1 != 0){
    stop(paste("The num_eigen parameter must be a non-negative integer",
               "between 1 and the number of columns in the target data."))
  } else if(class(contrasts) != "numeric" || length(contrasts) < 1 ||
            contrasts >= 0){
    stop("The contrasts parameter must be a non-negative numeric vector.")
  } else if(!is.null(start) && start <= 0){
    stop("The start parameter must be NULL or a positive real value.")
  } else if(!is.null(end) && end <= start){
    stop(paste("The end parameter must be NULL or a positive real value larger",
                "than the start parameter."))
  } else if(!is.null(num_contrasts) &&
            (num_contrasts < 2 || num_contrasts%%1 != 0)){
    stop("The num_contrasts parameter must be NULL or an integer larger than 1.")
  }

  # center both dataframes
  target <- scale(target, center = TRUE, scale = FALSE)
  background <- scale(background, center = TRUE, scale = FALSE)

  # calculate the empirical covariance matrices, correct scalling factor
  len_target <- nrow(target)
  c_target <- (len_target-1)/len_target*var(target)
  len_background <- nrow(background)
  c_background <- (len_background-1)/len_background*var(background)

  # determine the range of contrast parameters to use
  if(!is.null(start) && !is.null(end) && !is.null(num_contrasts)){
    contrasts <- exp(seq(log(start), log(end), length.out = num_contrasts))
  }

  # perform cPCA on the contrasted covariance matrices, get list of contrasts
  c_contrasts <- lapply(contrasts, function(x){c_target - x*c_background})

  # length of contrasts vector
  len_con <- length(c_contrasts)

  # for each contrasted covariance matrix, compute the eigenvectors
  loadings_mat <- lappy(1:len_con,
                        function(x){
                          RSpectra::eigs_sym(A = c_contrasts[[x]],
                                             k = num_eigen)$vectors
                        })

  # populate affinity matrix for spectral clustering using the principal angles
  aff_vect <- sapply(1:len_con,
                     function(i){
                       sapply((i+1):len_con,
                              function(j){
                                V_i <- loadings_mat[[i]]
                                V_j <- loadings_mat[[j]]
                                d <- svd(x = t(V_i)%*%V_j, nu = 2, nv = 2)$d
                                return(d[1]*d[2])
                              })
                     })
  aff_mat <- diag(x = 0.5, nrow = len_con)
  aff_mat[lower.tri(aff_mat, diag = FALSE)] <- aff_vect
  aff_mat <- t(aff_mat)
  # fix any computation errors
  aff_mat[is.nan(aff_mat)] <- 0
  aff_mat[is.na(aff_mat)] <- 0
  aff_mat[is.infinite(aff_mat)] <- 1
  aff_mat <- t(aff_mat) + aff_mat

  # perfrom spectral clustering using the affinity matrix
  spec_clust <- kernlab::specc(kernlab::as.kernelMatrix(aff_mat),
                               centers = num_medoids)

  # identify the medoids of the spectral clustering


  # compute the contrastive principal components for the medoids


  # return the alpha medoids with associated loadings and low-dim rep of target
  return(cpc_target)

}
