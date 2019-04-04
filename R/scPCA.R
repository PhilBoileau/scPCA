#' Sparse Constrastive Principal Component Analysis
#'
#' @description Given target and background dataframes or matrices, \code{scPCA}
#'   will perform the sparse contrastive principal component analysis of the
#'   target data for a given number of eigenvectors, a vector of real valued
#'   contrast parameters and a vector of penalty terms. For more information on
#'   the contrastice PCA method, which this method is an extension of, see
#'   \href{https://www.nature.com/articles/s41467-018-04608-8#ref-CR29}
#'   {Abid et al.}.
#'
#' @param target The target data. Either a numeric dataframe or a matrix with
#'   observations as rows and features as columns.
#' @param background The background data. Either a numeric dataframe or a matrix
#'   with observations as rows and features as columns. The number of features
#'   must match the number of features in the target data.
#' @param center Whether the target and background data should have their
#'   columns' centered. Defaults to \code{TRUE}.
#' @param scale Whether the target and background data should have their
#'   columns' scaled. Defaults to \code{TRUE}.
#' @param num_eigen The number of contrastive principal components to compute.
#'   Must be a non-negative integer between 1 and the number of columns in the
#'   target data. Default is 2.
#' @param contrasts The numeric vector of the contrastive parameters. Each
#'   element must be a unique non-negative real number. Defaults to 40
#'   logarithmically spaced values between 0.1 and 1000 when \code{start},
#'   \code{end} and \code{num_contrasts} parameters are \code{NULL}.
#'   Alternatively, see \code{start}, \code{end} and \code{num_contrasts}
#'   parameters.
#' @param penalties The numeric vector of penatly terms for the L1 pernalty on
#'   the loadings. Defaults to 11 equidistant values between  0 and 0.5.
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
#'   clustering of the constrastive subspaces. Defaults to
#'   \code{round(num_contrasts/5)} if there are at least 5 contrastive
#'   parameters, 1 otherwise.
#'
#' @return A list containing the vector of contrast and penalty parameter pairs
#'   that were selected as medoids by spectral clustering, the list of
#'   eigenvector matrices associated with each of these pairs and the list
#'   of reduced-dimension target data.
#'
#' @export
#'
#' @importFrom kernlab specc as.kernelMatrix
#' @importFrom elasticnet spca
#'
#' @author Philippe Boileau, \email{philippe_Boileau@@berkeley.edu}
#'
#' @examples
#' scPCA(target = toy_df[, 2:31],
#'      background = background_df)
scPCA <- function(target, background, center = TRUE, scale = TRUE,
                  num_eigen = 2, contrasts, penalties, start = NULL, end = NULL,
                  num_contrasts = NULL, num_medoids){

  # make sure that all parameters are input properly
  if(!(class(target) %in% c("tbl_df", "tbl", "spec_tbl_df",
                            "data.frame", "matrix"))){
    stop("Is your target data in the proper format? Check the documentation.")
  } else if(!(class(background) %in%
              c("tbl_df", "tbl", "spec_tbl_df", "data.frame", "matrix")) ||
            ncol(target) != ncol(background)){
    stop(paste("Is your background data in the proper format?",
               "Check the documentation."))
  } else if(class(num_eigen) != "numeric" || length(num_eigen) != 1 ||
            num_eigen < 1 || num_eigen > ncol(target) || num_eigen%%1 != 0){
    stop(paste("The num_eigen parameter must be a non-negative integer",
               "between 1 and the number of columns in the target data."))
  } else if(!missing(contrasts) &&
            (class(contrasts) != "numeric" || length(contrasts) < 1 ||
             contrasts <= 0)){
    stop("The contrasts parameter must be a non-negative numeric vector.")
  } else if(!is.null(start) && start <= 0){
    stop("The start parameter must be NULL or a positive real value.")
  } else if(!is.null(end) && end <= start){
    stop(paste("The end parameter must be NULL or a positive real value larger",
               "than the start parameter."))
  } else if(!is.null(num_contrasts) &&
            (num_contrasts < 2 || num_contrasts%%1 != 0)){
    stop("The num_contrasts parameter must be NULL or an integer larger than 1.")
  } else if(!missing(num_medoids) && num_medoids <= 0){
    stop(paste("The num_medoids parameter must be a positive integer that is",
               "smaller than the number of contrastive parameters."))
  } else if(center != TRUE && center != FALSE){
    stop("The center parameter should be set to TRUE or FALSE.")
  } else if(!missing(penalties) && (penalties <= 0)){
    stop("The penalties parameter must consist of positive numeric values.")
  }


  if(center && scale){
    target <- scale(target, center = TRUE, scale = TRUE)
    background <- scale(background, center = TRUE, scale = TRUE)
  } else if(center == TRUE && scale == FALSE){
    target <- scale(target, center = TRUE, scale = FALSE)
    background <- scale(background, center = TRUE, scale = FALSE)
  } else if(center == FALSE && scale == TRUE){
    target <- scale(target, center = FALSE, scale = TRUE)
    background <- scale(background, center = FALSE, scale = TRUE)
  }

  # calculate the empirical covariance matrices, correct scalling factor
  len_target <- nrow(target)
  c_target <- (len_target-1)/len_target*var(target)
  len_background <- nrow(background)
  c_background <- (len_background-1)/len_background*var(background)

  # determine the range of contrast parameters to use
  if(!is.null(start) && !is.null(end) && !is.null(num_contrasts)){
    contrasts <- exp(seq(log(start), log(end), length.out = num_contrasts))
  } else if(missing(contrasts)) {
    contrasts <- exp(seq(log(0.1), log(1000), length.out = 40))
  }

  # determine the range of penalty terms to use
  if(missing(penalties))
    penalties <- seq(0, 0.5, length.out = 11)

  # perform cPCA on the contrasted covariance matrices, get list of contrasts
  c_contrasts <- lapply(contrasts, function(x){c_target - x*c_background})

  # set length of contrasts and penalty vectors and number of medoids
  num_contrasts <- length(c_contrasts)
  num_penal <- length(penalties)
  if(missing(num_medoids) && num_contrasts >= 5)
    num_medoids <- round(num_contrasts / 5)
  else if(missing(num_medoids))
    num_medoids <- 1

  # create the grid of contrast and penalty paramters
  param_grid <- expand.grid(penalties, contrasts)
  colnames(param_grid) <- c("lambda", "alpha")

  # for each contrasted covariance matrix, compute the eigenvectors using
  # the penalization term
  loadings_mat <- lapply(1:num_contrasts,
                         function(x){
                           lapply(penalties,
                                  function(y){
                                    if(y == 0){
                                      eigen(c_contrasts[[x]],
                                        symmetric = TRUE)$vectors[, 1:num_eigen]
                                    } else {
                                      elasticnet::spca(c_contrasts[[x]],
                                                    K = num_eigen,
                                                    para = rep(y, num_eigen),
                                                    type = "Gram",
                                                    sparse = "penalty")$loadings
                                    }
                                  })
                         })

  # unlist the nested list into a single list
  loadings_mat <- unlist(loadings_mat, recursive = FALSE)

  # for each loadings matrix, project target onto constrastive subspace
  spaces <- lapply(1:(num_contrasts*num_penal),
                   function(x){
                     as.matrix(target) %*% loadings_mat[[x]]
                   })

  # remove all spaces projected to the 0 vector, update parameter grid, loadings
  param_grid <- param_grid[which(!duplicated(spaces)), ]
  loadings_mat <- loadings_mat[which(!duplicated(spaces))]
  spaces <- unique(spaces)

  # get the number of unique spaces
  num_spaces <- length(spaces)

  # produce the QR decomposition of these projections, extract Q
  qr_decomps <- lapply(1:num_spaces,
                       function(x){
                         qr.Q(qr(spaces[[x]]))
                       })

  # populate affinity matrix for spectral clustering using the principal angles
  aff_vect <- sapply(1:(num_spaces-1),
                     function(i){
                       sapply((i+1):num_spaces,
                              function(j){
                                Q_i <- qr_decomps[[i]]
                                Q_j <- qr_decomps[[j]]
                                d <- svd(x = t(Q_i)%*%Q_j, nu = 0, nv = 0)$d
                                return(d[1]*d[2])
                              })
                     })
  aff_mat <- diag(x = 0.5, nrow = num_spaces)
  aff_mat[lower.tri(aff_mat, diag = FALSE)] <- unlist(aff_vect)
  aff_mat <- t(aff_mat)
  # fix any computation errors, see numpy.nan_to_num
  aff_mat[is.nan(aff_mat)] <- 0
  aff_mat[is.na(aff_mat)] <- 0
  aff_mat[is.infinite(aff_mat)] <- 1
  aff_mat <- t(aff_mat) + aff_mat

  # perfrom spectral clustering using the affinity matrix
  spec_clust <- kernlab::specc(kernlab::as.kernelMatrix(aff_mat),
                               centers = num_medoids,
                               iterations = 10000)

  # identify the alpha medoids of the spectral clustering
  contrast_medoids <- sapply(1:num_medoids,
                             function(x){
                               sub_index <- which(spec_clust == x)
                               sub_aff_mat <- as.matrix(
                                 aff_mat[sub_index, sub_index])
                               aff_sums <- colSums(sub_aff_mat)
                               return(
                                 param_grid[sub_index[which.max(aff_sums)], ]
                               )
                             })

  # fix formating
  contrast_medoids <- matrix(unlist(contrast_medoids),
                             nrow = num_medoids,
                             byrow = TRUE)
  colnames(contrast_medoids) <- c("lambda", "alpha")

  # create the lists of contrastive parameter medoids, loadings and projections
  contrast_medoids <- contrast_medoids[order(contrast_medoids[, 2],
                                             contrast_medoids[, 1]), ]

  # get the index of the paramaters chosen as medoids
  combo <- rbind(param_grid, contrast_medoids)
  rownames(combo) <- seq(1, nrow(combo))
  med_index <- as.numeric(
    rownames(combo[duplicated(combo, fromLast = TRUE),, drop = TRUE])
  )
  med_loadings_mat <- loadings_mat[med_index]
  med_spaces <- spaces[med_index]

  # return the alpha medoids with associated loadings and low-dim rep of target
  return(list(contrast_medoids,
              med_loadings_mat,
              med_spaces))

}