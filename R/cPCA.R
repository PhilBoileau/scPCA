#' Constrastive Principal Component Analysis
#'
#' @description Given target and background dataframes or matrices, \code{cPCA}
#'   will perform the contrastive principal component analysis of the target
#'   data for a given number of eigenvectors and a vector of real valued
#'   contrast parameters. For more information on this method, see
#'   \href{https://www.nature.com/articles/s41467-018-04608-8#ref-CR29}{Abid et al.}.
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
#'   logarithmically spaced values between 0.1 and 1000.
#' @param num_medoids The number of medoids to select during the spectral
#'   clustering of the constrastive subspaces. Defaults to
#'   \code{round(length(contrasts)/5)} if there are at least 5 contrastive
#'   parameters, otherwise \code{length(contrasts)-1}. If there is a single
#'   contrastive parameter, then \code{num_medoids} should not be defined.
#'
#' @return A list containing the vector of contrast parameters that were
#'   selected as medoids by spectral clustering, the list of eigenvector
#'   matrices associated with each of these contrast parameters and the list
#'   of reduced-dimension target data associated with these contrast parameters.
#'
#' @export
#'
#' @importFrom kernlab specc as.kernelMatrix
#'
#' @author Philippe Boileau, \email{philippe_Boileau@@berkeley.edu}
#'
#' @examples
#' cPCA(
#'   target = toy_df[, 2:31],
#'   background = background_df
#' )
cPCA <- function(target, background, center = TRUE, scale = TRUE, num_eigen = 2,
                 contrasts = exp(seq(log(0.1), log(1000), length.out = 40)),
                 num_medoids) {

  # make sure that all parameters are input properly
  checkArgs(target, background, center, scale, num_eigen,
    contrasts,
    penalties = NULL, num_medoids
  )

  c_target <- covMat(target, center = center, scale = scale)
  c_background <- covMat(background, center = center, scale = scale)

  # determine the range of contrast parameters to use
  if (missing(contrasts)) {
    contrasts <- exp(seq(log(0.1), log(1000), length.out = 40))
  }

  # perform cPCA on the contrasted covariance matrices, get list of contrasts
  c_contrasts <- lapply(contrasts, function(x) {
    c_target - x * c_background
  })

  # set length of contrasts vector and number of medoids to consider.
  num_contrasts <- length(c_contrasts)
  if (missing(num_medoids) && num_contrasts >= 5) {
    num_medoids <- round(num_contrasts / 5)
  } else if (missing(num_medoids)) {
    num_medoids <- num_contrasts - 1
  }

  # for each contrasted covariance matrix, compute the eigenvectors
  loadings_mat <- lapply(
    1:num_contrasts,
    function(x) {
      eigen(c_contrasts[[x]],
        symmetric = TRUE
      )$vectors[, 1:num_eigen]
    }
  )

  # for each loadings matrix, project target onto constrastive subspace
  spaces <- lapply(
    1:num_contrasts,
    function(x) {
      as.matrix(target) %*% loadings_mat[[x]]
    }
  )


  # only perform this step if there is more than one contrastive parameter
  if (num_contrasts > 2) {
    # produce the QR decomposition of these projections, extract Q
    qr_decomps <- lapply(
      1:num_contrasts,
      function(x) {
        qr.Q(qr(spaces[[x]]))
      }
    )

    # populate affinity matrix for spectral clustering using the principal angles
    aff_vect <- sapply(
      1:(num_contrasts - 1),
      function(i) {
        sapply(
          (i + 1):num_contrasts,
          function(j) {
            Q_i <- qr_decomps[[i]]
            Q_j <- qr_decomps[[j]]
            d <- svd(x = t(Q_i) %*% Q_j, nu = 0, nv = 0)$d
            return(d[1] * d[2])
          }
        )
      }
    )
    aff_mat <- diag(x = 0.5, nrow = num_contrasts)
    aff_mat[lower.tri(aff_mat, diag = FALSE)] <- unlist(aff_vect)
    aff_mat <- t(aff_mat)
    # fix any computation errors, see numpy.nan_to_num
    aff_mat[is.nan(aff_mat)] <- 0
    aff_mat[is.na(aff_mat)] <- 0
    aff_mat[is.infinite(aff_mat)] <- 1000
    aff_mat <- t(aff_mat) + aff_mat

    # perfrom spectral clustering using the affinity matrix
    spec_clust <- kernlab::specc(kernlab::as.kernelMatrix(aff_mat),
      centers = num_medoids
    )

    # identify the alpha medoids of the spectral clustering
    contrast_medoids <- sapply(
      1:num_medoids,
      function(x) {
        sub_index <- which(spec_clust == x)
        sub_aff_mat <- as.matrix(
          aff_mat[sub_index, sub_index]
        )
        aff_sums <- colSums(sub_aff_mat)
        return(contrasts[sub_index[which.max(aff_sums)]])
      }
    )

    # create the lists of contrastive parameter medoids, loadings and projections
    contrast_medoids <- contrast_medoids[order(contrast_medoids)]
    med_index <- which(contrasts %in% contrast_medoids)
    med_loadings_mat <- loadings_mat[med_index]
    med_spaces <- spaces[med_index]
  } else {
    contrast_medoids <- contrasts
    med_loadings_mat <- loadings_mat
    med_spaces <- spaces
  }

  # return the alpha medoids with associated loadings and low-dim rep of target
  return(list(
    contrast_medoids,
    med_loadings_mat,
    med_spaces
  ))
}
