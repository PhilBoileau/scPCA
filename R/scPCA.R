#' Sparse Constrastive Principal Component Analysis
#'
#' @description Given target and background dataframes or matrices, \code{scPCA}
#'   will perform the sparse contrastive principal component analysis of the
#'   target data for a given number of eigenvectors, a vector of real valued
#'   contrast parameters and a vector of penalty terms. For more information on
#'   the contrastice PCA method, which this method is an extension of, see
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
#' @param penalties The numeric vector of penatly terms for the L1 pernalty on
#'   the loadings. Defaults to 11 equidistant values between  0 and 0.5.
#' @param num_medoids The number of medoids to select during the spectral
#'   clustering of the constrastive subspaces. Defaults to
#'   \code{round(length(contrasts)/5)} if there are at least 5 contrastive
#'   parameters, otherwise \code{length(contrasts)-1}. If there is a single
#'   contrastive parameter, then \code{num_medoids} should not be defined.
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
#' scPCA(
#'   target = toy_df[, 2:31],
#'   background = background_df
#' )
scPCA <- function(target, background, center = TRUE, scale = TRUE,
                  num_eigen = 2,
                  contrasts = exp(seq(log(0.1), log(1000), length.out = 40)),
                  penalties = seq(0, 0.5, length.out = 11),
                  num_medoids) {

  # make sure that all parameters are input properly
  checkArgs(func = "scPCA", target, background, center, scale, num_eigen,
            contrasts, penalties, num_medoids)

  # get the contrastive covariance matrices
  c_contrasts <- contrastiveCov(target, background, contrasts, center, scale)

  # set length of contrasts and penalty vectors and number of medoids
  num_contrasts <- length(c_contrasts)
  num_penal <- length(penalties)
  if (missing(num_medoids) && num_contrasts >= 5) {
    num_medoids <- round(num_contrasts / 5)
  } else if (missing(num_medoids)) {
    num_medoids <- num_contrasts
  }

  # create the grid of contrast and penalty paramters
  param_grid <- expand.grid(penalties, contrasts)
  colnames(param_grid) <- c("lambda", "alpha")

  # for each contrasted covariance matrix, compute the eigenvectors using
  # the penalization term
  loadings_mat <- lapply(
    seq_len(num_contrasts),
    function(x) {
      lapply(
        penalties,
        function(y) {
          if (y == 0) {
            eigen(c_contrasts[[x]],
              symmetric = TRUE
            )$vectors[, 1:num_eigen]
          } else {
            elasticnet::spca(c_contrasts[[x]],
              K = num_eigen,
              para = rep(y, num_eigen),
              type = "Gram",
              sparse = "penalty"
            )$loadings
          }
        }
      )
    }
  )

  # unlist the nested list into a single list
  loadings_mat <- unlist(loadings_mat, recursive = FALSE)

  # for each loadings matrix, project target onto constrastive subspace
  spaces <- lapply(
    seq_len(num_contrasts * num_penal),
    function(x) {
      as.matrix(target) %*% loadings_mat[[x]]
    }
  )

  # remove all spaces projected to the 0 vector, update parameter grid, loadings
  param_grid <- param_grid[which(!duplicated(spaces)), ]
  loadings_mat <- loadings_mat[which(!duplicated(spaces))]
  spaces <- unique(spaces)

  # get the number of unique spaces
  num_spaces <- length(spaces)

  # check if spectral clustering is necessary
  if (num_spaces > 2) {

    # produce the QR decomposition of these projections, extract Q
    qr_decomps <- lapply(
      seq_len(num_spaces),
      function(x) {
        qr.Q(qr(spaces[[x]]))
      }
    )

    # populate affinity matrix for spectral clustering using the principal angles
    aff_vect <- sapply(
      seq_len(num_spaces - 1),
      function(i) {
        sapply(
          (i + 1):num_spaces,
          function(j) {
            Q_i <- qr_decomps[[i]]
            Q_j <- qr_decomps[[j]]
            d <- svd(x = t(Q_i) %*% Q_j, nu = 0, nv = 0)$d
            return(d[1] * d[2])
          }
        )
      }
    )
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
      iterations = 10000
    )

    # identify the alpha medoids of the spectral clustering
    contrast_medoids <- sapply(
      seq_len(num_medoids),
      function(x) {
        sub_index <- which(spec_clust == x)
        sub_aff_mat <- as.matrix(
          aff_mat[sub_index, sub_index]
        )
        aff_sums <- colSums(sub_aff_mat)
        return(
          param_grid[sub_index[which.max(aff_sums)], ]
        )
      }
    )

    # fix formating
    contrast_medoids <- matrix(unlist(contrast_medoids),
      nrow = num_medoids,
      byrow = TRUE
    )
    colnames(contrast_medoids) <- c("lambda", "alpha")

    # create the lists of contrastive parameter medoids, loadings and projections
    contrast_medoids <- contrast_medoids[order(
      contrast_medoids[, 2],
      contrast_medoids[, 1]
    ), ]

    # get the index of the paramaters chosen as medoids
    combo <- rbind(param_grid, contrast_medoids)
    rownames(combo) <- seq(1, nrow(combo))
    med_index <- as.numeric(
      rownames(combo[duplicated(combo, fromLast = TRUE), , drop = TRUE])
    )
    med_loadings_mat <- loadings_mat[med_index]
    med_spaces <- spaces[med_index]
  } else {
    contrast_medoids <- param_grid
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
