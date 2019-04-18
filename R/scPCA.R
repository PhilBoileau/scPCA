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
  checkArgs(target, background, center, scale, num_eigen,
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

  # for each contrasted covariance matrix, compute components and projections
  c_proj <- projGridCP(target, center, scale, c_contrasts, contrasts,
                       penalties, num_eigen)

  num_spaces <- length(c_proj$spaces)

  # check if spectral clustering is necessary
  if (num_spaces > 2 && num_medoids > 1) {

    results <- specClustSelection(c_proj, num_medoids)

  } else {
    results <- list(
      medoids_params = c_proj$param_grid,
      med_loadings_mat = c_proj$loadings_mat,
      med_spaces = c_proj$spaces
    )
  }

  # return the alpha medoids with associated loadings and low-dim rep of target
  return(results)
}
