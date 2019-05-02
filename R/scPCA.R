#' Sparse Constrastive Principal Component Analysis
#'
#' @description Given target and background dataframes or matrices, \code{scPCA}
#'   will perform the sparse contrastive principal component analysis (scPCA) of
#'   the target data for a given number of eigenvectors, a vector of real valued
#'   contrast parameters and a vector of penalty terms. For more information on
#'   the contrastice PCA method, which this method is an extension of, consult
#'   \insertRef{abid2017contrastive}{scPCA}. Sparce PCA is performed using
#'   the method described in \insertRef{Zou2006}.
#'
#' @param target The target data. Either a numeric dataframe or a matrix with
#'   observations as rows and features as columns.
#' @param background The background data. Either a numeric dataframe or a matrix
#'   with observations as rows and features as columns. The number of features
#'   must match the number of features in the target data.
#' @param center Whether the target and background data should have their
#'   columns' centered. Defaults to \code{TRUE}.
#' @param scale Whether the target and background data should have their
#'   columns' scaled. Defaults to \code{FALSE}.
#' @param n_eigen The number of sparse contrastive components to compute.
#'   Default is 2.
#' @param contrasts The numeric vector of the contrastive parameters. Each
#'   element must be a unique non-negative real number. Defaults to 40
#'   logarithmically spaced values between 0.1 and 1000.
#' @param penalties The numeric vector of penatly terms for the L1 pernalty on
#'   the loadings. Defaults to 20 equidistant values between  0.05 and 1.
#' @param clust_method A \code{character} specifying the clustering method to
#'  use for choosing the optimal constrastive parameter. Currently, this is
#'  limited to either k-means or partitioning around medoids (PAM). The default
#'  is k-means.
#' @param n_centers The number of centers to use in the clustering algorithm.
#' @param max_iters The maximum number of iterations to use in k-means.
#'  clustering. Defaults to 10.
#'
#' @return A list containing the following components:
#'   \itemize{
#'     \item rotation - the matrix of variable loadings
#'     \item x - the rotated data, centred and scaled, if requested, data
#'     multiplied by the rotation matrix
#'     \item contrast - the optimal contrastive parameter used for cPCA
#'     \item penalty - the optimal L1 penalty term used for sparce PCA
#'     \item center - whether the target dataset was centered
#'     \item scale - whether the target dataset was scaled
#'   }
#'
#' @importFrom Rdpack reprompt
#'
#' @export
#'
#' @author Philippe Boileau, \email{philippe_Boileau@@berkeley.edu}
#'
#' @examples
#' scPCA(
#'   target = toy_df[, 1:30],
#'   background = background_df,
#'   contrasts = exp(seq(log(0.1), log(100), length.out = 10)),
#'   penalties = seq(0.1, 1, length.out = 9),
#'   n_centers = 4
#' )
#'
scPCA <- function(target, background, center = TRUE, scale = FALSE,
                  n_eigen = 2,
                  contrasts = exp(seq(log(0.1), log(1000), length.out = 40)),
                  penalties = seq(0.05, 1, length.out = 20),
                  clust_method = "kmeans", n_centers, max_iters = 10) {

  # make sure that all parameters are input properly
  checkArgs(target, background, center, scale, n_eigen,
            contrasts, penalties)

  # get the contrastive covariance matrices
  c_contrasts <- contrastiveCov(target, background, contrasts, center, scale)

  # find the optimal contrastive parameter and return its associated covariance
  # matrix, loading vector and rotation of the target data
  opt_params <- fitGrid(target, center, scale,
                        c_contrasts, contrasts, penalties, n_eigen,
                        clust_method = c("kmeans", "pam"), n_centers, max_iters)

  # create the list of results to output
  scpca <- list(
    rotation = opt_params$rotation,
    x = opt_params$x,
    contrast = opt_params$contrast,
    penalty = opt_params$penalty,
    center = center,
    scale = scale
  )

  return(scpca)
}
