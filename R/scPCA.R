#' Sparse Constrastive Principal Component Analysis
#'
#' @description Given target and background dataframes or matrices, \code{scPCA}
#'  will perform the sparse contrastive principal component analysis (scPCA) of
#'  the target data for a given number of eigenvectors, a vector of real valued
#'  contrast parameters and a vector of penalty terms. For more information on
#'  the contrastive PCA method, consult \insertRef{abid2018exploring}{scPCA}.
#'  Sparse PCA is performed via the method of \insertRef{zou2006sparse}{scPCA}.
#'
#' @param target The target (experimental) data set, in a standard format such
#'  as a \code{data.frame} or \code{matrix}.
#' @param background The background data set, in a standard format such as a
#'  \code{data.frame} or \code{matrix}. Note that the number of features must
#'  match the number of features in the target data.
#' @param center A \code{logical} indicating whether the target and background
#'  data sets should be centered to mean zero.
#' @param scale A \code{logical} indicating whether the target and background
#'  data sets should be scaled to unit variance.
#' @param n_eigen A \code{numeric} indicating the number of eigenvectors (or
#'  sparse contrastive components) to be computed. The default is to compute two
#'  such eigenvectors.
#' @param contrasts A \code{numeric} vector of the contrastive parameters. Each
#'  element must be a unique non-negative real number. The default is to use 40
#'  logarithmically spaced values between 0.1 and 1000.
#' @param penalties A \code{numeric} vector of the L1 penalty terms on the
#'  loadings. The default is to use 20 equidistant values between 0.05 and 1.
#' @param clust_method A \code{character} specifying the clustering method to
#'  use for choosing the optimal constrastive parameter. Currently, this is
#'  limited to either k-means or partitioning around medoids (PAM). The default
#'  is k-means clustering.
#' @param n_centers A \code{numeric} giving the number of centers to use in the
#'  clustering algorithm. If set to 1, cPCA, as first proposed by Abid et al.,
#'  is performed, regardless of what the \code{penalties} argument is set to.
#' @param max_iter A \code{numeric} giving the maximum number of iterations to
#'   be used in k-means clustering, defaulting to 10.
#' @param n_medoids A \code{numeric} indicating the number of medoids to
#'  consider if \code{n_centers} is set to 1. The default is 8 such medoids.
#' @param parallel A \code{logical} indicating whether to invoke parallel
#'  processing via the \code{BiocParallel} infrastructure. The default is
#'  \code{FALSE} for sequential evaluation.
#'
#' @return A list containing the following components:
#'   \itemize{
#'     \item rotation - the matrix of variable loadings
#'     \item x - the rotated data, centred and scaled if requested, multiplied
#'       by the rotation matrix
#'     \item contrast - the optimal contrastive parameter
#'     \item penalty - the optimal L1 penalty term
#'     \item center - whether the target dataset was centered
#'     \item scale - whether the target dataset was scaled
#'   }
#'
#' @importFrom Rdpack reprompt
#'
#' @export
#'
#' @examples
#' # perform cPCA on the simulated data set
#' scPCA(
#'   target = toy_df[, 1:30],
#'   background = background_df,
#'   contrasts = exp(seq(log(0.1), log(100), length.out = 5)),
#'   penalties = 0,
#'   n_centers = 4
#' )
#'
#' # perform scPCA on the simulated data set
#' scPCA(
#'   target = toy_df[, 1:30],
#'   background = background_df,
#'   contrasts = exp(seq(log(0.1), log(100), length.out = 5)),
#'   penalties = seq(0.1, 1, length.out = 3),
#'   n_centers = 4
#' )
#'
#' # cPCA as implemented in Abid et al.
#' scPCA(
#'   target = toy_df[, 1:30],
#'   background = background_df,
#'   contrasts = exp(seq(log(0.1), log(100), length.out = 10)),
#'   penalties = 0,
#'   n_centers = 1
#' )
scPCA <- function(target, background, center = TRUE, scale = FALSE,
                  n_eigen = 2,
                  contrasts = exp(seq(log(0.1), log(1000), length.out = 40)),
                  penalties = seq(0.05, 1, length.out = 20),
                  clust_method = "kmeans", n_centers, max_iter = 10,
                  n_medoids = 8, parallel = FALSE) {

  # check arguments to function
  checkArgs(
    target, background, center, scale, n_eigen,
    contrasts, penalties
  )

  # set target and background data sets to be matrices if from Matrix package
  target <- coerceMatrix(target)
  background <- coerceMatrix(background)

  # call parallelized function variants if so requested
  if (parallel == FALSE) {
    # create contrastive covariance matrices
    c_contrasts <- contrastiveCov(
      target, background, contrasts, center,
      scale
    )
    if (n_centers == 1) {
      opt_params <- fitCPCA(target, center, scale, c_contrasts, contrasts,
        n_eigen,
        n_medoids = 8
      )
    } else {
      opt_params <- fitGrid(target, center, scale, c_contrasts, contrasts,
        penalties, n_eigen,
        clust_method = c("kmeans", "pam"),
        n_centers, max_iter
      )
    }
  } else {
    # create contrastive covariance matrices
    c_contrasts <- bpContrastiveCov(
      target, background, contrasts,
      center, scale
    )
    if (n_centers == 1) {
      opt_params <- bpFitCPCA(target, center, scale, c_contrasts, contrasts,
        n_eigen,
        n_medoids = 8
      )
    } else {
      opt_params <- bpFitGrid(target, center, scale, c_contrasts, contrasts,
        penalties, n_eigen,
        clust_method = c("kmeans", "pam"), n_centers,
        max_iter
      )
    }
  }

  # create output object
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
