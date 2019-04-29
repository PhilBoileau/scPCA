#' Sparse Constrastive Principal Component Analysis
#'
#' @description Given target and background dataframes or matrices, \code{scPCA}
#'   will perform the sparse contrastive principal component analysis (PCA) of
#'   the target data for a given number of eigenvectors, a vector of real valued
#'   contrast parameters and a vector of penalty terms. For more information on
#'   the contrastice PCA method, which this method is an extension of, consult
#'   \insertRef{abid2017contrastive}{scPCA}. Sparce PCA is performed using
#'   \insertRed{WittenPMD2009}{scPCA}'s penalized matrix decomposition method.
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
#' @param num_eigen The number of contrastive principal components to compute.
#'   Must be a non-negative integer between 1 and the number of columns in the
#'   target data. Default is 2.
#' @param contrasts The numeric vector of the contrastive parameters. Each
#'   element must be a unique non-negative real number. Defaults to 40
#'   logarithmically spaced values between 0.1 and 1000.
#' @param penalties The numeric vector of penatly terms for the L1 pernalty on
#'   the loadings. Defaults to 11 equidistant values between  0 and 0.5.
#' @param clust_method A \code{character} specifying the clustering method to
#'  use for choosing the optimal constrastive parameter. Currently, this is
#'  limited to either k-means or partitioning around medoids (PAM). The default
#'  is k-means.
#' @param n_centers The number of centers to use in the clustering algorithm.
#' @param n_folds The number of folds in the cross-validation step used to
#'   find the optimal penalty term for sparse PCA. Defaults to 5.
#' @param sumabsvs_choice The choice of optimal L1 penalty term to use during
#'  the fitting of sparse PCA after performing the CV step. Either
#'  use the optimal penalty term (\code{"best"}) or the penalty inducing the
#'  most sparcity that is within 1 standard error of the smallest CV error
#'  (\code{"sparsest"}). Defaults to \code{"best"}.
#' @param n_iter The number of iterations of the SPC algorightm. Defaults to 20.
#' @param ... Additional arguments to pass to the clustering algorithm used to
#'   identify the optimal contrastive parameter
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
#' @importFrom PMA SPC SPC.cv
#'
#' @export
#'
#' @author Philippe Boileau, \email{philippe_Boileau@@berkeley.edu}
#'
#' @examples
#' scPCA(
#'   target = toy_df[, 1:30],
#'   background = background_df,
#'   n_centers = 4
#' )
#'
scPCA <- function(target, background, center = TRUE, scale = FALSE,
                  num_eigen = 2,
                  contrasts = exp(seq(log(0.1), log(1000), length.out = 40)),
                  penalties = seq(1.2, sqrt(ncol(target)), len=6),
                  clust_method = "kmeans", n_centers, n_folds = 5,
                  sumabsvs_choice = "best", n_iter = 20, ...) {

  # make sure that all parameters are input properly
  checkArgs(target, background, center, scale, num_eigen,
            contrasts, penalties)

  # get the contrastive covariance matrices
  c_contrasts <- contrastiveCov(target, background, contrasts, center, scale)

  # find the optimal contrastive parameter and return its associated covariance
  # matrix, loading vector and rotation of the target data
  opt_cont <- fitContrast(target, center, scale, c_contrasts,
                          contrasts, num_eigen, clust_method,
                          n_centers, ...)

  # find the optimal L1 penalty term based on the CV-MSE of the first loading
  v_init <- svd(opt_cont$c_cov, nu = 0, nv = 1)$v
  cv_out <- PMA::SPC.cv(x = opt_cont$c_cov, sumabsvs = penalties,
                        nfolds = n_folds, trace = FALSE, center = FALSE,
                        v = v_init)

  # determine which penalty value to use
  if(sumabsvs_choice == "best") {
    penalty <- cv_out$bestsumabsv
  } else {
    penalty <- cv_out$bestsumabsv1se
  }

  # perform sparce PCA using the selected penalty term
  scpca_out <- PMA::SPC(opt_cont$c_cov, sumabsv = penalty, K = num_eigen,
                        trace = FALSE, v = v_init, niter = n_iter,
                        center = FALSE, compute.pve = FALSE)

  # create the list of results to output
  scpca <- list(
    rotation = scpca_out$v,
    x = as.matrix(scale(target, center, scale)) %*% scpca_out$v,
    contrast = opt_cont$contrast,
    penalty = penalty,
    center = center,
    scale = scale
  )

  return(scpca)
}
