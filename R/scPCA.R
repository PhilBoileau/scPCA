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
#' @param cv A \code{numeric} indicating the number of cross-validation folds to
#'  use in choosing the optimal contrastive and penalization parameters from
#'  over the grids of \code{contrasts} and \code{penalties}. Cross-validation is
#'  expected to improve the robustness and generalization of the choice of these
#'  parameters; however, it increases the time the procedure costs, thus, the
#'  default is \code{NULL}, corresponding to no cross-validation.
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
#'  processing via the \pkg{BiocParallel} infrastructure. The default is
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
#' @importFrom origami cross_validate make_folds
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
                  n_eigen = 2, cv = NULL,
                  contrasts = exp(seq(log(0.1), log(1000), length.out = 40)),
                  penalties = seq(0.05, 1, length.out = 20),
                  clust_method = c("kmeans", "pam"), n_centers, max_iter = 10,
                  n_medoids = 8, parallel = FALSE) {
  # set defaults
  clust_method <- match.arg(clust_method)

  # check arguments to function
  checkArgs(
    target, background, center, scale, n_eigen,
    contrasts, penalties
  )

  # set target and background data sets to be matrices if from Matrix package
  target <- coerceMatrix(target)
  background <- coerceMatrix(background)

  if (is.null(cv)) {
    opt_params <- selectParams(
      target = target,
      background = background,
      center = center,
      scale = scale,
      n_eigen = n_eigen,
      contrasts = contrasts,
      penalties = penalties,
      clust_method = clust_method,
      n_centers = n_centers,
      max_iter = max_iter,
      n_medoids = n_medoids,
      parallel = parallel
    )
    if (n_centers > 1) {
      max_idx <- which.max(opt_params$ave_sil_widths)
      opt_params <- list(
        rotation = opt_params$rotation[[max_idx]],
        x = opt_params$x[[max_idx]],
        contrast = opt_params$contrast[max_idx],
        penalty = opt_params$penalty[max_idx]
      )
    }
  } else {
    # partition target and background data sets into CV-splits
    folds_target <- origami::make_folds(target, V = cv)
    folds_background <- origami::make_folds(background, V = cv)
    folds <- lapply(seq_len(cv), function(cv_fold) {
      out <- list(
        target = folds_target[[cv_fold]],
        background = folds_background[[cv_fold]]
      )
      return(out)
    })

    # select tuning parameters via cross-validation
    cv_opt_params <- origami::cross_validate(
      cv_fun = cvSelectParams,
      folds = folds,
      target = target,
      background = background,
      center = center,
      scale = scale,
      n_eigen = n_eigen,
      contrasts = contrasts,
      penalties = penalties,
      clust_method = clust_method,
      n_centers = n_centers,
      max_iter = max_iter,
      n_medoids = n_medoids,
      parallel = parallel,
      use_future = FALSE,
      .combine = FALSE
    )
    cv_sil_max <- do.call(c, (lapply(cv_opt_params$ave_sil_widths, max)))
    ave_sil_max <- mean(cv_sil_max)

    # re-fit on full data to get more stable estimates
    opt_params <- selectParams(
      target = target,
      background = background,
      center = center,
      scale = scale,
      n_eigen = n_eigen,
      contrasts = contrasts,
      penalties = penalties,
      clust_method = clust_method,
      n_centers = n_centers,
      max_iter = max_iter,
      n_medoids = n_medoids,
      parallel = parallel
    )
    if (n_centers > 1) {
      max_idx <- which.min(abs(opt_params$ave_sil_widths - ave_sil_max))
      opt_params <- list(
        rotation = opt_params$rotation[[max_idx]],
        x = opt_params$x[[max_idx]],
        contrast = opt_params$contrast[max_idx],
        penalty = opt_params$penalty[max_idx]
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
  class(scpca) <- "scpca"
  return(scpca)
}

################################################################################

#' Selection of Contrastive and Penalization Parameters
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
#'  processing via the \pkg{BiocParallel} infrastructure. The default is
#'  \code{FALSE} for sequential evaluation.
#'
#' @keywords internal
#'
selectParams <- function(target, background, center, scale, n_eigen,
                         contrasts, penalties, clust_method, n_centers,
                         max_iter, n_medoids, parallel) {
  # call parallelized function variants if so requested
  if (!parallel) {
    # create contrastive covariance matrices
    c_contrasts <- contrastiveCov(
      target = target, background = background, contrasts = contrasts,
      center = center, scale = scale
    )
    if (n_centers == 1) {
      opt_params <- fitCPCA(
        target = target, center = center, scale = scale,
        c_contrasts = c_contrasts, contrasts = contrasts,
        n_eigen = n_eigen, n_medoids = n_medoids
      )
    } else {
      opt_params <- fitGrid(
        target = target, center = center, scale = scale,
        c_contrasts = c_contrasts, contrasts = contrasts,
        penalties = penalties, n_eigen = n_eigen,
        clust_method = clust_method, n_centers = n_centers,
        max_iter = max_iter
      )
    }
  } else {
    # create contrastive covariance matrices
    c_contrasts <- bpContrastiveCov(
      target = target, background = background, contrasts = contrasts,
      center = center, scale = scale
    )
    if (n_centers == 1) {
      opt_params <- bpFitCPCA(
        target = target, center = center, scale = scale,
        c_contrasts = c_contrasts, contrasts = contrasts,
        n_eigen = n_eigen,
        n_medoids = n_medoids
      )
    } else {
      opt_params <- bpFitGrid(
        target = target, center = center, scale = scale,
        c_contrasts = c_contrasts, contrasts = contrasts,
        penalties = penalties, n_eigen = n_eigen,
        clust_method = clust_method, n_centers = n_centers,
        max_iter = max_iter
      )
    }
  }
  return(opt_params)
}

################################################################

#' Fold-Specific Selection of Contrastive and Penalization Parameters
#'
#' @param fold Object specifying cross-validation folds as generated by a call
#'  to \code{\link[origami]{make_folds}}.
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
#'  processing via the \pkg{BiocParallel} infrastructure. The default is
#'  \code{FALSE} for sequential evaluation.
#'
#' @importFrom origami training validation
#'
#' @keywords internal
#'
cvSelectParams <- function(fold, target, background, center, scale, n_eigen,
                           contrasts, penalties, clust_method, n_centers,
                           max_iter, n_medoids, parallel) {
  # make training and validation folds
  train_target <- origami::training(target, fold$target)
  valid_target <- origami::validation(target, fold$target)
  train_background <- origami::training(background, fold$background)

  # call parallelized function variants if so requested
  if (!parallel) {
    # create contrastive covariance matrices
    c_contrasts <- contrastiveCov(
      target = train_target, background = train_background,
      contrasts = contrasts, center = center, scale = scale
    )
    if (n_centers == 1) {
      opt_params <- fitCPCA(
        target = train_target, center = center,
        scale = scale, c_contrasts = c_contrasts,
        contrasts = contrasts, n_eigen = n_eigen,
        n_medoids = 8
      )
    } else {
      opt_params <- fitGrid(
        target = train_target,
        target_valid = valid_target,
        center = center, scale = scale,
        c_contrasts = c_contrasts,
        contrasts = contrasts,
        penalties = penalties,
        n_eigen = n_eigen,
        clust_method = clust_method,
        n_centers = n_centers,
        max_iter = max_iter
      )
    }
  } else {
    # create contrastive covariance matrices
    c_contrasts <- bpContrastiveCov(
      target = train_target, background = train_background,
      contrasts = contrasts, center = center, scale = scale
    )
    if (n_centers == 1) {
      opt_params <- bpFitCPCA(
        target = train_target, center = center,
        scale = scale, c_contrasts = c_contrasts,
        contrasts = contrasts, n_eigen = n_eigen,
        n_medoids = 8
      )
    } else {
      opt_params <- bpFitGrid(
        target = train_target, center = center,
        scale = scale, c_contrasts = c_contrasts,
        contrasts = contrasts, penalties = penalties,
        n_eigen = n_eigen,
        clust_method = clust_method,
        target_valid = valid_target,
        n_centers = n_centers,
        max_iter = max_iter
      )
    }
  }
  return(opt_params)
}
