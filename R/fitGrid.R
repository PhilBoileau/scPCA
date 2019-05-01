#' Identify the Optimal Contrastive Parameter
#'
#' @description This function is used to automatically select the optimal
#'   contrastive parameter based on k-means and average silhouette width.
#'
#' @param target The target data set.
#' @param center A \code{logical} indicating whether the data sets' columns
#'  should be centered so as to have mean zero.
#' @param scale A \code{logical} indicating whether the data sets' columns
#'  should be re-scaled to have unit variance.
#' @param c_contrasts List of of contrastive covariances.
#' @param contrasts Vector of contrastive parameter values used to compute the
#'   contrastive covariances,
#' @param penalties Vector of penalty parameters.
#' @param num_eigen The number of contrastive principal components to compute.
#' @param clust_method A \code{character} specifying the clustering method to
#'  use for choosing the optimal constrastive parameter. Currently, this is
#'  limited to either k-means or partitioning around medoids (PAM).
#' @param n_centers The number of n_centers to use in the clustering algorithm.
#' @param ... Additional arguments to pass to the clustering algorithm and the
#'   objective function.
#'
#' @return A list similar to that output by \code{\link[stats]{prcomp}}:
#'   \itemize{
#'     \item rotation - the matrix of variable loadings
#'     \item x - the rotated data, centred and scaled, if requested, data
#'     multiplied by the rotation matrix
#'     \item c_cov - the covariance matrix of the optimal contrastive parameter
#'     \item contrast - the optimal contrastive parameter
#'   }
#' Returns the optimal covariance matrix, contrastive parameter,
#'   loadings and projection of the target data based on the results of
#'   evaluating a clustering algorithm based on average silhouette width metric.
#'
#' @importFrom stats kmeans
#' @importFrom elasticnet spca
#' @importFrom cluster pam silhouette
#'
#' @author Philippe Boileau, \email{philippe_Boileau@@berkeley.edu}
#'
fitContrast <- function(target, center = TRUE, scale = TRUE,
                        c_contrasts, contrasts, penalties, num_eigen,
                        clust_method = c("kmeans", "pam"),
                        n_centers, ...){
  # preliminaries
  clust_method <- match.arg(clust_method)
  num_contrasts <- length(contrasts)
  num_penal <- length(penalties)

  # create the grid of contrast and penalty paramters
  param_grid <- expand.grid(penalties, contrasts)
  colnames(param_grid) <- c("lambda", "alpha")

  # create the loadings matrices
  loadings_mat <- lapply(
    seq_len(num_contrasts),
    function(x) {
      lapply(
        penalties,
        function(y) {
          if (y == 0) {
            res <- eigen(c_contrasts[[x]],
                         symmetric = TRUE
                         )$vectors[, 1:num_eigen]
          } else {
            res <- elasticnet::spca(c_contrasts[[x]],
                                    K = num_eigen,
                                    para = rep(y, num_eigen),
                                    type = "Gram",
                                    sparse = "penalty"
                                   )$loadings
          }
          colnames(res) <- paste0("V", as.character(seq(1, num_eigen)))
          return(res)
        }
      )
    }
  )

  # unlist the nested list into a single list
  loadings_mat <- unlist(loadings_mat, recursive = FALSE)

  # center and scale the target data
  target <- scale(target, center, scale)

  # for each loadings matrix, project target onto constrastive subspace
  subspaces <- lapply(
    seq_len(num_contrasts * num_penal),
    function(x) {
      as.matrix(target) %*% loadings_mat[[x]]
    }
  )

  # remove all duplicated spaces
  kernal_idx <- which(!duplicated(subspaces))
  param_grid <- param_grid[kernal_idx, ]
  loadings_mat <- loadings_mat[kernal_idx]
  subspaces <- unique(subspaces)

  # rescale all spaces to the unit hyperplane. now objective functions based
  # on metric spaces can be used
  norm_subspaces <- lapply(
    subspaces,
    function(subspace) {
      apply(subspace, 2,
            function(x){
              (x - min(x))/(max(x) - min(x))
            }
      )
    }
  )

  # get the objective function results for each space from clustering algorithm
  ave_sil_widths <- lapply(
    norm_subspaces,
    function(subspace) {
      if (clust_method == "pam") {
        clust_res <- cluster::pam(x = subspace, k = n_centers)
      } else if (clust_method == "kmeans") {
        clust_res <- stats::kmeans(x = subspace, centers = n_centers, list(...))
      }
      sil_width <- cluster::silhouette(clust_res$cluster, dist(subspace))[, 3]
      mean(sil_width)
    }
  )

  # select the best contrastive parameter, and return it's covariance matrix,
  # contrastive parameter, loadings and projection of the target data
  max_idx <- which.max(unlist(ave_sil_widths))
  return(list(
    rotation = loadings_mat[[max_idx]],
    x = subspaces[[max_idx]],
    contrast = param_grid[max_idx, 2],
    penalty = param_grid[max_idx, 1]
  ))
}
