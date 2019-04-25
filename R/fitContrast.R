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
#' @param num_eigen The number of contrastive principal components to compute.
#' @param clust_method A \code{character} specifying the clustering method to
#'  use for choosing the optimal constrastive parameter. Currently, this is
#'  limited to either k-means or partitioning around medoids (PAM).
#' @param n_centers The number of n_centers to use in the clustering algorithm.
#' @param ... Additional arguments to pass to the clustering algorithm and the
#'   objective function.
#'
#' @return Returns the optimal covariance matrix, contrastive parameter,
#'   loadings and projection of the target data based on the results of
#'   evaluating a clustering algorithm based on average silhouette width metric.
#'
#' @importFrom stats kmeans
#' @importFrom cluster pam silhouette
#'
fitContrast <- function(target, center = TRUE, scale = TRUE,
                        c_contrasts, contrasts, num_eigen,
                        clust_method = c("kmeans", "pam"),
                        n_centers, iter_max = 100, ...){
  # preliminaries
  clust_method <- match.arg(clust_method)
  num_contrasts <- length(contrasts)

  # get the loadings matrix of each contrastive covariance matrix
  loadings_mat <- lapply(
    seq_len(num_contrasts),
    function(x) {
        eigen(c_contrasts[[x]],
              symmetric = TRUE
        )$vectors[, seq_len(num_eigen)]
      }
    )

  # center and scale the target data
  target <- scale(target, center, scale)

  # for each loadings matrix, project target onto constrastive subspace
  subspaces <- lapply(
    seq_len(num_contrasts),
    function(x) {
      as.matrix(target) %*% loadings_mat[[x]]
    }
  )

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
    c_cov = c_contrasts[[max_idx]],
    contrast = contrasts[max_idx],
    loading_mat = loadings_mat[[max_idx]],
    proj = subspaces[[max_idx]]
  ))
}
