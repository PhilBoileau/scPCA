#' Identify the Optimal Contrastive Parameter
#'
#' @description This function is used to automatically select the optimal
#'   contrastive parameter based on k-means and average silhouette width.
#'
#' @param target The target data set.
#' @param center Whether the data sets' columns should be centered.
#' @param c_contrasts List of of contrastive covariances.
#' @param contrasts Vector of contrastive parameter values used to compute the
#'   contrastive covariances,
#' @param num_eigen The number of contrastive principal components to compute.
#' @param centers The number of centers to use in the clustering algorithm.
#' @param ... Additional arguments to pass to the cluster algorithm and the
#'   objective function.
#'
#' @return Returns the optimal covariance matrix, contrastive parameter,
#'   loadings and projection of the target data based on the kmeans results and
#'   the average silhouette width
#'
#' @imortFrom stats kmeans
#' @importFrom cluster silhouette

fitContrast <- function(target, center = TRUE, c_contrasts, contrasts,
                        num_eigen, centers, ...){

  num_contrasts <- length(contrasts)

  # get the loadings matrix of each contrastive covariance matrix
  loadings_mat <- lapply(
    seq_len(num_contrasts),
    function(x) {
        eigen(c_contrasts[[x]],
              symmetric = TRUE
        )$vectors[, 1:num_eigen]
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

  # get the objective function results for each space based on the clust alg
  ave_sil_widths <- lapply(
    norm_subspaces,
    function(subspace) {
      kmeans_res <- kmeans(subspace, centers = centers, ...)
      sil_width <- silhouette(kmeans_res$cluster, dist(subspace))$sil_width
      mean(sil_width)
    }
  )

  # select the best contrastive parameter, return it's covariance matrix,
  # contrastive parameter, loadings and projection of the target data
  min_idx <- which.min(unlist(sil_width))
  return(list(
    c_cov = c_contrasts[[min_idx]],
    contrast = contrasts[min_idx],
    loading_mat = loadings_mat[[min_idx]],
    subspaces[[min_idx]]
  ))
}
