#' Identify the Optimal Contrastive Parameter
#'
#' @param c_contrasts List of of contrastive covariances.
#' @param contrasts Vector of contrastive parameter values used to compute the
#'   contrastive covariances,
#' @param cluster_alg The clustering algorithm used to cluster each space.
#'   Defaults to k-means.
#' @param centers The number of centers to use in the clustering algorithm.
#' @param obj_fun The objective function to maximize. Defaults to average
#'   silhouette width.
#' @param ... Additional arguments to pass to the cluster algorithm and the
#'   objective function.
#'
#' @return Returns the contrastive projection that optimizes the specified
#'   objective function
#'
#' @imortFrom stats kmeans

fitContrast <- function(c_contrasts, contrasts, cluster_alg = "kmeans",
                        centers, obj_fun = "silhouette", ...){



}
