#' Identify the Optimal Contrastive and Penalty Parameters
#'
#' @description This function is used to automatically select the optimal
#'   contrastive parameter and L1 penalty term for scPCA based on a clustering
#'   algorithm and average silhouette width.
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
#' @param n_eigen The number of contrastive principal components to compute.
#' @param clust_method A \code{character} specifying the clustering method to
#'  use for choosing the optimal constrastive parameter. Currently, this is
#'  limited to either k-means or partitioning around medoids (PAM).
#' @param n_centers The number of n_centers to use in the clustering algorithm.
#' @param max_iters The maximum number of iterations to use in k-means
#'   clustering. Defaults to 10.
#'
#' @return A list similar to that output by \code{\link[stats]{prcomp}}:
#'   \itemize{
#'     \item rotation - the matrix of variable loadings
#'     \item x - the rotated data, centred and scaled, if requested, data
#'     multiplied by the rotation matrix
#'     \item contrast - the optimal contrastive parameter
#'     \item penalty - the optimal L1 penalty term
#'   }
#'
#' @importFrom elasticnet spca
#' @importFrom stats kmeans dist
#' @importFrom cluster pam silhouette
#'
#' @author Philippe Boileau, \email{philippe_Boileau@@berkeley.edu}
#'
fitGrid <- function(target, center, scale,
                    c_contrasts, contrasts, penalties, n_eigen,
                    clust_method = c("kmeans", "pam"),
                    n_centers, max_iters = 10){
  # preliminaries
  clust_method <- match.arg(clust_method)
  num_contrasts <- length(contrasts)
  num_penal <- length(penalties)

  # create the grid of contrast and penalty paramters
  param_grid <- expand.grid(penalties, contrasts)

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
                         )$vectors[, 1:n_eigen]
          } else {
            res <- elasticnet::spca(c_contrasts[[x]],
                                    K = n_eigen,
                                    para = rep(y, n_eigen),
                                    type = "Gram",
                                    sparse = "penalty"
                                   )$loadings
          }
          colnames(res) <- paste0("V", as.character(seq(1, n_eigen)))
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

  # remove all subspaces that had loading vectors consisting solely of zeros
  nz_load_idx <- which(sapply(norm_subspaces, function(s) sum(is.na(s))) == 0)
  norm_subspaces <- norm_subspaces[nz_load_idx]
  subspaces <- subspaces[nz_load_idx]
  param_grid <- param_grid[nz_load_idx, ]
  loadings_mat <- loadings_mat[nz_load_idx]

  # get the objective function results for each space from clustering algorithm
  ave_sil_widths <- sapply(
    norm_subspaces,
    function(subspace) {
      if (clust_method == "pam") {
        clust_res <- cluster::pam(x = subspace, k = n_centers)
      } else if (clust_method == "kmeans") {
        clust_res <- stats::kmeans(x = subspace, centers = n_centers,
                                   iter.max = max_iters)
      }
      sil_width <- cluster::silhouette(clust_res$cluster,
                                       stats::dist(subspace))[, 3]
      mean(sil_width)
    }
  )

  # select the best contrastive parameter, and return it's covariance matrix,
  # contrastive parameter, loadings and projection of the target data
  max_idx <- which.max(ave_sil_widths)
  return(
    list(
      rotation = loadings_mat[[max_idx]],
      x = subspaces[[max_idx]],
      contrast = param_grid[max_idx, 2],
      penalty = param_grid[max_idx, 1]
    )
  )
}
