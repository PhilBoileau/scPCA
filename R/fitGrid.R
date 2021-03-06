#' Identify the Optimal Contrastive and Penalty Parameters
#'
#' @description This function is used to automatically select the optimal
#'   contrastive parameter and L1 penalty term for scPCA based on a clustering
#'   algorithm and average silhouette width.
#'
#' @param target The target (experimental) data set, in a standard format such
#'  as a \code{data.frame} or \code{matrix}.
#' @param target_valid A holdout set of the target (experimental) data set, in
#'  a standard format such as a \code{data.frame} or \code{matrix}. \code{NULL}
#'  by default but used by \code{\link{cvSelectParams}} for cross-validated
#'  selection of the contrastive and penalization parameters.
#' @param center A \code{logical} indicating whether the target and background
#'  data sets should be centered to mean zero.
#' @param scale A \code{logical} indicating whether the target and background
#'  data sets should be scaled to unit variance.
#' @param c_contrasts A \code{list} of contrastive covariances.
#' @param contrasts A \code{numeric} vector of the contrastive parameters used
#'  to compute the contrastive covariances.
#' @param alg A \code{character} indicating the SPCA algorithm used to sparsify
#'   the contrastive loadings. Currently supports \code{iterative} for the
#'   \insertCite{zou2006sparse;textual}{scPCA} implemententation,
#'   \code{var_proj} for the non-randomized
#'   \insertCite{erichson2018sparse;textual}{scPCA} solution, and
#'   \code{rand_var_proj} for the randomized
#'   \insertCite{erichson2018sparse;textual}{scPCA} result.
#' @param penalties A \code{numeric} vector of the penalty terms.
#' @param n_eigen A \code{numeric} indicating the number of eigenvectors to be
#'  computed.
#' @param clust_method A \code{character} specifying the clustering method to
#'  use for choosing the optimal constrastive parameter. Currently, this is
#'  limited to either k-means, partitioning around medoids (PAM), and
#'  hierarchical clustering. The default is k-means clustering.
#' @param n_centers A \code{numeric} giving the number of centers to use in the
#'  clustering algorithm.
#' @param max_iter A \code{numeric} giving the maximum number of iterations to
#'   be used in k-means clustering, defaulting to 10.
#' @param linkage_method A \code{character} specifying the agglomerative
#'   linkage method to be used if \code{clust_method = "hclust"}. The options
#'   are \code{ward.D2}, \code{single}, \code{complete}, \code{average},
#'   \code{mcquitty}, \code{median}, and \code{centroid}. The default is
#'   \code{complete}.
#' @param clusters A \code{numeric} vector of cluster labels for observations in
#'   the \code{target} data. Defaults to \code{NULL}, but is otherwise used to
#'   identify the optimal set of hyperparameters when fitting the scPCA and the
#'   automated version of cPCA.
#' @param eigdecomp_tol A \code{numeric} providing the level of precision used by
#'   eigendecompositon calculations. Defaults to \code{1e-10}.
#' @param eigdecomp_iter A \code{numeric} indicating the maximum number of
#'   interations performed by eigendecompositon calculations. Defaults to
#'   \code{1000}.
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
#' @importFrom stats kmeans dist hclust cutree
#' @importFrom cluster pam silhouette
#' @importFrom RSpectra eigs_sym
#'
#' @references
#'   \insertAllCited{}
#'
#' @keywords internal
fitGrid <- function(target, target_valid = NULL, center, scale,
                    c_contrasts, contrasts, alg, penalties, n_eigen,
                    clust_method = c("kmeans", "pam", "hclust"),
                    n_centers, max_iter = 10,
                    linkage_method = "complete", clusters = NULL,
                    eigdecomp_tol = 1e-10, eigdecomp_iter = 1000) {
  # preliminaries
  num_contrasts <- length(contrasts)
  num_penal <- length(penalties)

  # create the grid of contrast and penalty parameters
  param_grid <- expand.grid(penalties, contrasts)

  # create the loadings matrices
  loadings_mat <- lapply(
    seq_len(num_contrasts),
    function(x) {
      lapply(
        penalties,
        function(y) {
          if (y == 0) {
            withCallingHandlers(
              res <- RSpectra::eigs_sym(c_contrasts[[x]],
                k = n_eigen,
                which = "LA",
                opts = list(tol = eigdecomp_tol, maxitr = eigdecomp_iter)
              )$vectors,
              warning = function(w) {
                warning(paste0(
                  "\nFor contrastive parameter = ",
                  round(contrasts[[x]], 3), ":\n")
                )
              }
            )
          } else {
            res <- spcaWrapper(
              alg = alg,
              contrast = contrasts[[x]],
              contrast_cov = c_contrasts[[x]],
              k = n_eigen,
              penalty = y,
              eigdecomp_tol = eigdecomp_tol,
              eigdecomp_iter = eigdecomp_iter
            )
          }
          colnames(res) <- paste0("V", as.character(seq_len(ncol(res))))
          res
        }
      )
    }
  )

  # unlist the nested list into a single list
  loadings_mat <- unlist(loadings_mat, recursive = FALSE)

  # center and scale the target data
  target <- safeColScale(target, center, scale)

  if (is.null(target_valid)) {
    # for each loadings matrix, project target onto constrastive subspace
    subspaces <- lapply(
      seq_len(num_contrasts * num_penal),
      function(x) {
        as.matrix(target %*% loadings_mat[[x]])
      }
    )
  } else {
    
    # center and scale the holdout set of target data
    target_valid <- safeColScale(target_valid, center, scale)

    # for each loadings matrix, project holdout target on constrastive subspace
    subspaces <- lapply(
      seq_len(num_contrasts * num_penal),
      function(x) {
        as.matrix(target_valid %*% loadings_mat[[x]])
      }
    )
  }

  if (is.null(n_centers)) {
    
    # remove rownames for spaces if not null (due to DelayedMatrix mult)
    if(!is.null(rownames(subspaces[[1]])[1]))
      rownames(subspaces[[1]]) <- NULL
    
    out <- list(
      rotation = loadings_mat[[1]],
      x = subspaces[[1]],
      contrast = contrasts[[1]],
      penalty = penalties[[1]]
    )
  } else {
  
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
        max_val <- max(subspace[, 1])
        min_val <- min(subspace[, 1])
        apply(
          subspace, 2,
          function(x) {
            x / (max_val - min_val)
          }
        )
      }
    )
  
    # remove all subspaces that had loading vectors consisting solely of zeros
    zero_subs <- do.call(c, lapply(subspaces, function(s) {
      any(apply(s, 2, function(l) all(l < 1e-6)))
    }))
    zero_subs_norm <- do.call(c, lapply(norm_subspaces, function(ns) {
      any(apply(ns, 2, function(l) all(l < 1e-6)))
    }))
    nz_load_idx <- which((zero_subs + zero_subs_norm) == 0)
    norm_subspaces <- norm_subspaces[nz_load_idx]
    subspaces <- subspaces[nz_load_idx]
    param_grid <- param_grid[nz_load_idx, ]
    loadings_mat <- loadings_mat[nz_load_idx]
  
    # get the objective function results for each space from clustering algorithm
    ave_sil_widths <- do.call(c, lapply(
      norm_subspaces, function(subspace) {
        if (!is.null(clusters)) {
          sil_width <- cluster::silhouette(
            clusters,
            stats::dist(subspace)
          )[, 3]
        } else if (clust_method == "pam") {
          clust_res <- cluster::pam(x = subspace, k = n_centers)
        } else if (clust_method == "kmeans") {
          clust_res <- stats::kmeans(
            x = subspace, centers = n_centers,
            iter.max = max_iter
          )
        } else if (clust_method == "hclust") {
          dist_matrix <- stats::dist(x = subspace, method = "euclidean")
          hclust_res <- stats::hclust(d = dist_matrix, method = linkage_method)
          clust_res <- stats::cutree(tree = hclust_res, k = n_centers)
          sil_width <- cluster::silhouette(
            clust_res,
            dist_matrix
          )[, 3]
        }
        if (clust_method %in% c("kmeans", "pam") && is.null(clusters)) {
          sil_width <- cluster::silhouette(
            clust_res$cluster,
            stats::dist(subspace)
          )[, 3]
        }
        mean(sil_width)
      }
    ))
  
    # remove rownames for spaces if not null (due to DelayedMatrix mult)
    if(!is.null(rownames(subspaces[[1]])[1])) {
      subspaces <- lapply(
        seq_len(length(subspaces)),
        function(id) {
          rownames(subspaces[[id]]) <- NULL
          subspaces[[id]]
        }
      )
    }
    
    # select the best contrastive parameter, and return its covariance matrix,
    # contrastive parameter, loadings and projection of the target data
    out <- list(
      rotation = loadings_mat,
      x = subspaces,
      contrast = param_grid[, 2],
      penalty = param_grid[, 1],
      ave_sil_widths = ave_sil_widths
    )
  }
  
  return(out)
}

################################################################################

#' Identify the Optimal Contrastive and Penalty Parameters in Parallel
#'
#' @description This function is used to automatically select the optimal
#'   contrastive parameter and L1 penalty term for scPCA based on a clustering
#'   algorithm and average silhouette width. Analogous to \code{\link{fitGrid}},
#'   but replaces all \code{lapply} calls by
#'   \code{\link[BiocParallel]{bplapply}}.
#'
#' @param target The target (experimental) data set, in a standard format such
#'  as a \code{data.frame} or \code{matrix}.
#' @param target_valid A holdout set of the target (experimental) data set, in a
#'  standard format such as a \code{data.frame} or \code{matrix}. \code{NULL} by
#'  default but used by \code{\link{cvSelectParams}} for cross-validated
#'  selection of the contrastive and penalization parameters.
#' @param center A \code{logical} indicating whether the target and background
#'  data sets should be centered to mean zero.
#' @param scale A \code{logical} indicating whether the target and background
#'  data sets should be scaled to unit variance.
#' @param c_contrasts A \code{list} of contrastive covariances.
#' @param contrasts A \code{numeric} vector of the contrastive parameters used
#'  to compute the contrastive covariances.
#' @param alg A \code{character} indicating the SPCA algorithm used to sparsify
#'  the contrastive loadings. Currently supports \code{iterative} for the
#'  \insertCite{zou2006sparse;textual}{scPCA} implemententation, \code{var_proj}
#'  for the non-randomized \insertCite{erichson2018sparse;textual}{scPCA}
#'  solution, and \code{rand_var_proj} fir the randomized
#'  \insertCite{erichson2018sparse;textual}{scPCA} result.
#' @param penalties A \code{numeric} vector of the penalty terms.
#' @param n_eigen A \code{numeric} indicating the number of eigenvectors to be
#'  computed.
#' @param clust_method A \code{character} specifying the clustering method to
#'  use for choosing the optimal constrastive parameter. Currently, this is
#'  limited to either k-means, partitioning around medoids (PAM), and
#'  hierarchical clustering. The default is k-means clustering.
#' @param n_centers A \code{numeric} giving the number of centers to use in the
#'  clustering algorithm.
#' @param max_iter A \code{numeric} giving the maximum number of iterations to
#'   be used in k-means clustering, defaulting to 10.
#' @param linkage_method A \code{character} specifying the agglomerative linkage
#'   method to be used if \code{clust_method = "hclust"}. The options are
#'   \code{ward.D2}, \code{single}, \code{complete},
#'   \code{average}, \code{mcquitty}, \code{median}, and \code{centroid}. The
#'   default is \code{complete}.
#' @param clusters A \code{numeric} vector of cluster labels for observations in
#'   the \code{target} data. Defaults to \code{NULL}, but is otherwise used to
#'   identify the optimal set of hyperparameters when fitting the scPCA and the
#'   automated version of cPCA.
#' @param eigdecomp_tol A \code{numeric} providing the level of precision used by
#'   eigendecompositon calculations. Defaults to \code{1e-10}.
#' @param eigdecomp_iter A \code{numeric} indicating the maximum number of
#'   interations performed by eigendecompositon calculations. Defaults to
#'   \code{1000}.
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
#' @importFrom stats kmeans dist hclust cutree
#' @importFrom cluster pam silhouette
#' @importFrom BiocParallel bplapply
#' @importFrom RSpectra eigs_sym
#'
#' @references
#'   \insertAllCited{}
#'
#' @keywords internal
bpFitGrid <- function(target, target_valid = NULL, center, scale,
                      c_contrasts, contrasts, penalties, n_eigen,
                      alg, clust_method = c("kmeans", "pam", "hclust"),
                      n_centers, max_iter = 10,
                      linkage_method = "complete", clusters = NULL,
                      eigdecomp_tol = 1e-10, eigdecomp_iter = 1000) {
  # preliminaries
  num_contrasts <- length(contrasts)
  num_penal <- length(penalties)

  # create the grid of contrast and penalty parameters
  param_grid <- expand.grid(penalties, contrasts)

  # create the loadings matrices
  loadings_mat <- BiocParallel::bplapply(
    seq_len(num_contrasts),
    function(x) {
      lapply(
        penalties,
        function(y) {
          if (y == 0) {
            withCallingHandlers(
              res <- RSpectra::eigs_sym(
                c_contrasts[[x]],
                k = n_eigen,
                which = "LA",
                opts = list(tol = eigdecomp_tol, maxitr = eigdecomp_iter)
              )$vectors,
              warning = function(w) {
                warning(paste0(
                  "\nFor contrastive parameter = ",
                  round(contrasts[[x]], 3), ":\n")
                )
              }
            )
          } else {
            res <- spcaWrapper(
              alg = alg,
              contrast = contrasts[[x]],
              contrast_cov = c_contrasts[[x]],
              k = n_eigen,
              penalty = y,
              eigdecomp_tol = eigdecomp_tol,
              eigdecomp_iter = eigdecomp_iter
            )
          }
          colnames(res) <- paste0("V", as.character(seq_len(ncol(res))))
          return(res)
        }
      )
    }
  )

  # unlist the nested list into a single list
  loadings_mat <- unlist(loadings_mat, recursive = FALSE)

  # center and scale the target data
  target <- safeColScale(target, center, scale)

  if (is.null(target_valid)) {
    # for each loadings matrix, project target onto constrastive subspace
    subspaces <- BiocParallel::bplapply(
      seq_len(num_contrasts * num_penal),
      function(x) {
        as.matrix(target %*% loadings_mat[[x]])
      }
    )
  } else {
    # center and scale the holdout set of target data
    target_valid <- safeColScale(target_valid, center, scale)

    # for each loadings matrix, project holdout target on constrastive subspace
    subspaces <- BiocParallel::bplapply(
      seq_len(num_contrasts * num_penal),
      function(x) {
        as.matrix(target_valid %*% loadings_mat[[x]])
      }
    )
  }

  # remove all duplicated spaces
  kernal_idx <- which(!duplicated(subspaces))
  param_grid <- param_grid[kernal_idx, ]
  loadings_mat <- loadings_mat[kernal_idx]
  subspaces <- unique(subspaces)

  # rescale all spaces to the unit hyperplane. now objective functions based
  # on metric spaces can be used
  norm_subspaces <- BiocParallel::bplapply(
    subspaces,
    function(subspace) {
      max_val <- max(subspace[, 1])
      min_val <- min(subspace[, 1])
      apply(
        subspace, 2,
        function(x) {
          x / (max_val - min_val)
        }
      )
    }
  )

  # remove all subspaces that had loading vectors consisting solely of zeros
  zero_subs <- do.call(c, lapply(
    subspaces,
    function(s) {
      any(apply(s, 2, function(l) all(l < 1e-6)))
    }
  ))
  zero_subs_norm <- do.call(c, lapply(
    norm_subspaces,
    function(ns) {
      any(apply(ns, 2, function(l) all(l < 1e-6)))
    }
  ))
  nz_load_idx <- which((zero_subs + zero_subs_norm) == 0)
  norm_subspaces <- norm_subspaces[nz_load_idx]
  subspaces <- subspaces[nz_load_idx]
  param_grid <- param_grid[nz_load_idx, ]
  loadings_mat <- loadings_mat[nz_load_idx]

  # get the objective function results for each space from clustering algorithm
  ave_sil_widths <- BiocParallel::bplapply(
    norm_subspaces,
    function(subspace) {
      if (!is.null(clusters)) {
        sil_width <- cluster::silhouette(
          clusters,
          stats::dist(subspace)
        )[, 3]
      } else if (clust_method == "pam") {
        clust_res <- cluster::pam(x = subspace, k = n_centers)
      } else if (clust_method == "kmeans") {
        clust_res <- stats::kmeans(
          x = subspace, centers = n_centers,
          iter.max = max_iter
        )
      } else if (clust_method == "hclust") {
        dist_matrix <- stats::dist(x = subspace, method = "euclidean")
        hclust_res <- stats::hclust(d = dist_matrix, method = linkage_method)
        clust_res <- stats::cutree(tree = hclust_res, k = n_centers)
        sil_width <- cluster::silhouette(
          clust_res,
          dist_matrix
        )[, 3]
      }
      if (clust_method %in% c("kmeans", "pam") && is.null(clusters)) {
        sil_width <- cluster::silhouette(
          clust_res$cluster,
          stats::dist(subspace)
        )[, 3]
      }
      mean(sil_width)
    }
  )
  ave_sil_widths <- unlist(ave_sil_widths)

  
  # remove rownames for spaces if not null (due to DelayedMatrix mult)
  if(!is.null(rownames(subspaces[[1]])[1])) {
    subspaces <- lapply(
      seq_len(length(subspaces)),
      function(id) {
        rownames(subspaces[[id]]) <- NULL
        subspaces[[id]]
      }
    )
  }
  
  # select the best contrastive parameter, and return its covariance matrix,
  # contrastive parameter, loadings and projection of the target data
  out <- list(
    rotation = loadings_mat,
    x = subspaces,
    contrast = param_grid[, 2],
    penalty = param_grid[, 1],
    ave_sil_widths = ave_sil_widths
  )
  return(out)
}
