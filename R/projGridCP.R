#' Project Subspaces for Contrasts, Penalty Parameters
#'
#' @description Project the target data into the contrastive subspaces defined
#'   by the unique combinations of the contrast parameters and the penalty
#'   terms. All loading matrices that are kernels are dropped.
#'
#' @param target The target data set
#' @param center Whether the data sets' columns should be centered
#' @param scale Whetehr the data sets' columns should be scaled to have variance
#' @param c_contrasts List of contrastve covariance matrices
#' @param contrasts Vector of contrast parameters
#' @param penalties Vector of penalty parameters
#' @param num_eigen The number of contrastive principal components to compute.
#'
#' @importFrom elasticnet spca
#'
#' @author Philippe Boileau, \email{philippe_Boileau@@berkeley.edu}
#'
#' @return A named list containing the loadings matrices, the projections of
#'   the target data into the spaces that they define and the grid of
#'   parameters.
projGridCP <- function(target, center, scale, c_contrasts, contrasts,
                       penalties, num_eigen) {
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
            eigen(c_contrasts[[x]],
              symmetric = TRUE
            )$vectors[, 1:num_eigen]
          } else {
            elasticnet::spca(c_contrasts[[x]],
              K = num_eigen,
              para = rep(y, num_eigen),
              type = "Gram",
              sparse = "penalty"
            )$loadings
          }
        }
      )
    }
  )

  # unlist the nested list into a single list
  loadings_mat <- unlist(loadings_mat, recursive = FALSE)

  # center and scale the target data set
  target <- scale(target, center, scale)

  # for each loadings matrix, project target onto constrastive subspace
  spaces <- lapply(
    seq_len(num_contrasts * num_penal),
    function(x) {
      as.matrix(target) %*% loadings_mat[[x]]
    }
  )

  # remove all spaces projected to the 0 vector, update parameter grid, loadings
  kernal_idx <- which(!duplicated(spaces))
  param_grid <- param_grid[kernal_idx, ]
  loadings_mat <- loadings_mat[kernal_idx]
  spaces <- unique(spaces)

  return(list(
    param_grid = param_grid,
    loadings_mat = loadings_mat,
    spaces = spaces
  ))
}
