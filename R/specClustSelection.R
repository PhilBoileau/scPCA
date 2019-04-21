#' Use Spectral Clustering to Select Contrastive Projections
#'
#' @param c_proj The contrastive projection list produced by
#'   \code{\link{projGridCP}}.
#' @param num_medoids The number of medoids to identify using spectral
#'   clustering.
#'
#' @return A named list containing the matrix of parameters, the loadings matrix
#'   and the contrastive projection of the target associated with each medoid.
#'
#' @importFrom kernlab specc as.kernelMatrix
#'
#' @author Philippe Boileau, \email{philippe_Boileau@@berkeley.edu}
#'
specClustSelection <- function(c_proj, num_medoids) {

  # get the number of spaces that are projected onto
  num_spaces <- length(c_proj$spaces)

  # produce the QR decomposition of these projections, extract Q
  qr_decomps <- lapply(
    seq_len(num_spaces),
    function(x) {
      qr.Q(qr(c_proj$spaces[[x]]))
    }
  )

  # populate affinity matrix for spectral clustering using the principal angles
  aff_vect <- sapply(
    seq_len(num_spaces - 1),
    function(i) {
      sapply(
        (i + 1):num_spaces,
        function(j) {
          Q_i <- qr_decomps[[i]]
          Q_j <- qr_decomps[[j]]
          d <- svd(x = t(Q_i) %*% Q_j, nu = 0, nv = 0)$d
          return(d[1] * d[2])
        }
      )
    }
  )
  aff_mat <- diag(x = 0.5, nrow = num_spaces)
  aff_mat[lower.tri(aff_mat, diag = FALSE)] <- unlist(aff_vect)
  aff_mat <- t(aff_mat)
  # fix any computation errors, see numpy.nan_to_num
  aff_mat[is.nan(aff_mat)] <- 0
  aff_mat[is.na(aff_mat)] <- 0
  aff_mat[is.infinite(aff_mat)] <- 1
  aff_mat <- t(aff_mat) + aff_mat

  # perfrom spectral clustering using the affinity matrix
  spec_clust <- kernlab::specc(kernlab::as.kernelMatrix(aff_mat),
    centers = num_medoids,
    iterations = 10000
  )

  # identify the alpha medoids of the spectral clustering
  contrast_medoids <- sapply(
    seq_len(num_medoids),
    function(x) {
      sub_index <- which(spec_clust == x)
      sub_aff_mat <- as.matrix(
        aff_mat[sub_index, sub_index]
      )
      aff_sums <- colSums(sub_aff_mat)
      return(
        c_proj$param_grid[sub_index[which.max(aff_sums)], ]
      )
    }
  )

  # fix formating
  contrast_medoids <- matrix(unlist(contrast_medoids),
    nrow = num_medoids,
    byrow = TRUE
  )
  colnames(contrast_medoids) <- c("lambda", "alpha")

  # create the lists of contrastive parameter medoids, loadings and projections
  contrast_medoids <- contrast_medoids[order(
    contrast_medoids[, 2],
    contrast_medoids[, 1]
  ), ]

  # get the index of the paramaters chosen as medoids
  combo <- rbind(c_proj$param_grid, contrast_medoids)
  rownames(combo) <- seq(1, nrow(combo))
  med_index <- as.numeric(
    rownames(combo[duplicated(combo, fromLast = TRUE), , drop = TRUE])
  )
  med_loadings_mat <- c_proj$loadings_mat[med_index]
  med_spaces <- c_proj$spaces[med_index]

  return(list(
    medoids_params = contrast_medoids,
    med_loadings_mat = med_loadings_mat,
    med_spaces = med_spaces
  ))
}
