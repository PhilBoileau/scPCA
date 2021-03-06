#' Contrastive Principal Component Analysis
#'
#' @description Given target and background dataframes or matrices, \code{cPCA}
#'   will perform contrastive principal component analysis (cPCA) of the target
#'   data for a given number of eigenvectors and a vector of real valued
#'   contrast parameters. This is identical to the implementation of cPCA
#'   method of \insertCite{abid2018exploring;textual}{scPCA}.
#'
#' @param target The target (experimental) data set, in a standard format such
#'  as a \code{data.frame} or \code{matrix}.
#' @param center A \code{logical} indicating whether the target and background
#'  data sets should be centered to mean zero.
#' @param scale A \code{logical} indicating whether the target and background
#'  data sets should be scaled to unit variance.
#' @param c_contrasts A \code{list} of contrastive covariances.
#' @param contrasts A \code{numeric} vector of the contrastive parameters used
#'  to compute the contrastive covariances.
#' @param n_eigen A \code{numeric} indicating the number of eigenvectors to be
#'  computed.
#' @param n_medoids A \code{numeric} indicating the number of medoids to
#'  consider. Not used if  \code{contrasts} is a single value.
#' @param eigdecomp_tol A \code{numeric} providing the level of precision used by
#'  eigendecompositon calculations. Defaults to \code{1e-10}.
#' @param eigdecomp_iter A \code{numeric} indicating the maximum number of
#'  interations performed by eigendecompositon calculations. Defaults to
#'  \code{1000}.
#'
#' @return A list of lists containing the cPCA results for each contrastive
#'   parameter deemed to be a medoid.
#'   \itemize{
#'     \item rotation - the list of matrices of variable loadings
#'     \item x - the list of rotated data, centred and scaled if requested,
#'     multiplied by the rotation matrix
#'     \item contrast - the list of contrastive parameters
#'     \item penalty - set to zero, since loadings are not penalized in cPCA
#'   }
#'
#' @references
#'   \insertAllCited{}
#'
#' @importFrom kernlab specc as.kernelMatrix
#' @importFrom Rdpack reprompt
#' @importFrom RSpectra eigs_sym
#'
#' @keywords internal
fitCPCA <- function(target, center, scale, c_contrasts, contrasts, n_eigen,
                    n_medoids, eigdecomp_tol, eigdecomp_iter) {
  # preliminaries
  num_contrasts <- length(contrasts)

  # for each contrasted covariance matrix, compute the eigenvectors
  loadings_mat <- lapply(
    seq_len(num_contrasts),
    function(x) {
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
    }
  )
  
  # center and scale the target data
  target <- safeColScale(target, center, scale)
  
  # for each loadings matrix, project target onto constrastive subspace
  spaces <- lapply(
    seq_len(num_contrasts),
    function(x) {
      as.matrix(target %*% loadings_mat[[x]])
    }
  )
  
  if (num_contrasts == 1) {
    
    # remove rownames for spaces if not null (due to DelayedMatrix mult)
    if(!is.null(rownames(spaces[[1]])[1]))
      rownames(spaces[[1]]) <- NULL
    
    out <- list(
      rotation = loadings_mat[[1]],
      x = spaces[[1]],
      contrast = contrasts,
      penalty = 0
    )
  
  } else {
  
    # produce the QR decomposition of these projections, extract Q
    qr_decomps <- lapply(
      seq_len(num_contrasts),
      function(x) {
        qr.Q(qr(spaces[[x]]))
      }
    )
  
    # populate affinity matrix for spectral clustering using the principal angles
    aff_vect <- lapply(
      seq(from = 1, to = num_contrasts - 1),
      function(i) {
        do.call(c, lapply(
          seq(from = i + 1, to = num_contrasts),
          function(j) {
            Q_i <- qr_decomps[[i]]
            Q_j <- qr_decomps[[j]]
            d <- svd(x = t(Q_i) %*% Q_j, nu = 0, nv = 0)$d
            d[1] * d[2]
          }
        ))
      }
    )
  
    aff_mat <- diag(x = 0.5, nrow = num_contrasts)
    aff_mat[lower.tri(aff_mat, diag = FALSE)] <- unlist(aff_vect)
    aff_mat <- t(aff_mat)
  
    # fix any computation errors, see numpy.nan_to_num
    aff_mat[is.nan(aff_mat)] <- 0
    aff_mat[is.na(aff_mat)] <- 0
    aff_mat[is.infinite(aff_mat)] <- 1000
    aff_mat <- t(aff_mat) + aff_mat
  
    # perfrom spectral clustering using the affinity matrix
    spec_clust <- kernlab::specc(kernlab::as.kernelMatrix(aff_mat),
      centers = n_medoids
    )
  
    # identify the alpha medoids of the spectral clustering
    contrast_medoids <- do.call(c, lapply(
      seq_len(n_medoids),
      function(x) {
        sub_index <- which(spec_clust == x)
        sub_aff_mat <- aff_mat[sub_index, sub_index]
        if (is.matrix(sub_aff_mat)) {
          aff_sums <- colSums(sub_aff_mat)
          contrasts[sub_index[which.max(aff_sums)]]
        } else {
          contrasts[sub_index]
        }
      }
    ))
  
    # create the lists of contrastive parameter medoids, loadings and projections
    med_index <- which(contrasts %in% contrast_medoids)
    
    # remove rownames for spaces if not null (due to DelayedMatrix mult)
    if(!is.null(rownames(spaces[[1]])[1])) {
      spaces <- lapply(
        med_index,
        function(id) {
          rownames(spaces[[id]]) <- NULL
          spaces[[id]]
        }
      )
    } else {
      spaces <- spaces[med_index]
    }
    
    out <- list(
      rotation = loadings_mat[med_index],
      x = spaces,
      contrast = contrasts[med_index],
      penalty = rep(0, length(med_index))
    )
  }
  
  return(out)
}

################################################################################

#' Contrastive Principal Component Analysis in Parallel
#'
#' @description Given target and background dataframes or matrices, \code{cPCA}
#'   will perform contrastive principal component analysis (cPCA) of the target
#'   data for a given number of eigenvectors and a vector of real valued
#'   contrast parameters. This is identical to the implementation of cPCA
#'   method by Abid et al. \insertCite{abid2018exploring;textual}{scPCA}.
#'   Analogous to \code{\link{fitCPCA}}, but replaces all \code{lapply} calls by
#'   \code{\link[BiocParallel]{bplapply}}.
#'
#' @param target The target (experimental) data set, in a standard format such
#'  as a \code{data.frame} or \code{matrix}.
#' @param center A \code{logical} indicating whether the target and background
#'  data sets should be centered to mean zero.
#' @param scale A \code{logical} indicating whether the target and background
#'  data sets should be scaled to unit variance.
#' @param c_contrasts A \code{list} of contrastive covariances.
#' @param contrasts A \code{numeric} vector of the contrastive parameters used
#'  to compute the contrastive covariances.
#' @param n_eigen A \code{numeric} indicating the number of eigenvectors to be
#'  computed.
#' @param n_medoids A \code{numeric} indicating the number of medoids to
#'  consider.
#' @param eigdecomp_tol A \code{numeric} providing the level of precision used by
#'  eigendecompositon calculations. Defaults to \code{1e-10}.
#' @param eigdecomp_iter A \code{numeric} indicating the maximum number of
#'  interations performed by eigendecompositon calculations. Defaults to
#'  \code{1000}.
#'
#' @return A list of lists containing the cPCA results for each contrastive
#'   parameter deemed to be a medoid.
#'   \itemize{
#'     \item rotation - the list of matrices of variable loadings
#'     \item x - the list of rotated data, centred and scaled if requested,
#'     multiplied by the rotation matrix
#'     \item contrast - the list of contrastive parameters
#'     \item penalty - set to zero, since loadings are not penalized in cPCA
#'   }
#'
#' @references
#'   \insertAllCited{}
#'
#' @importFrom kernlab specc as.kernelMatrix
#' @importFrom BiocParallel bplapply
#' @importFrom Rdpack reprompt
#' @importFrom RSpectra eigs_sym
#'
#' @keywords internal
bpFitCPCA <- function(target, center, scale, c_contrasts, contrasts, n_eigen,
                      n_medoids, eigdecomp_tol, eigdecomp_iter) {
  # preliminaries
  num_contrasts <- length(contrasts)

  # for each contrasted covariance matrix, compute the eigenvectors
  loadings_mat <- BiocParallel::bplapply(
    seq_len(num_contrasts),
    function(x) {
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
    }
  )

  # center and scale the target data
  target <- safeColScale(target, center, scale)

  # for each loadings matrix, project target onto constrastive subspace
  spaces <- BiocParallel::bplapply(
    seq_len(num_contrasts),
    function(x) {
      as.matrix(target %*% loadings_mat[[x]])
    }
  )

  # produce the QR decomposition of these projections, extract Q
  qr_decomps <- BiocParallel::bplapply(
    seq_len(num_contrasts),
    function(x) {
      qr.Q(qr(spaces[[x]]))
    }
  )

  # populate affinity matrix for spectral clustering using the principal angles
  aff_vect <- BiocParallel::bplapply(
    seq(from = 1, to = num_contrasts - 1),
    function(i) {
      do.call(c, lapply(
        seq(from = i + 1, to = num_contrasts),
        function(j) {
          Q_i <- qr_decomps[[i]]
          Q_j <- qr_decomps[[j]]
          d <- svd(x = t(Q_i) %*% Q_j, nu = 0, nv = 0)$d
          d[1] * d[2]
        }
      ))
    }
  )

  aff_mat <- diag(x = 0.5, nrow = num_contrasts)
  aff_mat[lower.tri(aff_mat, diag = FALSE)] <- unlist(aff_vect)
  aff_mat <- t(aff_mat)

  # fix any computation errors, see numpy.nan_to_num
  aff_mat[is.nan(aff_mat)] <- 0
  aff_mat[is.na(aff_mat)] <- 0
  aff_mat[is.infinite(aff_mat)] <- 1000
  aff_mat <- t(aff_mat) + aff_mat

  # perfrom spectral clustering using the affinity matrix
  spec_clust <- kernlab::specc(kernlab::as.kernelMatrix(aff_mat),
    centers = n_medoids
  )

  # identify the alpha medoids of the spectral clustering
  contrast_medoids <- BiocParallel::bplapply(
    seq_len(n_medoids),
    function(x) {
      sub_index <- which(spec_clust == x)
      sub_aff_mat <- aff_mat[sub_index, sub_index]
      if (is.matrix(sub_aff_mat)) {
        aff_sums <- colSums(sub_aff_mat)
        contrasts[sub_index[which.max(aff_sums)]]
      } else {
        contrasts[sub_index]
      }
    }
  )
  contrast_medoids <- unlist(contrast_medoids)

  # create the lists of contrastive parameter medoids, loadings and projections
  med_index <- which(contrasts %in% contrast_medoids)
  
  # remove rownames for spaces if not null (due to DelayedMatrix mult)
  if(!is.null(dimnames(spaces[[1]])[[1]])) {
    spaces <- lapply(
      med_index,
      function(id) {
        rownames(spaces[[id]]) <- NULL
        spaces[[id]]
      }
    )
  } else {
    spaces <- spaces[med_index]
  }
  
  out <- list(
    rotation = loadings_mat[med_index],
    x = spaces,
    contrast = contrasts[med_index],
    penalty = rep(0, length(med_index))
  )
  return(out)
}
