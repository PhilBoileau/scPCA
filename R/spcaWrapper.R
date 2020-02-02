#' Sparse PCA Wrapper
#'
#' @description This wrapper function specifies which implementation of sparse
#'   pricincipal component analysis (SPCA) is used to sparsify the loadings of
#'   the contrastive covariance matrix. Currently, the \code{scPCA} package
#'   supports the iterative algorithm detailed by
#'   \insertRef{zou2006sparse}{scPCA}, and
#'   \insertRef{erichson2018sparse}{scPCA}'s randomized and non-randomized
#'   versions of SPCA solved via variable projection. These methods are
#'   implemented in the \pkg{elasticnet} and \pkg{sparsepca} packages.
#'
#' @param alg A \code{character} indicating the SPCA algorithm used to sparsify
#'   the contrastive loadings. Currently supports \code{iterative} for the
#'   \insertRef{zou2006sparse}{scPCA} implemententation, \code{var_proj} for the
#'   non-randomized \insertRef{erichson2018sparse}{scPCA} solution, and
#'   \code{rand_var_proj} fir the randomized
#'   \insertRef{erichson2018sparse}{scPCA} result.
#' @param contrast_cov A contrastive covariance \code{matrix}.
#' @param k A \code{numeric} indicating the number of eigenvectors (or
#'   sparse contrastive components) to be computed.
#' @param penalty A \code{numeric} indicating the L1 penalty parameter applied
#'   to the loadings.
#'
#' @return A \code{p \times k} sparse loadings matrix, where \code{p} is the
#'   number of features, and \code{k} is the number of sparse contrastive
#'   components.
#'
#' @importFrom sparsepca spca rspca
#' @importFrom elasticnet spca
#'
#' @keywords internal
spcaWrapper <- function(alg, contrast_cov, k, penalty) {
  
  if (alg == "iterative") {
    
    loadings_mat <- elasticnet::spca(
      contrast_cov,
      K = k,
      para = rep(penalty, k),
      type = "Gram",
      sparse = "penalty"
    )$loadings
    
  } else {
    
    # get the rootmatrix, perform thresholding
    contrast_cov_eigen <- eigen(contrast_cov)
    eig_values <- (contrast_cov_eigen$values + abs(contrast_cov_eigen$values))/2
    right_sing_mat <- contrast_cov_eigen$vectors
    contrast_cov <- right_sing_mat%*%diag(sqrt(eig_values))%*%t(right_sing_mat)
    
    if (alg == "var_proj") {
    
      loadings_mat <- sparsepca::spca(
        contrast_cov,
        k = k,
        alpha = penalty,
        beta = 1e-6,
        center = FALSE,
        scale = FALSE,
        verbose = FALSE
      )$loadings
      
    } else if (alg == "rand_var_proj") {
    
      loadings_mat <- sparsepca::rspca(
        contrast_cov,
        k = k,
        alpha = penalty,
        beta = 1e-6,
        center = FALSE,
        scale = FALSE,
        verbose = FALSE
      )$loadings
    
    }
  }
  
  return(loadings_mat)
}