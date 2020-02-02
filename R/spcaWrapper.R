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
#' @keywords internal
spcaWrapper <- function(alg, contrast_cov, k, penalty) {
  
  if (alg == "iterative") {
    
  } else if (alg == "var_proj") {
    
  } else if (alg == "rand_var_proj") {
    
  } 
  
  return(loadings)
}