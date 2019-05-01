#' Cross-Validation Scheme for SPC on Target Data
#'
#' @details This function is a modified version of the \code{\link[PMA]{SPC.cv}}
#'   function. The cross-validation procedure is described in
#'   \insertRef{WittenPMD2009}{scPCA}.
#'
#' @param target The target data set.
#' @param c_proj The optimal contrastive projection of the target data.
#' @param cov_mat The optimal contrastive covariance matrix.
#' @param penalties The vector of penalties for the SPC.
#' @param n_eigen The number of components to calculate.
#' @param n_iter Number of iterations to be performed.
#' @param n_folds The number of folds used for CV.
#' @param v The first right singular vector of \code{cov_mat}.
#'
#' @importFrom PMA SPC
#'
#' @return The list of optimal penalty terms based on cross-validated mean
#'   squared error:
#'   \itemize{
#'     \item best - The optimal L1 penalty based of CV MSE
#'     \item best1se - The sparsest penalty within 1 standard error of CV MSE
#'   }
#'
cvSPC <- function(target, c_proj, cov_mat, penalties, n_eigen, n_iter,
                  n_folds, v){

  # identify the percent to remove
  percent_rm <- min(0.25, 1/n_folds)

  # matrix of random numbers, used to identify which entry to mask in each fold
  rand_mat <- matrix(runif(nrow(c_proj)*n_eigen), ncol = n_eigen)

  # initialize the matrix to contain the cv errors
  cv_err_mat <- matrix(NA, nrow = n_folds, ncol = length(penalties))

  # loop through each fold
  cv_err_ls <- lapply(
    seq_len(n_folds),
    function(x){
      rm_idx <- ((x - 1) * percent_rm < rand_mat) &
                  (rand_mat < x * percent_rm)
      c_proj_rm <- c_proj
      c_proj_rm[rm_idx] <- NA
      sapply(
        penalties,
        function(p){
          res <- SPC(cov_mat, sumabsv = p, niter = n_iter, v = v,
                     trace = FALSE, center = FALSE, K = n_eigen)
          c_proj_hat <- as.matrix(target) %*% res$v
          sum(((c_proj_hat - c_proj)[rm_idx])^2)
        })
    })

  # turn cv_err_ls into a matrix for easy searching
  cv_err_mat <- matrix(unlist(cv_err_ls),
                       nrow = n_folds,
                       byrow = TRUE)
  cv_risk <- apply(cv_err_mat, 2, mean)
  cv_sd <- apply(cv_err_mat, 2, sd)/sqrt(n_folds)

  # select the optimal penalty term (if there are ties, select smallest)
  best <- penalties[which.min(cv_risk)]
  # select the smallest penalty term 1 sd away from the optimal one
  best1se <- penalties[min(which(cv_risk < min(cv_risk) +
                                   cv_sd[which.min(cv_risk)]))]

  return(
    list(
      best = best,
      best1se = best1se
    )
  )
}
