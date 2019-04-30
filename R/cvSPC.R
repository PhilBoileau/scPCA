#' Cross-Validation Scheme for SPC on Target Data
#'
#' @param data The target data set.
#' @param n_eigen The number of components to calculate.
#' @param n_iter Number of iterations to be performed.
#' @param n_folds The number of folds used for CV.
#' @param penalties The vector of penalties for the SPC.
#' @param cov_mat The optimal contrastive covariance matrix.
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
cvSPC <- function(target, n_eigen, n_iter, n_folds, penalties, cov_mat, v){

  # identify the percent to remove
  percent_rm <- min(0.25, 1/n_folds)

  # matrix of random numbers, used to identify which entry to mask in each fold
  rand_mat <- matrix(runif(nrow(target)*ncol(target), ncol = ncol(target)))

  # initialize the matrix to contain the cv errors
  cv_err_mat <- matrix(NA, nrow = n_folds, ncol = length(penalties))

  # loop through each fold
  cv_err_ls <- lapply(
    seq_len(n_folds),
    function(x){
      rm_idx <- ((x - 1) * percent_rm < rand_mat) &
                  (rand_mat < x * percent_rm)
      target_rm <- target
      target_rm[rm_idx] <- NA
      sapply(
        penalties,
        function(p){
          res <- SPC(target_rm, sumabsvs = p, niter = n_iter, v = v,
                     trace = FALSE, center = FALSE, K = n_eigen)
          target_hat <- as.numeric(res$d) * res$u %*% t(res$v)
          sum(((target_hat - target)[target_rm])^2)
        })
    })
}
