cPCA <- function(target, center, scale, c_contrasts, contrasts, penalties,
                 n_eigen, num_medoids = 8){

  # preliminaries
  num_contrasts <- length(contrasts)

  # for each contrasted covariance matrix, compute the eigenvectors
  loadings_mat <- lappy(
    seq_len(num_contrasts),
    function(x) {
      res <- eigen(c_contrasts[[x]],
             symmetric = TRUE)$vectors[, seq_len(n_eigen)]
     }
   )

  # for each loadings matrix, project target onto constrastive subspace
  spaces <- lapply(
    seq_len(num_contrasts),
    function(x) {
      as.matrix(target) %*% loadings_mat[[x]]
    }
  )

  # produce the QR decomposition of these projections, extract Q
  qr_decomps <- lapply(
    seq_len(num_contrasts),
    function(x) {
      qr.Q(qr(spaces[[x]]))
    }
  )

  # populate affinity matrix for spectral clustering using the principal angles
  aff_vect <- sapply(
    seq_len(num_contrasts),
    function(i) {
      sapply(
        seq(from = i + 1, to = num_contrasts),
        function(j){
          Q_i <- qr_decomps[[i]]
          Q_j <- qr_decomps[[j]]
          d <- svd(x = t(Q_i)%*%Q_j, nu = 0, nv = 0)$d
          d[1]*d[2]
          }
      )
    }
  )

  aff_mat <- diag(x = 0.5, nrow = num_contrasts)
  aff_mat[lower.tri(aff_mat, diag = FALSE)] <- aff_vect
  aff_mat <- t(aff_mat)

  # fix any computation errors, see numpy.nan_to_num
  aff_mat[is.nan(aff_mat)] <- 0
  aff_mat[is.na(aff_mat)] <- 0
  aff_mat[is.infinite(aff_mat)] <- 1000
  aff_mat <- t(aff_mat) + aff_mat

  # perfrom spectral clustering using the affinity matrix
  spec_clust <- kernlab::specc(kernlab::as.kernelMatrix(aff_mat),
                               centers = num_medoids)

  # identify the alpha medoids of the spectral clustering
  contrast_medoids <- sapply(
    seq_len(num_medoids),
    function(x) {
      sub_index <- which(spec_clust == x)
      sub_aff_mat <- aff_mat[sub_index, sub_index]
      aff_sums <- colSums(sub_aff_mat)
      contrasts[sub_index[which.max(aff_sums)]]
    }
  )

  # create the lists of contrastive parameter medoids, loadings and projections
  med_index <- which(contrasts %in% contrast_medoids)
  list(
    rotation = loadings_mat[med_index],
    x = spaces[med_index],
    contrast = med_index,
    penalty = rep(0, length(med_index))
  )
}
