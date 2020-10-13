context("Test routine for fitting grid of contrastive/penalization parameters")
library(BiocParallel)
library(Matrix)
library(DelayedArray)

# initialization of inputs for testing function
set.seed(123)
target <- toy_df[, -31]
background <- background_df
center <- TRUE
scale <- TRUE
contrasts <- c(0, 1, 10, 100)
c_contrasts <- contrastiveCov(target, background, contrasts, TRUE, TRUE)
penalties <- c(0, 0.1, 1)
n_eigen <- 2
n_centers <- 4
fit_iterative <- fitGrid(
  target = target, center = center, scale = scale,
  c_contrasts = c_contrasts, contrasts = contrasts,
  penalties = penalties, n_eigen = n_eigen, n_centers = n_centers,
  clust_method = "kmeans", alg = "iterative"
)
fit_var_proj <- fitGrid(
  target = target, center = center, scale = scale,
  c_contrasts = c_contrasts, contrasts = contrasts,
  penalties = penalties, n_eigen = n_eigen, n_centers = n_centers,
  clust_method = "kmeans", alg = "var_proj"
)
set.seed(163)
fit_rand_var_proj <- fitGrid(
  target = target, center = center, scale = scale,
  c_contrasts = c_contrasts, contrasts = contrasts,
  penalties = penalties, n_eigen = n_eigen, n_centers = n_centers,
  clust_method = "kmeans", alg = "rand_var_proj"
)

# sparse versions
dgC_target <- as(as.matrix(toy_df[, -31]), "dgCMatrix")
dgC_background <- as(as.matrix(background_df), "dgCMatrix")
c_contrasts <- contrastiveCov(dgC_target, dgC_background, contrasts, TRUE, TRUE)
dgC_fit_iterative <- fitGrid(
  target = dgC_target, center = center, scale = scale,
  c_contrasts = c_contrasts, contrasts = contrasts,
  penalties = penalties, n_eigen = n_eigen, n_centers = n_centers,
  clust_method = "kmeans", alg = "iterative"
)
dgC_fit_var_proj <- fitGrid(
  target = dgC_target, center = center, scale = scale,
  c_contrasts = c_contrasts, contrasts = contrasts,
  penalties = penalties, n_eigen = n_eigen, n_centers = n_centers,
  clust_method = "kmeans", alg = "var_proj"
)
set.seed(163)
dgC_fit_rand_var_proj <- fitGrid(
  target = dgC_target, center = center, scale = scale,
  c_contrasts = c_contrasts, contrasts = contrasts,
  penalties = penalties, n_eigen = n_eigen, n_centers = n_centers,
  clust_method = "kmeans", alg = "rand_var_proj"
)

# DelayedMatrix versions
dm_target <- DelayedArray(toy_df[, -31])
dm_background <- DelayedArray(background_df)
c_contrasts <- contrastiveCov(dm_target, dm_background, contrasts, TRUE, TRUE)
dm_fit_iterative <- fitGrid(
  target = dm_target, center = center, scale = scale,
  c_contrasts = c_contrasts, contrasts = contrasts,
  penalties = penalties, n_eigen = n_eigen, n_centers = n_centers,
  clust_method = "kmeans", alg = "iterative"
)
dm_fit_var_proj <- fitGrid(
  target = dm_target, center = center, scale = scale,
  c_contrasts = c_contrasts, contrasts = contrasts,
  penalties = penalties, n_eigen = n_eigen, n_centers = n_centers,
  clust_method = "kmeans", alg = "var_proj"
)
set.seed(163)
dm_fit_rand_var_proj <- fitGrid(
  target = dm_target, center = center, scale = scale,
  c_contrasts = c_contrasts, contrasts = contrasts,
  penalties = penalties, n_eigen = n_eigen, n_centers = n_centers,
  clust_method = "kmeans", alg = "rand_var_proj"
)

test_that("The output is a list of length 4", {
  # df
  expect_equal(length(fit_iterative), 5)
  expect_equal(length(fit_var_proj), 5)
  expect_equal(length(fit_rand_var_proj), 5)
  
  # sparse
  expect_equal(length(dgC_fit_iterative), 5)
  expect_equal(length(dgC_fit_var_proj), 5)
  expect_equal(length(dgC_fit_rand_var_proj), 5)
  
  # DelayedMatrices
  expect_equal(length(dm_fit_iterative), 5)
  expect_equal(length(dm_fit_var_proj), 5)
  expect_equal(length(dm_fit_rand_var_proj), 5)
})

test_that("The rotation matrix is ncol(target) by n_eigen", {
  
  # df
  expect_equal(nrow(fit_iterative$rotation[[1]]), ncol(target))
  expect_equal(ncol(fit_iterative$rotation[[1]]), n_eigen)
  expect_equal(nrow(fit_var_proj$rotation[[1]]), ncol(target))
  expect_equal(ncol(fit_var_proj$rotation[[1]]), n_eigen)
  expect_equal(nrow(fit_rand_var_proj$rotation[[1]]), ncol(target))
  expect_equal(ncol(fit_rand_var_proj$rotation[[1]]), n_eigen)
  
  # sparse
  expect_equal(nrow(dgC_fit_iterative$rotation[[1]]), ncol(target))
  expect_equal(ncol(dgC_fit_iterative$rotation[[1]]), n_eigen)
  expect_equal(nrow(dgC_fit_var_proj$rotation[[1]]), ncol(target))
  expect_equal(ncol(dgC_fit_var_proj$rotation[[1]]), n_eigen)
  expect_equal(nrow(dgC_fit_rand_var_proj$rotation[[1]]), ncol(target))
  expect_equal(ncol(dgC_fit_rand_var_proj$rotation[[1]]), n_eigen)
  
  # DelayedMatrix
  expect_equal(nrow(dm_fit_iterative$rotation[[1]]), ncol(target))
  expect_equal(ncol(dm_fit_iterative$rotation[[1]]), n_eigen)
  expect_equal(nrow(dm_fit_var_proj$rotation[[1]]), ncol(target))
  expect_equal(ncol(dm_fit_var_proj$rotation[[1]]), n_eigen)
  expect_equal(nrow(dm_fit_rand_var_proj$rotation[[1]]), ncol(target))
  expect_equal(ncol(dm_fit_rand_var_proj$rotation[[1]]), n_eigen)
})

test_that("The projection matrix is nrow(target) by n_eigen", {
  
  # df
  expect_equal(nrow(fit_iterative$x[[2]]), nrow(target))
  expect_equal(ncol(fit_iterative$x[[2]]), n_eigen)
  expect_equal(nrow(fit_var_proj$x[[2]]), nrow(target))
  expect_equal(ncol(fit_var_proj$x[[2]]), n_eigen)
  expect_equal(nrow(fit_rand_var_proj$x[[2]]), nrow(target))
  expect_equal(ncol(fit_rand_var_proj$x[[2]]), n_eigen)
  
  # sparse
  expect_equal(nrow(dgC_fit_iterative$x[[2]]), nrow(target))
  expect_equal(ncol(dgC_fit_iterative$x[[2]]), n_eigen)
  expect_equal(nrow(dgC_fit_var_proj$x[[2]]), nrow(target))
  expect_equal(ncol(dgC_fit_var_proj$x[[2]]), n_eigen)
  expect_equal(nrow(dgC_fit_rand_var_proj$x[[2]]), nrow(target))
  expect_equal(ncol(dgC_fit_rand_var_proj$x[[2]]), n_eigen)
  
  # DelayedMatrix
  expect_equal(nrow(dm_fit_iterative$x[[2]]), nrow(target))
  expect_equal(ncol(dm_fit_iterative$x[[2]]), n_eigen)
  expect_equal(nrow(dm_fit_var_proj$x[[2]]), nrow(target))
  expect_equal(ncol(dm_fit_var_proj$x[[2]]), n_eigen)
  expect_equal(nrow(dm_fit_rand_var_proj$x[[2]]), nrow(target))
  expect_equal(ncol(dm_fit_rand_var_proj$x[[2]]), n_eigen)
})

test_that(paste(
  "The optimal contrast and penalty terms are real-valued",
  "and non-negative"
), {
  
  # df
  idx <- which.max(fit_iterative$ave_sil_widths)
  expect_gte(fit_iterative$contrast[[idx]], 0)
  expect_gte(fit_iterative$penalty[[idx]], 0)
  idx <- which.max(fit_var_proj$ave_sil_widths)
  expect_gte(fit_var_proj$contrast[[idx]], 0)
  expect_gte(fit_var_proj$penalty[[idx]], 0)
  idx <- which.max(fit_rand_var_proj$ave_sil_widths)
  expect_gte(fit_rand_var_proj$contrast[[idx]], 0)
  expect_gte(fit_rand_var_proj$penalty[[idx]], 0)
  
  # sparse
  idx <- which.max(dgC_fit_iterative$ave_sil_widths)
  expect_gte(dgC_fit_iterative$contrast[[idx]], 0)
  expect_gte(dgC_fit_iterative$penalty[[idx]], 0)
  idx <- which.max(dgC_fit_var_proj$ave_sil_widths)
  expect_gte(dgC_fit_var_proj$contrast[[idx]], 0)
  expect_gte(dgC_fit_var_proj$penalty[[idx]], 0)
  idx <- which.max(dgC_fit_rand_var_proj$ave_sil_widths)
  expect_gte(dgC_fit_rand_var_proj$contrast[[idx]], 0)
  expect_gte(dgC_fit_rand_var_proj$penalty[[idx]], 0)
  
  # DelayedMatrix
  idx <- which.max(dm_fit_iterative$ave_sil_widths)
  expect_gte(dm_fit_iterative$contrast[[idx]], 0)
  expect_gte(dm_fit_iterative$penalty[[idx]], 0)
  idx <- which.max(dm_fit_var_proj$ave_sil_widths)
  expect_gte(dm_fit_var_proj$contrast[[idx]], 0)
  expect_gte(dm_fit_var_proj$penalty[[idx]], 0)
  idx <- which.max(dm_fit_rand_var_proj$ave_sil_widths)
  expect_gte(dm_fit_rand_var_proj$contrast[[idx]], 0)
  expect_gte(dm_fit_rand_var_proj$penalty[[idx]], 0)
})

test_that("All cluster methods (kmeans, pam, hclust) run without errors", {
  expect_silent(fitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "kmeans",
    alg = "iterative"
  ))
  expect_silent(fitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "pam",
    alg = "iterative"
  ))
  expect_silent(fitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "hclust",
    linkage_method = "complete", alg = "iterative"
  ))
  expect_silent(fitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "hclust",
    linkage_method = "average", alg = "iterative"
  ))
  expect_silent(fitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "hclust",
    linkage_method = "single", alg = "iterative"
  ))
  expect_silent(fitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "hclust",
    linkage_method = "ward.D", alg = "iterative"
  ))
  expect_silent(fitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "hclust",
    linkage_method = "mcquitty", alg = "iterative"
  ))
  expect_silent(fitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "hclust",
    linkage_method = "median", alg = "iterative"
  ))
  expect_silent(fitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "hclust",
    linkage_method = "centroid", alg = "iterative"
  ))
  expect_silent(fitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "kmeans",
    alg = "var_proj"
  ))
  expect_silent(fitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "pam",
    alg = "var_proj"
  ))
  expect_silent(fitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "hclust",
    linkage_method = "complete", alg = "var_proj"
  ))
  expect_silent(fitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "hclust",
    linkage_method = "average", alg = "var_proj"
  ))
  expect_silent(fitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "hclust",
    linkage_method = "single", alg = "var_proj"
  ))
  expect_silent(fitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "hclust",
    linkage_method = "ward.D", alg = "var_proj"
  ))
  expect_silent(fitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "hclust",
    linkage_method = "mcquitty", alg = "var_proj"
  ))
  expect_silent(fitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "hclust",
    linkage_method = "median", alg = "var_proj"
  ))
  expect_silent(fitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "hclust",
    linkage_method = "centroid", alg = "var_proj"
  ))
  expect_silent(fitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "kmeans",
    alg = "rand_var_proj"
  ))
  expect_silent(fitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "pam",
    alg = "rand_var_proj"
  ))
  expect_silent(fitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "hclust",
    linkage_method = "complete", alg = "rand_var_proj"
  ))
  expect_silent(fitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "hclust",
    linkage_method = "average", alg = "rand_var_proj"
  ))
  expect_silent(fitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "hclust",
    linkage_method = "single", alg = "rand_var_proj"
  ))
  expect_silent(fitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "hclust",
    linkage_method = "ward.D", alg = "rand_var_proj"
  ))
  expect_silent(fitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "hclust",
    linkage_method = "mcquitty", alg = "rand_var_proj"
  ))
  expect_silent(fitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "hclust",
    linkage_method = "median", alg = "rand_var_proj"
  ))
  expect_silent(fitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "hclust",
    linkage_method = "centroid", alg = "rand_var_proj"
  ))
  register(SerialParam())
  expect_silent(bpFitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "kmeans",
    alg = "iterative"
  ))
  expect_silent(bpFitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "pam",
    alg = "iterative"
  ))
  expect_silent(bpFitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "hclust",
    linkage_method = "complete", alg = "iterative"
  ))
  expect_silent(bpFitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "hclust",
    linkage_method = "average", alg = "iterative"
  ))
  expect_silent(bpFitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "hclust",
    linkage_method = "single", alg = "iterative"
  ))
  expect_silent(bpFitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "hclust",
    linkage_method = "ward.D", alg = "iterative"
  ))
  expect_silent(bpFitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "hclust",
    linkage_method = "mcquitty", alg = "iterative"
  ))
  expect_silent(bpFitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "hclust",
    linkage_method = "median", alg = "iterative"
  ))
  expect_silent(bpFitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "hclust",
    linkage_method = "centroid", alg = "iterative"
  ))
  expect_silent(bpFitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "kmeans",
    alg = "var_proj"
  ))
  expect_silent(bpFitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "pam",
    alg = "var_proj"
  ))
  expect_silent(bpFitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "hclust",
    linkage_method = "complete", alg = "var_proj"
  ))
  expect_silent(bpFitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "hclust",
    linkage_method = "average", alg = "var_proj"
  ))
  expect_silent(bpFitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "hclust",
    linkage_method = "single", alg = "var_proj"
  ))
  expect_silent(bpFitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "hclust",
    linkage_method = "ward.D", alg = "var_proj"
  ))
  expect_silent(bpFitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "hclust",
    linkage_method = "mcquitty", alg = "var_proj"
  ))
  expect_silent(bpFitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "hclust",
    linkage_method = "median", alg = "var_proj"
  ))
  expect_silent(bpFitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "hclust",
    linkage_method = "centroid", alg = "var_proj"
  ))
  expect_silent(bpFitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "kmeans",
    alg = "rand_var_proj"
  ))
  expect_silent(bpFitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "pam",
    alg = "rand_var_proj"
  ))
  expect_silent(bpFitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "hclust",
    linkage_method = "complete", alg = "rand_var_proj"
  ))
  expect_silent(bpFitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "hclust",
    linkage_method = "average", alg = "rand_var_proj"
  ))
  expect_silent(bpFitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "hclust",
    linkage_method = "single", alg = "rand_var_proj"
  ))
  expect_silent(bpFitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "hclust",
    linkage_method = "ward.D", alg = "rand_var_proj"
  ))
  expect_silent(bpFitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "hclust",
    linkage_method = "mcquitty", alg = "rand_var_proj"
  ))
  expect_silent(bpFitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "hclust",
    linkage_method = "median", alg = "rand_var_proj"
  ))
  expect_silent(bpFitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "hclust",
    linkage_method = "centroid", alg = "rand_var_proj"
  ))
})

test_that("The following results are replicable/reproducible", {
  set.seed(12412)

  # using global example
  idx <- which.max(fit_iterative$ave_sil_widths)
  expect_equal(fit_iterative$contrast[[idx]], 1)
  expect_equal(fit_iterative$penalty[[idx]], 0)
  expect_equal(dgC_fit_iterative$contrast[[idx]], 1)
  expect_equal(dgC_fit_iterative$penalty[[idx]], 0)
  expect_equal(dm_fit_iterative$contrast[[idx]], 1)
  expect_equal(dm_fit_iterative$penalty[[idx]], 0)

  # using new example
  fit_pam <- fitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = c(1, 10, 100),
    penalties = c(0.1, 0.5, 1), n_eigen = n_eigen,
    n_centers = 4, clust_method = "pam", alg = "iterative"
  )
  dgC_fit_pam <- fitGrid(
    target = dgC_target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = c(1, 10, 100),
    penalties = c(0.1, 0.5, 1), n_eigen = n_eigen,
    n_centers = 4, clust_method = "pam", alg = "iterative"
  )
  dm_fit_pam <- fitGrid(
    target = dm_target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = c(1, 10, 100),
    penalties = c(0.1, 0.5, 1), n_eigen = n_eigen,
    n_centers = 4, clust_method = "pam", alg = "iterative"
  )

  idx <- which.max(fit_pam$ave_sil_widths)
  expect_equal(fit_pam$contrast[[idx]], 1)
  expect_equal(fit_pam$penalty[[idx]], 0.1)
  expect_equal(
    as.numeric(round(fit_pam$rotation[[idx]][11, ], 5)),
    c(-0.31483, 0)
  )
  expect_equal(
    as.numeric(round(fit_pam$x[[idx]][1, ], 5)),
    c(2.46554, 2.96500)
  )
  
  dgC_idx <- which.max(dgC_fit_pam$ave_sil_widths)
  expect_equal(dgC_fit_pam$contrast[[idx]], 1)
  expect_equal(dgC_fit_pam$penalty[[idx]], 0.1)
  expect_equal(
    as.numeric(round(dgC_fit_pam$rotation[[idx]][11, ], 5)),
    c(-0.31483, 0)
  )
  expect_equal(
    as.numeric(round(dgC_fit_pam$x[[idx]][1, ], 5)),
    c(2.46554, 2.96500)
  )
  
  dm_idx <- which.max(dm_fit_pam$ave_sil_widths)
  expect_equal(dm_fit_pam$contrast[[idx]], 1)
  expect_equal(dm_fit_pam$penalty[[idx]], 0.1)
  expect_equal(
    as.numeric(round(dm_fit_pam$rotation[[idx]][11, ], 5)),
    c(-0.31483, 0)
  )
  expect_equal(
    as.numeric(round(dm_fit_pam$x[[idx]][1, ], 5)),
    c(2.46554, 2.96500)
  )
})

test_that("The grid can be fit using user-defined clusterings", {
  expect_silent(fitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen, n_centers = n_centers,
    clust_method = "kmeans", alg = "iterative", clusters = toy_df[, 31]))
  expect_silent(bpFitGrid(
    target = dm_target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen, n_centers = n_centers,
    clust_method = "kmeans", alg = "iterative", clusters = toy_df[, 31]))
})

test_that("The parallelized analog matches the sequential variant exactly", {
  set.seed(123)
  fit <- fitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen, n_centers = n_centers,
    clust_method = "kmeans", alg = "iterative"
  )
  register(SerialParam())
  set.seed(123)
  bpfit <- bpFitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen, n_centers = n_centers,
    clust_method = "kmeans", alg = "iterative"
  )
  expect_equal(fit, bpfit)
  set.seed(432)
  fit <- fitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen, n_centers = n_centers,
    clust_method = "kmeans", alg = "var_proj"
  )
  set.seed(432)
  bpfit <- bpFitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen, n_centers = n_centers,
    clust_method = "kmeans", alg = "var_proj"
  )
  expect_equal(fit, bpfit)
  set.seed(432)
  fit <- fitGrid(
    target = dgC_target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen, n_centers = n_centers,
    clust_method = "kmeans", alg = "var_proj"
  )
  set.seed(432)
  bpfit <- bpFitGrid(
    target = dgC_target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen, n_centers = n_centers,
    clust_method = "kmeans", alg = "var_proj"
  )
  expect_equal(fit, bpfit)
  
  set.seed(432)
  fit <- fitGrid(
    target = dm_target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen, n_centers = n_centers,
    clust_method = "kmeans", alg = "var_proj"
  )
  set.seed(432)
  bpfit <- bpFitGrid(
    target = dm_target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen, n_centers = n_centers,
    clust_method = "kmeans", alg = "var_proj"
  )
  expect_equal(fit, bpfit)
})