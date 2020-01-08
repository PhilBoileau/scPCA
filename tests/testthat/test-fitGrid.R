context("Test routine for fitting grid of contrastive/penalization parameters")
library(BiocParallel)

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
fit <- fitGrid(
  target = target, center = center, scale = scale,
  c_contrasts = c_contrasts, contrasts = contrasts,
  penalties = penalties, n_eigen = n_eigen, n_centers = n_centers,
  clust_method = "kmeans"
)

test_that("The output is a list of length 4", {
  expect_equal(length(fit), 5)
})

test_that("The rotation matrix is ncol(target) by n_eigen", {
  expect_equal(nrow(fit$rotation[[1]]), ncol(target))
  expect_equal(ncol(fit$rotation[[1]]), n_eigen)
})

test_that("The projection matrix is nrow(target) by n_eigen", {
  expect_equal(nrow(fit$x[[2]]), nrow(target))
  expect_equal(ncol(fit$x[[2]]), n_eigen)
})

test_that(paste(
  "The optimal contrast and penalty terms are real-valued",
  "and non-negative"
), {
  idx <- which.max(fit$ave_sil_widths)
  expect_gte(fit$contrast[[idx]], 0)
  expect_gte(fit$penalty[[idx]], 0)
})

test_that("All cluster methods (kmeans, pam, hclust) run without errors", {
  expect_silent(fitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "kmeans"
  ))
  expect_silent(fitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "pam"
  ))
  expect_silent(fitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "hclust",
    linkage_method = "complete"
  ))
  expect_silent(fitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "hclust",
    linkage_method = "average"
  ))
  expect_silent(fitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "hclust",
    linkage_method = "single"
  ))
  expect_silent(fitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "hclust",
    linkage_method = "ward.D"
  ))
  expect_silent(fitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "hclust",
    linkage_method = "ward.D2"
  ))
  expect_silent(fitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "hclust",
    linkage_method = "mcquitty"
  ))
  expect_silent(fitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "hclust",
    linkage_method = "median"
  ))
  expect_silent(fitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen,
    n_centers = n_centers, clust_method = "hclust",
    linkage_method = "centroid"
  ))
})

test_that("The following results are replicable/reproducible", {
  set.seed(12412)

  # using global example
  idx <- which.max(fit$ave_sil_widths)
  expect_equal(fit$contrast[[idx]], 1)
  expect_equal(fit$penalty[[idx]], 0)

  # using new example
  fit_pam <- fitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = c(1, 10, 100),
    penalties = c(0.1, 0.5, 1), n_eigen = n_eigen,
    n_centers = 4, clust_method = "pam"
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
})

test_that("The parallelized analog matches the sequential variant exactly", {
  set.seed(123)
  fit <- fitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen, n_centers = n_centers,
    clust_method = "kmeans"
  )
  register(SerialParam())
  set.seed(123)
  bpfit <- bpFitGrid(
    target = target, center = center, scale = scale,
    c_contrasts = c_contrasts, contrasts = contrasts,
    penalties = penalties, n_eigen = n_eigen, n_centers = n_centers,
    clust_method = "kmeans"
  )
  expect_equal(fit, bpfit)
})
