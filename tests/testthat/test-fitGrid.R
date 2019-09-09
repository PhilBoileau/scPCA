context("Test fitGrid")
library(scPCA)

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
fit <- fitGrid(target, center, scale, c_contrasts, contrasts, penalties,
  n_eigen, n_centers,
  clust_method = "kmeans"
)

test_that("the output is a list of length 4", {
  expect_equal(length(fit), 4)
})

test_that("the rotation matrix is ncol(target) by n_eigen", {
  expect_equal(nrow(fit$rotation), ncol(target))
  expect_equal(ncol(fit$rotation), n_eigen)
})

test_that("the projection matrix is nrow(target) by n_eigen", {
  expect_equal(nrow(fit$x), nrow(target))
  expect_equal(ncol(fit$x), n_eigen)
})

test_that("the optimal contrast and penalty terms are real valued and non-neg", {
  expect_gte(fit$contrast, 0)
  expect_gte(fit$penalty, 0)
})

test_that("both cluster methods run without producing errors", {
  expect_silent(fitGrid(target, center, scale, c_contrasts, contrasts,
    penalties, n_eigen, n_centers,
    clust_method = "kmeans"
  ))
  expect_silent(fitGrid(target, center, scale, c_contrasts, contrasts,
    penalties, n_eigen, n_centers,
    clust_method = "pam"
  ))
})

test_that("the following results are repeatable", {
  set.seed(12412)

  # using global example
  expect_equal(fit$contrast, 1)
  expect_equal(fit$penalty, 0)

  # using new example
  fit_pam <- fitGrid(target, center, scale, c_contrasts, c(1, 10, 100),
    c(0.1, 0.5, 1), n_eigen,
    n_centers = 4,
    clust_method = "pam"
  )

  expect_equal(fit_pam$contrast, 1)
  expect_equal(fit_pam$penalty, 0.1)
  expect_equal(as.numeric(round(fit_pam$rotation[11, ], 5)), c(-0.31483, 0))
  expect_equal(as.numeric(round(fit_pam$x[1, ], 5)), c(2.46554, 2.96500))
})
