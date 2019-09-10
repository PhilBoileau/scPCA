context("Test coerceMatrix")
library(scPCA)
library(Matrix)

test_that("coerces sparse matrices to standard matrices", {
  expect_equal(
    class(coerceMatrix(data = as(as.matrix(background_df), "dgCMatrix"))),
    "matrix"
  )
  expect_equal(
    class(coerceMatrix(data = as(as.matrix(background_df), "dgeMatrix"))),
    "matrix"
  )
})
