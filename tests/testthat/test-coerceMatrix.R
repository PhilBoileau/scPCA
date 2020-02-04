context("Test routine for coercing inputs to matrix objects")
library(Matrix)

test_that("Coercion of sparse matrices to standard matrices works", {
  skip_on_travis("Bug in Travis CI")
  expect_equal(
    class(coerceMatrix(data = as(as.matrix(background_df), "dgCMatrix"))),
    "matrix"
  )
  expect_equal(
    class(coerceMatrix(data = as(as.matrix(background_df), "dgeMatrix"))),
    "matrix"
  )
})
