context("Test routine for coercing inputs to matrix objects")
library(Matrix)

skip_on_travis()
test_that("Coercion of sparse matrices to standard matrices works", {
  expect_equal(
    class(coerceMatrix(data = as(as.matrix(background_df), "dgCMatrix"))),
    "matrix"
  )
  expect_equal(
    class(coerceMatrix(data = as(as.matrix(background_df), "dgeMatrix"))),
    "matrix"
  )
})
