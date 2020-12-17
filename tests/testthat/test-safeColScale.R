context("Test routines for safe column scaling")

# initialization of inputs for testing function
set.seed(123)
target <- toy_df[, -31]
dgC_target <- as(as.matrix(toy_df[, -31]), "dgCMatrix")
dm_target <- DelayedArray::DelayedArray(toy_df[, -31])
rownames(dm_target) <- NULL
center <- TRUE
scale <- TRUE

# unit tests
test_that("`scale` and `safeColScale` produce the same output", {
  target_safeColScaled <- safeColScale(target, center, scale)
  dgc_target_safeColScaled <- safeColScale(dgC_target, center, scale)
  dm_target_safeColScaled <- safeColScale(dm_target, center, scale)
  sm_target_safeColScaled <- safeColScale(target, center, scale,
                                          scaled_matrix = TRUE)
  target_scaled <- scale(target, center, scale)

  # remove pesky attributed
  attr(target_scaled, "scaled:center") <- NULL
  attr(target_scaled, "scaled:scale") <- NULL

  expect_equal(target_scaled, target_safeColScaled)
  expect_equal(target_scaled, as.matrix(dgc_target_safeColScaled))
  expect_equal(target_scaled, as.matrix(dm_target_safeColScaled))
  expect_equal(round(target_scaled, 6),
               round(as.matrix(sm_target_safeColScaled), 6))
})

test_that("`safeColScale` avoid NA in its output even when `scale` fails to", {
  # make first column constant
  target[, 1] <- 1

  # check that scale fails to avoid NAs
  target_scaled <- scale(target, center, scale)
  expect_true(sum(colSums(is.na(target_scaled))) > 0)

  # check that safeColScale avoids NAs
  target_safeColScaled <- safeColScale(target, center, scale)
  expect_true(sum(colSums(is.na(target_safeColScaled))) == 0)

  dgc_target_safeColScaled <- safeColScale(dgC_target, center, scale)
  expect_true(sum(Matrix::colSums(is.na(dgc_target_safeColScaled))) == 0)

  dm_target_safeColScaled <- safeColScale(dm_target, center, scale)
  expect_true(sum(DelayedArray::colSums(is.na(dm_target_safeColScaled))) == 0)

  sm_target_safeColScaled <- safeColScale(target, center, scale,
                                          scaled_matrix = TRUE)

  expect_true(sum(colSums(is.na(sm_target_safeColScaled))) == 0)
})

# `safeColScale` is faster than `scale`
# library(ggplot2)
# library(microbenchmark)
# mb_scale <- microbenchmark(safeColScale(target), scale(target),
# times = 20, unit = "s")
# p_mb_scale <- ggplot(data = mb_scale, aes(y = time / 1e9, x = expr)) +
# geom_violin() + theme_grey(base_size = 20) +
# xlab("Method") + ylab("Time (seconds)")
# print(mb_scale)
# print(p_mb_scale)
