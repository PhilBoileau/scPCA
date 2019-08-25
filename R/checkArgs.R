#' Check Arguments passed to the scPCA Function
#'
#' @description Checks whether or not the all arguments in the scPCA functions
#'   are input properly.
#'
#' @param target Target data.
#' @param background Background data.
#' @param center Indicates whether target and background should be centered.
#' @param scale Indicates whether target and background should be scaled.
#' @param num_eigen The number of eigen vectors to compute.
#' @param contrasts The vector of contrasts parameters.
#' @param penalties The vector of penalty terms.
#'
#' @importFrom assertthat assert_that see_if is.count is.flag
#' @importFrom tibble is_tibble
#'
#' @keywords internal
#'
#' @author Philippe Boileau (\email{philippe_boileau@berkeley.edu}), with
#'  contributions from Nima Hejazi (\email{nh@nimahejazi.org}).
#'
#' @return Whether all argument conditions are satisfied
#'
checkArgs <- function(target, background, center, scale,
                      num_eigen, contrasts, penalties) {

  # assert that the target and background data frames are of the right class
  assertthat::assert_that(
    tibble::is_tibble(target) ||
      is.data.frame(target) ||
      is.matrix(target) ||
      is(background, "dgeMatrix") ||
      is(target, "dgCMatrix")
  )
  assertthat::assert_that(
    tibble::is_tibble(background) ||
      is.data.frame(background) ||
      is.matrix(background) ||
      is(background, "dgeMatrix") ||
      is(background, "dgCMatrix")
  )

  # assert that target and background have the same number of variables
  assertthat::assert_that(ncol(target) == ncol(background))

  # check the centering and scaling arguments
  assertthat::assert_that(assertthat::is.flag(center))
  assertthat::assert_that(assertthat::is.flag(scale))

  # check the number of eigenvectors to compute
  assertthat::assert_that(assertthat::is.count(num_eigen))
  assertthat::assert_that(num_eigen <= ncol(target))

  # check the contrastive parameters
  if (assertthat::see_if(!missing(contrasts))) {
    assertthat::assert_that(length(contrasts) > 0)
    assertthat::assert_that(is.numeric(contrasts))
    assertthat::assert_that(all(contrasts > 0))
  }

  # check penalty terms
  if (assertthat::see_if(!is.null(penalties))) {
    assertthat::assert_that(length(penalties) > 0)
    assertthat::assert_that(is.numeric(penalties))
    assertthat::assert_that(all(penalties >= 0))
  }
}
