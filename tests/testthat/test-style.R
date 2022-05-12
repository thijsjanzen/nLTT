test_that("Package style", {
  testthat::skip_if_not_installed("lintr")
  lintr::expect_lint_free()
})
