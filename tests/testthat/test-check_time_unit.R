test_that("use", {
  expect_silent(check_time_unit("since"))
  expect_silent(check_time_unit("ago"))
  expect_error(
    check_time_unit("nonsense"),
    "time_unit must be either 'since' or 'ago'"
  )
})
