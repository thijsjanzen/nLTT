test_that("use", {

  expect_silent(check_step_type("lower"))
  expect_silent(check_step_type("upper"))
  expect_error(check_step_type("nonsense"))
  expect_error(check_step_type(NA))
  expect_error(check_step_type(NULL))
  expect_error(check_step_type(Inf))
  expect_error(check_step_type(c()))
  expect_error(check_step_type(c("lower", "lower")))
})
