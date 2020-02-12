context("nltt_diff_exact_brts")
test_that("nltt_diff_exact_brts use", {
  n <- 10
  b_times <- seq(1, n)
  lineages <- b_times
  b_times2 <- b_times * b_times
  lineages2 <- b_times2

  measured <- expect_silent(
    nLTT::nltt_diff_exact_brts(
      b_times = b_times,
      lineages = lineages,
      b_times2 = b_times2,
      lineages2 = lineages2,
      distance_method = "abs",
      time_unit = "ago"
    )
  )
  expected <- 0.34649999999999986

  expect_equal(measured, expected, tolerance = 0.0001)
})

test_that("nltt_diff_exact_brts abuse", {
  n <- 10
  b_times <- seq(1, n)
  lineages <- b_times
  b_times2 <- b_times * b_times
  lineages2 <- b_times2

  expect_error(
    nLTT::nltt_diff_exact_brts(
      b_times = b_times,
      lineages = lineages,
      b_times2 = b_times2,
      lineages2 = lineages2,
      distance_method = "abs",
      time_unit = "since"
    ),
    regexp = "branching times must be negative,
              for example -3 time units since the present"
  )
  expect_error(
    nLTT::nltt_diff_exact_brts(
      b_times = - b_times,
      lineages = lineages,
      b_times2 = - b_times2,
      lineages2 = lineages2,
      distance_method = "abs",
      time_unit = "ago"
    ),
    regexp = "branching times must be positive, for example 3 time units ago"
  )

  expect_error(
    nLTT::nltt_diff_exact_brts(
      b_times = - b_times,
      lineages = lineages,
      b_times2 = - b_times2,
      lineages2 = lineages2,
      distance_method = "abs",
      time_unit = "nonsense"
    ),
    regexp = "time_unit must be either 'since' or 'ago'"
  )
})
