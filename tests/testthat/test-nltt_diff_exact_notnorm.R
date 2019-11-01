test_that("nltt_diff_exact_notnorm: use", {
  n <- 10
  b_times_n <- (seq(1, n) / n)
  lineages_n <- b_times_n
  b_times2_n <- b_times_n * b_times_n
  lineages2_n <- b_times2_n

  stored <- 0.3465

  expect_equal(
    measured <- nLTT::nltt_diff_exact_notnorm(
      event_times = b_times_n,
      species_number  = lineages_n,
      event_times2 = b_times2_n,
      species_number2 = lineages2_n,
      time_unit = "ago",
      distance_method = "abs"
    ),
    stored
  )
})


test_that("nltt_diff_exact_notnorm: abuse", {

  set.seed(42)
  p <- ape::rcoal(5)
  p$root.edge <- 0.1 # nolint ape variable name
  q <- ape::rcoal(5)
  p$root.edge <- 0.2 # nolint ape variable name

  n <- 10
  b_times_n <- (seq(1, n) / n)
  lineages_n <- b_times_n
  b_times2_n <- b_times_n * b_times_n
  lineages2_n <- b_times2_n

  # Stop if phylo with informative error
  expect_error(
    measured <- nLTT::nltt_diff_exact_notnorm(
      event_times = p,
      species_number  = lineages_n,
      event_times2 = q,
      species_number2 = lineages2_n,
      time_unit = "ago",
      distance_method = "abs"
    ),
    "event times should be a numeric vector of event times, not a
       phylogeny. Use nltt_diff_exact_norm_brts for phylo objects instead."
  )
})


