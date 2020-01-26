test_that("nltt_diff_exact_extinct: use don't normalize", {
  n <- 10
  b_times_n <- (seq(1, n) / n)
  lineages_n <- b_times_n
  b_times2_n <- b_times_n * b_times_n
  lineages2_n <- b_times2_n


  # Normalized data
  stored <- 0.3465
  expect_equal(
    measured <- nLTT::nltt_diff_exact_extinct(
      event_times = b_times_n,
      species_number  = lineages_n,
      event_times2 = b_times2_n,
      species_number2 = lineages2_n,
      time_unit = "ago",
      distance_method = "abs"
    ),
    stored
  )

  # Not normalized data
  set.seed(1)
  event_times <- sort(runif(15, 0, 10))
  event_times2 <- sort(runif(154, 0, 10))

  # Test of non-increasing lineages
  lineages1 <- 1:length(event_times)
  lineages2 <- 1:length(event_times2)
  lineages1[11:15] <- 9:5
  lineages2[133:154] <- 133:112

  # Distance method "abs"

  stored <- 677.391755469143
  expect_equal(
    measured <- nLTT::nltt_diff_exact_extinct(
      event_times = event_times,
      species_number  = lineages1,
      event_times2 = event_times2,
      species_number2 = lineages2,
      time_unit = "ago",
      distance_method = "abs",
      normalize = FALSE
    ),
    stored,
  )

  # Distance method "squ"
  stored <- 62770.448642959818
  expect_equal(
    measured <- nLTT::nltt_diff_exact_extinct(
      event_times = event_times,
      species_number  = lineages1,
      event_times2 = event_times2,
      species_number2 = lineages2,
      time_unit = "ago",
      distance_method = "squ",
      normalize = FALSE
    ),
    stored,
  )
})


test_that("nltt_diff_exact_extinct: use normalize", {
  n <- 10
  b_times_n <- (seq(1, n) / n)
  lineages_n <- b_times_n
  b_times2_n <- b_times_n * b_times_n
  lineages2_n <- b_times2_n


  # Normalized data
  stored <- 0.3465
  expect_equal(
    measured <- nLTT::nltt_diff_exact_extinct(
      event_times = b_times_n,
      species_number  = lineages_n,
      event_times2 = b_times2_n,
      species_number2 = lineages2_n,
      time_unit = "ago",
      distance_method = "abs"
    ),
    stored
  )

  # Not normalized data
  set.seed(1)
  event_times <- sort(runif(15, 0, 10))
  event_times2 <- sort(runif(154, 0, 10))

  # Test of non-increasing lineages
  lineages1 <- 1:length(event_times)
  lineages2 <- 1:length(event_times2)
  lineages1[11:15] <- 9:5
  lineages2[133:154] <- 133:112

  stored <- 0.21989107217693629
  expect_equal(
    measured <- nLTT::nltt_diff_exact_extinct(
      event_times = event_times,
      species_number  = lineages1,
      event_times2 = event_times2,
      species_number2 = lineages2,
      time_unit = "ago",
      distance_method = "abs",
      normalize = TRUE
    ),
    stored,
  )
})


test_that("nltt_diff_exact_extinct: abuse", {

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
    measured <- nLTT::nltt_diff_exact_extinct(
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

  # Stop if weird time unit with informative error
  expect_error(
    measured <- nLTT::nltt_diff_exact_extinct(
      event_times = b_times_n,
      species_number  = lineages_n,
      event_times2 = b_times2_n,
      species_number2 = lineages2_n,
      time_unit = "nonsense",
      distance_method = "abs"
    ),
    "time_unit must be either 'since' or 'ago'"
  )

  # Stop if time unit inconsistent with data and throw informative error
  expect_error(
    measured <- nLTT::nltt_diff_exact_extinct(
      event_times = b_times_n,
      species_number  = lineages_n,
      event_times2 = b_times2_n,
      species_number2 = lineages2_n,
      time_unit = "since",
      distance_method = "abs"
    ),
    "event times must be negative, for example -3 time units since the present"
  )

  expect_error(
    measured <- nLTT::nltt_diff_exact_extinct(
      event_times = b_times_n * -1,
      species_number  = lineages_n,
      event_times2 = b_times2_n * -1,
      species_number2 = lineages2_n,
      time_unit = "ago",
      distance_method = "abs"
    ),
    "event times must be positive, for example 3 time units ago"
  )
})

test_that("nltt_diff_exact_extinct: convergence with old function test", {

  # New function should give same resutls as old function when data is
  # normalized from the start

  n <- 10
  b_times_n <- (seq(1, n) / n)
  lineages_n <- b_times_n
  b_times2_n <- b_times_n * b_times_n
  lineages2_n <- b_times2_n

  expect_equal(
    measured_notnorm <- nLTT::nltt_diff_exact_extinct(
      event_times = b_times_n,
      species_number  = lineages_n,
      event_times2 = b_times2_n,
      species_number2 = lineages2_n,
      time_unit = "ago",
      distance_method = "abs"
    ),
    measured_norm <- nLTT::nltt_diff_exact_brts(
      b_times = b_times_n,
      lineages = lineages_n,
      b_times2 = b_times2_n,
      lineages2 = lineages2_n,
      time_unit = "ago",
      distance_method = "abs"
    )
  )
})
