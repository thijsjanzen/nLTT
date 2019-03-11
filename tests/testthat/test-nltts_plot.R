context("nltts_plot")

test_that("use", {
  # temporary fix to keep R-devel happy.
  # should be updated upon release of version 3.6

  suppressWarnings(RNGversion("3.5.0"))

  set.seed(42)
  n_tips <- 10
  phylos <- c(ape::rcoal(n = n_tips))
  testthat::expect_silent(nltts_plot(phylos))
})

test_that("abuse", {

  testthat::expect_error(
    nltts_plot(NULL),
    "there must be at least one phylogeny supplied"
  )
  testthat::expect_error(
    nltts_plot("nonsense"),
    "phylogenies must be of class 'multiPhylo' or 'list'"
  )

})
