context("nltts_plot")

test_that("use", {

  set.seed(42)
  phylos <- c(ape::rcoal(10))
  testthat::expect_silent(nltts_plot(phylos))
})

test_that("use, #38", {

  skip("Expose #38")

  # Works for 10 trees
  set.seed(42)
  phylos <- c(ape::rcoal(10, n = 10))
  testthat::expect_silent(nltts_plot(phylos))

  # Fails for 11 trees
  set.seed(42)
  phylos <- c(ape::rcoal(10, n = 11))
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
