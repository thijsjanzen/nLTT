context("nltts_plot")

test_that("use", {

  set.seed(42)
  n_tips <- 10
  phylos <- c(ape::rcoal(n = n_tips))
  testthat::expect_silent(nltts_plot(phylos))
})

test_that("use, #38", {

  # Works for 10 trees
  set.seed(42)
  n_tips <- 10
  n_trees <- 10
  phylos <- c(ape::rcoal(n = n_tips, N = n_trees))
  testthat::expect_silent(nltts_plot(phylos))

  # Fails for 11 trees
  set.seed(42)
  phylos <- c(ape::rcoal(n = n_tips, N = n_trees))
  testthat::expect_silent(nltts_plot(phylogenies = phylos))
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
