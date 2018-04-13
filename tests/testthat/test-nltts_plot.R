context("nltts_plot")

test_that("use", {

  phylos <- c(ape::rcoal(10))
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
