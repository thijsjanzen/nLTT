context("nltt_diff_exact")

test_that("multiplication works", {

  phylogeny_1 <- ape::read.tree(text = "(a:1,b:1):1;")
  phylogeny_2 <- ape::read.tree(
  text = "(((d:0.000000001,c:0.000000001):1,b:1):0.000000001,a:1.0):1;")


  testthat::expect_equal(
    nLTT::nltt_diff_exact(phylogeny_1, phylogeny_2),
    0.25, tolerance = 0.0001
  )

})
