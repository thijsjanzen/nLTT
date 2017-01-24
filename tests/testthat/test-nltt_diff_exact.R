context("nltt_diff_exact")

test_that("multiplication works", {

  tree1 <- ape::read.tree(text = "(a:1,b:1):1;")
  tree2 <- ape::read.tree(
    text = "(((d:0.000000001,c:0.000000001):1,b:1):0.000000001,a:1.000000001):1;")

  testthat::expect_equal(
    nLTT::nltt_diff(tree1 = tree1, tree2 = tree2),
    0.25, tolerance = 0.0001
  )


  testthat::expect_equal(
    nLTT::nltt_diff_exact(
      tree1 = tree1, tree2 = tree2, ignore_stem = FALSE),
    0.25, tolerance = 0.0001
  )

})
