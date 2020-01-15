context("get_norm_brts")
test_that("use", {
  phylogeny <- ape::read.tree(text = "((a:2,b:2):1,c:3);")
  phylogeny$root.edge <- 2 # nolint
  testthat::expect_true(
    all(nLTT::get_norm_brts(phylogeny) == c(0, 0.66666666666666663, 1.0, 1.0))
  )
  testthat::expect_silent(
    nLTT::get_norm_brts(phylogeny)
  )
})
