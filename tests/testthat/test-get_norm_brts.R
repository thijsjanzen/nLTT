context("get_norm_brts")
test_that("use", {
  phylogeny <- ape::read.tree(text = "((a:2,b:2):1,c:3);")
  phylogeny$root.edge <- 2 # nolint

  norm_brts <- as.vector(nLTT::get_norm_brts(phylogeny))
  ref_brts  <- as.vector(c(0, 0.66666666666666667, 1.0, 1.0))

  testthat::expect_true(all.equal(norm_brts, ref_brts))

  testthat::expect_silent(
    nLTT::get_norm_brts(phylogeny)
  )
})
