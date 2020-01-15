test_that("get_norm_n use", {
  phylogeny <- ape::read.tree(text = "((a:2,b:2):1,c:3);")
  phylogeny$root.edge <- 2 # nolint ape variable name
  expected <- c(0.0, 0.5, 1.0, 1.0)
  measured <- expect_silent(
    nLTT::get_norm_n(phylogeny = phylogeny)
  )
  expect_equal(expected, measured)
})
