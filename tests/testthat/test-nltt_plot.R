context("nltt_plot")

test_that("nltt_plot: use", {
  nltt_plot(ape::rcoal(10))
  nltt_lines(ape::rcoal(10))
})

test_that("nltt_plot: abuse", {
  expect_error(
    nltt_plot("not a phylogeny"),
    "nltt_plot: phylogeny must be of class 'phylo'"
  )
  expect_error(
    nltt_lines("not a phylogeny"),
    "nltt_lines: phylogeny must be of class 'phylo'"
  )

})
