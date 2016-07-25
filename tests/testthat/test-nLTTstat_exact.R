context("nLTTstat_exact")

test_that("nLTTstat_exact abuse", {

  phylo <- ape::rcoal(10)

  expect_error(
    nLTTstat_exact(tree1 = 42, tree2 = phylo, distance_method = "abs"),
    "nLTTstat_exact: tree1 must be of class 'phylo'"
  )
  expect_error(
    nLTTstat_exact(tree1 = phylo, tree2 = 42, distance_method = "abs"),
    "nLTTstat_exact: tree2 must be of class 'phylo'"
  )

  # Cannot test for incorrect distance_method due to Issue #6
})
