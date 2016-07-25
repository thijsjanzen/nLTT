context("nLTTstat")

test_that("nLTTstat use", {

  p <- ape::rcoal(10)
  q <- ape::rcoal(10)

  expect_equal(
    0.0, nLTTstat(tree1 = p, tree2 = p, distance_method = "abs"), # nolint nLTTstat should be all lowercase, left in for backwards compatibility
    tolerance = 0.0001
  )
  expect_equal(
    0.0, nLTTstat(tree1 = p, tree2 = p, distance_method = "squ"), # nolint nLTTstat should be all lowercase, left in for backwards compatibility
    tolerance = 0.0001
  )

  expect_true(
    0.0 < nLTTstat(tree1 = p, tree2 = q, distance_method = "abs"), # nolint nLTTstat should be all lowercase, left in for backwards compatibility
  )

  expect_equal(
    nLTTstat(tree1 = p, tree2 = q, distance_method = "abs"), # nolint nLTTstat should be all lowercase, left in for backwards compatibility
    nLTTstat(tree1 = p, tree2 = q, distance_method = "abs"),  # nolint nLTTstat should be all lowercase, left in for backwards compatibility
    tolerance = 0.0001
  )

  expect_equal(
    0.0, nLTTstat(tree1 = q, tree2 = q, distance_method = "abs"), # nolint nLTTstat should be all lowercase, left in for backwards compatibility
    tolerance = 0.0001
  )
  expect_equal(
    0.0, nLTTstat(tree1 = q, tree2 = q, distance_method = "squ"), # nolint nLTTstat should be all lowercase, left in for backwards compatibility
    tolerance = 0.0001
  )

  # Cannot test for incorrect distance_method due to Issue #6
})


test_that("nLTTstat abuse", {

  phylo <- ape::rcoal(10)

  expect_error(
    nLTTstat(tree1 = 42, tree2 = phylo, distance_method = "abs"), # nolint nLTTstat should be all lowercase, left in for backwards compatibility
    "nLTTstat: tree1 must be of class 'phylo'"
  )
  expect_error(
    nLTTstat(tree1 = phylo, tree2 = 42, distance_method = "abs"), # nolint nLTTstat should be all lowercase, left in for backwards compatibility
    "nLTTstat: tree2 must be of class 'phylo'"
  )

  # Cannot test for incorrect distance_method due to Issue #6
})
