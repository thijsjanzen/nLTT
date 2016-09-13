context("nLTTstat_exact")

test_that("nLTTstat_exact use", {

  p <- ape::rcoal(10)
  q <- ape::rcoal(10)

  expect_equal(
    0.0, nLTTstat_exact(tree1 = p, tree2 = p, distance_method = "abs"), # nolint nLTTstat_exact should be all lowercase, left in for backwards compatibility
    tolerance = 0.0001
  )
  expect_equal(
    0.0, nLTTstat_exact(tree1 = p, tree2 = p, distance_method = "squ"), # nolint nLTTstat_exact should be all lowercase, left in for backwards compatibility
    tolerance = 0.0001
  )

  expect_true(
    0.0 < nLTTstat_exact(tree1 = p, tree2 = q, distance_method = "abs"), # nolint nLTTstat_exact should be all lowercase, left in for backwards compatibility
  )

  expect_equal(
    nLTTstat_exact(tree1 = p, tree2 = q, distance_method = "abs"), # nolint nLTTstat_exact should be all lowercase, left in for backwards compatibility
    nLTTstat_exact(tree1 = p, tree2 = q, distance_method = "abs"),  # nolint nLTTstat_exact should be all lowercase, left in for backwards compatibility
    tolerance = 0.0001
  )

  expect_equal(
    0.0, nLTTstat_exact(tree1 = q, tree2 = q, distance_method = "abs"), # nolint nLTTstat_exact should be all lowercase, left in for backwards compatibility
    tolerance = 0.0001
  )
  expect_equal(
    0.0, nLTTstat_exact(tree1 = q, tree2 = q, distance_method = "squ"), # nolint nLTTstat_exact should be all lowercase, left in for backwards compatibility
    tolerance = 0.0001
  )

})


test_that("nLTTstat_exact abuse", {

  phylo <- ape::rcoal(10)

  expect_error(
    nLTTstat_exact(tree1 = 42, tree2 = phylo, distance_method = "abs"), # nolint nLTTstat_exact should be all lowercase, left in for backwards compatibility
    "nLTTstat_exact: tree1 must be of class 'phylo'"
  )
  expect_error(
    nLTTstat_exact(tree1 = phylo, tree2 = 42, distance_method = "abs"), # nolint nLTTstat_exact should be all lowercase, left in for backwards compatibility
    "nLTTstat_exact: tree2 must be of class 'phylo'"
  )
  expect_error(
    nLTTstat_exact(tree1 = phylo, tree2 = phylo, distance_method = "nonsense"), # nolint nLTTstat_exact should be all lowercase, left in for backwards compatibility
    "nLTTstat_exact: distance method unknown"
  )

  expect_error(
    nLTTstat_exact(tree1 = phylo, tree2 = phylo, ignore_stem = "no logical"), # nolint nLTTstat_exact should be all lowercase, left in for backwards compatibility
    "nLTTstat_exact: ignore_stem must be logical"
  )
})

test_that("nLTTstat_exact may ignore the stem", {

  tree1 <- ape::rcoal(4)
  tree2 <- tree1
  tree2$root.edge <- 5
  expect_true(nLTT::nLTTstat_exact(p, q, ignore_stem = TRUE) < 0.01)

  expect_true(nLTT::nLTTstat_exact(p, q, ignore_stem = FALSE) > 0.4)
})
