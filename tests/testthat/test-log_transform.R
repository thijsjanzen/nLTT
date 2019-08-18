context("log transformed")

test_that("Identical trees have an nLTTstat of zero,
          also when log transformed", {
  set.seed(314)
  p <- ape::rcoal(10)
  set.seed(314) # Same seed to generate the same tree
  q <- ape::rcoal(10)

  expect_equal(
    0.0, nLTTstat(tree1 = p, tree2 = p, distance_method = "abs",
                  log_transform = TRUE), # nolint nLTTstat should be all lowercase, left in for backwards compatibility
    tolerance = 0.0001
  )
  expect_equal(
    0.0, nLTTstat(tree1 = p, tree2 = p, distance_method = "squ",
                  log_transform = TRUE), # nolint nLTTstat should be all lowercase, left in for backwards compatibility
    tolerance = 0.0001
  )
})

test_that("compare with non-log transformed", {
  set.seed(42)
  p <- ape::rcoal(n = 100)
  q <- ape::rcoal(n = 100)

  stat_1 <- nLTTstat(tree1 = p, tree2 = q, distance_method = "abs",
                         log_transform = FALSE)

  stat_2 <- nLTTstat(tree1 = p, tree2 = q, distance_method = "abs",
                         log_transform = TRUE)

  stat_1
  stat_2

  # for some reason, log transforming richness leads to LOWER
  # nLTT values
  expect_lt(stat_1, stat_2)
})

# this function might benefit from additional tests.
