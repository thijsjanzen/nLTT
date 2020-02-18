test_that("use", {

  # Use 'c' to convert one phylo to multiPhylo
  expect_silent(
    check_phylogenies(c(ape::rcoal(2)))
  )
  # Use list of phylo's
  expect_silent(
    check_phylogenies(list(ape::rcoal(2), ape::rcoal(2)))
  )
  expect_error(
    check_phylogenies(c()),
    "there must be at least one phylogeny supplied"
  )
  expect_error(
    check_phylogenies(ape::rcoal(2)),
    "phylogenies must be of class 'multiPhylo' or 'list'"
  )
  expect_error(
    check_phylogenies("nonsense"),
    "phylogenies must be of class 'multiPhylo' or 'list'"
  )
  expect_error(check_phylogenies(NA))
  expect_error(check_phylogenies(Inf))
  expect_error(check_phylogenies(NULL))
})


