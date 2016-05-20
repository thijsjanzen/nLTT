context("get_average_nltt")

test_that("get_average_nltt: create some plots", {
  # The inner workings of get_average_nltt are done by get_average_nltt_matrix
  newick1 <- "((A:1,B:1):1,(C:1,D:1):1);"
  newick2 <- paste("((((XD:1,ZD:1):1,CE:2):1,(FE:2,EE:2):1):4,",
    "((AE:1,BE:1):1,(WD:1,YD:1):1):5);", sep = "")
  phylogeny1 <- ape::read.tree(text = newick1)
  phylogeny2 <- ape::read.tree(text = newick2)

  expect_silent(
    nLTT::get_average_nltt(c(phylogeny1, phylogeny2),
      dt = 0.20, plot_nltts = TRUE)
  )
})

test_that("get_average_nltt: check data types", {
  # Create a list or multiPhylo of phylogenies (of type phylo)
  # and run it through the get_average_nltt function

  n_trees <- 2
  n_tips <- 3
  set.seed(41)
  ape_phylogenies <- ape::rmtree(N = n_trees, n = n_tips)

  expect_silent(
    get_average_nltt(ape_phylogenies)
  )

  set.seed(41)
  treesim_phylogenies <- TreeSim::sim.bd.age(
    6, numbsim = n_trees, lambda = 0.4, mu = 0.0, complete = FALSE)

  expect_silent(
    get_average_nltt(treesim_phylogenies)
  )

  set.seed(41)
  combined_phylogenies <- c(ape::rcoal(10), ape::rcoal(20))

  expect_silent(
    get_average_nltt(combined_phylogenies)
  )

})



test_that("get_average_nltt: stop on incorrect input", {

  n_trees <- 2
  n_tips <- 3
  set.seed(41)
  ape_phylogenies <- ape::rmtree(N = n_trees, n = n_tips)

  single_phylogeny <- ape::rmtree(N = 1, n = n_tips)

  # must supply at least something
  expect_error(get_average_nltt(c()))

  #  dt must be from 0.0 to and including 1.0
  expect_error(get_average_nltt(ape_phylogenies, dt = -0.1))
  expect_error(get_average_nltt(ape_phylogenies, dt = 1.1))

  # must supply at least two trees
  expect_error(get_average_nltt(single_phylogeny))

  # must supply a phylogeny
  expect_error(get_average_nltt(c(1, 2, 3)))

  # must supply only phylogenies
  expect_error(get_average_nltt(list(c(1, 2), single_phylogeny)))
})
