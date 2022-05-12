context("mcmc_nltt")

test_that("mcmc_nltt use", {
  skip_on_cran() # These tests are very long
  testthat::skip_if_not_installed("TESS")

  set.seed(1)
  tree1 <- TESS::tess.sim.taxa(n = 1, nTaxa = 50,
                               max = 100, lambda = 1.0, mu = 0.0)[[1]]

  ll_bd <- function(params, phy) {
    lnl <- TESS::tess.likelihood(ape::branching.times(phy),
                                 lambda = params[1], mu = params[2],
                                 samplingProbability = 1, log = TRUE)
    prior1  <- dexp(params[1], rate = 10, log = TRUE)
    prior2  <- dexp(params[2], rate = 10, log = TRUE)
    return(lnl + prior1 + prior2)
  }

  tofit <- function(params) {
    if (params[1] <= 0) return(1e6)
    if (params[2] < 0) return(1e6)
    if (params[1] > 100) return(1e6)
    if (params[2] > 100) return(1e6)
    return(-1 * ll_bd(params, tree1))
  }

  max_lik <- optim(par = c(1, 0.001), fn = tofit)

  testthat::expect_output(
  mcmc_result <- mcmc_nltt(tree1, ll_bd, c(1, 0.001), c(TRUE, TRUE),
                           iterations = 10000, burnin = 1000,
                           thinning = 1, sigma = 1)
  )
  testthat::expect_equal(
    colMeans(mcmc_result)[[1]],
    max_lik$par[[1]],
    tolerance = 0.05
  )
testthat::expect_output(
  mcmc_result1 <- mcmc_nltt(tree1, ll_bd, c(1, 0.01), c(TRUE, TRUE),
                            iterations = 10000,
                            burnin = 1000, thinning = 1, sigma = 0.5)
)
testthat::expect_output(
  mcmc_result2 <- mcmc_nltt(tree1, ll_bd, c(1, 0.01), c(FALSE, FALSE),
                            iterations = 10000,
                            burnin = 1000, thinning = 1, sigma = 0.5)
)
  expect_equal(
    colMeans(mcmc_result1)[[1]],
    colMeans(mcmc_result2)[[1]],
    tolerance = 0.05
  )
})

test_that("mcmc_nltt abuse", {

  set.seed(1) #just to be safe
  tree1 <- ape::rphylo(n = 50, birth  = 1, death = 0)

  ll_bd <- function(params, phy) {
    lnl <- dexp(params[1]) * dexp(params[2]) # fake ll, not used here
    prior1  <- dexp(params[1], rate = 10, log = TRUE)
    prior2  <- dexp(params[2], rate = 10, log = TRUE)
    return(lnl + prior1 + prior2)
  }

  expect_output(
    expect_error(
      mcmc_nltt(tree1, ll_bd, c(1, 0.0), c(TRUE, TRUE),
                 iterations = 10000, burnin = 1000, thinning = 1, sigma = 0.5),
      "Cannot propose new value for a parameter with value 0.0."
    )
  )

  expect_error(
    mcmc_nltt(42, ll_bd, c(1, 0.01), c(TRUE, TRUE),
               iterations = 10000, burnin = 1000, thinning = 1, sigma = 0.5),
    "mcmc_nltt: phy must be of class 'phylo'"
  )

  expect_error(
    mcmc_nltt(tree1, ll_bd, c(1, -1), c(TRUE, TRUE),
               iterations = 10000, burnin = 1000, thinning = 1, sigma = 0.5),
    "mcmc_nltt: initial parameter values have to be above zero"
  )
})
