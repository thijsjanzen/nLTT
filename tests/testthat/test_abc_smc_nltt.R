context("abc_smc_nltt")

test_that("abc_smc_nltt use", {
  skip("Convergence seems to take extremely long,
       which makes it really hard to make good usage tests")
  treesim <- function(params) {
    t <- TESS::tess.sim.taxa(n = 1,
                             lambda = params[1], 
                             mu = params[2], nTaxa = 1000, max = 100000)[[1]]
    return (t)
  }

  prior_gen <- function() {
    return ( rexp(n = 2, rate = 1) )
  }

  prior_dens <- function(val) {
    if (val[2] > val[1]) {
      return(-1)
    }
    return( dexp( val[1], rate = 1) *
             dexp( val[2], rate = 1) )
  }

  set.seed(1)
  obs <- treesim(c(20,0))

  LL_BD <- function(params, phy) {
    lnl <- TESS::tess.likelihood(ape::branching.times(phy),
                                 lambda = params[1], mu = params[2],
                                 samplingProbability = 1, log = TRUE)
    prior1  <- dexp(params[1], rate = 1, log = TRUE)
    prior2  <- dexp(params[2], rate = 1, log = TRUE)
    return(lnl + prior1 + prior2)
  }

  tofit <- function(params) {
    if (params[1] <= 0) return(1e6)
    if (params[2] < 0) return(1e6)
    if (params[1] > 100) return(1e6)
    if (params[2] > 100) return(1e6)
    return(-1 * LL_BD(params, obs))
  }

  ML <- optim(par = c(1, 0.001), fn = tofit)

  statwrapper <- function(tree1) {
    return(nLTTstat(tree1, obs, "abs"))
  }

  v1 <-  abc_smc_nltt(
    obs, c(statwrapper), treesim, init_epsilon_values = 0.1,
    prior_generating_function = prior_gen,
    prior_density_function = prior_dens,
    number_of_particles = 100, sigma = 1, stop_rate = 0.1
  )

  expect_equal(
    colMeans( mcmc_nltt( tree1, LL_BD, c(1, 0.001), c(TRUE, TRUE),
                         iterations = 10000, burnin = 1000,
                         thinning = 1, sigma = 1))[[1]],
    ML$par[[1]],
    tolerance = 0.05
  )

})

test_that("abc_smc_nltt abuse", {

  treesim <- function(params) {
    t <- TESS::tess.sim.taxa(n = 1,
                             lambda = params[1],
                             mu = params[2], nTaxa = 1000, max = 100000)[[1]]
    return (t)
  }

  prior_gen <- function() {
    return( rexp(n = 2, rate = 1) )
  }

  prior_dens <- function(val) {
    if (val[2] > val[1]) {
      return(-1)
    }
    return ( dexp( val[1], rate = 1) *
               dexp( val[2], rate = 1) )
  }

  statwrapper <- function(tree1) {
    return( nLTTstat(tree1, obs, "abs"))
  }

  expect_error(
    abc_smc_nltt(
      42, c(statwrapper), treesim, init_epsilon_values = 0.1,
      prior_generating_function = prior_gen,
      prior_density_function = prior_dens,
      number_of_particles = 100, sigma = 1, stop_rate = 0.1),
    "abc_smc_nltt: tree must be of class 'phylo'"
  )

  obs <- treesim(c(1, 0))
  expect_error(
    abc_smc_nltt(obs, c(statwrapper), treeSim, init_epsilon_values = -0.5,
      prior_generating_function = prior_gen,
      prior_density_function = prior_dens,
      number_of_particles = 100, sigma = 1, stop_rate = 0.1),
    "abc_smc_nltt: epsilon values have to be positive"
  )

  expect_error(
    abc_smc_nltt(obs, statWrapper, treeSim, init_epsilon_values = -0.5,
                 prior_generating_function = prior_gen,
                 prior_density_function = prior_dens,
                 number_of_particles = 100, sigma = 1, stop_rate = 0.1),
    "abc_smc_nltt: the statistics function has to be given in vector style,"
  )
})

