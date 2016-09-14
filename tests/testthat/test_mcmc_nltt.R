context("mcmc_nltt")

test_that("mcmc_nltt use", {
  skip("@ThijsJanzen will fix this :-)")

  set.seed(1) #just to be safe
  tree1 <- TESS::tess.sim.taxa(n = 1, nTaxa = 100, max = 100, lambda = 1.0, mu = 0.0)[[1]]
  tree2 <- TESS::tess.sim.taxa(n = 1, nTaxa = 100, max = 100, lambda = 1.0, mu = 0.0)[[1]]

  LL_BD <- function(params, phy) {
    lnl <- TESS::tess.likelihood(ape::branching.times(phy), lambda = params[1], mu = params[2],
								  samplingProbability = 1, log = TRUE)
    prior1 <- dunif( params[1], 0, 100, log = TRUE)
    prior2 <- dunif( params[2], 0, 100, log = TRUE)
	return(lnl + prior1 + prior2)
  }

  #calculate Maximum Likelihood estimates
  #these should match our MCMC results
  ML1 <- DDD::bd_ML(ape::branching.times(tree1))
  ML2 <- DDD::bd_ML(ape::branching.times(tree2))

  expect_equal(
    colMeans( mcmc_nltt( tree1, LL_BD, c(1,0.01), c(TRUE,TRUE), 
      iterations = 10000, burnin=1000, thinning=1, sigma=0.05))[[1]],
    ML1[[1]],
    tolerance = 0.01
  )

  expect_equal(
    colMeans( mcmc_nltt( q, LL_BD, c(1,0.01), c(TRUE,TRUE),
	iterations = 10000, burnin = 1000, thinning = 1, sigma = 0.5))[[1]],
	ML2[[1]],
	tolerance = 0.001
  )

  expect_equal( #compare jumps in log space with jumps in normal space, should yield similar results
    colMeans( mcmc_nltt( q, LL_BD, c(1, 0.01), c(TRUE,TRUE), iterations = 10000, burnin = 1000, thinning = 1, sigma = 0.5))[[1]],
    colMeans( mcmc_nltt( q, LL_BD, c(1, 0.01), c(FALSE,FALSE), iterations = 10000, burnin = 1000, thinning = 1, sigma = 0.01))[[1]],
    tolerance = 0.001
  )
})

test_that("mcmc_nltt abuse", {
  skip("@ThijsJanzen will fix this :-)")
	set.seed(1) #just to be safe
	p <- TESS::tess.sim.taxa(n = 1, nTaxa = 50, max = 100, lambda = 1.0, mu = 0.0)[[1]]
	q <- TESS::tess.sim.taxa(n = 1, nTaxa = 50, max = 100, lambda = 1.0, mu = 0.0)[[1]]

	LL_BD <- function(params, phy) {
		lnl <- TESS::tess.likelihood(phy[[1]], lambda = params[1], mu = params[2],
	 							  samplingProbability = 1, log = TRUE)
	 	prior1 <- dunif( params[1], 0, 100, log = TRUE)
	 	prior2 <- dunif( params[2], 0, 100, log = TRUE)
	 	return(lnl + prior1 + prior2)
	}

	expect_error(
		mcmc_nltt( q, LL_BD, c(1, 0.0), c(TRUE,TRUE),
		  iterations = 10000, burnin = 1000, thinning = 1, sigma = 0.5),
		"Cannot propose new value for a parameter with value 0.0."
	)

	expect_error(
		mcmc_nltt( 42, LL_BD, c(1, 0.01), c(TRUE, TRUE), 
		  iterations = 10000, burnin = 1000, thinning = 1, sigma = 0.5),
		"mcmc_nltt: phy must be of class 'phylo'"
	)

	expect_error(
		mcmc_nltt( p, LL_BD,c(1, -1), c(TRUE, TRUE), 
		  iterations = 10000, burnin = 1000, thinning = 1, sigma = 0.5),
		"mcmc_nltt: initial parameter values have to be above zero"
	)
})
