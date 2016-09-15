################################################################################
#
# @brief calculate the weight of a parameter combination
#
# @date Last modified: 2016-03-07
# @author Thijs Janzen
# @since 2014-20-09, version 1.0
#
# @param    weights                 vector       Vector of weights
# @param    particles               list         List of parameter combinations
# @param    current                 vector       Current parameter combination for which we are determining the weight
# @param    sigma                   scalar       standard deviation of the perturbation
# @param    prior_density_function  function     Function to calculate the prior probability of a set of parameters
# @return                           scalar       Estimated weight
#
################################################################################

calculate_weight <- function(weights, particles,
                             current, sigma, prior_density_function) {
  vals <- c()
  for ( i in seq_along( particles)) {
    diff1 <- log(current[1]) - log(particles[[i]][1])
    diff2 <- log(current[2]) - log(particles[[i]][2])
    vals[i] <- weights[i] * stats::dnorm(diff1, mean = 0, sd = sigma) *
      stats::dnorm(diff2, mean = 0, sd = sigma)
  }

  numerator <- prior_density_function(current)

  return (numerator / sum( vals))
}

################################################################################
#
# @brief apply the ABC routine used in our Methods in Ecology and Evolution Paper
#
# @date Last modified: 2016-03-07
# @author Thijs Janzen
# @since 2014-20-09, version 1.0
#
# @param    tree                          phylo       Phylogenetic tree
# @param    statistics                    vector      A vector containing functions that take a tree as an argument and return a single scalar value (the statistic).
# @param    simulation_function           vector      A function that implements the diversification model and returns an object of class phylo
# @param    init_epsilon_values           vector      A vector containing the initial threshold values for the summary statistics from the vector statistics.
# @param    prior_generating_function     function    Function to generate parameters from the prior distribution of these parameters (e.g. a function returning lambda and mu in case of the birth-death model)
# @param    prior_density_function        function    Function to calculate the prior probability of a set of parameters (should match prior_generating_function in shape)
# @param    number_of_particles           scalar      Number of particles to be used per iteration of the ABC-SMC algorithm.
# @param    sigma                         scalar      Standard deviation of the perturbance distribution (perturbance distribution is a gaussian with mean 0).
# @param    stop_rate                     scalar      If the acceptance rate drops below \code{stopRate}, stop the ABC-SMC algorithm  and assume convergence.
# @return                                 matrix      A matrix with n columns, where n is the number of parameters you are trying to estimate.
#
################################################################################
#' A function to perform Approximate Bayesian Computation within an Sequential Markov Chain (ABC-SMC), for diversification analysis of phylogenetic trees.
#' @description This function performs ABC-SMC as described in Toni 2009 for given diversification model, provided a phylogenetic tree. ABC-SMC is not limited to only using the normalized LTT as statistic.
#' @usage
#'   abc_smc_nltt(
#'     tree, statistics, simulation_function, init_epsilon_values,
#'     prior_generating_function, prior_density_function,
#'     number_of_particles = 1000, sigma = 0.05, stop_rate = 1e-05
#'   )
#' @param tree an object of class \code{"phylo"}; the tree upon which we want to fit our diversification model
#' @param statistics A vector containing functions that take a tree as an argument and return a single scalar value (the statistic).
#' @param simulation_function A function that implements the diversification model and returns an object of class \code{"phylo"}.
#' @param init_epsilon_values A vector containing the initial threshold values for the summary statistics from the vector \code{statistics}.
#' @param prior_generating_function Function to generate parameters from the prior distribution of these parameters (e.g. a function returning lambda and mu in case of the birth-death model)
#' @param prior_density_function Function to calculate the prior probability of a set of parameters.
#' @param number_of_particles Number of particles to be used per iteration of the ABC-SMC algorithm.
#' @param sigma Standard deviation of the perturbance distribution (perturbance distribution is a gaussian with mean 0).
#' @param stop_rate   If the acceptance rate drops below \code{stopRate}, stop the ABC-SMC algorithm  and assume convergence.
#' @return A matrix with \code{n} columns, where \code{n} is the number of parameters you are trying to estimate.
#' @references  Toni, T., Welch, D., Strelkowa, N., Ipsen, A., & Stumpf, M.P.H. (2009). Approximate Bayesian computation scheme for parameter inference and model selection in dynamical systems. Journal of the Royal Society Interface, 6(31), 187-202.
#' @export
#' @author Thijs Janzen
#' @examples
#'   \dontrun{
#'
#'   prior_gen <- function() {
#'     return ( rexp(n=2, rate=0.1) )
#'   }
#'
#'   prior_dens <- function(val) {
#'     return ( dexp( val[1], rate = 0.1) * dexp( val[2], rate = 0.1) )
#'   }
#'
#'   require(TESS)
#'
#'   treeSim <- function(params) {
#'     t <- TESS.sim.age(n=1, lambda = params[1], mu = params[2], age = 10)[[1]]
#'     return (t)
#'   }
#'
#'   obs <- treeSim(c(0.5,0.1))
#'
#'   statWrapper <- function(tree1) {
#'     return ( nLTTstat_exact(tree1, obs, "abs"))
#'   }
#'
#'   stats <- c(statWrapper)
#'
#'   results <-  abc.smc.nltt(
#'     obs, stats, treeSim, init_epsilon_values = 0.2,
#'     prior_generating_function = prior_gen,
#'     prior_density_function = prior_dens,
#'     number_of_particles = 1000, sigma = 0.05, stop_rate = 1e-5
#'   )
#'
#'   } # end of dontrun
abc_smc_nltt <- function(tree,
                         statistics,
                         simulation_function,
                         init_epsilon_values,
                         prior_generating_function,
                         prior_density_function,
                         number_of_particles = 1000,
                         sigma = 0.05,
                         stop_rate = 1e-5) {

  if (!inherits(tree, "phylo")) {
    # Just checking
    stop("abc_smc_nltt: ",
         "tree must be of class 'phylo', ",
         "but was of type '", class(tree), "' instead")
  }

  #statistics has to be a vector of functions
  if (!inherits(statistics, "list")) {
    stop("abc_smc_nltt: ",
         "the statistics function has to be given in vector style, ",
         "e.g.: c(statisticsfunction), instead of statisticsfunction")
  }


  #just to get the number of parameters to be estimated.
  parameters <- prior_generating_function()

  # compute the observed statistics
  obs_statistics <- c()
  for (i in seq_along(statistics)) {
    obs_statistics[i] <- statistics[[i]](tree)
  }

  stats <- c()

  #generate a matrix with epsilon values
  #we assume that the SMC algorithm converges within 50 iterations
  epsilon <- matrix(nrow = 50, ncol = length(init_epsilon_values))
  for (j in seq_along(init_epsilon_values)) {
    if (init_epsilon_values[j] < 0) {
      stop("abc_smc_nltt: ",
           "epsilon values have to be positive,",
           "but were instead: ", init_epsilon_values[j])
    }

    for (i in 1:50) {
      epsilon[i, j] <- init_epsilon_values[j] * exp(-0.5 * (i - 1))
    }
  }

  #store weights
  new_weights <- c()
  new_params <- list( c( seq_along(parameters)))
  previous_weights <- c()
  previous_params  <- list( c( seq_along(parameters)))
  indices <- 1:number_of_particles

  #convergence is expected within 50 iterations
  #usually convergence occurs within 20 iterations
  for (i in 1:50 ) {
    cat("\nGenerating Particles for iteration\t", i, "\n")
    cat("0--------25--------50--------75--------100\n")
    cat("*")
    utils::flush.console()

    print_frequency <- 20
    tried <- 0
    number_accepted <- 0

    #replace all vectors
    if (i > 1) {
      #normalize the weights and store them as previous weights.
      previous_weights <- new_weights / sum(new_weights)
      new_weights <- c() #remove all currently stored weights
      previous_params <- new_params #store found params
      new_params <- list( c( seq_along(parameters))) #clear new params
    }

    while (number_accepted < number_of_particles) {
      #in this initial step, generate parameters from the prior
      if (i == 1) {
        parameters <- prior_generating_function()
      } else {
        #if not in the initial step, generate parameters
        #from the weighted previous distribution:
        index <- sample(x = indices, size = 1,
                        replace = TRUE, prob = previous_weights)

        for (p_index in seq_along(parameters)) {
          parameters[p_index] <- previous_params[[index]][p_index]
        }

        #only perturb one parameter, to avoid extremely
        #low acceptance rates due to simultaneous perturbation
        to_change <- sample(seq_along(parameters), 1)

        # perturb the parameter a little bit,
        #on log scale, so parameter doesn't go < 0
        eta <- log(parameters[to_change]) + stats::rnorm(1, 0, sigma)
        parameters[to_change] <- exp(eta)
      }

      #reject if outside the prior
      if (prior_density_function(parameters) > 0) {
        #simulate a new tree, given the proposed parameters
        new_tree <- simulation_function(parameters)
        accept <- TRUE

        #calculate the summary statistics for the simulated tree
        for (k in seq_along(statistics)) {
          stats[k] <- statistics[[k]](new_tree)
          if (is.na(stats[k])) stats[k] <- Inf
        }

        #check if the summary statistics are sufficiently
        #close to the observed summary statistics
        for (k in seq_along(statistics)) {
          if ( abs(stats[k] - obs_statistics[k]) > epsilon[i, k] ) {
            accept <- FALSE
            #the first step always accepts
            if (i == 1) accept <- TRUE
            break
          }
        }

        if ( accept ) {
          number_accepted <- number_accepted + 1
          new_params[[number_accepted]] <- parameters
          accepted_weight <- 1
          #calculate the weight
          if (i > 1) {
            accepted_weight <- calculate_weight(previous_weights,
                                                previous_params, parameters,
                                                sigma, prior_density_function)
          }
          new_weights[number_accepted] <- accepted_weight

          if ( (number_accepted) %%
               (number_of_particles / print_frequency) == 0) {
            cat("**")
            utils::flush.console()
          }
        }
      }

      #convergence if the acceptance rate gets too low
      tried <- tried + 1
      if (tried > (1 / stop_rate)) {
        if ( (number_accepted / tried) < stop_rate) {
          output <- c()
          for (k in seq_along(previous_params)) {
            add <- c()
            for (m in seq_along( parameters)) {
              add <- c( add, previous_params[[k]][m])
            }
            output <- rbind(output, add)
          }
          return (output)
        }
      }
    }
  }

  output <- c()
  for (k in seq_along(previous_params)) {
    add <- c()
    for (m in seq_along( parameters)) {
      add <- c( add, previous_params[[k]][m])
    }
    output <- rbind(output, add)
  }
  return (output)
}

################################################################################
#
# @brief Estimate the likelihood of a given tree, provided a likelihood function, using a Monte Carlo Markov Chain
#
# @date Last modified: 2014-20-09
# @author Thijs Janzen
# @since 2014-20-09, version 1.0
#
# @param    phy                   phylo       Vector of weights
# @param    likelihood_function   function    Function that calculates the likelihood of our diversification model, given the tree.
#                                             function should be of the format function(parameters, phy).
# @param    parameters            vector      Initial parameters to start the chain.
# @param    logtransforms         scalar      Whether to perform jumps on logtransformed parameters (TRUE) or not (FALSE)
# @param    iterations            scalar      Length of the chain
# @param    burnin                scalar      Length of the burnin, default is 30% of iterations
# @param    thinning              scalar      Size of thinning, default = 1
# @param    sigma                 scalar      Standard deviation of the jumping distribution, which is N(0, sigma).
# @return                         mcmc        An MCMC object, as used by the package "coda".
#
################################################################################

mcmc_nltt <- function(phy, likelihood_function,
                      parameters, logtransforms, iterations,
                      burnin = round(iterations / 3), thinning = 1, sigma = 1) {

  #check data type of phy
  if (!inherits(phy, "phylo")) {
    # Just checking
    stop("mcmc_nltt: ",
         "phy must be of class 'phylo', ",
         "but was of type '", class(phy), "' instead")
  }

  # create a list for the samples & reserve memory for the chain
  chain <- array(dim = c( floor( iterations / thinning) + 1,
                          length(parameters)))

  if (parameters[2] < 0) {
    #Just checking
    stop("mcmc_nltt: ",
         "initial parameter values have to be above zero\n",
         "but mu was ", parameters[2], " instead")
  }

  # pre-compute current posterior probability
  pp <- likelihood_function(parameters, phy)

  cat("\nGenerating Chain\n")
  cat("0--------25--------50--------75--------100\n")
  cat("*")
  utils::flush.console()
  print_frequency <- 20

  for (i in seq_len(burnin + iterations)) {
    #propose new values
    for ( j in seq_along(parameters) ) {
      if ( logtransforms[j] == TRUE ) {
        if ( parameters[j] == 0) {
          stop("Cannot propose new value for a parameter with value 0.0.")
        }

        eta           <- log(parameters[j])
        new_eta       <- eta + stats::rnorm(1, 0, sigma)
        new_val       <- exp(new_eta)
        # calculate the Hastings ratio
        hr            <- log(new_val / parameters[j])
        parameters[j] <- new_val
        new_pp        <- likelihood_function(parameters, phy)

        #accept or reject
        if ( is.finite(new_pp) &&
             is.finite(hr) &&
             new_pp - pp + hr > log(stats::runif(1, 0, 1)) ) {
          pp <- new_pp
        } else {
          parameters[j] <- exp(eta)
        }
      } else {

        eta           <- parameters[j]
        new_val       <- eta + stats::rnorm(1, 0, sigma)
        #calculate the Hastings ratio
        hr            <- 0.0
        parameters[j] <- new_val

        if (parameters[j] >= 0 & parameters[1] > 0) {
          new_pp        <- likelihood_function(parameters, phy)

          #accept or reject
          if ( is.finite(new_pp) &&
               is.finite(hr) &&
               new_pp - pp + hr > log(stats::runif(1, 0, 1)) ) {
            pp <- new_pp
          } else {
            parameters[j] <- eta
          }
        } else {
          parameters[j] <- eta
        }
      }
    }

    # sample the parameter
    if (i >= burnin) {
      if ( (i) %% ( (iterations - burnin) / print_frequency) == 0) {
        cat("**")
        utils::flush.console()
      }
      if ( (i - burnin) %% thinning == 0 ) {
        chain[ (i - burnin) / thinning + 1, ] <- parameters
      }
    }
  }
  cat("\nFinished MCMC.\n")
  #return a mcmc object, used by coda to plot
  return( coda::as.mcmc(chain) )
}
