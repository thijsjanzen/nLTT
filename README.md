# nLTT

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/nLTT)](https://cran.r-project.org/package=nLTT)
[![Build Status](https://travis-ci.org/richelbilderbeek/nLTT.svg?branch=master)](https://travis-ci.org/richelbilderbeek/nLTT)
[![codecov](https://codecov.io/gh/richelbilderbeek/nLTT/branch/master/graph/badge.svg)](https://codecov.io/gh/richelbilderbeek/nLTT)
[![](http://cranlogs.r-pkg.org/badges/grand-total/nLTT)](http://cran.rstudio.com/web/packages/nLTT/index.html)
[![](http://cranlogs.r-pkg.org/badges/nLTT)](http://cran.rstudio.com/web/packages/nLTT/index.html)


Repository for the R nLTT package

## What is the nLTT statistic?
The nLTT statistic is a likelihood free summary statistic to compare the similarity between two phylogenetic trees.  It calculates the distance between the lineage through time curves of the two trees, after normalizing the lineage through time curves with respect to the maximum number of lineages obtained in each tree, and with respect to the total time between the root and the tips of the tree (see also the wiki).

A more detailed description, and a detailed analysis of the performance of the nLTT statistic can be found in the following paper:

Janzen, Thijs, Sebastian Höhna, and Rampal S. Etienne. "Approximate Bayesian computation of diversification rates from molecular phylogenies: introducing a new efficient summary statistic, the nLTT." Methods in Ecology and Evolution 6.5 (2015): 566-575. [link](http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12350/full)

## What else does the package do?
Apart from providing functions that calculate the nLTT statistic, the nLTT package for R also provides functions to:
- plot the normalized Lineages-Through-Time plot for a single, or for multiple, trees
- calculate and plot the average Lineages-Through-Time plot for multiple trees
- estimate the parameters of a model for which the likelihood is known (for comparison), using MCMC (as used in the paper)
- estimate the parameters of a model for which no likelihood is available, using ABC-SMC (as used in the paper)
