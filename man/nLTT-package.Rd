\name{nLTT-package}
\alias{nLTT-package}
\alias{nLTT}
\docType{package}
\title{
	Package providing functions to visualize the normalized Lineage-Through-Time statistic, and calculate the difference between two nLTT curves
}
\description{
	This package provides a function to visualize the normalized Lineage-Through-Time (nLTT) statistic, where the number of lineages relative to the maximum number of lineages in a phylogenetic tree is plotted against the relative time between the most common recent ancestor and the present. Furthermore the package provides a function to calculate the difference between two nLTT curves, including two different distance measurements. 
\cr
\cr
Updates:\cr
Version 1.2: added an "exact" nLTT function. This function is faster for small trees, and provides an exact measurement of the nLTT function. Comparison between "old" and "exact" estimates show that these are highly correlated, although the "exact" values are slightly higher than the "old" values. The "exact" function should generally be preferred, unless dealing with extremely large trees (500+ tips) in which case the old function is much faster. \cr\cr 
Version 1.2: updated the example for the ABC_SMC_nLTT function, prior generating and prior density functions are now more realistic
\cr\cr
Version 1.1.1: fixed a minor bug in the ABC_SMC_nLTT function \cr
Version 1.1.1: removed some intermediate output in ABC_SMC_nLTT function \cr
Version 1.1: Made a universal nLTT function called "nLTTstat", with argument "distanceMethod", this serves as a more elegant wrapper for the functions "normLTTdiffABS" and "normLTTdiffSQ" \cr
Version 1.1: Updated references in the manual \cr
	
}
\details{
\tabular{ll}{
Package: \tab nLTT\cr
Type: \tab Package\cr
Version: \tab 1.2\cr
Date: \tab 2016-04-25\cr
License: GPL 2.0
}

}

\references{
 Janzen,T. Hoehna,S., Etienne,R.S. (2015) Approximate Bayesian Computation of diversification rates from molecular phylogenies: introducing a new efficient summary statistic, the nLTT. Methods in Ecology and Evolution. doi: 10.1111/2041-210X.12350
}

\author{
Thijs Janzen

Maintainer: Thijs Janzen <thijsjanzen@gmail.com>
}
