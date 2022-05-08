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
Updates:
Version 1.4.6: Fixed testing, to comply _R_CHECK_LENGTH_1_CONDITION_. \cr
Version 1.4.4: Added support for phylogenies with extinct lineages. \cr
Version 1.4.3: Added support for log transformation before normalization. \cr
Version 1.4: Added the following four functions: get_branching_times, get_n_lineages, get_norm_brts and get_norm_n. Furthermore, vignette building has improved, and the underlying code base has been polished up as well. \cr
Version 1.3.1: Added walkthrough vignette, and updated several typos in the manual \cr
Version 1.3: Version 1.3 adds a lot of extended functionality: firstly, we have added functions to calculate, and plot, the average nLTT across a number of phylogenies. Furthermore, we have added vignettes, and we have added a GitHub repository. On the GitHub repository the vignettes are separately accessible through the wiki. Lastly we have added an extra option to the nLTT functions, where the user can specify if the used trees are rooted, or not. Under the hood, some changes have been made as well, the majority of the code is now conforming to the lintR code conventions, and we have written formalized tests that check correctness of all code (code coverage 100%).
\cr

Version 1.2.1: updated comments and coding style to adhere to the general coding rules. Backwards compatibility has been favoured for the nLTT stat functions. ABC related functions are no longer backwards compatible (variable names have been changed to adhere to coding style). \cr
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
Version: \tab 1.4.4 \cr
Date: \tab 2020-02-13\cr
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
