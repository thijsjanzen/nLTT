\name{nltt_plot}
\alias{nltt_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Normalized version of the ape function ltt.plot
}
\description{
    This function uses a modified version of the ltt.plot function from \code{"ape"} to plot the normalized number of lineages through normalized time, where the number of lineages is normalized by dividing by the number of tips of the tree, and the time is normalized by the total time between the most common recent ancestor and the present, such that t(MRCA) = 0 & t(present) = 1.
}
\usage{
nltt_plot(phy, xlab = "Normalized Time", ylab = "Normalized Lineages", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{
	an object of class \code{"phylo"};
}
  \item{xlab}{
a character string (or a variable of mode character)
    giving the label for the \eqn{x}-axis (default is "Normalized Time").
}
  \item{ylab}{
a character string (or a variable of mode character)
    giving the label for the \eqn{y}-axis (default is "Normalized Lineages").
}
  \item{\dots}{further graphical arguments that can be passed to \code{plot()}.
}
}
\author{
 Thijs Janzen
}
\examples{
	data(exampleTrees);
	nltt_plot(exampleTrees[[1]]);
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
