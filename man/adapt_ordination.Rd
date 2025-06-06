% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internal_builders.R
\name{adapt_ordination}
\alias{adapt_ordination}
\title{Adapt format of foreign multivariate analyses}
\usage{
adapt_ordination(ordination)
}
\arguments{
\item{ordination}{An object of class \code{"gm.prcomp"}, \code{"bgPCA"},
\code{"PCA"}, \code{"pls2B"}, \code{"pls"}, \code{"phyl.pca"} or \code{"mvgls.pca"}.}
}
\value{
An object of equivalent class with scores, eigenvectors, eigenvalues/
standard deviations, and original centroid arranged so they can be used
in downstream \code{morphospace} operations.
}
\description{
Reorganize results from \code{\link[geomorph]{gm.prcomp}},
\code{\link[Morpho]{groupPCA}}, \code{\link[Morpho]{pls2B}},
\code{\link[phytools]{phyl.pca}} and \code{\link[mvMORPH]{mvgls.pca}} for
their use in the mspace workflow. Used internally.
}
\examples{
#load data
data("tails")

#perform PLS
pls <- Morpho::pls2B(x = tails$sizes, y = geomorph::two.d.array(tails$shapes))

#adapt format
pls_adapted <- adapt_ordination(pls)

#compare contents
names(pls)
names(pls_adapted)
}
\keyword{internal}
