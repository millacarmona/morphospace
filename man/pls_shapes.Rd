% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multivariate_analyses.R
\name{pls_shapes}
\alias{pls_shapes}
\title{2B Partial Least Squares for shape data}
\usage{
pls_shapes(
  X,
  shapes,
  tree = NULL,
  evmodel = "BM",
  LOOCV = FALSE,
  recompute = FALSE
)
}
\arguments{
\item{X}{A matrix with variables as columns and observations as rows,
representing the external variables that supervise ordination
(corresponding to the first block).}

\item{shapes}{Shape data (corresponding to the second block of
\code{pls2b}.}

\item{tree}{A \code{"phylo"} object containing a phylogenetic tree. Tip
labels should match the row number and names from \code{x} and \code{y}.}

\item{evmodel}{Character, specifying an evolutionary model to perform
ancestral character reconstruction; options are "BM" (Brownian motion),
"EB" (Early burst) and "lambda" (Pagel's lambda transformation) (see
\code{\link[mvMORPH]{mvgls}} for more details).}

\item{LOOCV}{Logical; whether to apply leave-one-out cross-validation.}

\item{recompute}{Logical; whether to re-compute rotation matrix using the
scores resulting from LOOCV.}
}
\value{
A \code{"pls_shapes"} or \code{"phy_pls_shapes"} object formatted
following the \code{"prcomp"} class:
\itemize{
\item \code{$x:} the scores from the shapes block.
\item \code{$x2:} the scores from the supervising block (i.e., the
\code{x} scores).
\item \code{$sdev:} the standard deviations of the PLS axes (i.e., the
standard deviations calculated from the scores of each axis in the
shape block).
\item \code{$rotation:} a matrix of vector coefficients for the shape
block.
\item \code{$center:} the mean values of the original variables from
the shape block.
\item \code{$totvar:} the sum of the variances of the original variables
in the shape block.
}
}
\description{
A wrapper for \code{\link[=pls2b]{pls2b()}} specifically aimed at summarizing
covariation between shape data and other external, non-shape variable(s).
}
\details{
This function finds the linear combination maximizing covariation
between a block of variables assumed to be shape variables formatted as
a 2-margins matrix and another block which can be either another set of
shape variables or one or more non-shape variables for supervising the
analysis. If a phylogenetic tree is supplied, observations from \code{x}
and \code{shapes} blocks are treated as data measured for its tips. In this
case, the resulting axes maximize the residual covariation among blocks
left after removing covariation among blocks accounted by phylogenetic
history (assuming a Brownian model of evolution and 100\% phylogenetic
signal).

It has been reported that PLS scores from the different blocks can show
spurious covariation when the number of variables exceeds the number of
observations (a common situation in geometric morphometrics analyses). This
problem can be alleviated by carrying out a leave-one-out cross-validation
(LOOCV; i.e., each observation is excluded from the calculation of PLS axes
before its projection in the resulting ordination as a way to calculate its
score).

The number of axes in both blocks will be equal to the number of variables
in the smallest block (i.e., the one with less variables), unless the
number of observations is lower (in which case this will be equal to the
number PLS axes).
}
\examples{
#load data
data("tails")

#extract shapes and sizes, compute mean shapes and sizes for all species,
#extract tree
shapes <- tails$shapes
sizes <- log(tails$sizes)
sp_shapes <- expected_shapes(shapes, tails$data$species)
sp_sizes <- tapply(X = log(sizes), INDEX = tails$data$species, FUN = mean)
tree <- tails$tree

#perform PLS between shape and size
pls <- pls_shapes(shapes = shapes, X = sizes)

#inspect results
names(pls) #the contents of the resulting object
exp_var(pls) #variance explained by each axis
plot(pls$x2, pls$x) #ordination

#perform PLS between shape and size with leave-one-out CV
pls_cv <- pls_shapes(shapes = shapes, X = sizes, LOOCV = TRUE)

#inspect results
names(pls_cv) #the contents of the resulting object
exp_var(pls_cv) #variance explained by each axis
plot(pls_cv$x2, pls_cv$x) #ordination

#perform phylogenetic PLS between shape and size (for species means)
ppls <- pls_shapes(shapes = sp_shapes, X = sp_sizes, tree = tree)

#inspect results
names(ppls) #the contents of the resulting object
exp_var(ppls) #variance explained by each axis
plot(ppls$x2, ppls$x) #ordination

#perform phylogenetic PLS between shape and size with leave-one-out CV
ppls <- pls_shapes(shapes = sp_shapes, X = sp_sizes, tree = tree,
                   LOOCV = TRUE)

#inspect results
names(ppls) #the contents of the resulting object
exp_var(ppls) #variance explained by each axis
plot(ppls$x2, ppls$x) #ordination
}
\references{
Rohlf, F. J., & Corti, M. (2000). \emph{Use of two-block partial
least-squares to study covariation in shape}. Systematic Biology, 49(4),
740-753.

Zelditch, M. L., Swiderski, D. L., & Sheets, H. D. (2012). \emph{Geometric
morphometrics for biologists: A primer}. Second ed. Academic Press.

Bookstein, F. L. (2019). \emph{Pathologies of between-groups principal
components analysis in geometric morphometrics}. Evolutionary Biology,
46(4), 271-302.
}
\seealso{
\code{\link{pls2b}}
}
