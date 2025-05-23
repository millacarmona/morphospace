% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mspace_workflow.R
\name{proj_groups}
\alias{proj_groups}
\title{Delimit groups in morphospace}
\usage{
proj_groups(
  mspace,
  shapes = NULL,
  groups = NULL,
  ellipse = FALSE,
  conflev = 0.95,
  density = TRUE,
  labels = NULL,
  pipe = TRUE,
  ...
)
}
\arguments{
\item{mspace}{An \code{"mspace"} object.}

\item{shapes}{Optional shape data (if \code{NULL}, the function will look for
them first in \code{mspace$projected$scores} and then in
\code{mspace$ordination$x}).}

\item{groups}{Factor; classification of observations into groups. Its length
must be the same as the number of shapes provided in \code{shapes}.}

\item{ellipse}{Logical; whether to plot confidence ellipses (if \code{FALSE},
convex hulls will be used instead).}

\item{conflev}{Numeric; confidence level used for confidence ellipse(s).}

\item{density}{Logical; whether to add density distribution for groups
(univariate ordinations only).}

\item{labels}{Either logical, indicating whether to include all group labels,
or a character string with the exact names of the groups whose labels
should be included.}

\item{pipe}{Logical; is the function being included in a pipe?}

\item{...}{Further arguments passed to \code{\link[=hulls_by_group_2D]{hulls_by_group_2D()}} or
\code{\link[=ellipses_by_group_2D]{ellipses_by_group_2D()}}.}
}
\value{
If a plot device with a morphospace is open, convex hulls or
confidence ellipses enclosing the scores corresponding to \code{groups}
are projected into morphospace. If \code{pipe = TRUE} the supplied
\code{"mspace"} object will be modified by appending a \code{$gr_class}
slot to \code{$projected}, as well as by adding some graphical parameters
(stored into the \code{$plotinfo} slot), and returned invisibly.
}
\description{
Project convex hulls or confidence ellipses enclosing
\emph{a priori} groups into an existing morphospace.
}
\details{
The goal of this function is to add a classification for shapes
populating the morphospace to \code{"mspace"} objects, as well as to
facilitate group visualization. Other than that, it is just a wrapper for
\code{hulls_by_group_2D} and \code{ellipses_by_group_2D}.
}
\examples{
#load and extract relevant data, packages and information
library(magrittr)
data("shells")
shapes <- shells$shapes
species <- shells$data$species
sp_shapes <- expected_shapes(shapes, species)

#generate basic morphospace, add sampled shapes and convex hulls for species
mspace(shapes, links = links, mag = 0.7, axes = c(1,2),
       bg.model = "light gray") \%>\%
  proj_shapes(shapes = shapes, col = c(1:4)[species], pch = 1) \%>\%
  proj_groups(groups = species, col = 1:4, lwd = 2, alpha = 0.5)

#generate basic morphospace, add sampled shapes and 95\% confidence for
#species
mspace(shapes, links = links, mag = 0.7, axes = c(1,2),
       bg.model = "light gray") \%>\%
  proj_shapes(shapes = shapes, col = c(1:4)[species], pch = 1) \%>\%
  proj_groups(groups = species, col = 1:4, lwd = 2, ellipse = TRUE,
              conflev = 0.95, alpha = 0.1)

#add labels for all the groups
mspace(shapes, links = links, mag = 0.7, axes = c(1,2),
       bg.model = "light gray") \%>\%
  proj_shapes(shapes = shapes, col = c(1:4)[species], pch = 1) \%>\%
  proj_groups(groups = species, col = 1:4, lwd = 2, ellipse = TRUE,
              conflev = 0.95, alpha = 0.1, labels = TRUE)

#add labels for selected the groups
msp <- mspace(shapes, links = links, mag = 0.7, axes = c(1,2),
              bg.model = "light gray") \%>\%
  proj_shapes(shapes = shapes, col = c(1:4)[species], pch = 1) \%>\%
  proj_groups(groups = species, col = 1:4, lwd = 2, ellipse = TRUE,
              conflev = 0.95, alpha = 0.1, labels = c("coihuicoensis", "esbelta"))

#add legend using plot_mspace()
plot_mspace(msp, legend = TRUE, cex.legend = 1, pch.groups = 16)
}
\seealso{
\code{\link{hulls_by_group_2D}}, \code{\link{ellipses_by_group_2D}}
}
