% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mspace_workflow.R
\name{plot_mspace}
\alias{plot_mspace}
\title{Plot morphospaces and combine them with other variables}
\usage{
plot_mspace(
  mspace,
  axes = NULL,
  links = NULL,
  template = NULL,
  x = NULL,
  y = NULL,
  nh,
  nv,
  mag,
  rescale,
  invax = NULL,
  adj_frame = c(1, 1),
  points = TRUE,
  models = TRUE,
  groups = TRUE,
  phylo = TRUE,
  shapeax = TRUE,
  landsc = TRUE,
  legend = FALSE,
  scalebar = FALSE,
  cex.legend = 1,
  asp = 1,
  xlab,
  ylab,
  xlim = NULL,
  ylim = NULL,
  size.models,
  asp.models,
  rot.models,
  col.models,
  bg.models,
  lwd.models,
  alpha.models,
  cex.ldm,
  col.ldm,
  pch.points = 1,
  col.points = 1,
  bg.points = 1,
  cex.points = 1,
  density.points = TRUE,
  labels.points,
  col.groups = 1,
  bg.groups = 1,
  pch.groups = 16,
  cex.groups = 1,
  ellipse.groups = mspace$plotinfo$ellipse.groups,
  conflev.groups = 0.95,
  lwd.groups = 1,
  lty.groups = 1,
  alpha.groups = 0,
  boxplot.groups = FALSE,
  density.groups = TRUE,
  labels.groups,
  col.phylo = 1,
  lwd.phylo = 1,
  lty.phylo = 1,
  labels.nodes,
  col.nodes,
  pch.nodes,
  bg.nodes,
  cex.nodes,
  col.tips,
  pch.tips,
  bg.tips,
  labels.tips,
  cex.tips,
  type.axis,
  lwd.axis = 1,
  lty.axis = 1,
  col.axis = 1,
  palette.landsc = grDevices::heat.colors,
  display.landsc,
  alpha.landsc,
  nlevels.landsc = 50,
  drawlabels.landsc = FALSE,
  lty.landsc = 1,
  lwd.landsc = 1
)
}
\arguments{
\item{mspace}{An \code{"mspace"} object created using the
\code{\link{mspace}} \%>\% \code{proj_*} pipeline.}

\item{axes}{Numeric of length 1 or 2, indicating the axes to be plotted.}

\item{links}{A list with the indices of the coordinates defining a wireframe
connecting landmarks (following the format used in \code{Morpho}), or a
2-columns matrix indicating the pairs of landmarks that should be linked
(following the format used in \code{geomorph}).}

\item{template}{Either a 2-column matrix with landmarks/semilandmarks and
template curves coordinates (for 2D shape data) or a \code{"mesh3d"}
object representing the mean shape of the sample (for 3D shape data).}

\item{x, y}{Optional vector with a non-morphometric variable (numerical or
categorical) to be plotted in the x or y axis against an ordination axis.
Alternatively, a \code{"phylo"} object can be provided.}

\item{nh, nv}{Positive integers; the number of shape models along the x
(\code{nh}) and the y (\code{nv}) axes.}

\item{mag}{Numeric; magnifying factor for shape models.}

\item{rescale}{Logical; whether to re-scale background shape models so shape
variation is shown more clearly.}

\item{invax}{Optional numeric indicating which of the axes provided in
\code{axes} needs to be inverted (options are \code{1}, \code{2} or
\code{c(1,2)}).}

\item{adj_frame}{Numeric of length 2, providing \emph{a posteriori} scaling
factors for the width and height of the frame, respectively.}

\item{points}{Logical; whether to plot the scatter points corresponding
to the sampled shapes stored in \code{mspace$projected$scores}.}

\item{models}{Logical; whether to plot background shape models (stored
in \code{mspace$projected$shapemodels}).}

\item{groups}{Logical; whether to plot the convex hulls/confidence
ellipses enclosing the groups stored in
\code{mspace$projected$gr_class}.}

\item{phylo}{Logical; whether to plot phylogenetic relationships stored
in \code{mspace$projected$phylo}.}

\item{shapeax}{Logical; whether to plot morphometric axes stored in
\code{mspace$projected$shape_axis}.}

\item{landsc}{Logical; whether to plot landscape surface stored in
\code{mspace$projected$landsc}.}

\item{legend}{Logical; whether to show legend for groups.}

\item{scalebar}{Logical; whether to show scalebar for landscapes.}

\item{cex.legend}{Numeric; size of legend labels/symbols.}

\item{xlim, ylim, xlab, ylab, asp}{Standard arguments passed to the generic
\code{\link[graphics:plot.default]{graphics::plot()}} function.}

\item{size.models}{Numeric; size factor for shape models.}

\item{asp.models}{Numeric; y/x aspect ratio of shape models.}

\item{rot.models}{Numeric; angle (in degrees) to rotate shape models.}

\item{col.models}{Color for wireframes/outlines of shape models.}

\item{bg.models}{Background color for outlines/meshes of shape models.}

\item{lwd.models}{Integer; width of border lines in wireframes/outlines
of shape models.}

\item{alpha.models}{Numeric; transparency factor for background models
(3D only).}

\item{cex.ldm}{Numeric; size of landmarks/semilandmarks in the background
models.}

\item{col.ldm}{Color of landmarks/semilandmarks in the background
models.}

\item{pch.points}{Symbol of the scatterpoints.}

\item{col.points}{Color of the scatterpoints.}

\item{bg.points}{Background color of the scatterpoints.}

\item{cex.points}{Numeric; size of the scatterpoints}

\item{density.points}{Logical; whether to add density distribution for
points (univariate ordinations only). Overriden by
\code{density.groups = TRUE}}

\item{labels.points}{Either logical, indicating whether to include all point
labels, or a character string with the exact names of the points whose
labels should be included.}

\item{col.groups}{Color of the hulls/ellipses and/or scatterpoints
corresponding to groups' mean shapes.}

\item{bg.groups}{Background color of the scatterpoints corresponding to
groups' mean shapes.}

\item{pch.groups}{Symbol of the scatterpoints corresponding to
groups' mean shapes.}

\item{cex.groups}{Numeric; size of the scatterpoints corresponding to
groups' mean shapes.}

\item{ellipse.groups}{Logical; whether to use confidence ellipses to
delimit groups (if \code{FALSE} convex hulls are used instead).}

\item{conflev.groups}{Numeric; confidence level used for confidence
ellipse(s).}

\item{lwd.groups}{Integer; width of the lines in groups' ellipses/hulls.}

\item{lty.groups}{Integer; type of the lines in groups' ellipses/hulls.}

\item{alpha.groups}{Numeric; transparency factor for groups'
ellipses/hulls/density distributions.}

\item{boxplot.groups}{Logical; whether to plot represent categorical data
using boxplots. If \code{FALSE}, violin plots are displayed instead.
Used only if \code{x} or \code{y} are provided.}

\item{density.groups}{Logical; whether to add density distribution for
groups (univariate ordinations only).}

\item{labels.groups}{Either logical, indicating whether to include all group
labels, or a character string with the exact names of the groups whose
labels should be included.}

\item{col.phylo}{Color of the lines depicting phylogenetic
branches.}

\item{lwd.phylo}{Integer; width of the lines depicting phylogenetic
branches.}

\item{lty.phylo}{Integer; type of the lines depicting phylogenetic
branches.}

\item{labels.nodes}{Either logical, indicating whether to include node
labels, or a character string with the exact names of the nodes whose
labels should be included (with the form of "node_n", e.g., "node_14"
corresponds to the root of a tree with 13 tips).}

\item{col.nodes}{Color of the scatterpoints representing the nodes of the
phylogeny.}

\item{pch.nodes}{Symbol of the scatterpoints representing the nodes of
the phylogeny.}

\item{bg.nodes}{Background color of the scatterpoints representing the
nodes of the phylogeny.}

\item{cex.nodes}{Numeric; size of the scatterpoints representing the
nodes of the phylogeny.}

\item{col.tips}{Color of the scatterpoints representing the tips of the
phylogeny.}

\item{pch.tips}{Symbol of the scatterpoints representing the tips of the
phylogeny.}

\item{bg.tips}{Background color of the scatterpoints representing the
tips of the phylogeny.}

\item{labels.tips}{Either logical, indicating whether to include tip labels,
or a character string with the exact names of the tips whose labels
should be included.}

\item{cex.tips}{Numeric; size of the scatterpoints representing the tips
of the phylogeny.}

\item{type.axis}{Integer; type of arrows (\code{0} = no arrow;
\code{1} = pointing towards the maximum; \code{2} = pointing towards
the maximum, \code{3} = pointing in both directions).}

\item{lwd.axis}{Integer; width of the lines depicting a morphometric
axis.}

\item{lty.axis}{Integer; type of the lines depicting a morphometric
axis.}

\item{col.axis}{Color of the lines depicting a morphometric axis.}

\item{palette.landsc}{A function defining a color palette to use for
landscape representation.}

\item{display.landsc}{How to display landscape representation; options
are \code{"contour"} and \code{"filled.contour"}. For bivariate
landscapes only.}

\item{alpha.landsc}{Numeric; transparency factor for filled contours
depicting landscapes.}

\item{nlevels.landsc}{Number of levels (i.e., contours) to use in
landscape representation.}

\item{drawlabels.landsc}{Logical; should the labels indicating the value
of each surface contour be plotted?}

\item{lty.landsc}{Integer; type of the contour lines depicting
landscapes.}

\item{lwd.landsc}{Integer; width of the contour lines depicting
landscapes.}
}
\value{
None
}
\description{
Flexible representation of morphospaces, including their
combination with either other variables or a phylogeny.
}
\details{
This function allows to regenerate/tweak morphospaces contained in
\code{"mspace"} objects already in existence. By default, \link{plot_mspace}
regenerates the morphospace with all its projected elements, preserving
graphical parameters used originally during the \code{\link{mspace}} +
\code{proj_*} pipeline (and stored in \code{mspace$plotinfo}). However,
all the graphical parameters can be modified to customize  representation.
Also, \link{plot_mspace} can be used to add a legend and/or a scalebar to aid
identification of groups and interpretation of landscapes, respectively.

In addition, this function expands the range of graphical options available
beyond 'pure' morphospaces. If a non-shape variable (assumed to be measured
for the same specimens in \code{mspace$projected$scores}) is fed to one of
the \code{x} or \code{y} arguments, a 'hybrid' morphospace is produced
(i.e. the bivariate plot will be constructed from the combination of
\code{x} or \code{y} and a morphometric axis; background shape models
will represent variation only for the latter)If a \code{"phylo"} object
(assumed to describe the phylogenetic relationships among tips scores
stored in \code{mspace$projected$phylo_scores}) is provided instead for
either \code{x} or \code{y}, a vertical or horizontal phenogram will be
deployed (the function assumes the phylogenetic tree is time-calibrated).

\emph{Note}: when regenerating landscapes, it's important to keep in mind
that this surface has been calculated for a specific set of shapes and
ordination axes. Hence, its regeneration using \link{plot_mspace} is only
warranted if the axes being depicted are the same than those used when the
surface landscape was originally computed using \code{\link{mspace}} +
\code{\link{proj_landscape}}. The only exception is when one of the
original axes (i.e., those specified with the \code{axes} argument in
\code{\link{mspace}}) is dropped (i.e., not specified with the \code{axes}
argument of \code{\link{plot_mspace}}). This will result in the collapse of
the 3D landscape projected into a bivariate morphospace into a 2D landscape
projected into a univariate one.
}
\examples{
\dontrun{
#load and extract relevant data, packages and information
library(magrittr)
data("tails")
shapes <- tails$shapes
species <- tails$data$species
sizes <- tails$sizes
sp_shapes <- expected_shapes(shapes, species)
sp_sizes <- cbind(tapply(sizes, species, mean))
tree <- tails$tree
links <- tails$links

#generate basic morphospace, add sampled shapes, species classification,
#phylogenetic structure and performance landscape
msp <- mspace(shapes, axes = c(1,2), plot = FALSE) \%>\%
 proj_shapes(shapes = shapes) \%>\%
 proj_groups(groups = species) \%>\%
 proj_phylogeny(shapes = sp_shapes, tree = tree) \%>\%
 proj_landscape(FUN = morphospace:::computeLD)


##Regenerating/modifying 'pure' morphospaces:

#plot mspace object as it is
plot_mspace(msp)

#remove landscape
plot_mspace(msp, landsc = FALSE)

#add colors for points, by species
plot_mspace(msp, col.points = species, landsc = FALSE,
            col.groups = 1:nlevels(species))

#add links for landmark configurations
plot_mspace(msp, links = links, landsc = FALSE,
            col.points = species, col.groups = 1:nlevels(species))

#change number and sizes of shape models in the background
plot_mspace(msp, nh = 2, nv = 2, links = links,
            size.models = 0.5, landsc = FALSE,
            col.points = species, col.groups = 1:nlevels(species))

#magnify deformation and highlight landmarks
plot_mspace(msp, mag = 1.5, nh = 2, nv = 2, links = links,
            size.models = 0.5, col.points = species, landsc = FALSE,
            col.groups = 1:nlevels(species), cex.ldm = 5, col.ldm = "red")

#change axes 1,2 for 1,3
plot_mspace(msp, axes = c(1,3), mag = 1.5, nh = 2, nv = 2, links = links,
            size.models = 0.5, col.points = species,
            col.groups = 1:nlevels(species), cex.ldm = 5, col.ldm = "red")

#change colors for as tips and nodes of the phylogeny
plot_mspace(msp, axes = c(1,3), mag = 1.5, nh = 2, nv = 2, links = links,
            size.models = 0.5, cex.ldm = 5, col.ldm = "red",
            col.tips = "red", col.nodes = "blue")

#plot only first PC axis, with general distribution of specimens
plot_mspace(msp, axes = 1, mag = 1.5, nh = 2, nv = 2, links = links,
            size.models = 0.5, cex.ldm = 5, density.groups = FALSE,
            col.ldm = "red", col.tips = "red", col.nodes = "blue")

#plot only first PC axis, with distribution of specimens by group
plot_mspace(msp, axes = 1, mag = 1.5, nh = 2, nv = 2, links = links,
            alpha.groups = 0.5, size.models = 0.5, cex.ldm = 5,
            density.groups = TRUE, col.ldm = "red", col.tips = "red",
            col.nodes = "blue")

#add legend
plot_mspace(msp, axes = 1, mag = 1.5, nh = 2, nv = 2, links = links,
            alpha.groups = 0.5, size.models = 0.5, cex.ldm = 5,
            density.groups = TRUE, col.ldm = "red", col.tips = "red",
            col.nodes = "blue", legend = TRUE, cex.legend = 1.3,
            pch.groups = 16)


#remove all the other elements to see only shape variation captured by the
#first two PCs
plot_mspace(msp, axes = c(1,2), links = links,
            col.points = species, col.groups = 1:nlevels(species),
            points = FALSE, groups = FALSE, phylo = FALSE)


##Plotting 'hybrid' morphospaces:

#plot size against first PC
plot_mspace(msp, x = sizes,  axes = 1, links = links, col.points = species,
            col.groups = 1:nlevels(species), pch.points = 16,
            xlab = "Centroid size")

#plot species agains first PC (violin plot, horizontal)
plot_mspace(msp, y = species,  axes = 1, links = links, col.points = species,
            col.groups = 1:nlevels(species), pch.points = 16)

#plot species agains first PC (boxplot, vertical)
plot_mspace(msp, x = species,  axes = 1, links = links, col.points = species,
            col.groups = 1:nlevels(species), pch.points = 16)


##Plotting phenograms:

#plot vertical phenogram against PC1
plot_mspace(msp, y = tree,  axes = 1, links = links, col.points = species,
            col.groups = 1:nlevels(species), pch.points = 16,
            xlab = "Branch lengths", pch.tips = 21, bg.tips = "red",
            pch.nodes = 21, bg.nodes = "blue")

#plot horizontal phenogram against PC2
plot_mspace(msp, x = tree,  axes = 2, links = links, col.points = species,
            col.groups = 1:nlevels(species), pch.points = 16,
            ylab = "Branch lengths", pch.tips = 21, bg.tips = "red",
            pch.nodes = 21, bg.nodes = "blue")
}
}
