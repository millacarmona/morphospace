% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shapes_operations.R
\name{extract_shapes}
\alias{extract_shapes}
\title{Retrieve / compute shapes from an existing morphospace}
\usage{
extract_shapes(
  mspace,
  axis = NULL,
  nshapes = NULL,
  scores = NULL,
  range = NULL,
  keep.template = TRUE,
  mag = 1
)
}
\arguments{
\item{mspace}{An \code{"mspace"} object.}

\item{axis}{Optional integer indicating an axis along which shapes should
be sampled.}

\item{nshapes}{Optional integer indicating the number of shapes the user
wishes to extract.}

\item{scores}{An optional vector of length 2 or 2-column matrix indicating
the (x,y) coordinates in the morphospace that the user wishes to extract
as shapes. If \code{NULL}, a new device will open and the user will be
asked choose the coordinates interactively. Ignored when \code{axis} is
provided.}

\item{range}{Optional integer vector of length 2, indicating the range of
values the axis should be sampled over.}

\item{keep.template}{Logical; should warped templates be returned as well?}

\item{mag}{Optional numeric indicating a magnifying factor for representing
shape transformation.}
}
\value{
A list containing sampled shapes (\code{$shapes}), as well as their
associated templates (\code{$templates}), when warranted.
}
\description{
Extracts shapes from \code{"mspace"} objects in various ways
(background shape models, specific axes, or particular coordinates - either
provided or chosen interactively by the user).
}
\details{
This function provides the user with an easy way to extract
theoretical shapes from an existing morphospace. If only an \code{"mspace"}
object is provided, the set of background shape models (optionally
amplified by a factor of \code{mag}) will be returned.

If \code{axis} is provided, a sample of \code{nshapes} shapes computed at
regular intervals along the specified ordination axis (either over its
empirical range -optionally amplified by a factor of \code{mag}- or, if
provided, between the extremes of \code{range}) will be returned.

If \code{axis = NULL}, this function will let the user to select arbitrary
coordinates in the morphospace to be back-transformed into shapes. There
are two alternatives: 1) if \code{scores = NULL} (the default option) the
user will be asked to interactively select the location(s) of
\code{nshapes} points in a new graphical device. 2) Otherwise, the set of
shapes represented by the (x,y) coordinates provided in \code{scores} will
be returned.
}
\examples{
#load all the relevant data and packages
library(Morpho)
library(geomorph)

data("tails")
shapes <- tails$shapes
sizes <- tails$sizes
species <- tails$data$species
type <- tails$data$type
links <- tails$links
sp_shapes <- expected_shapes(shapes, species)
tree <- tails$tree

#build phylomorphospace
phylomsp <- mspace(shapes, links = links) \%>\%
  proj_phylogeny(sp_shapes, tree = tree)


##Extracting background shape models

#extract background shape models
background_shapes <- extract_shapes(phylomsp)

#pile shapes and visualise the corresponding coordinates sampled in the
#morphospace
pile_shapes(background_shapes$shapes, links = links)

plot_mspace(phylomsp, phylo = FALSE)
background_scores <- proj_eigen(two.d.array(background_shapes$shapes),
                                phylomsp$ord$rotation, phylomsp$ord$center)
points(background_scores, pch = 21, bg = "red")


##Sampling a particular ordination axis

#extract shapes along PC2
PC2_shapes <- extract_shapes(phylomsp, axis = 2, nshapes = 8)

#pile shapes and visualise the corresponding coordinates sampled in the
#morphospace
pile_shapes(PC2_shapes$shapes, links = links, mshape = FALSE)

plot_mspace(phylomsp, phylo = FALSE)
PC2_scores <- proj_eigen(two.d.array(PC2_shapes$shapes),
                         phylomsp$ord$rotation, phylomsp$ord$center)
points(PC2_scores, pch = 21, bg = "blue")


##Sampling particular (x,y) locations

#1. Interactively

if (interactive()) {

#select 1 shape in the new window
arbitrary_shape <- extract_shapes(phylomsp, nshapes = 1)

#plot shape and visualise the corresponding coordinates sampled in the
#morphospace
plot(arbitrary_shape$shapes[,,1], axes = FALSE, xlab = "", ylab = "")
lineplot(arbitrary_shape$shapes[,,1], links)

plot_mspace(phylomsp, phylo = FALSE)
arbitrary_scores <- proj_eigen(two.d.array(arbitrary_shape$shapes),
                                phylomsp$ord$rotation, phylomsp$ord$center)
points(arbitrary_scores, pch = 21, bg = "magenta")

}


#2. Specifying coordinates

#get scores of the nodes of the phylogeny for the first two PCs
nodes_scores0 <- phylomsp$projected$phylo_scores[14:25,1:2]

#extract shapes from morphospace
nodes_shapes <- extract_shapes(phylomsp, scores = nodes_scores0)

plot_mspace(phylomsp, phylo = TRUE)
nodes_scores <- proj_eigen(two.d.array(nodes_shapes$shapes),
                         phylomsp$ord$rotation, phylomsp$ord$center)
points(nodes_scores, pch = 21, bg = "green")
}
