% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shapes_operations.R
\name{rotate_coords}
\alias{rotate_coords}
\title{Rotate x,y coordinates}
\usage{
rotate_coords(xy, degrees)
}
\arguments{
\item{xy}{(x,y) coordinates}

\item{degrees}{Numeric; angle to rotate \code{(x,y)}}
}
\description{
Rotate x,y coordinates by an arbitrary amount of degrees
}
\examples{
#load data
data(wings)

shapes <- wings$shapes
links <- wings$links

# rotate first shape from the data set
rot_shape <- rotate_coords(shapes[,,1], 180)

#plot and compare
plot(shapes[,,1], pch = 1)
Morpho::lineplot(shapes[,,1], links)

points(rot_shape, pch = 16, col = "red")
Morpho::lineplot(rot_shape, links, col = "red")
}
