
#############################################################################################

#' Pile shapes
#'
#' @param shapes A \code{k x p x n} array of superimposed landmarks.
#' @param links An optional list containing the indices of the coordinates
#'   defining the wireframe in Morpho format.
#' @param mshape Logical; whether to plot the mean configuration.
#' @param ... Additional arguments passed to [plot()].
#'
#' @return
#' @export
#'
#' @examples
stack <- function(shapes, links = NULL, mshape = TRUE, ...){

  longlist <- c()
  for(i in 1:dim(shapes)[3]) longlist <- rbind(longlist, shapes[,,i])

  plot(longlist, col = "#708095", axes = FALSE, xlab = "", ylab = "", ...)
  if(mshape == TRUE) {
    points(consensus(shapes), pch = 16, ...)
    if(!is.null(links)) { for(l in 1:length(links)) lines(consensus(shapes)[links[[l]],]) }
  }
}


#############################################################################################

#' Plot 2D convex hulls for a series of groups
#'
#' @description Plot convex hulls for different groups in 2D scatterplots.
#'
#' @param xy Coordinates of the scatterplot.
#' @param fac A factor grouping data points.
#' @param col A vector (either character or numeric) indicating the colors used
#'   for each group.
#' @param ... Further arguments passed to [polygon()].
#'
#' @return
#' @export
#'
#' @examples
hulls_by_group_2D <- function(xy, fac, col = 1:nlevels(fac), ...) {

  if(length(col) == 1) col <- rep(col, nlevels(fac))

  for(i in 1:nlevels(fac)) {
    x <- xy[fac == levels(fac)[i], 1]
    y <- xy[fac == levels(fac)[i], 2]
    hullp <- grDevices::chull(x = x, y = y)
    graphics::polygon(x[hullp], y[hullp], border = col[i], ...)
    }
}


#############################################################################################

#' Plot 3D convex hulls for a series of groups
#'
#' @param xyz Coordinates for the scatterplot
#' @param fac A factor grouping data points.
#' @param col A vector (either character or numeric) indicating the colors used
#'   for each group.
#' @param ... Further arguments passed to [rgl.triangles()] (e.g.
#'   \code{specular}, \code{alpha}).
#'
#' @return
#' @export
#'
#' @examples
hulls_by_group_3D<-function(xyz, fac, col = 1:nlevels(fac), ...) {

  if(length(col) == 1) col <- rep(col, nlevels(fac))

  for(i in 1:nlevels(fac)) {
    matsp <- xyz[fac == levels(fac)[i], 1:3]
    surf <- t(geometry::convhulln(matsp))
    convex <- rgl::rgl.triangles(matsp[surf, 1], matsp[surf, 2], matsp[surf, 3], col = col[i], ...)
    }
}






