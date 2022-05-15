
#############################################################################################

#' Pile shapes
#'
#' @description Superimpose all the shapes in the sample.
#'
#' @param shapes Shape data.
#' @param links An optional list with the indices of the coordinates defining
#'   the wireframe (following the format used in \code{Morpho}).
#' @param mshape Logical; whether to plot the mean configuration.
#' @param ... Additional arguments passed to [plot()].
#'
#' @export
#'
#' @examples
#' #load landmark data
#' data("tails")
#' shapes <- tails$shapes
#' links <- tails$links
#'
#' #pile shapes
#' pile_shapes(shapes, mshape = FALSE) #bare
#' pile_shapes(shapes, mshape = FALSE, links = links) #with links
#' pile_shapes(shapes, mshape = TRUE, links = links) #with links and mean shape
#'
#' #load outline data
#' data("shells")
#' shapes <- shells$shapes
#' links <- shells$links
#'
#' #pile shapes
#' pile_shapes(shapes, mshape = FALSE) #bare
#' pile_shapes(shapes, mshape = TRUE, links = links) #with mean shape
pile_shapes <- function(shapes, links = NULL, mshape = TRUE, ...) {

  dat <- shapes_mat(shapes)

  if(dat$datype == "fcoef") {
    coords_l <- lapply(1:nrow(dat$data2d), function(i){
      inv_efourier(coe = dat$data2d[i,], nb.p = 300)
    })
    shapes_coord <- abind::abind(coords_l, along = 3)
  } else {
    shapes_coord <- shapes
  }

  longlist <- c()
  for(i in 1:dim(shapes_coord)[3]) longlist <- rbind(longlist, shapes_coord[,,i])

  if(dat$datype == "landm") {
    plot(longlist, type = "n", axes = FALSE, xlab = "", ylab = "", ...)

    if(!is.null(links)) {
      for(i in 1:dim(shapes_coord)[3]){
        for(l in 1:length(links)) lines(shapes_coord[links[[l]],,i], col = "gray")
      }
    }

    points(longlist, col = "#708095")

    if(mshape == TRUE) {
      points(consensus(shapes_coord), pch = 16, ...)
      if(!is.null(links)) {
        for(l in 1:length(links)) lines(consensus(shapes_coord)[links[[l]],])
      }
    }
  } else {
    plot(longlist, axes = FALSE, xlab = "", ylab = "", type = "n")

    for(i in 1:dim(shapes_coord)[3]) lines(shapes_coord[,,i], col = "#708095")

    if(mshape == TRUE) {
      lines(consensus(shapes_coord), lwd = 4)
    }
  }

}


#############################################################################################

#' Plot 2D convex hulls for a series of groups
#'
#' @description Plot convex hulls for different groups in 2D scatterplots
#'   created using the generic [plot()] function.
#'
#' @param xy Coordinates of the scatterplot.
#' @param fac A factor grouping data points.
#' @param col A vector (either character or numeric) indicating the colors used
#'   for each group.
#' @param ... Further arguments passed to [polygon()].
#'
#' @export
#'
#' @examples
#' #load landmark data
#' data("tails")
#' shapes <- tails$shapes
#' species <- tails$data$species
#'
#' #perform PCA
#' pca <- prcomp(geomorph::two.d.array(shapes))
#'
#' #plot and add convex hulls
#' plot(pca$x)
#' hulls_by_group_2D(pca$x, fac = species, col = "black")
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
#' @description Plot convex hulls for different groups in 3D scatterplots
#'   created using \code{rgl}.
#'
#' @param xyz Coordinates for the scatterplot
#' @param fac A factor grouping data points.
#' @param col A vector (either character or numeric) indicating the colors used
#'   for each group.
#' @param ... Further arguments passed to [rgl.triangles()] (e.g.
#'   \code{specular}, \code{alpha}).
#'
#' @export
#'
#' @examples
#' #load landmark data
#' data("tails")
#' shapes <- tails$shapes
#' species <- tails$data$species
#'
#' #perform PCA
#' pca <- prcomp(geomorph::two.d.array(shapes))
#'
#' #plot and add convex hulls
#' \dontrun{
#'
#' plot3d(pca$x)
#' hulls_by_group_3D(pca$x, fac = species, col = "gray")
#'
#' }
hulls_by_group_3D<-function(xyz, fac, col = 1:nlevels(fac), ...) {

  if(length(col) == 1) col <- rep(col, nlevels(fac))

  for(i in 1:nlevels(fac)) {
    matsp <- xyz[fac == levels(fac)[i], 1:3]
    surf <- t(geometry::convhulln(matsp))
    convex <- rgl::rgl.triangles(matsp[surf, 1], matsp[surf, 2], matsp[surf, 3], col = col[i], ...)
    }
}






