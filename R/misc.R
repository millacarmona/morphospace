
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
      inv_efourier(coe = dat$data2d[i,], nb.pts = 300)
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
        for(l in 1:length(links)) graphics::lines(shapes_coord[links[[l]],,i], col = "gray")
      }
    }

    graphics::points(longlist, col = "#708095")

    if(mshape == TRUE) {
      graphics::points(consensus(shapes_coord), pch = 16, ...)
      if(!is.null(links)) {
        for(l in 1:length(links)) graphics::lines(consensus(shapes_coord)[links[[l]],])
      }
    }
  } else {
    plot(longlist, axes = FALSE, xlab = "", ylab = "", type = "n")

    for(i in 1:dim(shapes_coord)[3]) graphics::lines(shapes_coord[,,i], col = "#708095")

    if(mshape == TRUE) {
      graphics::lines(consensus(shapes_coord), lwd = 4)
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
#' #load landmark data and nencessary packages
#' library(geomorph)
#' data("tails")
#' shapes <- tails$shapes
#' species <- tails$data$species
#'
#' #perform PCA
#' pca <- prcomp(two.d.array(shapes))
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
#' @param ... Further arguments passed to [rgl::rgl.triangles()] (e.g.
#'   \code{specular}, \code{alpha}).
#'
#' @export
#'
#' @examples
#' #load landmark data and necessary packages
#' library(geomorph)
#' data("tails")
#' shapes <- tails$shapes
#' species <- tails$data$species
#'
#' #perform PCA
#' pca <- prcomp(two.d.array(shapes))
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



######################################################################################

#' Build template for 2D shape data
#'
#' @description Create a template (i.e. a set of curves describing the structure
#'   the landmarks are placed upon) to aid morphospace visualization, interactively.
#'
#' @param image Character; the path to an image of the structure, in png format.
#' @param nlands Numeric; the number of landmarks to be placed.
#' @param ncurves Numeric; the number of curves to be drawn.
#'
#' @details This functions let the user create a template interactively. The user
#'   will be first asked to place the landmarks (the same place and order than the
#'   shape data of interest). Once \code{nland} landmarks have been placed, the
#'   user will be asked to place the number of curves specified in \code{ncurves};
#'   the number of coordinates used to drawn each is arbitrary, and is up the user
#'   to decide when the curve is ready (press the 'Finish' button in the top-right
#'   corner of the Plots pane).
#'
#' @return A 2-column matrix of the landmarks followed by the coordinates defining
#'   the curves drawn, separated by \code{NA}s.
#'
#' @export
#'
#' @examples
#' #generate template interactively
#'
#' \dontrun{
#' temp <- build_template2d(image, nlands = 9, ncurves = 9)
#' }
#'
#' plot(temp, type = "n", asp = 1)
#' points(temp[c(1:9),], col = "red", pch = 16)
#' points(temp[-c(1:9),], type = "l)
build_template2d <- function(image, nlands, ncurves) {

  im <- png::png(image)

  plot(c(1, dim(im)[2]), c(1, dim(im)[1]), type = "n", xlab = "", ylab = "",
       asp = 1, axes = FALSE)
  rasterImage(im, 1, 1, dim(im)[2], dim(im)[1])

  lands <- locator(nlands, type = "p", pch = 8, col = "white")
  cat(paste0("Place the ", nlands, " landmarks"))

  curves <- lapply(1:ncurves, function(i) {
    cat(paste0("When the ", i,
               " curve is ready click Finish (top-right corner of Plots pane) or enter <Esc> in the console"))
    locator(type = "l", pch = 8, col = i)
  })

  template <- cbind(c(lands$x, unlist(lapply(curves, function(x) {c(NA, x$x)}))),
                    c(lands$y, unlist(lapply(curves, function(x) {c(NA, x$y)}))))

  return(template)
}


