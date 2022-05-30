
#############################################################################################

#' Pile shapes
#'
#' @description Superimpose all the shapes in the sample.
#'
#' @param shapes Shape data.
#' @param links An optional list with the indices of the coordinates defining
#'   the wireframe (following the format used in \code{Morpho}).
#' @param mshape Logical; whether to plot the mean configuration.
#' @param ... Additional arguments passed to [plot()] or [rgl::points3d()].
#'
#' @export
#'
#' @examples
#' #load 2D landmark data
#' data("tails")
#' shapes <- tails$shapes
#' links <- tails$links
#'
#' #pile shapes
#' pile_shapes(shapes, mshape = FALSE) #bare
#' pile_shapes(shapes, mshape = FALSE, links = links) #with links
#' pile_shapes(shapes, mshape = TRUE, links = links) #with links and mean shape
#'
#'
#' #load 3D landmark data
#' data("shells3D")
#' shapes <- shells3D$shapes
#'
#' \dontrun{
#' #pile shapes
#' pile_shapes(shapes, mshape = FALSE) #bare
#' pile_shapes(shapes, mshape = FALSE, links = list(1:10)) #with false links (just as an example)
#' pile_shapes(shapes, mshape = TRUE, links = list(1:10)) #with false links and mean shape
#' }
#'
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

    if(ncol(shapes_coord) == 2) {
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
      rgl::plot3d(longlist, aspect = FALSE, axes = FALSE, col = "#708095",
                  xlab = "", ylab = "", zlab = "", size = 5, ...)

      if(!is.null(links)) {
        for(i in 1:dim(shapes_coord)[3]){
          for(l in 1:length(links)) rgl::lines3d(shapes_coord[links[[l]],,i], col = "gray", ...)
        }
      }

      if(mshape == TRUE) {
        rgl::points3d(consensus(shapes_coord), size = 10)
        if(!is.null(links)) {
          for(l in 1:length(links)) rgl::lines3d(consensus(shapes_coord)[links[[l]],])
        }
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
#'   created using the generic [plot()] function. Used internally (mostly).
#'
#' @param xy Coordinates of the scatterplot.
#' @param fac A factor grouping data points.
#' @param col A vector (either character or numeric) indicating the colors used
#'   for each group.
#' @param lty A vector (either character or numeric) indicating the type of line
#'   used to draw hulls.
#'
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
hulls_by_group_2D <- function(xy, fac, col = 1:nlevels(fac), lty = 1, ...) {

  if(length(col) == 1) col <- rep(col, nlevels(fac))
  if(length(lty) == 1) lty <- rep(lty, nlevels(fac))

  for(i in 1:nlevels(fac)) {
    x <- xy[fac == levels(fac)[i], 1]
    y <- xy[fac == levels(fac)[i], 2]
    hullp <- grDevices::chull(x = x, y = y)
    graphics::polygon(x[hullp], y[hullp], border = col[i], lty = lty[i], ...)
    }
}


#############################################################################################

#' Plot 3D convex hulls for a series of groups
#'
#' @description Plot convex hulls for different groups in 3D scatterplots
#'   created using \code{rgl}. Used internally (mostly).
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
#' @references
#' Claude, J. (2008). \emph{Morphometrics with R}. Springer Science & Business Media,
#' 316.
#'
#' @examples
#' #generate template interactively
#' \dontrun{
#' temp <- build_template2d(image = "extdata/sample_wing.jpg", nlands = 9, ncurves = 9)
#'
#' plot(temp, type = "n", asp = 1)
#' points(temp[c(1:9),], col = "red", pch = 16)
#' points(temp[-c(1:9),], type = "l")
#' }
build_template2d <- function(image, nlands, ncurves) {

  if(any(grepl(x = image, pattern = ".jpg"),
         grepl(x = image, pattern = ".jpeg"),
         grepl(x = image, pattern = ".JPG"),
         grepl(x = image, pattern = ".JPEG"))) {ras <- jpeg::readJPEG(image)}

  if(any(grepl(x = image, pattern = ".png"),
         grepl(x = image, pattern = ".PNG"))) {ras <- png::readPNG(image)}


  plot(c(1, dim(ras)[2]), c(1, dim(ras)[1]), type = "n", xlab = "", ylab = "",
       asp = 1, axes = FALSE)
  graphics::rasterImage(ras, 1, 1, dim(ras)[2], dim(ras)[1])

  lands <- graphics::locator(nlands, type = "p", pch = 8, col = "white")
  cat(paste0("Place the ", nlands, " landmarks"))

  curves <- lapply(1:ncurves, function(i) {
    cat(paste0("\nWhen the ", i,
               " curve is ready click Finish (top-right corner of Plots pane) or enter <Esc> in the console"))
    graphics::locator(type = "l", pch = 8, col = i)
  })

  template <- cbind(c(lands$x, unlist(lapply(curves, function(x) {c(NA, x$x)}))),
                    c(lands$y, unlist(lapply(curves, function(x) {c(NA, x$y)}))))

  return(template)
}
