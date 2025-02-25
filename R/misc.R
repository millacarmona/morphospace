################################################################################

#' Pile shapes
#'
#' @description Superimpose all the shapes in the sample.
#'
#' @param shapes Shape data.
#' @param links An optional list with the indices of the coordinates defining
#'   the wireframe (following the format used in \code{Morpho}).
#' @param mshape Logical; whether to plot the mean configuration.
#' @param ... Additional arguments passed to \code{\link[graphics]{lines}} or
#'   \code{\link[rgl]{lines3d}}.
#'
#' @export
#'
#' @return None
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
#' pile_shapes(shapes, mshape = FALSE, links = list(1:10)) #with false links
#' pile_shapes(shapes, mshape = TRUE, links = list(1:10)) #false links + mshape
#' }
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
    coords_l <- lapply(seq_len(nrow(dat$data2d)), function(i) {
      inv_efourier(coe = dat$data2d[i,], nb.pts = 300)
    })
    shapes_coord <- abind::abind(coords_l, along = 3)
  } else {
    shapes_coord <- shapes
  }

  longlist <- NULL
  for(i in seq_len(dim(shapes_coord)[3])) longlist <- rbind(longlist, shapes_coord[,,i])

  if(dat$datype == "landm") {

    if(ncol(shapes_coord) == 2) {
      plot(longlist, type = "n", axes = FALSE, xlab = "", ylab = "")

      if(!is.null(links)) {
        for(i in seq_len(dim(shapes_coord)[3])) {
          for(l in seq_len(length(links))) graphics::lines(shapes_coord[links[[l]],,i],
                                                           col = "gray", ...)
        }
      }

      graphics::points(longlist, col = "#708095")

      if(mshape == TRUE) {
        graphics::points(expected_shapes(shapes_coord), pch = 16)
        if(!is.null(links)) {
          for(l in seq_len(length(links))) graphics::lines(expected_shapes(shapes_coord)[links[[l]],])
        }
      }

    } else {
      rgl::plot3d(longlist, aspect = FALSE, axes = FALSE, col = "#708095",
                  xlab = "", ylab = "", zlab = "", size = 5)

      if(!is.null(links)) {
        for(i in seq_len(dim(shapes_coord)[3])) {
          for(l in seq_len(length(links))) rgl::lines3d(shapes_coord[links[[l]],,i], col = "gray", ...)
        }
      }

      if(mshape == TRUE) {
        rgl::points3d(expected_shapes(shapes_coord), size = 10)
        if(!is.null(links)) {
          for(l in seq_len(length(links))) rgl::lines3d(expected_shapes(shapes_coord)[links[[l]],])
        }
      }

    }

  } else {
    plot(longlist, axes = FALSE, xlab = "", ylab = "", type = "n")

    for(i in seq_len(dim(shapes_coord)[3])) graphics::lines(shapes_coord[,,i], col = "#708095", ...)

    if(mshape == TRUE) {
      graphics::lines(expected_shapes(shapes_coord), lwd = 4)
    }
  }
}


################################################################################

#' Build template for 2D shape data
#'
#' @description Create a template (i.e., a set of curves describing the
#'   structure the landmarks are placed upon) to aid morphospace visualization,
#'   interactively.
#'
#' @param image Character; the path to an image of the structure, in png format.
#' @param nlands Integer; the number of landmarks to be placed.
#' @param ncurves Integer; the number of curves to be drawn.
#'
#' @details This functions let the user create a template interactively. The
#'   user will be first asked to place the landmarks (the same place and order
#'   than the shape data of interest). Once \code{nland} landmarks have been
#'   placed, the user will be asked to place the number of curves specified in
#'   \code{ncurves}; the number of coordinates used to drawn each is arbitrary,
#'   and is up the user to decide when the curve is ready (press the 'Finish'
#'   button in the top-right corner of the Plots pane).
#'
#' @return A 2-column matrix with the landmark configuration (standardized for
#'   scale and position) followed by the coordinates defining the curves drawn,
#'   separated by \code{NA}s.
#'
#' @export
#'
#' @references
#' Claude, J. (2008). \emph{Morphometrics with R}. Springer Science & Business
#'   Media, 316.
#'
#' @examples
#' #generate template interactively
#' if (interactive()) {
#' temp <- build_template2d(image = "extdata/sample_wing.jpg", nlands = 9,
#'                          ncurves = 9)
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

  curves <- lapply(seq_len(ncurves), function(i) {
    cat(paste0("\nWhen the ", i,
               " curve is ready click Finish (top-right corner of Plots pane) or enter <Esc> in the console"))
    graphics::locator(type = "l", pch = 8, col = i)
  })

  template <- cbind(c(lands$x, unlist(lapply(curves, function(x) {c(NA, x$x)}))),
                    c(lands$y, unlist(lapply(curves, function(x) {c(NA, x$y)}))))
  template_trans <- t(t(rbind(template)) - colMeans(template, na.rm = TRUE))
  template_trans_scald <- template_trans / Morpho::cSize(template_trans[seq_len(nlands),])

  return(template_trans_scald)
}


################################################################################

#' Plot 3D convex hulls for a series of groups
#'
#' @description Plot convex hulls for different groups in 3D scatterplots
#'   created using \code{rgl}. Used internally (mostly).
#'
#' @param xyz Coordinates for the scatterplot
#' @param fac A factor grouping data points.
#' @param col A vector (either character or numeric) indicating the colors used
#'   for each group.
#' @param ... Further arguments passed to [rgl::triangles3d()] (e.g.
#'   \code{specular}, \code{alpha}).
#'
#' @return None
#'
#' @export
#' @keywords internal
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
#' if (interactive()) {
#'
#' plot3d(pca$x)
#' hulls_by_group_3D(pca$x, fac = species, col = "gray")
#'
#' }
hulls_by_group_3D <- function(xyz, fac, col = seq_len(nlevels(fac)), ...) {

  if(length(col) == 1) col <- rep(col, nlevels(fac))

  for(i in seq_len(nlevels(fac))) {
    matsp <- xyz[fac == levels(fac)[i], 1:3]
    surf <- t(geometry::convhulln(matsp))
    convex <- rgl::triangles3d(matsp[surf, 1], matsp[surf, 2], matsp[surf, 3], col = col[i], ...)
  }
}


################################################################################

#' Express colors as hexadecimal coding
#'
#' @description Little function intended for internal use; will transform colors
#'   (expressed either as numerical or characters) into hexadecimal code.
#'
#' @param col Either numeric or character, specifying the color(s) to be
#'    transformed.
#'
#' @return The color(s) used as input, but expressed as hexadecimal code(s).
#'
#' @export
#' @keywords internal
#'
#' @examples
#' plot(rnorm(n = 100), pch = 16, col = "red")
#' plot(rnorm(n = 100), pch = 16, col = col2hex("red"))
#'
#' plot(rnorm(n = 100), pch = 16, col = 2)
#' plot(rnorm(n = 100), pch = 16, col = col2hex(2))
col2hex <- function(col) {

  rgbcols <- grDevices::col2rgb(col)
  hexcols <- NULL
  for(i in seq_len(ncol(rgbcols))) {
    hexcols[i] <- grDevices::rgb(rgbcols[1,i], rgbcols[2,i], rgbcols[3,i],
                                 maxColorValue = 255)
  }
  return(hexcols)
}


################################################################################

#' Compute lift-to-drag ratio for shapes from the 'tails' data set
#'
#' @description Compute lift-to-drag ratio (a proxy for aerodynamic performance
#'   on shapes of tails of Tyrannus species (landmark data). To be used as a
#'   complement for the 'tails' data set in tests and examples.
#'
#' @param model landmark shape data from the "tails" data set.
#' @param metric Character, specifying whether "LD" (lift/drag ratio) or "MD"
#'   (moment/drag ratio) should be computed.
#' @param center Logical, whether to center center tail configuration
#'
#' @noRd
flight_performance <- function(model, metric = "LD", center = TRUE) {

  #1. prepare tail 'model': format as landmark configuration, then rotate
  model <- scale(matrix(c(model), ncol = 2, byrow = TRUE),
                 scale = FALSE, center = center)

  ang <- compute_angle(colMeans(model[1:4,]), model[7,]) #get angle bw center of the base and central rectrix
  rotmodel <- rotate_coords(model, -ang) #rotate to make it horizontal

  #if the base is not at the left (i.e., negative x values) flip around x-axis
  base <- rotmodel[1:4,]
  if(sign(mean(base)) != -1) rotmodel <- cbind(rotmodel[,1] *-1, rotmodel[,2])
  mcsx <- min(rotmodel[c(5:9),1]) #get position of MCS along x-axis

  #2. get different relevant landmarks:
  #2.1 - points at the outer margin marking the MCS, and the distance bw them
  mcs_points <- rbind(find_intersection(cbind(mcsx, range(rotmodel[,2])), rotmodel[c(4,9),]),
                      find_intersection(cbind(mcsx, range(rotmodel[,2])), rotmodel[c(3,5),]))
  mcs_dist <- stats::dist(mcs_points)

  #2.2 - get chordwise points ( )
  cw_points <- rbind(find_intersection(mcs_points,
                                       rbind(colMeans(rotmodel[1:4,]), rotmodel[7,])),
                     colMeans(rotmodel[1:4,]))
  cw_dist <- stats::dist(cw_points)


  #3. calculate aerodynamic parameters (From Thomas & Balmford 1995)

  #constants
  p <- 1
  a <- 1
  U <- 1
  Cdf <- 1

  #tail area
  S <- Momocs::coo_area(rotmodel[c(3,5,6,7,8,9,4),])

  #lift generated by the tail
  lift <- ((pi/4) * a * (U^2) * (mcs_dist^2))# / model_area

  #induced drag of the tail (produced by generating lift)
  drag.i <- (1/2) * lift * a

  #profile drag of the tail (produced by the friction of tail against air)
  drag.f <- (1/2) * p * (U^2) * S * Cdf

  #moment from tail apex to tails aerodynamic center
  moment.apex <- (lift * (2/3) * cw_dist)# /  model_area

  if(metric == "LD") value <- lift / (drag.i + drag.f)
  if(metric == "MD") value <- moment.apex / (drag.i + drag.f)

  return(value)

}
computeLD <- function(model) {flight_performance(model, metric = "LD")}
computeMD <- function(model) {flight_performance(model, metric = "MD")}


################################################################################

#' Calculate angle between two vectors
#'
#' @description Calculate angle between two vectors in a bivariate space,
#'   in degrees.
#'
#' @param p1,p2 Vectors of length 2 containing the vectors in 2D space.
#'
#' @noRd
compute_angle <- function(p1, p2) {
  dx <- p2[1] - p1[1]
  dy <- p2[2] - p1[2]
  rad <- atan2(dy, dx)
  angle <- rad * (180 / pi)

  if(sign(angle) == -1) angle <- 180 + diff(c(-180, angle))

  return(angle)
}

################################################################################

#' Find intersection between two lines
#'
#' @description Find point of intersection between two non-parallel lines.
#'
#' @param line1,line2 Matrices defining (x,y) values for points at the
#'   extreme of a straight line.
#'
#' @noRd
find_intersection <- function(line1, line2) {

  p1 <- line1[1,]
  p2 <- line1[2,]
  p3 <- line2[1,]
  p4 <- line2[2,]

  # Check if first line is vertical
  if (p1[1] == p2[1]) {
    x_vert1 <- p1[1]
    vertical1 <- TRUE
  } else {
    m1 <- (p2[2] - p1[2]) / (p2[1] - p1[1])
    b1 <- p1[2] - m1 * p1[1]
    vertical1 <- FALSE
  }

  # Check if second line is vertical
  if (p3[1] == p4[1]) {
    x_vert2 <- p3[1]
    vertical2 <- TRUE
  } else {
    m2 <- (p4[2] - p3[2]) / (p4[1] - p3[1])
    b2 <- p3[2] - m2 * p3[1]
    vertical2 <- FALSE
  }

  # Handle cases where one or both lines are vertical
  if (vertical1 && vertical2) {
    if (x_vert1 == x_vert2) {
      return("The lines are coincident and do not intersect at a unique point.")
    } else {
      return("The lines are parallel and do not intersect.")
    }
  } else if (vertical1) {
    x <- x_vert1
    y <- m2 * x + b2
  } else if (vertical2) {
    x <- x_vert2
    y <- m1 * x + b1
  } else {
    # Check if the lines are parallel (no intersection)
    if (m1 == m2) {
      if (b1 == b2) {
        return("The lines are coincident and do not intersect at a unique point.")
      } else {
        return("The lines are parallel and do not intersect.")
      }
    }
    # Calculate intersection point
    x <- (b2 - b1) / (m1 - m2)
    y <- m1 * x + b1
  }

  return(c(x, y))

}
