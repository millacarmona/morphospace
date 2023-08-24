################################################################################

#' Pile shapes
#'
#' @description Superimpose all the shapes in the sample.
#'
#' @param shapes Shape data.
#' @param links An optional list with the indices of the coordinates defining
#'   the wireframe (following the format used in \code{Morpho}).
#' @param mshape Logical; whether to plot the mean configuration.
#' @param ... Additional arguments passed to [graphics::lines()] or
#'   [rgl::lines3d()].
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
#' \dontrun{
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

#' Plot 2D convex hulls for a series of groups
#'
#' @description Plot convex hulls for different groups in 2D scatterplots
#'   created using the generic [graphics::plot()] function. Used internally
#'   (mostly).
#'
#' @param xy Coordinates of the scatterplot.
#' @param fac A factor grouping data points.
#' @param col A vector (either character or numeric) indicating the colors used
#'   for each group.
#' @param lty A numeric vector indicating the type of line used to draw hulls.
#' @param alpha Numeric; transparency factor for hulls.
#' @param ... Further arguments passed to [graphics::polygon()].
#'
#' @seealso \code{\link{ellipses_by_group_2D}}, \code{\link{hulls_by_group_3D}}
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
#' plot(pca$x)
#' hulls_by_group_2D(pca$x, fac = species, col = "black")
hulls_by_group_2D <- function(xy, fac, col = seq_len(nlevels(fac)),
                              lty = 1, alpha = 0, ...) {

  col <- if(length(col) == 1) rep(col, nlevels(fac)) else {
    if(nlevels(fac) > length(col)) {
      c(rep(NA, nlevels(fac) - length(unique(col))),
        unique(col))[order(c((1:nlevels(fac))[-which(levels(fac) %in% unique(fac))],
                             which(levels(fac) %in% unique(fac))))]
    } else col
  }

  lty <- if(length(lty) == 1) rep(lty, nlevels(fac)) else if(nlevels(fac) > length(lty)) {
    c(rep(0, nlevels(fac) - length(lty)),
      lty)[order(c((1:nlevels(fac))[-which(levels(fac) %in% unique(fac))],
                   which(levels(fac) %in% unique(fac))))]
  } else lty


  for(i in seq_len(nlevels(fac))) {
    x <- xy[fac == levels(fac)[i], 1]
    y <- xy[fac == levels(fac)[i], 2]
    hullp <- grDevices::chull(x = x, y = y)
    graphics::polygon(x[hullp], y[hullp], border = col[i],
                      col = grDevices::adjustcolor(col[i], alpha.f = alpha), lty = lty[i], ...)
  }
}


################################################################################

#' Plot 2D confidence ellipses for a series of groups
#'
#' @description Plot confidence ellipses for different groups in 2D scatterplots
#'   created using the generic [graphics::plot()] function. Used internally
#'   (mostly).
#'
#' @param xy Coordinates of the scatterplot.
#' @param fac A factor grouping data points.
#' @param col A vector (either character or numeric) indicating the colors used
#'   for each group.
#' @param lty A numeric vector indicating the type of line used to draw
#'   ellipses.
#' @param alpha Numeric; transparency factor for ellipses.
#' @param conflev Numeric, specifying the confidence level for drawing ellipses.
#' @param ... Further arguments passed to [graphics::polygon()].
#'
#' @seealso \code{\link{hulls_by_group_2D}}, \code{\link[car]{ellipse}}
#'
#' @export
#'
#' @examples
#' #load landmark data and necessary packages
#' library(geomorph)
#' data("wings")
#' shapes <- wings$shapes
#' species <- wings$data$species
#'
#' #perform PCA
#' pca <- prcomp(two.d.array(shapes))
#'
#' #plot and add 95% confidence ellipses
#' plot(pca$x)
#' ellipses_by_group_2D(pca$x, fac = species, col = "black", conflev = 0.95)
ellipses_by_group_2D <- function(xy, fac, col = seq_len(nlevels(fac)),
                                 lty = 1, conflev = 0.95, alpha = 0, ...) {

  col <- if(length(col) == 1) rep(col, nlevels(fac)) else {
    if(nlevels(fac) > length(col)) {
      c(rep(NA, nlevels(fac) - length(unique(col))),
        unique(col))[order(c((1:nlevels(fac))[-which(levels(fac) %in% unique(fac))],
                             which(levels(fac) %in% unique(fac))))]
    } else col
  }

  lty <- if(length(lty) == 1) rep(lty, nlevels(fac)) else if(nlevels(fac) > length(lty)) {
    c(rep(0, nlevels(fac) - length(lty)),
      lty)[order(c((1:nlevels(fac))[-which(levels(fac) %in% unique(fac))],
                   which(levels(fac) %in% unique(fac))))]
  } else lty

  for(i in seq_len(nlevels(fac))) {
    cent <- colMeans(xy[fac == levels(fac)[i], 1:2])
    vcv <- stats::var(xy[fac == levels(fac)[i], 1:2])
    if(any(!is.na(vcv))) {
      ell <- car::ellipse(center = cent, shape = vcv,
                          radius = stats::qnorm((1 - conflev) / 2, lower.tail = F),
                          draw = FALSE)
      graphics::polygon(ell, border = col[i],
                        col = grDevices::adjustcolor(col[i], alpha.f = alpha), lty = lty[i], ...)
    }
  }
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
hulls_by_group_3D <- function(xyz, fac, col = seq_len(nlevels(fac)), ...) {

  if(length(col) == 1) col <- rep(col, nlevels(fac))

  for(i in seq_len(nlevels(fac))) {
    matsp <- xyz[fac == levels(fac)[i], 1:3]
    surf <- t(geometry::convhulln(matsp))
    convex <- rgl::triangles3d(matsp[surf, 1], matsp[surf, 2], matsp[surf, 3], col = col[i], ...)
  }
}


################################################################################

#' Plot univariate density distributions for a series of groups
#'
#' @description Plot density distribution for different groups in "univariate"
#'   scatterplots. Used internally.
#'
#' @param xy Coordinates of the scatterplot.
#' @param fac A factor grouping data points.
#' @param ax the axis of \code{xy} corresponding to the active variable.
#' @param lty A vector indicating the type of line (integer) used to draw density
#'   distributions.
#' @param lwd A vector indicating the width of line (integer) used to draw
#'   density distributions.
#' @param alpha Numeric; transparency factor for density distributions.
#' @param col A vector (either character or numeric) indicating the colors used
#'   for each group.
#'
#' @export
#'
#' @examples
#' #load Fourier data and necessary packages
#' library(geomorph)
#' data("shells")
#' shapes <- shells$shapes$coe
#' species <- shells$data$species
#'
#' #perform PCA
#' pca <- prcomp(shapes)
#'
#' #bind 1st axis with a column of 0s, plot and add density distributions
#' xy <- cbind(pca$x[,1], 0)
#' plot(xy, ylim = c(0,1))
#' density_by_group_2D(xy, fac = species, ax = 1)
density_by_group_2D <- function(xy, fac, ax, alpha = 0.2, lwd = 1, lty = 1,
                                col = seq_len(nlevels(fac))) {

  if(length(col) == 1) col <- rep(col, nlevels(fac))
  lty <- if(length(lty) == 1) rep(lty, nlevels(fac)) else if(nlevels(fac) > length(lty)) {
    c(rep(0, nlevels(fac) - length(lty)),
      lty)[order(c((1:nlevels(fac))[-which(levels(fac) %in% unique(fac))],
                   which(levels(fac) %in% unique(fac))))]
  } else lty

  dens <- lapply(seq_len(nlevels(fac)), function(i) {
    if(sum(fac == levels(fac)[i]) > 2) {
      subdens <- stats::density(xy[fac == levels(fac)[i], ax])
      list(x = subdens$x, y = subdens$y)
    }
  })
  ymax <- max(unlist(lapply(dens, function(x) {x$y})))


  for(i in seq_len(nlevels(fac))) {
    graphics::polygon(dens[[i]]$x, dens[[i]]$y / ymax, lwd = lwd, border = col[i],
                      lty = lty[i], col = grDevices::adjustcolor(col[i], alpha.f = alpha))
  }
}


################################################################################

#' Plot scatterpoints into univariate morphospace
#'
#' @description An \code{ad hoc} wrapper for [stats::density()] +
#'   [graphics::polygon()] / [graphics::points()]. Used internally for
#'   projecting scatterpoints into univariate morphospaces.
#'
#' @param scores a matrix of (x,y) coordinates to plot.
#' @param density Logical; whether to plot density curve.
#' @param col Color of the scatterpoints.
#' @param pch Symbol of the scatterpoints.
#' @param bg Background color for the scatterpoints.
#' @param cex Numeric; size of the scatterpoints.
#' @param ... Further arguments passed to [graphics::points()].
plot_univ_scatter <- function(scores, density, col = 1, bg = 1, pch = 1, cex = 1, ...) {

  if(is.null(dim(scores))) scores <- rbind(scores)

  if(density & nrow(scores) > 1) {
    dens <- stats::density(scores)
    graphics::polygon(dens$x, dens$y / max(dens$y), lwd = 2,
                      col = grDevices::adjustcolor(1, alpha.f = 0.5))
  }

  graphics::abline(h = 0)

  if(any(pch %in% c(21:25))) {
    graphics::points(cbind(scores, 0), pch = pch, bg = bg, cex = cex, ...)
  } else {
    graphics::points(cbind(scores, 0), pch = pch, col = col, cex = cex, ...)
  }
}


################################################################################

#' Plot scatterpoints into bivariate morphospace
#'
#' @description An \code{ad hoc} wrapper for [graphics::points()]. Used
#'   internally for projecting scatterpoints into bivariate morphospaces.
#'
#' @param scores a matrix of (x,y) coordinates to plot.
#' @param col Color of the scatterpoints.
#' @param pch Symbol of the scatterpoints.
#' @param bg Background color for the scatterpoints.
#' @param cex Numeric; size of the scatterpoints.
#' @param ... Further arguments passed to [graphics::points()].
plot_biv_scatter <- function(scores, col = 1, bg = 1, pch = 1, cex = 1, ...) {

  if(is.null(dim(scores))) scores <- rbind(scores)

  if(any(pch %in% c(21:25))) {
    graphics::points(scores, pch = pch, bg = bg, cex = cex, ...)
  } else {
    graphics::points(scores, pch = pch, col = col, cex = cex, ...)
  }
}


################################################################################

#' Plot bivariate landscape into morphospace
#'
#' @description An \code{ad hoc} wrapper for [graphics::contour()] /
#'   [graphics::.filled.contour()] /  [graphics::image()]. Used internally for
#'   projecting landscapes into bivariate morphospaces.
#'
#' @param landscape A list containing the results of [akima::interp()].
#' @param display Either \code{"contour"} or \code{"filled.contour"}.
#' @param type Type of landscape to be plotted. Options are \code{"theoretical"}
#'   or \code{"empirical"}.
#' @param levels Levels to be used to create the contours.
#' @param lty Integer; type of the lines depicting contours.
#' @param lwd Integer; width of the lines depicting contours.
#' @param col Colors used to represent the landscape contours.
#' @param drawlabels Logical; should the labels indicating the value of each
#'   surface contour be plotted?
#' @param alpha Transparency factor for filled contours.
plot_biv_landscape <- function(landscape, display, type, levels, lwd, lty, col, drawlabels, alpha) {
  if(display == "contour") {
    graphics::contour(landscape$x, landscape$y, landscape$z, levels = levels,
                      lwd = lwd, lty = lty, col = col, labels = round(levels, digits = 3),
                      drawlabels = drawlabels, add = TRUE)
    graphics::box()
  }

  if(display == "filled.contour") {
    if(type == "theoretical") {
      graphics::.filled.contour(landscape$x, landscape$y, landscape$z, levels = levels,
                                col = grDevices::adjustcolor(col = col, alpha = alpha))
    }
    if(type == "empirical") {
      graphics::image(landscape$x, landscape$y, landscape$z, add = TRUE,
                      col = grDevices::adjustcolor(col = col, alpha = alpha))
    }
  }
}


################################################################################

#' Plot univariate landscape into morphospace
#'
#' @description Plot a curve representing a "univariate" landscape. Used
#'   internally for projecting landscapes into univariate morphospaces.
#'
#' @param landscape A list containing the results of [akima::interp()].
#' @param drawlabels Logical; should the labels indicating the value of each
#'   surface contour be plotted?
#' @param col Colors used to represent the landscape curve.
#' @param lwd Integer; width of the lines depicting the landscape curve.
plot_univ_landscape <- function(landscape, drawlabels, col, lwd) {
  if(drawlabels) {
    w.transp <-round(stats::quantile(x = 1:length(landscape$z), probs = c(0.25, 0.5, 0.75)))
    w.transp <- sort(c(w.transp - 1, w.transp, w.transp + 1))

    label_col <- col[w.transp[c(2,5,8)]]
    col[w.transp] <- NA
  }

  for(i in 1:(length(landscape$x) - 1)) {
    graphics::lines(rbind(c(landscape$x[i], landscape$z[i]),
                          c(landscape$x[i + 1], landscape$z[i + 1])),
                    col = col[i], lwd = lwd)
  }
  graphics::box()

  if(drawlabels == TRUE) {
    x_text <- colMeans(matrix(landscape$x[w.transp + 1], nrow = 3))
    y_text <- colMeans(matrix(landscape$z[w.transp + 1], nrow = 3))
    labels <- round(colMeans(matrix(landscape$z[w.transp], nrow = 3)), 2)
    graphics::text(cbind(x_text, y_text), labels = labels, cex = 0.7, col = label_col)
  }
}


################################################################################

# #' Compute lift-to-drag ratio for shapes from the 'tails' data set
# #'
# #' @description Compute lift-to-drag ratio (a proxy for aerodynamic performance
# #'   on shapes of tails of Tyrannus species (landmark data). To be used as a
# #'   complement for the 'tails' data set in tests and examples.
# #'
# #' @param shapes landmark shape data from the "tails" data set.
# #' @param MSC Logical; whether to use the maximum continuous span of the tail
# #'   as a proxy for the amount of lift produced, instead of the lifting area.

computeLD <- function(model, MCS = FALSE) {

  tail <- matrix(model, ncol = 2, byrow = TRUE)
  Ax <- tail[9, 1] ; Ay <- tail[9, 2]
  Bx <- tail[5, 1] ; By <- tail[5, 2]
  X <- tail[7, 1] ; Y <- tail[7, 2]

  # determine position of the tip of the inner rectrix relative to the tips
  # of the outermost rectrices (positive: anterior to the ORT; negative:
  # posterior to the ORT)
  tip_pos1 <- sign((Bx - Ax) * (Y - Ay) - (By - Ay) * (X - Ax))

  # do the same but for the inner outer rectrices
  Ax <- tail[8, 1] ; Ay <- tail[8, 2]
  Bx <- tail[6, 1] ; By <- tail[6, 2]
  X <- tail[7, 1] ; Y <- tail[7, 2]

  # determine position of the tip of the inner rectrix relative to the tips
  # of the outermost rectrices (positive: anterior to the ORT; negative:
  # posterior to the ORT)
  tip_pos2 <- sign((Bx - Ax) * (Y - Ay) - (By - Ay) * (X - Ax))

  # compute the area of the polygon representing the whole tail
  tail_area <- Momocs::coo_area(tail[c(1, 3, 5, 6, 7, 8, 9, 4, 2, 1), ])


  # if posterior, compute MCS as the distance between the ORTs, and the
  # lifting area as the polygon enclosed by the ORTs and tail base
  if(tip_pos1 == -1) {
    mcs <- stats::dist(rbind(tail[5, ], tail[9, ]))
    lifting_area <- Momocs::coo_area(tail[c(1, 3, 5, 9, 4, 2, 1),])
  }

  # if anterior, find the maximum continuous span as the line perpendicular
  # to the saggital axis that intersects the ORTs and the IRT
  if(tip_pos1 == 1 & tip_pos2 == 1) {

    # center tail around the IRT and set MCS as a vertical line passing
    # through the origin
    tail_cent <- cbind(tail[, 1] - tail[7, 1],
                       tail[, 2] - tail[7, 2])
    mcs_vec <- c(0,10e10)

    ort1 <- rbind(tail_cent[4, ], tail_cent[9, ])
    ort1_vec <- stats::lm(ort1[, 2] ~ ort1[, 1])$coef

    ort2 <- rbind(tail_cent[3, ], tail_cent[5, ])
    ort2_vec <- stats::lm(ort2[, 2] ~ ort2[, 1])$coef

    # find tips of MCS
    A <- matrix(c(mcs_vec[2], -1, ort1_vec[2], -1), byrow = TRUE, nrow = 2)
    b <- c(-mcs_vec[1], -ort1_vec[1])
    p1 <- solve(A, b)

    A <- matrix(c(mcs_vec[2], -1, ort2_vec[2], -1), byrow = TRUE, nrow = 2)
    b <- c(-mcs_vec[1], -ort2_vec[1])
    p2 <- solve(A, b)

    # define "new tail base", enclosing the polygon using the tips of the MCS
    newbase <- rbind(tail_cent[c(2, 4), ],
                     p1,
                     tail_cent[7, ],
                     p2,
                     tail_cent[c(3, 1), ])

    # compute MCS and lifting area
    lifting_area <- Momocs::coo_area(newbase[c(1, 2, 3, 4, 5, 6, 7, 1),])
    mcs <- stats::dist(rbind(p1, p2))

  }

  # if anterior, find the maximum continuous span as the line perpendicular
  # to the saggital axis that intersects the ORTs and the IRT
  if(tip_pos1 == 1 & tip_pos2 == -1) {

    # center tail around the IRT and set MCS as a vertical line passing
    # through the origin
    tail_cent <- cbind(tail[, 1] - mean(tail[6, 1], tail[8, 1]),
                       tail[, 2] - mean(tail[6, 2], tail[8, 2]))
    mcs_vec <- c(0,10e10)

    ort1 <- rbind(tail_cent[4, ], tail_cent[9, ])
    ort1_vec <- stats::lm(ort1[, 2] ~ ort1[, 1])$coef

    ort2 <- rbind(tail_cent[3, ], tail_cent[5, ])
    ort2_vec <- stats::lm(ort2[, 2] ~ ort2[, 1])$coef

    # find tips of MCS
    A <- matrix(c(mcs_vec[2], -1, ort1_vec[2], -1), byrow = TRUE, nrow = 2)
    b <- c(-mcs_vec[1], -ort1_vec[1])
    p1 <- solve(A, b)

    A <- matrix(c(mcs_vec[2], -1, ort2_vec[2], -1), byrow = TRUE, nrow = 2)
    b <- c(-mcs_vec[1], -ort2_vec[1])
    p2 <- solve(A, b)

    # define "new tail base", enclosing the polygon using the tips of the MCS
    newbase <- rbind(tail_cent[c(2, 4), ],
                     p1,
                     tail_cent[7, ],
                     p2,
                     tail_cent[c(3, 1), ])

    # compute MCS and lifting area
    lifting_area <- Momocs::coo_area(newbase[c(1, 2, 3, 4, 5, 6, 7, 1),])
    mcs <- stats::dist(rbind(p1, p2))

  }

  #compute and return lift/drag ratio
  if(!MCS) LD_ratio <- lifting_area / tail_area
  if(MCS)  LD_ratio <- (mcs ^ 2) / tail_area

  return(as.numeric(LD_ratio))

}

