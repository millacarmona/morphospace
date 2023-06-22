################################################################################

#' Pile shapes
#'
#' @description Superimpose all the shapes in the sample.
#'
#' @param shapes Shape data.
#' @param links An optional list with the indices of the coordinates defining
#'   the wireframe (following the format used in \code{Morpho}).
#' @param mshape Logical; whether to plot the mean configuration.
#' @param ... Additional arguments passed to [graphics::plot()] or
#'   [rgl::points3d()].
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
    # coords_l <- lapply(1:nrow(dat$data2d), function(i) {
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
      plot(longlist, type = "n", axes = FALSE, xlab = "", ylab = "", ...)

      if(!is.null(links)) {
        for(i in seq_len(dim(shapes_coord)[3])) {
          for(l in seq_len(length(links))) graphics::lines(shapes_coord[links[[l]],,i], col = "gray")
        }
      }

      graphics::points(longlist, col = "#708095")

      if(mshape == TRUE) {
        graphics::points(expected_shapes(shapes_coord), pch = 16, ...)
        if(!is.null(links)) {
          for(l in seq_len(length(links))) graphics::lines(expected_shapes(shapes_coord)[links[[l]],])
        }
      }

    } else {
      rgl::plot3d(longlist, aspect = FALSE, axes = FALSE, col = "#708095",
                  xlab = "", ylab = "", zlab = "", size = 5, ...)

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

    for(i in seq_len(dim(shapes_coord)[3])) graphics::lines(shapes_coord[,,i], col = "#708095")

    if(mshape == TRUE) {
      graphics::lines(expected_shapes(shapes_coord), lwd = 4)
    }
  }

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

  if(length(col) == 1) col <- rep(col, nlevels(fac))
  if(length(lty) == 1) lty <- rep(lty, nlevels(fac))

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

  if(length(col) == 1) col <- rep(col, nlevels(fac))
  if(length(lty) == 1) lty <- rep(lty, nlevels(fac))

  for(i in seq_len(nlevels(fac))) {
    cent <- colMeans(xy[fac == levels(fac)[i], 1:2])
    vcv <- stats::var(xy[fac == levels(fac)[i], 1:2])
    ell <- car::ellipse(center = cent, shape = vcv,
                        radius = stats::qnorm((1 - conflev) / 2, lower.tail = F),
                        draw = FALSE)
    graphics::polygon(ell, border = col[i],
                      col = grDevices::adjustcolor(col[i], alpha.f = alpha), lty = lty[i], ...)
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
hulls_by_group_3D<-function(xyz, fac, col = seq_len(nlevels(fac)), ...) {

  if(length(col) == 1) col <- rep(col, nlevels(fac))

  for(i in seq_len(nlevels(fac))) {
    matsp <- xyz[fac == levels(fac)[i], 1:3]
    surf <- t(geometry::convhulln(matsp))
    convex <- rgl::triangles3d(matsp[surf, 1], matsp[surf, 2], matsp[surf, 3], col = col[i], ...)
    }
}


################################################################################

#' Build template for 2D shape data
#'
#' @description Create a template (i.e. a set of curves describing the structure
#'   the landmarks are placed upon) to aid morphospace visualization,
#'   interactively.
#'
#' @param image Character; the path to an image of the structure, in png format.
#' @param nlands Numeric; the number of landmarks to be placed.
#' @param ncurves Numeric; the number of curves to be drawn.
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
#' Media, 316.
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

#' Plot scatterpoints into univariate morphospace
#'
#' @description Used internally.
#'
#' @export
plot_univ_scatter <- function(scores, density, col = 1, bg = 1, pch = 1, cex = 1, ...) {
  if(density) {
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
#' @description Used internally.
#'
#' @export
plot_biv_scatter <- function(scores, col = 1, bg = 1, pch = 1, cex = 1, ...) {
  if(any(pch %in% c(21:25))) {
    graphics::points(scores, pch = pch, bg = bg, cex = cex, ...)
  } else {
    graphics::points(scores, pch = pch, col = col, cex = cex, ...)
  }
}


################################################################################

#' Plot univariate density distributions for a series of groups
#'
#' @description Used internally.
#'
#' @export
density_by_group_2D <- function(xy, fac, ax, alpha = 0.2, lwd = 1, lty = 1,
                                col = seq_len(nlevels(fac))) {

  if(length(col) == 1) col <- rep(col, nlevels(fac))
  if(length(lty) == 1) lty <- rep(lty, nlevels(fac))

  dens <- lapply(seq_len(nlevels(fac)), function(i) {
    subdens <- stats::density(xy[fac == levels(fac)[i], ax])
    list(x = subdens$x, y = subdens$y)
  })
  ymax <- max(unlist(lapply(dens, function(x) {x$y})))

  #graphics::abline(h = 0)
  for(i in seq_len(nlevels(fac))) {
    graphics::polygon(dens[[i]]$x, dens[[i]]$y / ymax, lwd = lwd, border = col[i],
                      lty = lty[i], col = grDevices::adjustcolor(col[i], alpha.f = alpha))
  }
}


################################################################################

#' Plot bivariate landscape into morphospace
#'
#' @description Used internally.
#'
#' @export
plot_biv_landscape <- function(landscape, display, type, levels, lwd, lty, col, drawlabels, alpha) {
  if(display == "contour") {
    graphics::contour(landscape$x, landscape$y, landscape$z, levels = levels,
                      lwd = lwd, lty = lty, col = col, labels = round(levels, digits = 3),
                      drawlabels = drawlabels, add = TRUE)
    box()
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

#' Plot bivariate landscape into morphospace
#'
#' @description Used internally.
#'
#' @export
plot_univ_landscape <- function(landscape, drawlabels, col, lwd) {
  if(drawlabels) {
    w.transp <-round(quantile(x = 1:length(landscape$z), probs = c(0.25, 0.5, 0.75)))
    w.transp <- sort(c(w.transp - 1, w.transp, w.transp + 1))

    label_col <- col[w.transp[c(2,5,8)]]
    col[w.transp] <- NA
  }

  for(i in 1:(length(landscape$x) - 1)) {
    lines(rbind(c(landscape$x[i], landscape$z[i]),
                c(landscape$x[i + 1], landscape$z[i + 1])),
          col = col[i], lwd = lwd)
  }
  box()

  if(drawlabels == TRUE) {
    x_text <- colMeans(matrix(landscape$x[w.transp + 1], nrow = 3))
    y_text <- colMeans(matrix(landscape$z[w.transp + 1], nrow = 3))
    labels <- round(colMeans(matrix(landscape$z[w.transp], nrow = 3)), 2)
    text(cbind(x_text, y_text), labels = labels, cex = 0.7, col = label_col)
  }
}

