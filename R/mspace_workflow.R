
###########################################################################

#' Create morphospace
#'
#' @description Create a morphospace as an ordination comprising set of axes
#'   synthesizing shape variation. Allows a variety of multivariate methods for
#'   building the ordination.
#'
#' @param shapes Shapes data.
#' @param axes Numeric of length 2, indicating the axes to be plotted.
#' @param links A list with the indices of the coordinates defining the
#'   wireframe (following the format used in \code{Morpho}).
#' @param template A 2-column matrix containing 1) the actual
#'   landmarks/semilandmarks being analized, followed by 2) the (x,y) cartesian
#'   coordinates defining a curve or set of curves that will be warped using the
#'   deformation interpolated from the changes between landmarks/semilandmarks
#'   (the actual positions of which must be marked with a row of NA, see the
#'   tails dataset).
#' @param p Numeric, indicating the number of landmarks/semilandmarks used (for
#'   landmark data only).
#' @param k Numeric, indicating the number of cartesian dimensions of
#'   landmarks/semilandmarks (for landmark data only).
#' @param FUN The function to be used for synthesizing geometric morphometric
#'   variation. The options include \code{prcomp}, \code{bg_prcomp} and
#'   \code{phy_prcomp}.
#' @param nh Numeric; the number of shape models along the x axis.
#' @param nv Numeric; the number of shape models along the y axis.
#' @param mag Numeric; magnifying factor for shape models.
#' @param rot.models  Numeric; angle (in degrees) to rotate shape models.
#' @param size.models Numeric; size factor for shape models.
#' @param asp.models Numeric; the y/x aspect ratio of shape models.
#' @param col.models The color for wireframes/outlines.
#' @param bg.models Background color for outlines.
#' @param lwd.models Numeric; the width of the lines in wireframes/outlines.
#' @param points Logical; whether to plot the scatter points.
#' @param cex.ldm Numeric; size of landmarks/semilandmarks in the background
#'   models.
#' @param col.ldm The color of landmarks/semilandmarks in the background models.
#' @param plot Logical; whether to plot morphospace.
#' @param xlim,ylim,xlab,ylab,asp Standard arguments passed to the generic plot
#'   function
#' @param ... Further arguments passed to [FUN].
#'
#' @return
#'
#' @seealso \code{\link{proj_shapes}}, \code{\link{proj_consensus}},
#'   \code{\link{proj_groups}}, \code{\link{proj_phylogeny}},
#'   \code{\link{proj_axis}}, \code{\link{plot_consensus}}
#'
#' @export
#'
#' @examples
mspace <- function(shapes,
                   axes = c(1,2),
                   links = NULL,
                   template = NULL,
                   p = NULL,
                   k = NULL,
                   FUN = prcomp,
                   nh = 5,
                   nv = 4,
                   mag = 1,
                   asp = NA,
                   xlim = NULL,
                   ylim = NULL,
                   xlab = NULL,
                   ylab = NULL,
                   rot.models = 0,
                   size.models = 1,
                   asp.models = 1,
                   col.models = "#708095",
                   bg.models = NULL,
                   lwd.models = 1,
                   points = FALSE,
                   cex.ldm = 1,
                   col.ldm = "black",
                   plot = TRUE,
                   ...) {

  dat <- shapes_mat(shapes)
  data2d <- dat$data2d
  datype <- dat$datype

  if(datype == "landm") {
    if(is.null(p) | is.null(k)) {
      if(length(dim(shapes)) == 3){
        p <- nrow(shapes)
        k <- ncol(shapes)
      } else {stop("Provide values for p and k, or provide shapes as an array")}
    }
  } else {
    links <- NULL
    p <- 300
    k <- 2
  }

  ordtype <- match.call()$FUN
  if(is.null(ordtype)) ordtype <- "prcomp"

  FUN <- match.fun(FUN)
  ordination <- FUN(data2d, ...)

  shapemodels <- morphogrid(ordination = ordination, axes = axes, datype = datype, template = template,
                            x = NULL, y = NULL, p = p, k = k, nh = nh, nv = nv, mag = mag,
                            asp = asp, xlim = xlim, ylim = ylim, rot.models = rot.models,
                            size.models = size.models, asp.models = asp.models)
  models_mat <- shapemodels$models_mat
  models_arr <- shapemodels$models_arr

  if(!is.null(xlim)) xlim <- range(c(models_mat[,1]))
  if(!is.null(ylim)) ylim <- range(c(models_mat[,2]))

  if(is.null(xlab)) xlab <- paste0("PC", axes[1])
  if(is.null(ylab)) ylab <- paste0("PC", axes[2])

  if(ordtype == "bg_prcomp") {
    xlab <- paste0("bg", xlab)
    ylab <- paste0("bg", ylab)
  }
  if(ordtype == "phy_prcomp") {
    xlab <- paste0("phy", xlab)
    ylab <- paste0("phy", ylab)
  }

  if(plot == TRUE) {

    plot(models_mat, type = "n", asp = asp,
         xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab)

    for(i in 1:dim(models_arr)[3]) {
      if(datype == "landm") {
        points(models_arr[,,i], pch = 16, cex = cex.ldm * 0.1, col = col.ldm)

        if(!is.null(template)) {
          lines(models_arr[,,i], col = col.models, lwd = lwd.models)
        } else {
          for(l in 1:length(links)) lines(models_arr[,,i][links[[l]],],
                                          col = col.models, lwd = lwd.models)
        }

      } else {
        graphics::polygon(models_arr[,,i], pch = 16,
                          col = bg.models, border = col.models, lwd = lwd.models)
      }
    }

    if(points == TRUE) points(ordination$x)

  }


  plotinfo <- list(p = p, k = k, links = links, template = template, axes = axes, nh = nh, nv = nv, mag = mag,
                   asp = asp, asp.models = asp.models, rot.models = rot.models, size.models = size.models,
                   lwd.models = lwd.models, bg.models = bg.models, col.models = col.models,
                   cex.ldm = cex.ldm, col.ldm = col.ldm)

  results <- list(x = ordination$x, rotation = ordination$rotation, center = ordination$center,
                  datype = datype, ordtype = ordtype, plotinfo = plotinfo)
  class(results) <- "mspace"

  return(invisible(results))

}


###############################################################################################

#' Project shapes into morphospace
#'
#' @description Project a set of shapes into an existing morphospace.
#'
#' @param shapes Shape data.
#' @param mspace An \code{"mspace"} object created using [mspace()].
#' @param pipe Logical; is the function being included in a pipe?
#' @param ... Further arguments passed to [points()].
#'
#' @return
#'
#' @export
#'
#' @examples
proj_shapes <- function(shapes, mspace, pipe = TRUE, ...) {

  dat <- shapes_mat(shapes)
  data2d <- dat$data2d
  datype <- dat$datype

  if(mspace$datype != datype) stop("shapes and mspace types are not compatible")

  scores <- proj_eigen(x = data2d, vectors = mspace$rotation,
                       center = mspace$center)

  if(.Device != "null device") points(scores[, mspace$plotinfo$axes], ...)

  if(pipe == FALSE) return(invisible(scores))
  if(pipe == TRUE) return(invisible(mspace))

}

###############################################################################################

#' Project consensus shape(s) into morphospace
#'
#' @description Project one or more mean shapes into an existing morphospace.
#'
#' @param shapes Shapes data.
#' @param mspace An \code{"mspace"} object created using [mspace()].
#' @param pipe Logical; is the function being included in a pipe?
#' @param ... Further arguments passed to [points()].
#'
#' @return
#'
#' @seealso \code{\link{consensus}}
#'
#' @export
#'
#' @examples
proj_consensus <- function(shapes, mspace, pipe = TRUE, ...) {

  dat <- shapes_mat(shapes)
  data2d <- dat$data2d
  datype <- dat$datype

  if(mspace$datype != datype) stop("shapes and mspace types are not compatible")

  gr_centroids <- proj_eigen(x = data2d, vectors = mspace$rotation,
                             center = mspace$center)

  if(.Device != "null device") points(gr_centroids[, mspace$plotinfo$axes], ...)

  mspace$gr_centroids <- gr_centroids

  if(pipe == FALSE) return(invisible(gr_centroids))
  if(pipe == TRUE) return(invisible(mspace))

}


###############################################################################################

#' Delimit groups in morphospace
#'
#' @description Project convex hulls ennclosing a priori groups in an existing
#'   morphospace.
#'
#' @param mspace An \code{"mspace"} object created using [mspace()].
#' @param shapes Optional shapes data.
#' @param groups Factor; classification of observations into groups.
#' @param pipe Logical; is the function being included in a pipe?
#' @param ... Further arguments passed to [hulls_by_group_2D()].
#'
#' @return
#'
#' @seealso \code{\link{hulls_by_group2D}}
#'
#' @export
#'
#' @examples
proj_groups <- function(mspace, shapes = NULL, groups, pipe = TRUE, ...) {

  if(is.null(shapes)) {
    data2d <- mspace$x
  } else {
    dat <- shapes_mat(shapes)
    datype <- dat$datype
    data2d <- proj_eigen(x = dat$data2d, vectors = mspace$rotation,
                         center = mspace$center)

    if(mspace$datype != datype) stop("shapes and mspace types are not compatible")

  }

  hulls_by_group_2D(data2d[, mspace$plotinfo$axes], fac = groups, ...)

  mspace$gr_class <- groups

  if(pipe == TRUE) return(invisible(mspace))
  #if(pipe == FALSE) return(  volume_of_hulls  )

}


###############################################################################################

#' Project morphometric axis into morphospace
#'
#' @description Project one or more morphometric axes (i.e., linear combinations
#'   of shape variables) into an existing morphospace.
#'
#' @param obj An object containing either a multivariate ordination of class
#'   \code{"prcomp", "bg_prcomp", "phy_prcomp"} or \code{"pls_shape"} or a
#'   \code{"mlm"} object fitted using [lm()].
#' @param mspace An \code{"mspace"} object created using [mspace()].
#' @param axis An optional vector indicating the axis from \code{obj} to be
#'   projected.
#' @param mag Numeric; magnifying factor for representing shape transformation.
#' @param pipe Logical; is the function being included in a pipe?
#' @param ... Further arguments passed to [lines()].
#'
#' @return
#'
#' @export
#'
#' @seealso \code{\link{ax_transformation}}
#'
#' @examples
proj_axis <- function(obj, mspace, axis = 1, mag = 1, pipe = TRUE, ...) {

  ext_shapes2d <- ax_transformation(obj = obj, axis = axis, mag = mag)
  ext_scores <- proj_eigen(x = ext_shapes2d, vectors = mspace$rotation,
                           center = mspace$center)

  if(.Device != "null device") lines(ext_scores[, mspace$plotinfo$axes], ...)

  if(pipe == FALSE) return(invisible(ext_scores))
  if(pipe == TRUE) return(invisible(mspace))

}


###############################################################################################

#' Project phylogenetic structure into morphospace
#'
#' @description Project phylogenetic relationships among a set of shapes
#'   (representing the consensuses of phhylogenetic terminals) into an existing
#'   morphospace.
#'
#' @param tree A \code{"phy"} object containing a phylogenetic tree.
#' @param mspace An \code{"mspace"} object created using [mspace()].
#' @param pipe Logical; is the function being included in a pipe?
#' @param ... Further arguments passed to [lines()].
#'
#' @return
#'
#' @export
#'
#' @examples
proj_phylogeny <- function(tree, mspace, pipe = TRUE, ...) {

  if(is.null(mspace$gr_centroids)) stop("Group centroids have not been provided; add proj_consensus() before")

  nodes_scores <- apply(mspace$gr_centroids[tree$tip.label,], 2, phytools::fastAnc, tree = tree)
  phylo_scores <- rbind(mspace$gr_centroids[tree$tip.label,], nodes_scores)


  if(.Device != "null device") {
    for(i in 1:nrow(tree$edge)) {
      lines(rbind(phylo_scores[tree$edge[i, 1], mspace$plotinfo$axes],
                  phylo_scores[tree$edge[i, 2], mspace$plotinfo$axes]), ...)
    }
  }

  mspace$phylo_scores <- phylo_scores
  mspace$phylo <- tree

  if(pipe == FALSE) {
    phylo_shapes <- rev_eigen(scores = phylo_scores, vectors = mspace$rotation,
                              center = mspace$center)
    return(invisible(phylo_shapes))
  } else {
    return(invisible(mspace))
  }

}



#########################################################################################

#' Plot morphospaces
#'
#' @description Flexible deployment of morphospaces, including their combination
#'   with other variables or a phylogeny.
#'
#'
#' @param mspace An \code{"mspace"} object created using [mspace()].
#' @param axes Numeric of length 1 or 2, indicating the axes to be plotted.
#' @param links A list with the indices of the coordinates defining the
#'   wireframe (following the format used in \code{Morpho}).
#' @param template A 2-column matrix containing 1) the actual
#'   landmarks/semilandmarks being analized, followed by 2) the (x,y) cartesian
#'   coordinates defining a curve or set of curves that will be warped using the
#'   deformation interpolated from the changes between landmarks/semilandmarks
#'   (the actual positions of which must be marked with a row of NA, see the
#'   tails dataset).
#' @param x Optional vector with a non-morphometric variable to be plotted in
#'   the x axis. Alternatively, a \code{"phy"} object can be provided.
#' @param y Optional vector with a non-morphometric variable to be plotted in
#'   the y axis. Alternatively, a \code{"phy"} object can be provided.
#' @param nh Numeric; the number of shape models along the x axis.
#' @param nv Numeric; the number of shape models along the y axis.
#' @param mag Numeric; magnifying factor for shape models.
#' @param points Logical; whether to plot the scatter points corresponding to
#'   the sampled shapes stored in \code{mspace$x}.
#' @param mshapes Logical; whether to plot the scatter points corresponding to
#'   groups' mean shapes stored in \code{mspace$gr_centroids}.
#' @param groups Logical; whether to plot the convex hulls enclosing the groups
#'   stored in \code{mspace$gr_class}
#' @param phylo  Logical; whether to plot the phylogenetic relationships stored
#'   in \code{mspace$phylo}
#' @param size.models Numeric; size factor for shape models.
#' @param asp.models Numeric; the y/x aspect ratio of shape models.
#' @param rot.models  Numeric; angle (in degrees) to rotate shape models.
#' @param col.models The color for wireframes/outlines.
#' @param bg.models Background color for outlines.
#' @param lwd.models Numeric; the width of the lines in wireframes/outlines.
#' @param cex.ldm Numeric; size of landmarks/semilandmarks in the background
#'   models.
#' @param col.ldm The color of landmarks/semilandmarks in the background models.
#' @param pch.points Numeric; the symbol of the scatter points corresponding to
#'   sampled shapes.
#' @param col.points The color of the scatter points corresponding to sampled
#'   shapes.
#' @param cex.points Numeric; the size of the scatter points corresponding to
#'   sampled shapes.
#' @param col.groups The color of the scatter points corresponding to groups
#'   mean shapes.
#' @param pch.groups Numeric; the symbol of the scatter points corresponding to
#'   groups mean shapes.
#' @param cex.groups Numeric; the size of the scatter points corresponding to
#'   groups mean shapes.
#' @param lwd.branches Numeric; the width of the lines depicting phylogenetic
#'   branches.
#' @param xlim,ylim,xlab,ylab,asp Standard arguments passed to the generic plot
#'   function
#'
#' @return
#' @export
#'
#' @examples
plot_mspace <- function(mspace,
                        axes,
                        links = NULL,
                        template = NULL,
                        x = NULL,
                        y = NULL,
                        nh,
                        nv,
                        mag,
                        points = TRUE,
                        mshapes = TRUE,
                        groups = TRUE,
                        phylo = TRUE,
                        asp = NA,
                        xlab,
                        ylab,
                        xlim = NULL,
                        ylim = NULL,
                        size.models,
                        asp.models,
                        rot.models,
                        col.models ,
                        bg.models,
                        lwd.models,
                        cex.ldm,
                        col.ldm,
                        pch.points = 1,
                        col.points = 1,
                        cex.points = 1,
                        col.groups = 1,
                        pch.groups = 16,
                        cex.groups = 1,
                        lwd.branches = 1) {


  supplied <- names(as.list(match.call()))[-1]
  new_args <- lapply(1:length(supplied), function(i) {get(supplied[i])})
  names(new_args) <- supplied
  inh_args <- mspace$plotinfo
  merged_args <- utils::modifyList(inh_args, new_args)
  args <- merged_args


  ordination <- list(x = mspace$x, rotation = mspace$rotation, center = mspace$center)

  if(!is.null(x) & !is.null(y)) stop("Only one of x or y can be specified")



  if(is.null(x) & is.null(y)) { #if neither x nor y have been provided, plot pure morphospace

    shapemodels <- morphogrid(ordination = ordination, axes = args$axes, datype = mspace$datype,
                              template = args$template, x = NULL, y = NULL, p = mspace$plotinfo$p,
                              k = mspace$plotinfo$k, nh = args$nh, nv = args$nv, mag = args$mag,
                              asp = args$asp, xlim = args$xlim, ylim = args$ylim, rot.models = args$rot.models,
                              size.models = args$size.models, asp.models = args$asp.models)

    models_mat <- shapemodels$models_mat
    models_arr <- shapemodels$models_arr

    if(!is.null(xlim)) xlim <- range(c(models_mat[,1]))
    if(!is.null(ylim)) ylim <- range(c(models_mat[,2]))


    if(is.null(args$xlab)) xlab <- paste0("PC", args$axes[1])
    if(is.null(args$ylab)) ylab <- paste0("PC", args$axes[2])

    if(mspace$ordtype == "bg_prcomp") {
      xlab <- paste0("bg", xlab)
      ylab <- paste0("bg", ylab)
    }
    if(mspace$ordtype == "phy_prcomp") {
      xlab <- paste0("phy", xlab)
      ylab <- paste0("phy", ylab)
    }

    plot(models_mat, type = "n", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab)
    for(i in 1:dim(models_arr)[3]) {
      if(mspace$datype == "landm") {
        points(models_arr[,,i], pch = 16, cex = args$cex.ldm * 0.1, col = args$col.ldm)

        if(!is.null(args$template)) {
          lines(models_arr[,,i], col = args$col.models, lwd = args$lwd.models)
        } else {
          for(l in 1:length(args$links)) lines(models_arr[,,i][args$links[[l]],],
                                               col = args$col.models, lwd = args$lwd.models)
        }

      } else {
        graphics::polygon(models_arr[,,i], pch = 16,
                          col = args$bg.models, border = args$col.models, lwd = args$lwd.models)
      }
    }


    #add points, hulls, phylogeny, and/or consensus
    if(points == TRUE) points(mspace$x[, args$axes],
                              pch = pch.points, col = col.points, cex = cex.points)
    if(groups == TRUE) {
      if(is.null(mspace$gr_class)) {
        stop("groups classification has not been added to mspace object")
      } else {
        hulls_by_group_2D(mspace$x[, args$axes], fac = mspace$gr_class,
                          col = col.groups)
      }
    }
    if(phylo == TRUE) {
      if(is.null(mspace$phylo)) {
        stop("phylogenetic relationships have not been added to mspace object")
      } else {
        for(i in 1:nrow(mspace$phylo$edge)) {
          lines(rbind(mspace$phylo_scores[mspace$phylo$edge[i, 1], args$axes],
                      mspace$phylo_scores[mspace$phylo$edge[i, 2], args$axes]),
                lwd = lwd.branches)
        }
      }
    }
    if(mshapes == TRUE) {
      if(is.null(mspace$gr_centroids)) {
        stop("groups centroids have not been added to mspace object")
      } else {
        points(mspace$gr_centroids[, args$axes],
               col = col.groups, pch = pch.groups, cex = cex.groups)
      }
    }



  } else { #if x or y have been provided, show hybrid morphospace

    #if x/y is a phy object, prepare the ground for a phenogram
    if(any(class(x) == "phylo", class(y) == "phylo")) {
      phenogr <- TRUE
    } else {
      phenogr <- FALSE
    }

    if(!is.null(x)) {
      if(class(x) == "phylo") {
        tree <- x
        heights <- phytools::nodeHeights(tree)
        x <- c(rep(max(heights), length(tree$tip.label)),
               unique(heights[,1]))
        args$xlim <- range(x)
      }
    } else {
      if(class(y) == "phylo") {
        tree <- y
        heights <- phytools::nodeHeights(tree)
        y <- c(rep(max(heights), length(tree$tip.label)),
               unique(heights[,1]))
        args$ylim <- range(y)
      }
    }


    shapemodels <- morphogrid(ordination = ordination, x = x, y = y, axes = args$axes, template = args$template,
                              datype = mspace$datype, p = mspace$plotinfo$p, k = mspace$plotinfo$k, nh = args$nh,
                              nv = args$nv, mag = args$mag, asp = args$asp, xlim = args$xlim, ylim = args$ylim,
                              rot.models = args$rot.models, size.models = args$size.models,
                              asp.models = args$asp.models)

    models_mat <- shapemodels$models_mat
    models_arr <- shapemodels$models_arr

    if(!is.null(xlim)) xlim <- range(c(models_mat[,1]))
    if(!is.null(ylim)) ylim <- range(c(models_mat[,2]))

    if(is.null(args$xlab)) {
      if(!is.null(x)) {
        xlab <- "x"
      } else {
        xlab <- paste0("PC", args$axes[1])
        if(mspace$ordtype == "bg_prcomp") {
          xlab <- paste0("bg", xlab)
        }
        if(mspace$ordtype == "phy_prcomp") {
          xlab <- paste0("phy", xlab)
        }
      }
    }
    if(is.null(args$ylab)) {
      if(!is.null(y)) {
        ylab <- "y"
      } else {
        ylab <- paste0("PC", args$axes[1])
        if(mspace$ordtype == "bg_prcomp") {
          ylab <- paste0("bg", ylab)
        }
        if(mspace$ordtype == "phy_prcomp") {
          ylab <- paste0("phy", ylab)
        }
      }
    }


    plot(models_mat, type = "n", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab)

    for(i in 1:dim(models_arr)[3]) {
      if(mspace$datype == "landm") {
        points(models_arr[,,i], pch = 16, cex = args$cex.ldm * 0.1, col = args$col.ldm)

        if(!is.null(args$template)) {
          lines(models_arr[,,i], col = args$col.models, lwd = args$lwd.models)
        } else {
          for(l in 1:length(args$links)) lines(models_arr[,,i][args$links[[l]],],
                                          col = args$col.models, lwd = args$lwd.models)
        }

      } else {
        graphics::polygon(models_arr[,,i], pch = 16,
                          col = args$bg.models, border = args$col.models, lwd = args$lwd.models)
      }
    }


    if(phenogr == TRUE) { #if x/y is a phy object, plot a phenogram
      for(i in 1:nrow(tree$edge)) {
        phyloxy <- cbind(x, mspace$phylo_scores[,args$axes[1]], y)
        lines(rbind(phyloxy[tree$edge[i, 1],],
                    phyloxy[tree$edge[i, 2],]), lwd = lwd.branches)
      }
      if(points == TRUE) {
        points(phyloxy[-c(1:length(tree$tip.label)),], pch = 16)
        points(phyloxy[c(1:length(tree$tip.label)),][rownames(mspace$gr_centroids),],
               bg = col.groups, pch = 21, cex = cex.groups)

        }
    } else { #else, go for a generic hybrid morphospace

      xy <- cbind(x, mspace$x[,args$axes[1]], y)

      #add points, hulls, phylogeny, and/or consensus
      if(points == TRUE) points(xy, pch = pch.points,
                                col = col.points, cex = cex.points)
      if(groups == TRUE) {
        if(is.null(mspace$gr_class)) {
          stop("groups classification has not been added to mspace object")
        } else {
          hulls_by_group_2D(xy, fac = mspace$gr_class, col = col.groups)
        }
      }

      if(phylo == TRUE) {
        if(is.null(mspace$phylo)) {
          stop("phylogenetic relationships have not been added to mspace object")
        } else {
          if(is.null(mspace$gr_class)) {
            stop("groups classification has not been added to mspace object")
          } else {
            if(!is.null(x)) {
              meanx <- tapply(x, mspace$gr_class, mean)
            } else {
              meanx <- NULL
            }
            if(!is.null(y)) {
              meany <- tapply(y, mspace$gr_class, mean)
            } else {
              meany <- NULL
            }

            meanxy <- cbind(meanx, mspace$gr_centroids[,args$axes[1]], meany)
            nodesxy <- apply(meanxy[tree$tip.label,], 2, phytools::fastAnc, tree = tree)
            phyloxy <- rbind(meanxy, nodesxy)

          }

          for(i in 1:nrow(mspace$phylo$edge)) {
            lines(rbind(phyloxy[mspace$phylo$edge[i, 1],],
                        phyloxy[mspace$phylo$edge[i, 2],]),
                  lwd = lwd.branches)
          }
        }
      }
      if(mshapes == TRUE) {
        if(is.null(mspace$gr_class)) {
          stop("groups centroids have not been added to mspace object")
        } else {

          if(!is.null(x)) {
            meanx <- tapply(x, mspace$gr_class, mean)
          } else {
            meanx <- NULL
          }
          if(!is.null(y)) {
            meany <- tapply(y, mspace$gr_class, mean)
          } else {
            meany <- NULL
          }

          meanxy <- cbind(meanx, mspace$gr_centroids[,args$axes[1]], meany)

          points(meanxy, col = col.groups,
                 pch = pch.groups, cex = cex.groups)
        }
      }
    }
  }

}

