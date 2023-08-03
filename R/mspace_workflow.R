################################################################################

#' Generate morphospace
#'
#' @description Build morphospaces using a variety of multivariate methods,
#'    and depict shape variation represented by the resulting ordination axes.
#'
#' @param shapes Shape data.
#' @param axes Numeric of length 1 (univariate morphospace) or 2 (bivariate
#'    morphospace), indicating the axes to be plotted.
#' @param links A list with the indices of the coordinates defining the
#'    wireframe (following the format used in \code{Morpho}).
#' @param template Either a 2-column matrix with landmarks/semilandmarks and
#'    template curves coordinates (for 2D shape data) or a \code{"mesh3d"}
#'    object \strong{representing the mean shape of the sample} (for 3D shape
#'    data). See details below.
#' @param FUN The function (method) to be used for ordination of shape
#'    variation. Supported alternatives include \code{\link[stats]{prcomp}},
#'    \code{\link{phy_prcomp}}, \code{\link{bg_prcomp}} and
#'    \code{\link{pls_shapes}}.
#' @param nh Numeric; number of shape models along the x axis.
#' @param nv Numeric; number of shape models along the y axis.
#' @param mag Numeric; magnifying factor for shape models.
#' @param invax Optional numeric indicating which of the axes provided in
#'    \code{axes} needs to be inverted (options are \code{1}, \code{2} or
#'    \code{c(1,2)}).
#' @param rescale Logical; whether to re-scale background shape models so shape
#'    variation is shown more clearly.
#' @param adj_frame Numeric of length 2, providing \emph{a posteriori} scaling
#'    factors for the width and height of the frame, respectively.
#' @param rot.models  Numeric; angle (in degrees) to rotate shape models.
#' @param size.models Numeric; size factor for shape models.
#' @param asp.models Numeric; the y/x aspect ratio of shape models.
#' @param col.models Color for wireframes/outlines.
#' @param bg.models Background color for outlines/meshes.
#' @param lwd.models Numeric; width of the lines in wireframes/outlines.
#' @param alpha.models Numeric; transparency factor for background models
#'    (3D only).
#' @param points Logical; whether to plot the scatter points.
#' @param cex.ldm Numeric; size of landmarks/semilandmarks in the background
#'    models.
#' @param col.ldm Color of landmarks/semilandmarks in the background models.
#' @param plot Logical; whether to plot morphospace.
#' @param models Logical; whether to plot background shape models.
#' @param xlim,ylim,xlab,ylab,asp Standard arguments passed to the generic plot
#'    function.
#' @param ... Further arguments passed to \code{FUN}.
#'
#' @details This function is the heart of the \code{morphospace} workflow. It
#'    computes a new ordination space from a sample of normalized shapes using
#'    eigenanalysis-based multivariate methods (supported alternatives include
#'    PCA, between groups PCA, phylogenetic PCA, two-block PLS); it will also
#'    generate a series of shape models depicting the range of variation. The
#'    resulting \code{"mspace"} object stores all the information necessary to
#'    project and/or retrieve new shapes into the ordination space. The output
#'    of \code{mspace} can be expanded using the \code{proj_*} family of
#'    functions and the \code{%>%} operator from \code{magrittr}.
#'
#'    For landmark data, representation of shape variation can be further
#'    aided by providing a template, which contains additional geometric
#'    features from the structure the landmark/semilandmarks were placed upon.
#'    Templates will be warped using TPS interpolation to produce the set of
#'    background shape models. For 2D landmark data, templates must be provided
#'    as a 2-column matrix containing the actual landmarks/semilandmarks,
#'    followed by the coordinates defining a curve or set of curves, separated
#'    from the former and from each other by a row of\code{NA}s
#'    (see \code{\link{build_template2d}}). For 3D landmark data, the template
#'    must be a \code{"mesh3d"} object corresponding to the \strong{actual mean
#'    shape of the sample} (which can be obtained using
#'    \code{\link{expected_shapes}} + [Morpho::tps3d()]; see examples below).
#'
#'
#' @return An object of class \code{"mspace"}, which is a list containing:
#'   \itemize{
#'   \item{$ordination:} { a list with the output from the ordination method
#'      used, styled in the \code{\link{prcomp}} format (typically containing at
#'      least an \code{$x}, \code{$rotation}, and \code{$center} slots. Also
#'      contains the type of data (\code{$datype}) and ordination method
#'      (\code{$ordtype}) used.}
#'   \item{$projected:} { a list containing the elements that have been
#'      projected into the morphospace stored in the \code{"mspace"} object.
#'      Initially includes only the background shape models
#'      (\code{$shapemodels}), but more elements can be added using the
#'      \code{proj_*} family of functions.
#'   \item{$plotinfo:} { a list containing the graphical parameters used to
#'      create the plot. Passed to \code{\link{plot_mspace()}} to regenerate
#'      morphospaces.}
#'   }
#'
#' @seealso \code{\link{{print.mspace}}, \code{\link{proj_shapes}},
#'   \code{\link{proj_consensus}}, \code{\link{proj_groups}},
#'   \code{\link{proj_phylogeny}}, \code{\link{proj_axis},
#'   \code{\link{proj_landscape}}, \code{\link{extract_shapes},
#'   \code{\link{prcomp}}, \code{\link{bg_prcomp}}, \code{\link{phy_prcomp}},
#'   \code{\link{pls_shapes}},
#'
#' @export
#'
#' @examples
#' ##2D Landmark data
#'
#' #load and extract relevant data and information
#' data("tails")
#' shapes <- tails$shapes
#' species <- tails$data$species
#' links <- tails$links
#' tree <- tails$tree
#'
#' #generate morphospace using the basic sample of shapes, PCA as ordination
#' #method and the links between landmarks provided for backround models
#' mspace(shapes, links = links, mag = 0.7, axes = c(1,2), points = TRUE)
#'
#' #increase magnification factor x2:
#' mspace(shapes, links = links, mag = 1.5, axes = c(1,2), points = TRUE)
#'
#' #plot PCs 1 and 3
#' mspace(shapes, links = links, mag = 1.5, axes = c(1,3), points = TRUE)
#'
#' #generate morphospace using the basic sample of shapes, bgPCA as ordination
#' #method and links between landmarks for backround models
#' mspace(shapes, links = links, FUN = bg_prcomp, groups = species, mag = 0.7,
#'        axes = c(1,2), invax = 1, points = TRUE)
#'
#' #generate morphospace using species' consensus shapes, phylogenetic PCA as
#' #ordination method, and links between landmarks for background models
#' sp_shapes <- expected_shapes(shapes, species)
#' mspace(sp_shapes, links = links, FUN = phy_prcomp, tree = tree, mag = 0.7,
#'        axes = c(1,2), points = TRUE)
#'
#' #just create a morphospace without plotting, save into an object, and inspect
#' morphosp <- mspace(shapes, links = links, mag = 0.7, axes = c(1,2),
#'                    plot = FALSE)
#' morphosp #general info about the object
#' names(morphosp) #slots
#'
#'
#' #load wing data for a quick demo with templates
#' data("wings")
#' shapes <- wings$shapes
#' links <- wings$links
#' template <- wings$template
#'
#' #generate morphospace using links
#' mspace(shapes, links = links, mag = 3, axes = c(1,2), points = TRUE)
#'
#' #generate morphospace using template
#' mspace(shapes, template = template, mag = 3, axes = c(1,2), points = TRUE)
#'
#'
#'
#' ##3D Landmark data
#'
#' \dontrun{
#' #load data and packages and extract relevant data and information
#' library(Morpho)
#' library(geomorph)
#' data("shells3D")
#' shapes <- shells3D$shapes
#' mesh_meanspec <- shells3D$mesh_meanspec
#'
#' #generate morphospace. This is interactive, you need to rotate the shape by
#' #yourself and then press enter into the console.
#' mspace(shapes, mag = 1, axes = c(1,2), col.ldm = "black", cex.ldm = 2,
#'        points = TRUE)
#'
#' #generate morphospace using a mesh template that improves visualization:
#' #first, get shape corresponding to shells3D$mesh_meanspec using
#' #geomorph::findMeanSpec
#' meanspec_id <- findMeanSpec(shapes)
#' meanspec_shape <- shapes[,,meanspec_id]
#'
#' #then get the consensus shape and warp the sample mesh to get the mesh
#' #corresponding to the consensus using Morpho::tps3d
#' meanshape <- expected_shapes(shapes)
#' meanmesh <- tps3d(x = mesh_meanspec , refmat = meanspec_shape,
#'                   tarmat = meanshape)
#'
#' #finally, generate morphospace providing template (this function used the
#' #mesh warped to the mean shape of the entire sample, hence the previous
#' #lines)
#' morphosp_3d <- mspace(shapes, mag = 1, axes = c(1,2), template = meanmesh,
#'                       bg.models = "gray", nh = 4, nv = 4, cex.ldm = 0,
#'                       points = TRUE)
#'
#' #inspect the contents of the object
#' morphosp_3d
#' }
#'
#'
#' ##Outline data
#'
#' #load and extract relevant data and information
#' data("shells")
#' shapes <- shells$shapes$coe
#'
#' #generate morphospace using all the raw variation
#' mspace(shapes, mag = 1, axes = c(1,2), nh = 5, nv = 4, size.models = 1,
#'        rescale = F, asp.models = 1, bg.model = "light gray")
#'
#' #shapes in the background can be rescaled to avhieve a (slightly) better
#' #visualization. Also, save the ordination into an object.
#' morphosp_F <- mspace(shapes, mag = 1, axes = c(1,2), nh = 5, nv = 4,
#'                      size.models = 1, asp.models = 1,
#'                      bg.model = "light gray")
#'
#' #inspect the contents of the object
#' morphosp_F
mspace <- function(shapes,
                   axes = c(1,2),
                   links = NULL,
                   template = NULL,
                   FUN = stats::prcomp,
                   nh = 5,
                   nv = 4,
                   mag = 1,
                   invax = NULL,
                   asp = NA,
                   rescale = TRUE,
                   xlim = NULL,
                   ylim = NULL,
                   xlab = NULL,
                   ylab = NULL,
                   adj_frame = c(1, 1),
                   rot.models = 0,
                   size.models = 1,
                   asp.models = 1,
                   col.models = "#708095",
                   bg.models = NULL,
                   lwd.models = 1,
                   alpha.models = 1,
                   points = FALSE,
                   cex.ldm = 1,
                   col.ldm = "black",
                   plot = TRUE,
                   models = TRUE,
                   ...) {

  dat <- shapes_mat(shapes)
  data2d <- dat$data2d
  datype <- dat$datype

  if(datype == "landm") {
    p <- nrow(shapes)
    k <- ncol(shapes)
  } else {
    links <- NULL
    p <- 300
    k <- 2
  }

  FUN <- match.fun(FUN)
  ordination <- FUN(data2d, ...)
  ordination$ordtype <- class(ordination)
  ordination$datype <- datype

  if(!is.null(invax)) {
    ordination$x[,axes[invax]] <- ordination$x[,axes[invax]] * -1
    ordination$rotation[,axes[invax]] <- ordination$rotation[,axes[invax]] * -1
  }

  if(ncol(ordination$x) == 1 | length(axes) == 1) {
    y <- rep(1, nrow(ordination$x))
    ylim <- c(0, 1)
    ylab <- "relative density"
  } else {
    y <- NULL
  }

  if(ncol(ordination$x) == 1) axes <- 1

  rotation_matrix <- NULL

  if(k == 3 & datype == "landm") {

    shapemodels <- morphogrid(ordination = ordination, axes = axes, datype = datype, template = NULL,
                              x = NULL, y = y, p = p, k = k, nh = nh, nv = nv, mag = mag,
                              asp = asp, xlim = xlim, ylim = ylim, rot.models = rot.models,
                              size.models = size.models, asp.models = 1, rescale = rescale)

    refshape <- expected_shapes(shapes)
    xlim <- range(ordination$x[,axes[1]])
    if(ncol(ordination$x) > 1 | length(axes) > 1) ylim <- range(ordination$x[,axes[2]])

    plot_morphogrid3d(x = NULL, y = y, morphogrid = shapemodels, refshape = refshape,
                      template = template, links = links, ordtype = ordination$ordtype,
                      adj_frame = adj_frame, axes = axes, xlim = xlim, ylim = ylim, xlab = xlab,
                      ylab = ylab, cex.ldm = cex.ldm, col.ldm = col.ldm, col.models = col.models,
                      lwd.models = lwd.models, bg.models = bg.models, size.models = size.models,
                      asp.models = asp.models, alpha.models = alpha.models, plot = plot, models = models)

    rotation_matrix <- rgl::par3d()$userMatrix

  } else {

    shapemodels <- morphogrid(ordination = ordination, axes = axes, datype = datype, template = template,
                              x = NULL, y = y, p = p, k = k, nh = nh, nv = nv, mag = mag,
                              asp = asp, xlim = xlim, ylim = ylim, rot.models = rot.models,
                              size.models = size.models, asp.models = asp.models, rescale = rescale)

    plot_morphogrid2d(x = NULL, y = y, morphogrid = shapemodels, template = template,
                      links = links, datype = datype, ordtype = ordination$ordtype, axes = axes,
                      adj_frame = adj_frame, p = p, xlab = xlab, ylab = ylab, cex.ldm = cex.ldm,
                      col.ldm = col.ldm, col.models = col.models, lwd.models = lwd.models,
                      bg.models = bg.models, plot = plot, models = models)

  }

  if(points == TRUE) graphics::points(ordination$x[,axes])

  if(length(axes) == 1) {
    if(datype == "landm") shapemodels <- shapemodels$shapemodels[,,1:nh]
    if(datype == "fcoef") shapemodels <- shapemodels$shapemodels[1:nh,]
  } else {
    shapemodels <- shapemodels$shapemodels
  }


  plotinfo <- list(p = p, k = k, links = links, template = template, axes = axes, nh = nh, nv = nv, mag = mag,
                   xlim = xlim, ylim = ylim, asp = asp, rescale = rescale, adj_frame = adj_frame,
                   asp.models = asp.models, rot.models = rot.models, size.models = size.models,
                   lwd.models = lwd.models, bg.models = bg.models, col.models = col.models,
                   alpha.models = alpha.models, cex.ldm = cex.ldm, col.ldm = col.ldm,
                   models = models, rotation_matrix = rotation_matrix)

  results <- list(ordination = ordination,
                  projected = list(shapemodels = shapemodels),
                  plotinfo = plotinfo)
  class(results) <- "mspace"

  return(invisible(results))

}


################################################################################

#' Project shapes into morphospace
#'
#' @description Project a set of shapes as scatterpoints into an existing
#'    morphospace.
#'
#' @param mspace An \code{"mspace"} object.
#' @param shapes Shape data.
#' @param density Logical; whether to add density distribution for points
#'   (univariate ordinations only).
#' @param pipe Logical; is the function being included in a pipe?
#' @param ... Further arguments passed to [graphics::points()].
#'
#' @details The purpose of this function is to maintain morphospace building
#'   and sample representation as separated steps, as well as to add flexibility
#'   to graphical representation of scatterpoints.
#'
#' @return If a plot device with a morphospace is open, shapes fed to
#'   \code{shapes} are projected into morphospace. If \code{pipe = FALSE}
#'   those scores are returned invisibly. If \code{pipe = TRUE} the supplied
#'   \code{"mspace"} object will be modified by appending a \code{$scores} slot
#'   to \code{$projected} and adding some graphical parameters (stored into the
#'   \code{$plotinfo} slot), and returned invisibly.
#'
#' @importFrom magrittr %>%
#'
#' @export
#'
#' @examples
#' #load and extract relevant data, packages and information
#' library(magrittr)
#' data("wings")
#' shapes <- wings$shapes
#' cactus <- wings$data$sex
#' species <- wings$data$species
#' template <- wings$template
#'
#' #generate basic morphospace, add sampled shapes
#' mspace(shapes, template = template, mag = 0.7, axes = c(1,2)) %>%
#'   proj_shapes(shapes = shapes)
#'
#' #change colors, symbols, etc for the scatter
#' mspace(shapes, template = template, mag = 0.7, axes = c(1,2)) %>%
#'   proj_shapes(shapes = shapes[,,species == "Db"],  col = c("green"),
#'               pch = c(1,16)[cactus[species == "Db"]]) %>%
#'   proj_shapes(shapes = shapes[,,species == "Dk"],  col = c("blue"),
#'               pch = c(1,16)[cactus[species == "Dk"]])
proj_shapes <- function(mspace, shapes, density = TRUE, pipe = TRUE, ...) {

  args <- c(as.list(environment()), list(...))

  dat <- shapes_mat(shapes)
  data2d <- dat$data2d
  datype <- dat$datype

  if(mspace$ordination$datype != datype) stop("shapes and mspace types are not compatible")

  scores <- proj_eigen(x = data2d, vectors = mspace$ordination$rotation,
                       center = mspace$ordination$center)

  if(.Device != "null device") {

    if(ncol(mspace$ordination$x) > 1) {
      if(length(mspace$plotinfo$axes) > 1) {
        plot_biv_scatter(scores = scores[, mspace$plotinfo$axes], ...)
      } else {
        plot_univ_scatter(scores = cbind(scores[, mspace$plotinfo$axes[1]]),
                          density = density, ...)
      }
    } else {
      plot_univ_scatter(scores = scores, density = density, ...)
    }
  }

  mspace$projected$scores <- rbind(cbind(mspace$projected$scores), cbind(scores))

  if(is.null(args$pch)) args$pch <- 1
  if(length(args$pch) == 1) args$pch <- rep(args$pch, nrow(scores))
  mspace$plotinfo$pch.points <- c(mspace$plotinfo$pch.points, args$pch)

  if(is.null(args$col)) args$col <- 1
  if(length(args$col) == 1) args$col <- rep(args$col, nrow(scores))
  if(is.factor(args$col)) args$col <- as.numeric(args$col)
  mspace$plotinfo$col.points <- c(mspace$plotinfo$col.points, col2hex(args$col))

  if(is.null(args$bg)) args$bg <- 1
  if(length(args$bg) == 1) args$bg <- rep(args$bg, nrow(scores))
  if(is.factor(args$bg)) args$bg <- as.numeric(args$bg)
  mspace$plotinfo$bg.points <- c(mspace$plotinfo$bg.points, col2hex(args$bg))

  if(is.null(args$cex)) args$cex <- 1
  if(length(args$cex) == 1) args$cex <- rep(args$cex, nrow(scores))
  mspace$plotinfo$cex.points <- c(mspace$plotinfo$cex.points, args$cex)

  mspace$plotinfo$density.points <- args$density

  if(pipe == FALSE) return(invisible(scores))
  if(pipe == TRUE) return(invisible(mspace))

}


################################################################################

#' Delimit groups in morphospace
#'
#' @description Project convex hulls or confidence ellipses enclosing
#'   \emph{a priori} groups into an existing morphospace.
#'
#' @param mspace An \code{"mspace"} object.
#' @param shapes Optional shape data (if \code{NULL}, the function will look for
#'    them first in \code{mspace$projected$scores} and then in
#'    \code{mspace$ordination$x}).
#' @param groups Factor; classification of observations into groups. Its length
#'    must be the same as the number of shapes provided in \code{shapes}.
#' @param ellipse Logical; whether to plot confidence ellipses (if \code{FALSE},
#'   convex hulls will be used instead).
#' @param conflev Numeric; confidence level used for confidence ellipse(s).
#' @param density Logical; whether to add density distribution for groups
#'   (univariate ordinations only).
#' @param pipe Logical; is the function being included in a pipe?
#' @param ... Further arguments passed to [hulls_by_group_2D()] or
#'   [ellipses_by_group_2D()].
#'
#' @details The goal of this function is to add a classification for shapes
#'   populating the morphospace to \code{"mspace"} objects, as well as to
#'   facilitate group visualization. Other than that, it is just a wrapper for
#'   \code{hulls_by_group_2D} and \code{ellipses_by_group_2D}.
#'
#' @return If a plot device with a morphospace is open, convex hulls or
#'   confidence ellipses enclosing the scores corresponding to \code{groups}
#'   are projected into morphospace. If \code{pipe = TRUE} the supplied
#'   \code{"mspace"} object will be modified by appending a \code{$gr_class}
#'   slot to \code{$projected}, as well as by adding some graphical parameters
#'   (stored into the \code{$plotinfo} slot), and returned invisibly.
#'
#' @seealso \code{\link{hulls_by_group_2D}}, \code{\link{ellipses_by_group_2D}}
#'
#' @export
#'
#' @examples
#' #load and extract relevant data, packages and information
#' library(magrittr)
#' data("shells")
#' shapes <- shells$shapes
#' species <- shells$data$species
#' sp_shapes <- expected_shapes(shapes, species)
#'
#' #generate basic morphospace, add sampled shapes and convex hulls for species
#' mspace(shapes, links = links, mag = 0.7, axes = c(1,2),
#'        bg.model = "light gray") %>%
#'   proj_shapes(shapes = shapes, col = c(1:4)[species], pch = 1) %>%
#'   proj_groups(groups = species, col = 1:4, lwd = 2, alpha = 0.5)
#'
#' #generate basic morphospace, add sampled shapes and 95% confidence for
#' #species
#' msp <- mspace(shapes, links = links, mag = 0.7, axes = c(1,2),
#'               bg.model = "light gray") %>%
#'   proj_shapes(shapes = shapes, col = c(1:4)[species], pch = 1) %>%
#'   proj_groups(groups = species, col = 1:4, lwd = 2, ellipse = TRUE,
#'               conflev = 0.95, alpha = 0.1)
#'
#' #add legend using plot_mspace()
#' plot_mspace(msp, legend = TRUE, cex.legend = 2, pch.groups = 16)
proj_groups <- function(mspace, shapes = NULL, groups, ellipse = FALSE,
                        conflev = 0.95, density = TRUE, pipe = TRUE, ...) {

  args <- c(as.list(environment()), list(...))

  if(is.null(shapes)) {
    if(!is.null(mspace$projected$scores)) {
      data2d <- mspace$projected$scores
    } else {
      data2d <- mspace$ordination$x
    }
  } else {
    dat <- shapes_mat(shapes)
    datype <- dat$datype
    data2d <- proj_eigen(x = dat$data2d,
                         vectors = mspace$ordination$rotation,
                         center = mspace$ordination$center)

    if(mspace$ordination$datype != datype) stop("shapes and mspace types are not compatible")

  }

  if(.Device != "null device") {

    if(ncol(mspace$ordination$x) > 1) {
      if(length(mspace$plotinfo$axes) > 1) {

        if(ellipse == FALSE) {
          hulls_by_group_2D(data2d[, mspace$plotinfo$axes], fac = groups, ...)
        } else {
          ellipses_by_group_2D(data2d[, mspace$plotinfo$axes], fac = groups,
                                conflev = conflev, ...)
        }
      } else {
        if(density == TRUE) {
          density_by_group_2D(data2d, groups, ax = mspace$plotinfo$axes[1], ...)
        }
      }

    } else {
      if(density == TRUE) {
        density_by_group_2D(data2d, groups, ax = 1, ...)
      }
    }
  }


  mspace$projected$gr_class <- if(is.null(mspace$projected$gr_class)) groups else {
    if(all(levels(groups) %in% levels(mspace$projected$gr_class))) {
      c(factor(paste0(mspace$projected$gr_class,"_bis."),
               levels = paste0(levels(mspace$projected$gr_class),"_bis.")), groups)
    } else c(mspace$projected$gr_class, groups)
  }
  mspace$projected$gr_scores <- rbind(cbind(mspace$projected$gr_scores), cbind(data2d))


  if(is.null(args$col)) args$col <- sort(unique(as.numeric(groups)))
  args$col <- if(length(args$col) == 1) rep(args$col, nlevels(groups)) else {
    if(nlevels(groups) > length(args$col)) {
      c(rep(NA, nlevels(groups) - length(unique(args$col))),
        unique(args$col))[order(c((1:nlevels(groups))[-which(levels(groups) %in% unique(groups))],
                                  which(levels(groups) %in% unique(groups))))]
    } else args$col
  }
  if(is.factor(args$col)) args$col <- as.numeric(args$col)

  if(is.null(args$lty)) args$lty <- 1
  args$lty <- if(length(args$lty) == 1) rep(args$lty, nlevels(groups)) else
    if(nlevels(groups) > length(args$lty)) {
      c(rep(0, nlevels(groups) - length(args$lty)),
        args$lty)[order(c((1:nlevels(groups))[-which(levels(groups) %in% unique(groups))],
                          which(levels(groups) %in% unique(groups))))]
    } else args$lty

  args$alpha <- if(is.null(args$alpha))
    if(length(mspace$plotinfo$axes) == 1 | ncol(mspace$ordination$x) == 1) 0.2 else 0
  else args$alpha


  mspace$plotinfo$col.groups <- c(mspace$plotinfo$col.groups, col2hex(args$col))
  mspace$plotinfo$bg.groups <- c(mspace$plotinfo$bg.groups, col2hex(args$col))
  mspace$plotinfo$lty.groups <- c(mspace$plotinfo$lty.groups, args$lty)
  mspace$plotinfo$density.groups <- args$density
  mspace$plotinfo$ellipse.groups <- args$ellipse
  mspace$plotinfo$conflev.groups <- args$conflev
  mspace$plotinfo$lwd.groups <- args$lwd
  mspace$plotinfo$alpha.groups <- args$alpha

  if(pipe == TRUE) return(invisible(mspace))

}


################################################################################

#' Project a morphometric axis into morphospace
#'
#' @description Project one or more morphometric axes (i.e., linear combinations
#'   of shape variables) as vectors into an existing bivariate morphospace.
#'
#' @param mspace An \code{"mspace"} object.
#' @param obj An object containing either a multivariate ordination of class
#'   \code{"prcomp", "bg_prcomp", "phy_prcomp"} or \code{"pls_shape"}, or a
#'   \code{"mlm"} object fitted using [stats::lm()].
#' @param axis Optional; which axis from \code{obj} should to be projected?
#' @param mag Numeric; magnifying factor for representing shape transformation.
#' @param type Integer; type of arrows ([0] = no arrow; [1] = pointing towards
#'   the maximum; [2] = pointing towards the maximum, [3] = pointing in both
#'   directions).
#' @param pipe Logical; is the function being included in a pipe?
#' @param ... Further arguments passed to [graphics::arrows()].
#'
#' @details This function is primarily aimed at graphically representing
#'   morphometric axes (estimated using either linear models or multivariate
#'   ordination methods) into an existing morphospace for heuristic exploration
#'   of patterns in the data. It can also be used to extract theoretical shapes
#'   at the extremes of those axes, although [ax_transformation()] does the
#'   same thing in a more flexible and straightforward way.
#'
#'   Axes computed by fitting linear models to shape data can differ in
#'   extension from axes obtained through supervised ordination using the same
#'   supervising variable (i.e., shape transformations will be either stretched
#'   or truncated) due to the former assuming that the explanatory variable has
#'   been measured without error. Also, the former will not be necessarily
#'   centered.
#'
#'   For statistical analysis of axes (e.g., trajectory analysis) their vector
#'   coefficients can be extracted directly from slope coefficients stored in
#'   \code{"mlm"} objects or eigenvector coefficients stored in the
#'   \code{$rotation} slot returned by multivariate ordination methods.
#'
#' @return If a plot device with a morphospace is open, a straight line marking
#'   the scores representing shapes at the extremes of the morphometric axis is
#'   projected into morphospace. If \code{pipe = FALSE} those scores are
#'   returned invisibly. If \code{pipe = TRUE} the supplied \code{"mspace"}
#'   object will be modified by appending a \code{$shape_axis} slot to
#'   \code{$projected}, as well as by adding some graphical parameters (stored
#'   into the \code{$plotinfo} slot), and returned invisibly.
#'
#' @export
#'
#' @seealso \code{\link{ax_transformation}}
#'
#' @examples
#' #load and extract relevant data, packages and information
#' library(magrittr)
#' library(geomorph)
#' data("tails")
#' shapes <- tails$shapes
#' logsizes <- log(tails$sizes)
#' species <- tails$data$species
#' sp_shapes <- expected_shapes(shapes, species)
#' tree <- tails$tree
#' links <- tails$links
#'
#' ##Compare (orientations) axes resulting from different version of PCA
#'
#' #first perform the different variants of PCA on tail shape data
#' pca <- prcomp(two.d.array(shapes))
#' bgpca <- bg_prcomp(two.d.array(shapes), groups = species)
#' phypca <- phy_prcomp(two.d.array(sp_shapes), tree = tree)
#'
#' #then project the first 2 axes from each into morphospace
#' mspace(shapes, links = links, mag = 0.7, axes = c(1,2)) %>%
#'   proj_axis(obj = pca, axis = 1, col = "red", lwd = 2) %>%
#'   proj_axis(obj = pca, axis = 2, col = "red", lwd = 2) %>%
#'   proj_axis(obj = bgpca, axis = 1, col = "blue", lwd = 2) %>%
#'   proj_axis(obj = bgpca, axis = 2, col = "blue", lwd = 2) %>%
#'   proj_axis(obj = phypca, axis = 1, col = "black", lwd = 2) %>%
#'   proj_axis(obj = phypca, axis = 2, col = "black", lwd = 2)
#'
#'
#' ##Linear models vs ordination methods
#'
#' #compute intraspecific allometric axis detrend_shapes, using lm and
#' #pls_shapes
#' detr_shapes <- arrayspecs(
#'   detrend_shapes(lm(two.d.array(shapes) ~ species)),
#'                           p = 9, k = 2)
#' intrasp_allo_mod <- lm(two.d.array(detr_shapes) ~ logsizes)
#' intrasp_allo_pls <- pls_shapes(shapes = two.d.array(detr_shapes),
#'                                X = logsizes)
#'
#' #compute intraspecific allometric axis using tapply, lm and pls_shapes
#' sp_logsizes <- tapply(logsizes, species, max)
#' intersp_allo_mod <- lm(two.d.array(sp_shapes) ~ sp_logsizes)
#' intersp_allo_pls <- pls_shapes(shapes = two.d.array(sp_shapes),
#'                                X = sp_logsizes)
#'
#' #generate basic morphospace, add intraspecific (red) and interspecific (blue)
#' #axes for lm models
#' mspace(shapes, links = links, mag = 0.7, axes = c(1,2)) %>%
#'   proj_axis(obj = intrasp_allo_mod, col = "red", lwd = 2, type = 2) %>%
#'   proj_axis(obj = intersp_allo_mod, col = "blue", lwd = 2, type = 2)
#' #for pls ordination
#' mspace(shapes, links = links, mag = 0.7, axes = c(1,2)) %>%
#'   proj_axis(obj = intrasp_allo_pls, col = "red", lwd = 2, type = 2) %>%
#'   proj_axis(obj = intersp_allo_pls, col = "blue", lwd = 2, type = 2)
proj_axis <- function(mspace, obj, axis = 1, mag = 1, pipe = TRUE, type = 3, ...) {

  args <- c(as.list(environment()), list(...))

  ext_shapes2d <- ax_transformation(obj = obj, axis = axis, mag = mag)
  ext_scores <- proj_eigen(x = ext_shapes2d, vectors = mspace$ordination$rotation,
                           center = mspace$ordination$center)

  if(.Device != "null device") {
    if(length(mspace$plotinfo$axes) > 1) {
      graphics::arrows(x0 = ext_scores[1, mspace$plotinfo$axes[1]],
                       y0 = ext_scores[1, mspace$plotinfo$axes[2]],
                       x1 = ext_scores[2, mspace$plotinfo$axes[1]],
                       y1 = ext_scores[2, mspace$plotinfo$axes[2]],
                       code = type, length = 0.1, ...)
    } else {
      print("\nphenotypic change vectors are omitted from univariate morphospaces")
    }
  }

  mspace$projected$shape_axis <- c(mspace$projected$shape_axis, list(ext_scores))

  mspace$plotinfo$type.axis <- c(mspace$plotinfo$type.axis, args$type)
  mspace$plotinfo$lwd.axis <- c(mspace$plotinfo$lwd.axis, args$lwd)
  mspace$plotinfo$lty.axis <- c(mspace$plotinfo$lty.axis, args$lty)
  mspace$plotinfo$col.axis <- c(mspace$plotinfo$col.axis, args$col)

  if(pipe == FALSE) return(invisible(ext_scores))
  if(pipe == TRUE) return(invisible(mspace))

}



################################################################################

#' Project phylogenetic structure into morphospace
#'
#' @description Project phylogenetic relationships among a set of shapes
#'   (representing the tips of a phylogenetic tree) into an existing bivariate
#'   morphospace.
#'
#' @param mspace An \code{"mspace"} object.
#' @param shapes Shape data, with 3rd margin names matching tip labels
#'   from \code{tree}.
#' @param tree A \code{"phylo"} object containing a phylogenetic tree.
#' @param pch.tips Numeric; symbol of the scatter points corresponding to the
#'   tips of the phylogeny.
#' @param col.tips Color of the hulls/ellipses and/or scatter points
#'   corresponding to the tips of the phylogeny.
#' @param bg.tips Background color of the scatter points corresponding to the
#'   tips of the phylogeny.
#' @param cex.tips Numeric; size of the scatter points corresponding to the
#'   tips of the phylogeny.
#' @param pch.nodes Numeric; symbol of the scatter points corresponding to the
#'   nodes of the phylogeny.
#' @param col.nodes Color of the hulls/ellipses and/or scatter points
#'   corresponding to the nodes of the phylogeny.
#' @param bg.nodes Background color of the scatter points corresponding to the
#'   nodes of the phylogeny.
#' @param cex.nodes Numeric; size of the scatter points corresponding to the
#'   nodes of the phylogeny.
#' @param pipe Logical; is the function being included in a pipe?
#' @param ... Further arguments passed to [graphics::lines()] (commonly, [pch],
#'   [col]/[bg] and [cex].
#'
#' @details The purpose of this function is twofold. First, it is meant to
#'   transform a morphospace into a phylomorphospace by projecting node shapes
#'   and phylogenetic relationships. To this end, a set of named shapes must be
#'   provided; \code{dim(shapes)[3]} must match \code{tree$tip.labels}. Second,
#'   this function can be used to retrieve the scores corresponding to nodes of
#'   the phylogenetic tree (\code{$projected$phylo_scores}, which can then be
#'   used to compute the node shapes using \code{\link{extract_shapes}}. The
#'   position of these shapes in morphospace is estimated using the
#'   squared-changes parsimony algorithm as performed by [phytools::fastAnc()].
#'
#' @return If a plot device with a morphospace is open, shapes representing the
#'   tips and nodes of the phylogenetic tree, as well as the lines connecting
#'   them, are projected into morphospace. If \code{pipe = FALSE} scores for
#'   nodes and tips of the phylogeny are returned invisibly.
#'   If \code{pipe = TRUE} the supplied \code{"mspace"} object will be modified
#'   by appending a \code{$phylo_Scores} and a \code{$phylo} slots to
#'   \code{$projected}, as well as by adding some graphical parameters (stored
#'   into the \code{$plotinfo} slot), and returned invisibly.
#'
#' @export
#'
#' @examples
#' #load and extract relevant data, packages and information
#' library(magrittr)
#' data("tails")
#' shapes <- tails$shapes
#' species <- tails$data$species
#' sp_shapes <- expected_shapes(shapes, species)
#' tree <- tails$tree
#' links <- tails$links
#'
#' #generate basic morphospace, add sampled shapes, species mean shapes, and
#' #phylogenetic structure
#' mspace(shapes, links = links, mag = 0.7, axes = c(1,2), cex.ldm = 0) %>%
#'   proj_shapes(shapes = shapes, col = c(1:13)[species], pch = 1,
#'               cex = 0.7) %>%
#'   proj_phylogeny(shapes = sp_shapes, tree = tree)
proj_phylogeny <- function(mspace, shapes = NULL, tree, pipe = TRUE,
                           pch.nodes = 16, col.nodes = "gray", bg.nodes = 1, cex.nodes = 0.8,
                           pch.tips = 16, bg.tips = 1, col.tips = 1, cex.tips = 1, ...) {

  args <- as.list(environment())

  dat <- shapes_mat(shapes)
  data2d <- dat$data2d
  datype <- dat$datype

  if(mspace$ordination$datype != datype) stop("shapes and mspace types are not compatible")

  if(is.null(rownames(data2d))) {
    stop("Shapes must be given names, which must match tree$tip.labels")
  } else {
    if(length(tree$tip.label) != nrow(data2d)) {
      stop("Number of tips in the tree does not match the number of shapes provided")
    } else {
      if(all(tree$tip.label %in% rownames(data2d))) {
        data2d <- data2d[tree$tip.label,]
      } else {
        stop("Names in phylogenetic tree does not match shape names")
      }
    }
  }

  tips_scores <- proj_eigen(x = data2d,
                            vectors = mspace$ordination$rotation,
                            center = mspace$ordination$center)

  nodes_scores <- apply(tips_scores, 2, phytools::fastAnc, tree = tree)
  phylo_scores <- rbind(tips_scores, nodes_scores)

  if(.Device != "null device") {

    if(length(mspace$plotinfo$axes) > 1) {
      for(i in seq_len(nrow(tree$edge))) {
        graphics::lines(rbind(phylo_scores[tree$edge[i, 1], mspace$plotinfo$axes],
                              phylo_scores[tree$edge[i, 2], mspace$plotinfo$axes]), ...)
      }

      tips <- seq_len(length(tree$tip.label))
      plot_biv_scatter(scores = phylo_scores[-tips,], pch = pch.nodes, bg = bg.nodes,
                       col = col.nodes, cex = cex.nodes)
      plot_biv_scatter(scores = phylo_scores[tips,], pch = pch.tips, bg = bg.tips,
                       col = col.tips, cex = cex.tips)

    } else {
      warning("phylogenetic relationships are not projected into univariate morphospaces")
    }
  }

  mspace$projected$phylo_scores <- phylo_scores
  mspace$projected$phylo <- tree

  if(is.null(args$col)) args$col <- 1

  mspace$plotinfo$col.phylo <- args$col
  mspace$plotinfo$lwd.phylo <- args$lwd
  mspace$plotinfo$lty.phylo <- args$lty
  mspace$plotinfo$col.nodes <- args$col.nodes
  mspace$plotinfo$bg.nodes <- args$bg.nodes
  mspace$plotinfo$pch.nodes <- args$pch.nodes
  mspace$plotinfo$cex.nodes <- args$cex.nodes
  mspace$plotinfo$col.tips <- args$col.tips
  mspace$plotinfo$bg.tips <- args$bg.tips
  mspace$plotinfo$pch.tips <- args$pch.tips
  mspace$plotinfo$cex.tips <- args$cex.tips


  if(pipe == FALSE) {
    return(invisible(phylo_scores))
  } else {
    return(invisible(mspace))
  }

}


################################################################################

#' Project landscape into morphospace
#'
#' @description Compute and project a landscape surface as a contour map over an
#'   existing morphospace.
#'
#' @param mspace An \code{"mspace"} object.
#' @param shapes Optional shape data. If provided, a landscape will be computed
#'   for the region of the morphospace encompassing that sample of shapes
#'   ("empirical landscape"). If \code{NULL}, the landscape will be computed for
#'   the set of background shape models ("theoretical landscape").
#' @param FUN An optional \emph{ad hoc} function to be applied to a set of
#'   shapes, stored in "two-dimensional" format, along its first margin (i.e.
#'   rows), and returning a single numeric value from each.
#' @param X An optional vector containing the values assigned to each shape
#'   (vector length and order must match those from the shapes provided in
#'   \code{shapes} or from the background shape models, depending on whether or
#'   not \code{shapes} have been provided).
#' @param linear Logical; whether to use linear interpolation (if \code{FALSE},
#'   a cubic spline interpolation is used instead. See [akima::interp()].
#' @param resolution Numeric; the resolution used for interpolation.
#' @param expand Numeric; Magnification factor to extend (adjust) the reach of
#'   the landscape, attained by extrapolating the x, y, and z values. Only
#'   available for \code{linear = FALSE}.
#' @param display Either \code{"contour"} or \code{"filled.contour"}. For
#'   bivariate landscapes only.
#' @param alpha Numeric; transparency factor for filled contours.
#' @param palette A function defining a color palette to use for landscape
#'   representation.
#' @param nlevels Number of levels (i.e., contours) to use for landscape
#'   representation.
#' @param drawlabels Logical; should the labels indicating the value of each
#'   surface contour be plotted?
#' @param lty Numeric; type of the lines depicting contours.
#' @param lwd Numeric; width of the lines depicting contours.
#' @param spar Numeric; smoothing parameter used to smooth landscape outline in
#'   univariate representations.
#' @param pipe Logical; is the function being included in a pipe?
#' @param ... Further arguments passed to \code{FUN}.
#'
#' @details The purpose of this function is to generate and depict a 2- (for
#'   univariate morphospaces) or 3-dimensonal (for bivariate morphospaces)
#'   surface (i.e., a landscape), interpolated from values assigned to the set
#'   of shapes projected into an existing morphospace. These can be a sample of
#'   shapes specified by the user, producing a surface for a specific region of
#'   the morphospace ("empirical landscapes"). Alternatively, the set of
#'   background shape models can be used to generate a surface for the entire
#'   observed morphospace ("theoretical morphospace"). Generally, the values
#'   that are interpolated will represent a variable measuring functional
#'   performance (although it can be any kind of continuous variable), and can
#'   be either provided directly through the \code{X} argument or computed
#'   automatically using an \emph{ad hoc} function through the \code{FUN}
#'   argument.
#'
#'   If the \code{FUN} argument is used, the function supplied must include a
#'   \code{model} argument feeding the \emph{ad hoc} function with a single
#'   shape (stored as a vector of shape descriptors), and return a single
#'   numeric value computed for or from that shape. If the \code{X} argument is
#'   used instead, values should be in the same order than the shapes provided
#'   in \code{shapes}, or than the background shape models (i.e., from left to
#'   right and from bottom to top; see \code{\link{morphogrid}},
#'   \code{\link{plot_morphogrid2d}} and \code{\link{plot_morphogrid3d}}) if
#'   \code{shapes = NULL}. Otherwise, landscape topography will be dissociated
#'   from the shapes represented in the morphospace. See examples below.
#'
#'   Directly providing the values to be interpolated using \code{X} can be
#'   useful when importing variables obtained using other software or pipelines.
#'   If the user desires to compute those variables for the background shape
#'   models, they can be extracted using \code{\link{extract_shapes}}) prior to
#'   exporting or feeding them to their preferred analytical protocol.
#'
#' @return If a plot device with a morphospace is open, the landscape surface is
#'   projected into it as a contour map using [akima::interp()]. If
#'   \code{pipe = FALSE}, a list containing the x, y and z values used to plot
#'   the landscape (x and z for univariate morphospaces) is returned invisibly.
#'   If \code{pipe = TRUE} the supplied \code{"mspace"} object will be modified
#'   by appending a \code{$landsc} slot to \code{$projected}, as well as by
#'   adding some graphical parameters (stored into the \code{$plotinfo} slot),
#'   and returned invisibly.
#'
#' @seealso \code{\link{morphogrid}}, \code{\link{mspace}},
#'   \code{\link{plot_morphogrid2d}}, \code{\link{plot_morphogrid3d}},
#'   \code{\link{extract_shapes}}
#'
#' @export
#'
#' @examples
#' #load data and packages
#' library(geomorph)
#' library(Morpho)
#' library(Momocs)
#'
#' data("tails")
#' shapes <- tails$shapes
#' links <- tails$links
#' type <- tails$data$type
#'
#'
#' #Compute and plot adaptive landscape for wing tail shape:
#'
#' ##Using the FUN argument
#'
#' #a function to compute lift/drag ratio on tail shape:
#' computeLD <- function(model, MCS = FALSE) {
#'
#'   tail <- matrix(model, ncol = 2, byrow = TRUE)
#'   Ax <- tail[9, 1] ; Ay <- tail[9, 2]
#'   Bx <- tail[5, 1] ; By <- tail[5, 2]
#'   X <- tail[7, 1] ; Y <- tail[7, 2]
#'
#'   # determine position of the tip of the inner rectrix relative to the tips
#'   # of the outermost rectrices (positive: anterior to the ORT; negative:
#'   # posterior to the ORT)
#'   tip_pos1 <- sign((Bx - Ax) * (Y - Ay) - (By - Ay) * (X - Ax))
#'
#'   # do the same but for the inner outer rectrices
#'   Ax <- tail[8, 1] ; Ay <- tail[8, 2]
#'   Bx <- tail[6, 1] ; By <- tail[6, 2]
#'   X <- tail[7, 1] ; Y <- tail[7, 2]
#'
#'   # determine position of the tip of the inner rectrix relative to the tips
#'   # of the outermost rectrices (positive: anterior to the ORT; negative:
#'   # posterior to the ORT)
#'   tip_pos2 <- sign((Bx - Ax) * (Y - Ay) - (By - Ay) * (X - Ax))
#'
#'   # compute the area of the polygon representing the whole tail
#'   tail_area <- coo_area(tail[c(1, 3, 5, 6, 7, 8, 9, 4, 2, 1), ])
#'
#'
#'   # if posterior, compute MCS as the distance between the ORTs, and the
#'   # lifting area as the polygon enclosed by the ORTs and tail base
#'   if(tip_pos1 == -1) {
#'     mcs <- dist(rbind(tail[5, ], tail[9, ]))
#'     lifting_area <- coo_area(tail[c(1, 3, 5, 9, 4, 2, 1),])
#'   }
#'
#'   # if anterior, find the maximum continuous span as the line perpendicular
#'   # to the saggital axis that intersects the ORTs and the IRT
#'   if(tip_pos1 == 1 & tip_pos2 == 1) {
#'
#'     # center tail around the IRT and set MCS as a vertical line passing
#'     # through the origin
#'     tail_cent <- cbind(tail[, 1] - tail[7, 1],
#'                        tail[, 2] - tail[7, 2])
#'     mcs_vec <- c(0,10e10)
#'
#'     ort1 <- rbind(tail_cent[4, ], tail_cent[9, ])
#'     ort1_vec <- lm(ort1[, 2] ~ ort1[, 1])$coef
#'
#'     ort2 <- rbind(tail_cent[3, ], tail_cent[5, ])
#'     ort2_vec <- lm(ort2[, 2] ~ ort2[, 1])$coef
#'
#'     # find tips of MCS
#'     A <- matrix(c(mcs_vec[2], -1, ort1_vec[2], -1), byrow = TRUE, nrow = 2)
#'     b <- c(-mcs_vec[1], -ort1_vec[1])
#'     p1 <- solve(A, b)
#'
#'     A <- matrix(c(mcs_vec[2], -1, ort2_vec[2], -1), byrow = TRUE, nrow = 2)
#'     b <- c(-mcs_vec[1], -ort2_vec[1])
#'     p2 <- solve(A, b)
#'
#'     # define "new tail base", enclosing the polygon using the tips of the MCS
#'     newbase <- rbind(tail_cent[c(2, 4), ],
#'                      p1,
#'                      tail_cent[7, ],
#'                      p2,
#'                      tail_cent[c(3, 1), ])
#'
#'     # compute MCS and lifting area
#'     lifting_area <- coo_area(newbase[c(1, 2, 3, 4, 5, 6, 7, 1),])
#'     mcs <- dist(rbind(p1, p2))
#'
#'   }
#'
#'   # if anterior, find the maximum continuous span as the line perpendicular
#'   # to the saggital axis that intersects the ORTs and the IRT
#'   if(tip_pos1 == 1 & tip_pos2 == -1) {
#'
#'     # center tail around the IRT and set MCS as a vertical line passing
#'     # through the origin
#'     tail_cent <- cbind(tail[, 1] - mean(tail[6, 1], tail[8, 1]),
#'                        tail[, 2] - mean(tail[6, 2], tail[8, 2]))
#'     mcs_vec <- c(0,10e10)
#'
#'     ort1 <- rbind(tail_cent[4, ], tail_cent[9, ])
#'     ort1_vec <- lm(ort1[, 2] ~ ort1[, 1])$coef
#'
#'     ort2 <- rbind(tail_cent[3, ], tail_cent[5, ])
#'     ort2_vec <- lm(ort2[, 2] ~ ort2[, 1])$coef
#'
#'     # find tips of MCS
#'     A <- matrix(c(mcs_vec[2], -1, ort1_vec[2], -1), byrow = TRUE, nrow = 2)
#'     b <- c(-mcs_vec[1], -ort1_vec[1])
#'     p1 <- solve(A, b)
#'
#'     A <- matrix(c(mcs_vec[2], -1, ort2_vec[2], -1), byrow = TRUE, nrow = 2)
#'     b <- c(-mcs_vec[1], -ort2_vec[1])
#'     p2 <- solve(A, b)
#'
#'     # define "new tail base", enclosing the polygon using the tips of the MCS
#'     newbase <- rbind(tail_cent[c(2, 4), ],
#'                      p1,
#'                      tail_cent[7, ],
#'                      p2,
#'                      tail_cent[c(3, 1), ])
#'
#'     # compute MCS and lifting area
#'     lifting_area <- coo_area(newbase[c(1, 2, 3, 4, 5, 6, 7, 1),])
#'     mcs <- dist(rbind(p1, p2))
#'
#'   }
#'
#'   #compute and return lift/drag ratio
#'   if(!MCS) LD_ratio <- lifting_area / tail_area
#'   if(MCS)  LD_ratio <- (mcs ^ 2) / tail_area
#'
#'   return(as.numeric(LD_ratio))
#'
#' }
#'
#' ###"Theoretical" landscapes
#'
#' #plot morphospace with its associated adaptive landscape
#' mspace(shapes, links = links, nh = 8, nv = 8, size.model = 1.5,
#'        cex.ldm = 0) %>%
#'   proj_shapes(shapes, pch = 16) %>%
#'   proj_landscape(nlevels = 60, FUN = computeLD, expand = 1.2, lwd = 2)
#'
#' ##Using the X argument
#'
#' #first, create morphospace and extract background shapes without plotting
#' msp <- mspace(shapes, links = links, nh = 8, nv = 8, size.model = 1.5,
#'               cex.ldm = 0, plot = FALSE)
#' shapemodels2d <- two.d.array(extract_shapes(msp)$shapes)
#'
#' #run computeLD (defined above) through the "two-dimensional" matrix of shapes
#' #models (this is the same thing the proj_landscape function is doing
#' #internally when FUN is used, but this vector could be replaced with some
#' #other variable obtained in a different way)
#' LDs <- apply(X = shapemodels2d, FUN = computeLD, MARGIN = 1)
#'
#' #second, plot morphospace with its associated adaptive landscape
#' msp <- mspace(shapes, links = links, nh = 8, nv = 8, size.model = 1.5,
#'               cex.ldm = 0) %>%
#'   proj_shapes(shapes, pch = 16) %>%
#'   proj_landscape(nlevels = 20, X = LDs, expand = 1.2,
#'                  palette = terrain.colors, display = "filled.contour")
#'
#' #add scalebar using plot_mspace()
#' plot_mspace(msp, scalebar = TRUE)
#'
#' #be careful to use the same morphospaces in the second and first steps.
#' #For example, retaining the LD values computed above but changing the axes
#' #represented by mspace will result in the same surface wrongly being
#' #projected over a different morphospace!
#' mspace(shapes, links = links, nh = 8, nv = 8, size.model = 1.5,
#'        cex.ldm = 0, axes = c(2,3)) %>%
#'   proj_shapes(shapes, pch = 16) %>%
#'   proj_landscape(nlevels = 20, X = LDs, expand = 1.2,
#'                  palette = terrain.colors, display = "filled.contour")
#'
#' ###"Empirical" landscapes
#'
#' #it is essentially the same, but providing a set of shapes with the shapes
#' #argument of proj_landscape. Let's compute it it only for Tyrannus species with
#' #non-deep forked tail shapes:
#' mspace(shapes, links = links, nh = 8, nv = 8, size.model = 1.5,
#'        cex.ldm = 0) %>%
#'   proj_shapes(shapes, pch = 16) %>%
#'   proj_landscape(shapes = shapes[,,type == "NDF"], nlevels = 20, linear = TRUE,
#'                  FUN = computeLD, expand = 1.2, display = "filled.contour")
#'
#' #in this case, the resolution of the projected surface can be improved using
#' #the argument resolution (however, be aware this can be computationally intense!
#' #especially if a theoretical landscape is being computed)
#' mspace(shapes, links = links, nh = 8, nv = 8, size.model = 1.5,
#'        cex.ldm = 0) %>%
#'   proj_shapes(shapes, pch = 16) %>%
#'   proj_landscape(shapes = shapes[,,type == "NDF"], nlevels = 20, linear = TRUE,
#'                  resolution = 200, FUN = computeLD, expand = 1.2,
#'                  display = "filled.contour")
proj_landscape <- function(mspace, shapes = NULL, FUN = NULL, X = NULL, linear = FALSE,
                           resolution = 50, expand = 1, display = "contour", nlevels = 50,
                           palette = grDevices::heat.colors, alpha = 0.5, lwd = 1, lty = 1,
                           drawlabels = FALSE, spar = 0.5, pipe = TRUE, ...) {

  args <- c(as.list(environment()), list(...))

  if(is.null(shapes)) {
    extrap <- TRUE
    type <- "theoretical"

    data2d <- shapes_mat(mspace$projected$shapemodels)$data2d

  } else {
    extrap <- FALSE
    type <- "empirical"

    dat <- shapes_mat(shapes)
    data2d <- dat$data2d
    datype <- dat$datype

    if(mspace$ordination$datype != datype) stop("shapes and mspace types are not compatible")
  }

  gridcoords <- proj_eigen(x = data2d, vectors = mspace$ordination$rotation[, mspace$plotinfo$axes],
                           center = mspace$ordination$center)

  if(is.null(X)) X <- apply(X = data2d, FUN = FUN, MARGIN = 1, ...)

  if(type == "theoretical") {
    frame <- list(xlim = graphics::par("usr")[1:2], ylim = graphics::par("usr")[3:4])
  } else {
    frame <- list(xlim = range(gridcoords[,1]))
    if(length(mspace$plotinfo$axes) > 1) frame$ylim <- range(gridcoords[,2])
  }

  xlim <- frame$xlim * expand
  xo <- seq(from = xlim[1], to = xlim[2], length.out = resolution)

  if(length(mspace$plotinfo$axes) == 1) {
    gridcoords <- expand.grid(gridcoords, 0:(mspace$plotinfo$nh - 1))
    X <- rep(X, times = mspace$plotinfo$nh)

    yo <- 0
  } else {
    ylim <- frame$ylim * expand
    yo <- seq(from = ylim[1], to = ylim[2], length.out = resolution)

    xo[which.x0 <- which.min(abs(xo))] <- 0
    yo[which.y0 <- which.min(abs(yo))] <- 0
  }


  landscape <- suppressWarnings({akima::interp(x = gridcoords[, 1], y = gridcoords[, 2],
                                               z = X, extrap = extrap, linear = linear,
                                               xo = xo, yo = yo)})


  if(type == "empirical" & length(mspace$plotinfo$axes) == 1) {
    spline <- stats::smooth.spline(x = landscape$x[!is.na(landscape$z)],
                                   y = landscape$z[!is.na(landscape$z)], spar = spar)
    smoothed_landsc <- Morpho::equidistantCurve(x = cbind(spline$x, spline$y), n = resolution)

    landscape$x <- smoothed_landsc[,1]
    landscape$z <- smoothed_landsc[,2]
  }

  if(length(mspace$plotinfo$axes) > 1) {
    landscape$which.y0 <- which.y0
    landscape$which.x0 <- which.x0
  }


  if(.Device != "null device") {

    if(length(mspace$plotinfo$axes) > 1) {

      cols <- palette(n = nlevels)
      levs <- seq(from = min(landscape$z, na.rm = TRUE),
                  to = max(landscape$z, na.rm = TRUE),
                  length.out = nlevels)

      plot_biv_landscape(landscape = landscape, display = display, levels = levs, col = cols,
                         lwd = lwd, lty = lty, drawlabels = drawlabels, alpha = alpha, type = type)
    } else {

      cols <- palette(n = length(landscape$z))[order(order(landscape$z))]

      plot_univ_landscape(landscape = landscape, drawlabels = drawlabels,
                          col = cols, lwd = lwd)
    }
  }

  mspace$projected$landsc <- landscape
  mspace$projected$landsc$type <- type

  mspace$plotinfo$palette.landsc <- args$palette
  mspace$plotinfo$display.landsc <- args$display
  mspace$plotinfo$nlevels.landsc <- args$nlevels
  mspace$plotinfo$resolution.landsc <- args$resolution
  mspace$plotinfo$lty.landsc <- args$lty
  mspace$plotinfo$lwd.landsc <- args$lwd
  mspace$plotinfo$alpha.landsc <- args$alpha
  mspace$plotinfo$drawlabels.landsc <- args$drawlabels


  if(pipe == FALSE) {
    return(invisible(landscape))
  } else {
    return(invisible(mspace))
  }
}


################################################################################

#' Print \code{"mspace"} objects
#'
#' @description Print contents of \code{"mspace"} objects in a friendly way :)
#'
#' @param mspace An \code{"mspace"} object.
#'
#' @details Will return information about how the ordination was obtained
#'   (method and data) as well as the elements projected into it.
#'
#' @export
print.mspace <- function(mspace){

  if(mspace$ordination$ordtype == "prcomp") ordtype <- "Principal Component Analysis"
  if(mspace$ordination$ordtype == "bg_prcomp") ordtype <- "Between-Groups Principal Component Analysis"
  if(mspace$ordination$ordtype == "phy_prcomp") ordtype <- "Phylogenetic Principal Component Analysis"
  if(mspace$ordination$ordtype == "pls_shapes") ordtype <- "Two-Block Partial Least Squares"

  datmssg <- paste0("data from ", nrow(mspace$ordination$x), " shapes (quantified using ")
  if(mspace$ordination$datype == "landm") {
    datype <- paste0(datmssg, mspace$plotinfo$p, " (", mspace$plotinfo$k, "D) landmark coordinates)")
    nbs <- dim(mspace$projected$shapemodels)[3]
  }
  if(mspace$ordination$datype == "fcoef") {
    datype <- paste0(datmssg, dim(mspace$projected$shapemodels)[2], " Fourier coefficients)")
    nbs <- dim(mspace$projected$shapemodels)[1]
  }

  if(ncol(mspace$ordination$rotation) > 1) nax <- " axes" else nax <- " axis"

  cat("an mspace object containing")
  cat(paste0("\n* an ordination space made out of ", ncol(mspace$ordination$rotation),
             nax, " built using:"))
  cat(paste0("\n   - ", ordtype))
  cat(paste0("\n   - ", datype))

  cat(paste0("\n* ", length(mspace$projected), " elements projected:"))
  cat(paste0("\n   - a set of ", paste0(nbs, " background shape models")))
  if(!is.null(mspace$projected$scores)) {
    cat(paste0("\n   - a set of ", nrow(mspace$projected$scores), " scores"))
  }
  if(!is.null(mspace$projected$gr_class)) {
    cat(paste0("\n   - a classification into ", nlevels(mspace$projected$gr_class), " groups"))
  }
  if(!is.null(mspace$projected$shape_axis)) {
    if(length(mspace$projected$shape_axis) > 1) nax <- " axes" else nax <- " axis"
    cat(paste0("\n   - ", length(mspace$projected$shape_axis), " shape", nax))
  }
  if(!is.null(mspace$projected$phylo_scores)) {
    cat(paste0("\n   - a phylogenetic structure for ", mspace$projected$phylo$Nnode + 1,
               " tips and ", mspace$projected$phylo$Nnode, " nodes"))
  }
  if(!is.null(mspace$projected$landsc)) {
    if(length(mspace$plotinfo$axes) == 1) {
      wax <- paste0("axis ", mspace$plotinfo$axes[1])
    } else {
      wax <- paste0("axes ", mspace$plotinfo$axes[1], " and ", mspace$plotinfo$axes[2])
    }
    cat(paste0("\n   - a ", mspace$projected$landsc$type, " landscape computed for ", wax))
  }

}


################################################################################

#' Plot morphospaces and combine them with other variables
#'
#' @description Flexible representation of morphospaces, including their
#'   combination with either other variables or a phylogeny.
#'
#' @param mspace An \code{"mspace"} object created using the
#'   \code{\link{mspace()}} %>% \code{proj_*} pipeline.
#' @param axes Numeric of length 1 or 2, indicating the axes to be plotted.
#' @param links A list with the indices of the coordinates defining the
#'   wireframe (following the format used in \code{Morpho}).
#' @param template Either a 2-column matrix with landmarks/semilandmarks and
#'    template curves coordinates (for 2D shape data) or a \code{"mesh3d"}
#'    object representing the mean shape of the sample (for 3D shape data).
#' @param x Optional vector with a non-morphometric variable to be plotted in
#'   the x axis. Alternatively, a \code{"phylo"} object can be provided.
#' @param y Optional vector with a non-morphometric variable to be plotted in
#'   the y axis. Alternatively, a \code{"phylo"} object can be provided.
#' @param nh Numeric; the number of shape models along the x axis.
#' @param nv Numeric; the number of shape models along the y axis.
#' @param mag Numeric; magnifying factor for shape models.
#' @param rescale Logical; whether to re-scale background shape models so shape
#'    variation is shown more clearly.
#' @param invax Optional numeric indicating which of the axes provided in
#'   \code{axes} needs to be inverted (options are \code{1}, \code{2} or
#'   \code{c(1,2)}).
#' @param adj_frame Numeric of length 2, providing \emph{a posteriori} scaling
#'   factors for the width and height of the frame, respectively.
#' @param points Logical; whether to plot the scatter points corresponding to
#'   the sampled shapes stored in \code{mspace$projected$scores}.
#' @param models Logical; whether to plot background shape models.
#' @param mshapes Logical; whether to plot the scatter points corresponding to
#'   groups' mean shapes stored in \code{mspace$projected$gr_centroids}.
#' @param groups Logical; whether to plot the convex hulls/confidence ellipses
#'   enclosing the groups stored in \code{mspace$projected$gr_class}.
#' @param phylo Logical; whether to plot phylogenetic relationships stored
#'   in \code{mspace$projected$phylo}.
#' @param shapeax Logical; whether to plot morphometric axes stored in
#'   \code{mspace$projected$shape_axis}.
#' @param landsc Logical; whether to plot landscape surface stored in
#'   \code{mspace$projected$landsc}.
#' @param legend Logical; whether to show legend for groups.
#' @param scalebar Logical; whether to show scalebar for landscapes.
#' @param cex.legend Numeric; size of legend labels/symbols.
#' @param size.models Numeric; size factor for shape models.
#' @param asp.models Numeric; y/x aspect ratio of shape models.
#' @param rot.models Numeric; angle (in degrees) to rotate shape models.
#' @param col.models Color for wireframes/outlines.
#' @param bg.models Background color for outlines/meshes.
#' @param lwd.models Numeric; width of the lines in wireframes/outlines.
#' @param alpha.models Numeric; transparency factor for background models (3D
#'   only).
#' @param cex.ldm Numeric; size of landmarks/semilandmarks in the background
#'   models.
#' @param col.ldm Color of landmarks/semilandmarks in the background models.
#' @param pch.points Numeric; symbol of the scatter points corresponding to
#'   sampled shapes.
#' @param col.points Color of the scatter points corresponding to sampled
#'   shapes.
#' @param bg.points Background color of the scatter points corresponding
#'   to sampled shapes.
#' @param cex.points Numeric; size of the scatter points corresponding to
#'   sampled shapes.
#' @param density.points Logical; whether to add density distribution for points
#'   (univariate ordinations only). Overriden by \code{density.groups = TRUE}
#' @param col.groups Color of the hulls/ellipses and/or scatter points
#'   corresponding to groups mean shapes.
#' @param bg.groups Background color of the scatter points corresponding
#'   to groups mean shapes.
#' @param pch.groups Numeric; symbol of the scatter points corresponding to
#'   groups mean shapes.
#' @param cex.groups Numeric; size of the scatter points corresponding to
#'   groups mean shapes.
#' @param ellipse.groups Logical; whether to plot confidence ellipses
#'   (if \code{FALSE} convex hulls will be used instead).
#' @param conflev.groups Numeric; confidence level used for confidence
#'   ellipse(s).
#' @param lwd.groups Numeric; width of the lines in groups' ellipses/hulls.
#' @param lty.groups Numeric; type of the lines in groups' ellipses/hulls.
#' @param alpha.groups Numeric; transparency factor for groups'
#'   ellipses/hulls/density distributions.
#' @param density.groups Logical; whether to add density distribution for groups
#'   (univariate ordinations only).
#' @param pch.tips Numeric; symbol of the scatter points corresponding to the
#'   tips of the phylogeny.
#' @param col.tips Color of the hulls/ellipses and/or scatter points
#'   corresponding to the tips of the phylogeny.
#' @param bg.tips Background color of the scatter points corresponding to the
#'   tips of the phylogeny.
#' @param cex.tips Numeric; size of the scatter points corresponding to the
#'   tips of the phylogeny.
#' @param pch.nodes Numeric; symbol of the scatter points corresponding to the
#'   nodes of the phylogeny.
#' @param col.nodes Color of the hulls/ellipses and/or scatter points
#'   corresponding to the nodes of the phylogeny.
#' @param bg.nodes Background color of the scatter points corresponding to the
#'   nodes of the phylogeny.
#' @param cex.nodes Numeric; size of the scatter points corresponding to the
#'   nodes of the phylogeny.
#' @param lwd.phylo Numeric; width of the lines depicting phylogenetic
#'   branches.
#' @param lty.phylo Numeric; type of the lines depicting phylogenetic branches.
#' @param col.phylo Numeric; color of the lines depicting phylogenetic branches.
#' @param type.axis Integer; type of arrows ([0] = no arrow; [1] = pointing
#'   towards the maximum; [2] = pointing towards the maximum, [3] = pointing in
#'   both directions).
#' @param lwd.axis Numeric; width of the lines depicting a morphometric axis.
#' @param lty.axis Numeric; type of the lines depicting  a morphometric axis.
#' @param col.axis Numeric; color of the lines depicting a morphometric axis.
#' @param display.landsc How to display landscape representation; options are
#'   \code{"contour"} and \code{"filled.contour"}. For bivariate landscapes
#'   only.
#' @param alpha.landsc Numeric; transparency factor for filled contours
#'   depicting landscapes.
#' @param palette.landsc A function defining a color palette to use for
#'   landscape representation.
#' @param nlevels.landsc Number of levels (i.e., contours) to use in landscape
#'   representation.
#' @param drawlabels.landsc Logical; should the labels indicating the value of
#'   each surface contour be plotted?
#' @param lty.landsc Numeric; type of the contour lines depicting landscapes.
#' @param lwd.landsc Numeric; width of the contour lines depicting landscapes.
#' @param xlim,ylim,xlab,ylab,asp Standard arguments passed to the generic
#'   [graphics::plot()] function.
#'
#' @details This function allows to regenerate/tweak morphospaces contained in
#'   \code{"mspace"} objects already in existence. By default, [plot_mspace]
#'   regenerates the morphospace with all its projected elements preserving the
#'   graphical parameters used originally during the \code{\link{mspace}} +
#'   \code{proj_*} pipeline (and stored in \code{mspace$plotinfo}). However,
#'   all the graphical parameters can be modified to customize morphospaces
#'   graphical representations. Also, [plot_mspace] can be used to add a legend
#'   and/or a scalebar to aid identification of groups and interpretation of
#'   landscapes, respectively.
#'
#'   In addition, this function expands the range of graphical options available
#'   beyond 'pure' morphospaces. If a numeric non-shape variable (assumed to be
#'   measured for the same specimens in \code{mspace$projected$scores}) is fed
#'   to one of \code{x} or \code{y}, a 'hybrid' morphospace is produced (i.e.
#'   the bivariate plot will be constructed from the combination of \code{x} or
#'   \code{y} and a morphometric axis; background shape models will represent
#'   variation only for the latter). If instead a \code{"phylo"} object (assumed
#'   to describe the phylogenetic relationships among tips scores stored in
#'   \code{mspace$projected$phylo_scores}) is provided for either \code{x} or
#'   \code{y}, a vertical or horizontal phenogram will be deployed (the x/y axis
#'   range will correspond to branch lengths, so caution should be exercised
#'   when interpreting the output).
#'
#'   \emph{Note}: when regenerating landscapes, it's important to keep in mind
#'   that this surface has been calculated for a specific set of shapes and
#'   ordination axes. Hence, its regeneration using [plot_mspace] is only
#'   warranted if the axes being depicted are the same than those used when the
#'   surface landscape was originally computed using \code{\link{mspace}} +
#'   \code{\link{proj_landscape}}. The only exception is when one of the
#'   original axes (i.e., those specified with the \code{axes} argument in
#'   \code{\link{mspace}}) is dropped (i.e., not specified with the \code{axes}
#'   argument of \code{\link{plot_mspace}}). This will result in the collapse of
#'   the 3D landscape projected into a bivariate morphospace into a 2D landscape
#'   projected into a univariate one.
#'
#' @export
#'
#' @examples
#' #load and extract relevant data, packages and information
#' library(magrittr)
#' data("tails")
#' shapes <- tails$shapes
#' species <- tails$data$species
#' sizes <- tails$sizes
#' sp_shapes <- expected_shapes(shapes, species)
#' sp_sizes <- cbind(tapply(sizes, species, mean))
#' tree <- tails$tree
#' links <- tails$links
#'
#' #generate basic morphospace, add sampled shapes, species mean shapes, species
#' #classification, and phylogenetic structure
#' msp <- mspace(shapes, axes = c(1,2), plot = FALSE) %>%
#'  proj_shapes(shapes = shapes) %>%
#'  proj_consensus(shapes = sp_shapes) %>%
#'  proj_groups(groups = species) %>%
#'  proj_phylogeny(tree = tree)
#'
#' ##Plotting 'pure' morphospaces:
#'
#' #plot mspace object as it is
#' plot_mspace(msp, axes = c(1,2),
#'             points = TRUE, groups = FALSE, mshapes = FALSE, phylo = FALSE)
#'
#' #add colors for points, by species
#' plot_mspace(msp, axes = c(1,2), col.points = species, col.groups = 1:nlevels(species),
#'             points = TRUE, groups = FALSE, mshapes = FALSE, phylo = FALSE)
#'
#' #add links
#' plot_mspace(msp, axes = c(1,2), links = links,
#'             col.points = species, col.groups = 1:nlevels(species),
#'             points = TRUE, groups = FALSE, mshapes = FALSE, phylo = FALSE)
#'
#' #change number and sizes of shape models in the background
#' plot_mspace(msp, axes = c(1,2), nh = 2, nv = 2, links = links, size.models = 0.5,
#'             col.points = species, col.groups = 1:nlevels(species),
#'             points = TRUE, groups = FALSE, mshapes = FALSE, phylo = FALSE)
#'
#' #magnify deformation and highlight landmarks
#' plot_mspace(msp, axes = c(1,2), mag = 1.5, nh = 2, nv = 2, links = links, size.models = 0.5,
#'             col.points = species, col.groups = 1:nlevels(species), cex.ldm = 5, col.ldm = "red",
#'             points = TRUE, groups = FALSE, mshapes = FALSE, phylo = FALSE)
#'
#' #change axes 1,2 for 1,3
#' plot_mspace(msp, axes = c(1,3), mag = 1.5, nh = 2, nv = 2, links = links, size.models = 0.5,
#'             col.points = species, col.groups = 1:nlevels(species), cex.ldm = 5, col.ldm = "red",
#'             points = TRUE, groups = FALSE, mshapes = FALSE, phylo = FALSE)
#'
#' #represent ranges and centroids of groups
#' plot_mspace(msp, axes = c(1,3), mag = 1.5, nh = 2, nv = 2, links = links, size.models = 0.5,
#'             col.points = species, col.groups = 1:nlevels(species), cex.ldm = 5, col.ldm = "red",
#'             points = TRUE, groups = TRUE, mshapes = TRUE, phylo = FALSE)
#'
#' #add phylogeny
#' plot_mspace(msp, axes = c(1,3), mag = 1.5, nh = 2, nv = 2, links = links, size.models = 0.5,
#'             col.points = species, col.groups = 1:nlevels(species), cex.ldm = 5, col.ldm = "red",
#'             points = TRUE, groups = TRUE, mshapes = TRUE, phylo = TRUE)
#'
#'
#'
#' ##Plotting 'hybrid' morphospaces:
#'
#' #plot size against first PC
#' plot_mspace(msp, x = sizes,  axes = 1, links = links, col.points = species,
#'             col.groups = 1:nlevels(species), pch.points = 16,
#'             points = TRUE, groups = FALSE, mshapes = FALSE, phylo = FALSE, xlab = "Centroid size")
#'
#'
#'
#' ##Plotting phenograms:
#'
#' #plot horizontal phenogram against PC1
#' plot_mspace(msp, x = tree,  axes = 1, links = links, col.points = species,
#'             col.groups = 1:nlevels(species), pch.points = 16,
#'             points = TRUE, groups = FALSE, mshapes = FALSE, phylo = FALSE, xlab = "Branch lengths")
#'
#' #plot horizontal phenogram against PC1
#' plot_mspace(msp, y = tree,  axes = 2, links = links, col.points = species,
#'             col.groups = 1:nlevels(species), pch.points = 16,
#'             points = TRUE, groups = FALSE, mshapes = FALSE, phylo = FALSE, ylab = "Branch lengths")

# plot_mspace <- function(mspace,
#                         axes = NULL,
#                         links = NULL,
#                         template = NULL,
#                         x = NULL,
#                         y = NULL,
#                         nh,
#                         nv,
#                         mag,
#                         invax = NULL,
#                         adj_frame = c(1,1),
#                         points = TRUE,
#                         models = TRUE,
#                         mshapes = TRUE,
#                         groups = TRUE,
#                         phylo = TRUE,
#                         shapeax = TRUE,
#                         landsc = TRUE,
#                         legend = FALSE,
#                         scalebar = FALSE,
#                         cex.legend = 1,
#                         asp = NA,
#                         xlab,
#                         ylab,
#                         xlim = NULL,
#                         ylim = NULL,
#                         size.models,
#                         asp.models,
#                         rot.models,
#                         col.models,
#                         bg.models,
#                         lwd.models,
#                         alpha.models,
#                         cex.ldm,
#                         col.ldm,
#                         pch.points = 1,
#                         col.points = 1,
#                         bg.points = 1,
#                         cex.points = 1,
#                         density.points = TRUE,
#                         col.groups = 1,
#                         bg.groups = 1,
#                         pch.groups = 16,
#                         cex.groups = 1,
#                         ellipse.groups = mspace$plotinfo$ellipse.groups,
#                         conflev.groups = 0.95,
#                         lwd.groups = 1,
#                         lty.groups = 1,
#                         alpha.groups = 0,
#                         density.groups = TRUE,
#                         lwd.phylo = 1,
#                         lty.phylo = 1,
#                         col.phylo = 1,
#                         lwd.axis = 1,
#                         lty.axis = 1,
#                         col.axis = 1,
#                         palette.landsc = heat.colors,
#                         display.landsc,
#                         alpha.landsc,
#                         nlevels.landsc = 50,
#                         drawlabels.landsc = FALSE,
#                         lty.landsc = 1,
#                         lwd.landsc = 1) {
#
#
#   supplied <- names(as.list(match.call()))[-1]
#   new_args <- lapply(seq_len(length(supplied)), function(i) {get(supplied[i])})
#   names(new_args) <- supplied
#   inh_args <- mspace$plotinfo
#   merged_args <- utils::modifyList(inh_args, new_args)
#   args <- merged_args
#
#
#   if(!is.null(invax)) {
#     mspace$ordination$x[,axes[invax]] <- mspace$ordination$x[,axes[invax]]  * -1
#     mspace$ordination$rotation[,axes[invax]] <- mspace$ordination$rotation[,axes[invax]] * -1
#     if(!is.null(mspace$projected$gr_centroids)) {
#       mspace$projected$gr_centroids[,axes[invax]] <- mspace$projected$gr_centroids[,axes[invax]] * -1
#     }
#     if(!is.null(mspace$projected$phylo_scores)) {
#       mspace$projected$phylo_scores[,axes[invax]] <- mspace$projected$phylo_scores[,axes[invax]] * -1
#     }
#   }
#
#   ordination <- mspace$ordination
#
#   if(!is.null(x) & !is.null(y)) stop("Only one of x or y can be specified")
#
#   if(any(legend, scalebar)) {
#
#     layout_mat <- matrix(c(rep(1, 4), 2), nrow = 1)
#     if(all(legend, scalebar)) layout_mat <- rbind(layout_mat, c(rep(1, 4), 3))
#
#     orig_par <- graphics::par(names(graphics::par())[-c(13, 19, 21:23, 54)])
#     on.exit(graphics::par(orig_par))
#
#     graphics::layout(layout_mat)
#     graphics::par(mar = c(5.1, 4.1, 4.1, 2.1))
#
#   }
#
#   if(is.null(x) & is.null(y)) { #if neither x nor y have been provided, plot pure morphospace
#
#     if(ncol(ordination$x) == 1 | length(args$axes) == 1) {
#       y <- rep(1, nrow(ordination$x))
#       args$ylim <- c(0, 1)
#       args$ylab <- "relative density"
#     } else {
#       y <- NULL
#     }
#
#     if(mspace$plotinfo$k == 3 & mspace$ordination$datype == "landm") {
#
#       shapemodels <- morphogrid(ordination = ordination, axes = args$axes, datype = mspace$ordination$datype,
#                                 template = NULL, x = NULL, y = y, p = mspace$plotinfo$p,
#                                 k = mspace$plotinfo$k, nh = args$nh, nv = args$nv, mag = args$mag,
#                                 asp = args$asp, xlim = args$xlim, ylim = args$ylim, rot.models = args$rot.models,
#                                 size.models = args$size.models, asp.models = args$asp.models)
#
#       refshape <- matrix(rev_eigen(0,
#                                    ordination$rotation[, 1],
#                                    ordination$center),
#                          nrow = mspace$plotinfo$p, ncol = mspace$plotinfo$k, byrow = TRUE)
#
#       args$xlim <- range(ordination$x[,args$axes[1]])
#       if(ncol(ordination$x) > 1 | length(args$axes) > 1) args$ylim <- range(ordination$x[, args$axes[2]])
#
#       plot_morphogrid3d(x = NULL, y = y, morphogrid = shapemodels, refshape = refshape,
#                         template = args$template, links = args$links, ordtype = mspace$ordination$ordtype,
#                         axes = args$axes, xlim = args$xlim, ylim = args$ylim, adj_frame = args$adj_frame,
#                         xlab = args$xlab, ylab = args$ylab, cex.ldm = args$cex.ldm,
#                         col.ldm = args$col.ldm, col.models = args$col.models, lwd.models = args$lwd.models,
#                         bg.models = args$bg.models,  size.models = args$size.models,
#                         asp.models = args$asp.models, alpha.models = args$alpha.models, models = args$models)
#     } else {
#
#       shapemodels <- morphogrid(ordination = ordination, axes = args$axes, datype = mspace$ordination$datype,
#                                 template = args$template, x = NULL, y = y, p = mspace$plotinfo$p,
#                                 k = mspace$plotinfo$k, nh = args$nh, nv = args$nv, mag = args$mag,
#                                 asp = args$asp, xlim = args$xlim, ylim = args$ylim, rot.models = args$rot.models,
#                                 size.models = args$size.models, asp.models = args$asp.models)
#
#       plot_morphogrid2d(x = x, y = y, morphogrid = shapemodels, template = args$template,
#                         links = args$links, datype = mspace$ordination$datype, ordtype = mspace$ordination$ordtype,
#                         axes = args$axes, adj_frame = args$adj_frame, p = mspace$plotinfo$p,
#                         xlab = args$xlab, ylab = args$ylab, cex.ldm = args$cex.ldm, col.ldm = args$col.ldm,
#                         col.models = args$col.models, lwd.models = args$lwd.models,
#                         bg.models = args$bg.models, models = args$models)
#     }
#
#     #add points, hulls, phylogeny, and/or consensus
#     if(points & !is.null(mspace$projected$scores)) {
#       if(args$density.groups) density.points <- FALSE else density.points <- TRUE
#       if(ncol(mspace$projected$scores) > 1) {
#         if(length(args$axes) > 1) {
#           if(any(args$pch.points %in% c(21:25))) {
#             graphics::points(mspace$projected$scores[, args$axes], pch = args$pch.points,
#                              bg = args$bg.points, cex = args$cex.points)
#           } else {
#             graphics::points(mspace$projected$scores[, args$axes], pch = args$pch.points,
#                              col = args$col.points, cex = args$cex.points)
#           }
#         } else {
#           if(density.points) {
#             dens <- stats::density(mspace$projected$scores[,args$axes[1]])
#             graphics::polygon(dens$x, dens$y / max(dens$y), lwd = 2,
#                               col = grDevices::adjustcolor(1, alpha.f = 0.5))
#           }
#           graphics::abline(h = 0)
#           if(any(args$pch.points %in% c(21:25))) {
#             graphics::points(cbind(mspace$projected$scores[, args$axes[1]], 0), pch = args$pch.points,
#                              bg = args$bg.points, cex = args$cex.points)
#           } else {
#             graphics::points(cbind(mspace$projected$scores[, args$axes[1]], 0), pch = args$pch.points,
#                              col = args$col.points, cex = args$cex.points)
#           }
#         }
#       } else {
#         if(density.points) {
#           dens <- stats::density(mspace$projected$scores[, 1])
#           graphics::polygon(dens$x, dens$y/max(dens$y), lwd = 2,
#                             col = grDevices::adjustcolor(1, alpha.f = 0.5))
#         }
#         graphics::abline(h = 0)
#         if(any(args$pch.points %in% c(21:25))) {
#           graphics::points(cbind(mspace$projected$scores[, 1], 0), pch = args$pch.points,
#                            bg = args$bg.points, cex = args$cex.points)
#         } else {
#           graphics::points(cbind(mspace$projected$scores[, 1], 0), pch = args$pch.points,
#                            col = args$col.points, cex = args$cex.points)
#         }
#       }
#     }
#
#     if(groups & !is.null(mspace$projected$gr_class)) {
#       if(ncol(mspace$projected$scores) > 1) {
#         if(length(args$axes) > 1) {
#           if(!args$ellipse.groups) {
#             hulls_by_group_2D(mspace$projected$scores[, args$axes], fac = mspace$projected$gr_class,
#                               col = args$col.groups, lty = args$lty.groups, lwd = args$lwd.groups,
#                               alpha = args$alpha.groups)
#           } else {
#             ellipses_by_group_2D(mspace$projected$scores[, args$axes], fac = mspace$projected$gr_class,
#                                  conflev = args$conflev.groups, col = args$col.groups,
#                                  lty = args$lty.groups, lwd = args$lwd.groups,
#                                  alpha = args$alpha.groups)
#           }
#         } else {
#           if(args$density.groups) {
#             dens <- lapply(seq_len(nlevels(mspace$projected$gr_class)), function(i) {
#               subdens <- stats::density(mspace$projected$scores[mspace$projected$gr_class == levels(mspace$projected$gr_class)[i],
#                                                  args$axes[1]])
#               list(x = subdens$x, y = subdens$y)
#             })
#             ymax <- max(unlist(lapply(dens, function(x) {x$y})))
#
#             graphics::abline(h = 0)
#             for(i in seq_len(nlevels(mspace$projected$gr_class))) {
#               graphics::polygon(dens[[i]]$x, dens[[i]]$y / ymax, lwd = 2,
#                                 col = grDevices::adjustcolor(i, alpha.f = alpha.groups))
#             }
#           }
#         }
#
#       } else {
#         if(args$density.groups) {
#           dens <- lapply(seq_len(nlevels(mspace$projected$gr_class)), function(i) {
#             subdens <- stats::density(mspace$projected$scores[mspace$projected$gr_class == levels(mspace$projected$gr_class)[i], 1])
#             list(x = subdens$x, y = subdens$y)
#           })
#           ymax <- max(unlist(lapply(dens, function(x) {x$y})))
#
#           graphics::abline(h=0)
#           for(i in seq_len(nlevels(mspace$projected$gr_class))) {
#             graphics::polygon(dens[[i]]$x, dens[[i]]$y/ymax, lwd = 2,
#                               col = grDevices::adjustcolor(i, alpha.f = alpha.groups))
#           }
#         }
#       }
#     }
#
#     if(phylo & !is.null(mspace$projected$phylo)) {
#       if(length(args$axes) > 1) {
#         for(i in seq_len(nrow(mspace$projected$phylo$edge))) {
#           graphics::lines(rbind(mspace$projected$phylo_scores[mspace$projected$phylo$edge[i, 1], args$axes],
#                                 mspace$projected$phylo_scores[mspace$projected$phylo$edge[i, 2], args$axes]),
#                           lwd = args$lwd.phylo, lty = args$lty.phylo, col = args$col.phylo)
#         }
#       } else {
#         cat("\nphylogenetic relationships are omitted from univariate morphospaces")
#       }
#     }
#
#     if(mshapes & !is.null(mspace$projected$gr_centroids)) {
#       if(ncol(mspace$projected$scores) > 1) {
#         if(length(args$axes) > 1) {
#           if(any(args$pch.groups %in% c(21:25))) {
#             graphics::points(mspace$projected$gr_centroids[, args$axes], bg = args$bg.groups,
#                              pch = args$pch.groups, cex = args$cex.groups)
#           } else {
#             graphics::points(mspace$projected$gr_centroids[, args$axes], col = args$col.groups,
#                              pch = args$pch.groups, cex = args$cex.groups)
#           }
#         } else {
#           graphics::abline(h = 0)
#           if(any(args$pch.groups %in% c(21:25))) {
#             graphics::points(cbind(mspace$projected$gr_centroids[, args$axes[1]], 0), bg = args$bg.groups,
#                              pch = args$pch.groups, cex = args$cex.groups)
#           } else {
#             graphics::points(cbind(mspace$projected$gr_centroids[, args$axes[1]], 0), col = args$col.groups,
#                              pch = args$pch.groups, cex = args$cex.groups)
#           }
#         }
#       } else {
#         graphics::abline(h = 0)
#         if(any(args$pch.groups %in% c(21:25))) {
#           graphics::points(cbind(mspace$projected$gr_centroids[, 1], 0), bg = args$bg.groups,
#                            pch = args$pch.groups, cex = args$cex.groups)
#         } else {
#           graphics::points(cbind(mspace$projected$gr_centroids[, 1], 0), col = args$col.groups,
#                            pch = args$pch.groups, cex = args$cex.groups)
#         }
#       }
#     }
#
#     if(shapeax & !is.null(mspace$projected$shape_axis)) {
#       if(length(args$axes) > 1) {
#
#         if(length(args$lwd.axis) != length(mspace$projected$shape_axis)) {
#           args$lwd.axis <- rep(1, length(mspace$projected$shape_axis))
#         }
#         if(length(args$lty.axis) != length(mspace$projected$shape_axis)) {
#           args$lty.axis <- rep(1, length(mspace$projected$shape_axis))
#         }
#         if(length(args$col.axis) != length(mspace$projected$shape_axis)) {
#           args$col.axis <- rep(1, length(mspace$projected$shape_axis))
#         }
#
#         for(i in seq_len(length(mspace$projected$shape_axis))) {
#           graphics::lines(mspace$projected$shape_axis[[i]][, args$axes],
#                           col = args$col.axis[i], lwd = args$lwd.axis[i], lty = args$lty.axis[i])
#         }
#       } else {
#         cat("\nmorphometric axes are excluded from univariate morphospaces")
#       }
#     }
#
#     if(landsc & !is.null(mspace$projected$landsc)) {
#       if(all(args$axes %in% mspace$plotinfo$axes)) { #all specified vectors must be in the original mspace object
#
#         if(length(args$axes) > 1) { # if number of specified axes is 2...
#           if(all(args$axes == mspace$plotinfo$axes)) { # ...and axes are in the same order
#             cols.landsc <- args$palette.landsc(n = args$nlevels.landsc)
#             levs.landsc <- seq(from = min(mspace$projected$landsc$z, na.rm = TRUE),
#                                to = max(mspace$projected$landsc$z, na.rm = TRUE),
#                                length.out = args$nlevels.landsc)
#
#             if(args$display.landsc == "contour") {
#               graphics::contour(mspace$projected$landsc$x, mspace$projected$landsc$y,
#                                 mspace$projected$landsc$z, col = cols.landsc, levels = levs.landsc,
#                                 lwd = args$lwd.landsc, lty = args$lty.landsc,
#                                 drawlabels = args$drawlabels.landsc, labels = round(levs.landsc, digits = 3),
#                                 add = TRUE)
#               box()
#             }
#
#             if(args$display.landsc == "filled.contour") {
#               if(mspace$projected$landsc$type == "empirical") {
#                 graphics::image(mspace$projected$landsc$x, mspace$projected$landsc$y,
#                                 mspace$projected$landsc$z, add = TRUE,
#                                 col = grDevices::adjustcolor(cols.landsc, alpha = args$alpha.landsc))
#               }
#               if(mspace$projected$landsc$type == "theoretical") {
#                 graphics::.filled.contour(mspace$projected$landsc$x, mspace$projected$landsc$y,
#                                           mspace$projected$landsc$z,
#                                           levels = levs.landsc,
#                                           col = grDevices::adjustcolor(cols.landsc,
#                                                                        alpha = args$alpha.landsc))
#               }
#
#             }
#
#
#           } else { # ...but axes are NOT in the same order
#             cat("\naxes used to generate landscape are in the wrong order; won't be regenerated")
#           }
#
#         } else { # if the number of specified axes is 1....
#
#           if(length(mspace$plotinfo$axes) > 1) { # ....and the number of original axes is 2
#             if(args$axes == mspace$plotinfo$axes[1]) {
#               landsc.z <- mspace$projected$landsc$z[, mspace$projected$landsc$which.y0]
#               landsc.x <- mspace$projected$landsc$x
#             } else {
#               landsc.z <- mspace$projected$landsc$z[mspace$projected$landsc$which.x0,]
#               landsc.x <- mspace$projected$landsc$y
#             }
#           } else { # ...and the number of original axes is 1
#             landsc.z <- mspace$projected$landsc$z
#             landsc.x <- mspace$projected$landsc$x
#           }
#
#           cols.landsc <- args$palette.landsc(n = args$resolution.landsc)[order(order(landsc.z))]
#
#           if(args$drawlabels.landsc == TRUE) {
#             w.transp <-round(quantile(x = 1:length(landsc.z), probs = c(0.25, 0.5, 0.75)))
#             w.transp <- sort(c(w.transp - 1, w.transp, w.transp + 1))
#
#             label_cols <- cols.landsc[w.transp[c(2,5,8)]]
#             cols.landsc[w.transp] <- NA
#           }
#
#           for(i in 1:(length(landsc.x) - 1)) {
#             lines(rbind(c(landsc.x[i], landsc.z[i]),
#                         c(landsc.x[i + 1], landsc.z[i + 1])),
#                   col = cols.landsc[i], lwd = args$lwd.landsc)
#           }
#           box()
#
#           if(args$drawlabels.landsc == TRUE) {
#             x_text <- colMeans(matrix(landsc.x[w.transp + 1], nrow = 3))
#             y_text <- colMeans(matrix(landsc.z[w.transp + 1], nrow = 3))
#             labels <- round(colMeans(matrix(landsc.z[w.transp], nrow = 3)), 2)
#             text(cbind(x_text, y_text), labels = labels, cex = 0.7, col = label_cols)
#           }
#
#         }
#       } else {
#         cat("\nlandscape was originally generated for a different set of axes; won't be regenerated")
#       }
#     }
#
#     ##############################################################################################
#
#   } else { #if either x or y have been provided, show hybrid morphospace
#
#     if(is.null(axes)) {
#       args$axes <- rep(args$axes[1], 2)
#     }
#
#     #if x/y is a phy object, prepare the ground for a phenogram
#     if(any(any(class(x) == "phylo"), any(class(y) == "phylo"))) {
#       phenogr <- TRUE
#     } else {
#       phenogr <- FALSE
#     }
#
#     if(!is.null(x)) {
#       if(any(class(x) == "phylo")) {
#         tree <- x
#         heights <- phytools::nodeHeights(tree)
#         x <- c(rep(max(heights), length(tree$tip.label)),
#                unique(heights[, 1]))
#         args$xlim <- range(x)
#       }
#     } else {
#       if(any(class(y) == "phylo")) {
#         tree <- y
#         heights <- phytools::nodeHeights(tree)
#         y <- c(rep(max(heights), length(tree$tip.label)),
#                unique(heights[, 1]))
#         args$ylim <- range(y)
#       }
#     }
#
#
#
#     if(mspace$plotinfo$k == 3 & mspace$ordination$datype == "landm") {
#
#       shapemodels <- morphogrid(ordination = ordination, axes = args$axes, datype = mspace$ordination$datype,
#                                 template = NULL, x = x, y = y, p = mspace$plotinfo$p,
#                                 k = mspace$plotinfo$k, nh = args$nh, nv = args$nv, mag = args$mag,
#                                 asp = args$asp, xlim = args$xlim, ylim = args$ylim, rot.models = args$rot.models,
#                                 size.models = args$size.models, asp.models = args$asp.models)
#
#       refshape <- matrix(rev_eigen(0,
#                                    ordination$rotation[, 1],
#                                    ordination$center),
#                          nrow = mspace$plotinfo$p, ncol = mspace$plotinfo$k, byrow = TRUE)
#
#       if(is.null(x)) xlim <- range(ordination$x[,args$axes[1]])
#       if(is.null(y)) ylim <- range(ordination$x[,args$axes[1]])
#
#       plot_morphogrid3d(x = x, y = y, morphogrid = shapemodels, refshape = refshape,
#                         template = args$template, links = args$links, ordtype = mspace$ordination$ordtype,
#                         axes = args$axes, xlim = xlim, ylim = ylim, adj_frame = args$adj_frame,
#                         xlab = args$xlab, ylab = args$ylab, cex.ldm = args$cex.ldm,
#                         col.ldm = args$col.ldm, col.models = args$col.models, lwd.models = args$lwd.models,
#                         bg.models = args$bg.models,  size.models = args$size.models,
#                         asp.models = args$asp.models, alpha.models = args$alpha.models, models = args$models)
#     } else {
#
#       shapemodels <- morphogrid(ordination = ordination, axes = args$axes, datype = mspace$ordination$datype,
#                                 template = args$template, x = x, y = y, p = mspace$plotinfo$p,
#                                 k = mspace$plotinfo$k, nh = args$nh, nv = args$nv, mag = args$mag,
#                                 asp = args$asp, xlim = args$xlim, ylim = args$ylim, rot.models = args$rot.models,
#                                 size.models = args$size.models, asp.models = args$asp.models)
#
#       plot_morphogrid2d(x = x, y = y, morphogrid = shapemodels, template = args$template,
#                         links = args$links, datype = mspace$ordination$datype, ordtype = mspace$ordination$ordtype,
#                         axes = args$axes, adj_frame = args$adj_frame, p = mspace$plotinfo$p,
#                         xlab = args$xlab, ylab = args$ylab, cex.ldm = args$cex.ldm, col.ldm = args$col.ldm,
#                         col.models = args$col.models, lwd.models = args$lwd.models,
#                         bg.models = args$bg.models, models = args$models)
#     }
#
#
#     if(phenogr == TRUE) { #if x/y is a phy object, plot a phenogram
#
#       if(length(args$cex.groups) == 1) args$cex.groups <- rep(args$cex.groups, nlevels(mspace$projected$gr_class))
#       if(length(args$pch.groups) == 1) args$pch.groups <- rep(args$pch.groups, nlevels(mspace$projected$gr_class))
#       if(length(args$col.groups) == 1) args$col.groups <- rep(args$col.groups, nrow(mspace$projected$gr_centroids))
#
#       maptips <-  order(match(rownames(mspace$projected$gr_centroids), tree$tip.label))
#       plot_phenogram(x = x, y = y, tree = tree, axis = args$axes, points = points,
#                      pch.groups = args$pch.groups[maptips], phylo_scores = mspace$projected$phylo_scores,
#                      cex.groups = args$cex.groups[maptips], col.groups = args$col.groups[maptips],
#                      bg.groups = args$bg.groups[maptips], lwd.phylo = args$lwd.phylo,
#                      lty.phylo = args$lty.phylo, col.phylo = args$col.phylo)
#
#     } else { #else, go for a generic hybrid morphospace
#
#       xy <- cbind(x, mspace$projected$scores[,args$axes[1]], y)
#
#       #add points, hulls, phylogeny, and/or consensus
#       if(points) {
#         if(any(args$pch.points %in% c(21:25))) {
#           graphics::points(xy, pch = args$pch.points, bg = args$bg.points,
#                            cex = args$cex.points)
#         } else {
#           graphics::points(xy, pch = args$pch.points, col = args$col.points,
#                            cex = args$cex.points)
#         }
#       }
#
#       if(groups & !is.null(mspace$projected$gr_class)) {
#         if(!ellipse.groups) {
#           hulls_by_group_2D(xy, fac = mspace$projected$gr_class, col = args$col.groups,
#                             lty = args$lty.groups, lwd = args$lwd.groups, alpha = args$alpha.groups)
#         } else {
#           ellipses_by_group_2D(xy, fac = mspace$projected$gr_class, conflev = args$conflev.groups,
#                                col = args$col.groups, lty = args$lty.groups,
#                                lwd = args$lwd.groups, alpha = args$alpha.groups)
#         }
#       }
#
#       if(phylo & !is.null(mspace$projected$phylo)) {
#         if(is.null(mspace$projected$gr_class)) {
#           stop("groups classification has not been added to mspace object")
#         } else {
#           if(!is.null(x)) {
#             meanx <- tapply(x, mspace$projected$gr_class, mean)
#           } else {
#             meanx <- NULL
#           }
#           if(!is.null(y)) {
#             meany <- tapply(y, mspace$projected$gr_class, mean)
#           } else {
#             meany <- NULL
#           }
#
#           meanxy <- cbind(meanx, mspace$projected$gr_centroids[,args$axes[1]], meany)[mspace$projected$phylo$tip.label,]
#           nodesxy <- apply(meanxy[mspace$projected$phylo$tip.label,], 2, phytools::fastAnc, tree = mspace$projected$phylo)
#           phyloxy <- rbind(meanxy, nodesxy)
#
#         }
#
#         for(i in seq_len(nrow(mspace$projected$phylo$edge))) {
#           graphics::lines(rbind(phyloxy[mspace$projected$phylo$edge[i, 1],],
#                                 phyloxy[mspace$projected$phylo$edge[i, 2],]),
#                           lwd = args$lwd.phylo, lty = args$lty.phylo, col = args$col.phylo)
#         }
#       }
#
#       if(mshapes & !is.null(mspace$projected$gr_centroids)) {
#         if(!is.null(x)) {
#           meanx <- tapply(x, mspace$projected$gr_class, mean)
#         } else {
#           meanx <- NULL
#         }
#         if(!is.null(y)) {
#           meany <- tapply(y, mspace$projected$gr_class, mean)
#         } else {
#           meany <- NULL
#         }
#
#         meanxy <- cbind(meanx, mspace$projected$gr_centroids[,args$axes[1]], meany)
#
#         if(any(args$pch.groups %in% c(21:25))) {
#           graphics::points(meanxy, bg = args$bg.groups, pch = args$pch.groups,
#                            cex = args$cex.groups)
#         } else {
#           graphics::points(meanxy, col = args$col.groups, pch = args$pch.groups,
#                            cex = args$cex.groups)
#         }
#       }
#     }
#   }
#
#   if(any(legend, scalebar)) {
#
#     if(legend) {
#       if(is.null(mspace$projected$gr_class)) {
#         stop("Groups ($gr_class) are necessary to generate legend labels")
#       } else {
#
#         if(is.null(args$pch.groups)) args$pch.groups <- args$pch.points
#
#         graphics::par(mar = c(5.1, 1, 4.1, 1))
#         plot(0, type = "n", axes = FALSE, xlab = "", ylab = "",
#              ylim = c(0,1), xlim = c(0,1), xaxs = "i", yaxs = "i")
#
#         if(any(args$pch.groups %in% c(21:25))) {
#           graphics::legend("topleft", legend = levels(mspace$projected$gr_class),
#                            cex = cex.legend, pch = args$pch.groups,
#                            pt.bg = args$bg.groups)
#         } else {
#           graphics::legend("topleft", legend = levels(mspace$projected$gr_class),
#                            cex = cex.legend, pch = args$pch.groups,
#                            col = args$col.groups)
#         }
#       }
#     }
#
#     if(scalebar) {
#       if(is.null(mspace$projected$landsc)) {
#         stop("A landscape ($landsc) is necessary to generate a scalebar")
#       } else {
#
#         graphics::par(mar = c(5.1, 1, 4.1, 10))
#         plot(0, type = "n", axes = FALSE, xlab = "", ylab = "",
#              ylim = range(levs.landsc), xlim = c(0,1), xaxs = "i", yaxs = "i")
#
#         if(args$display.landsc == "filled.contour") {
#           alpha.landsc <- args$alpha.landsc
#         } else {
#           alpha.landsc <- 1
#         }
#
#         rect(0, levs.landsc[-length(levs.landsc)], 1, levs.landsc[-1L],
#              col = adjustcolor(cols.landsc, alpha = alpha.landsc), border = "#00000000")
#         ticks <- quantile(levs.landsc, probs = seq(0, 0.95, length.out = 4))
#         axis(side = 4, at = round(ticks, 1))
#         box()
#         mtext("Landscape", side = 4, line = 3)
#
#       }
#     }
#
#     graphics::layout(matrix(1, nrow=1))
#
#   }
#
# }


################################################################################

#' Plot morphospaces and combine them with other variables
#'
#' @description Flexible representation of morphospaces, including their
#'   combination with either other variables or a phylogeny.
#'
#' @param mspace An \code{"mspace"} object created using the
#'   \code{\link{mspace()}} %>% \code{proj_*} pipeline.
#' @param axes Numeric of length 1 or 2, indicating the axes to be plotted.
#' @param links A list with the indices of the coordinates defining the
#'   wireframe (following the format used in \code{Morpho}).
#' @param template Either a 2-column matrix with landmarks/semilandmarks and
#'    template curves coordinates (for 2D shape data) or a \code{"mesh3d"}
#'    object representing the mean shape of the sample (for 3D shape data).
#' @param x,y Optional vector with a non-morphometric variable to be plotted in
#'   the x or y axis. Alternatively, a \code{"phylo"} object can be provided.
#' @param nh,nv Numeric; the number of shape models along the x (\code{nh}) and
#'   the y (\code{nv}) axes.
#' @param mag Numeric; magnifying factor for shape models.
#' @param rescale Logical; whether to re-scale background shape models so shape
#'    variation is shown more clearly.
#' @param invax Optional numeric indicating which of the axes provided in
#'   \code{axes} needs to be inverted (options are \code{1}, \code{2} or
#'   \code{c(1,2)}).
#' @param adj_frame Numeric of length 2, providing \emph{a posteriori} scaling
#'   factors for the width and height of the frame, respectively.
#' @param _.models,_.ldm Control how background shape models are displayed.
#'   \itemize{
#'   \item{models} { Logical; whether to plot background shape models (stored
#'      in \code{mspace$projected$shapemodels}).}
#'   \item{size.models} { Numeric; size factor for shape models.}
#'   \item{col.models} { Color for wireframes/outlines of shape models.}
#'   \item{bg.models} { Background color for outlines/meshes of shape models.}
#'   \item{asp.models} { Numeric; y/x aspect ratio of shape models.}
#'   \item{rot.models} { Numeric; angle (in degrees) to rotate shape models.}
#'   \item{lwd.models} { Numeric; width of border lines in wireframes/outlines
#'      of shape models.}
#'   \item{alpha.models} { Numeric; transparency factor for background models
#'      (3D only).}
#'   \item{cex.ldm} { Numeric; size of landmarks/semilandmarks in the background
#'      models.}
#'   \item{col.ldm} { Color of landmarks/semilandmarks in the background
#'      models.}
#'   }
#' @param _.points Control how scatterpoints representing projected shapes are
#'   represented.
#'   \itemize{
#'   \item{points} { Logical; whether to plot the scatter points corresponding
#'      to the sampled shapes stored in \code{mspace$projected$scores}.}
#'   \item{pch.points} { Numeric; symbol of the scatterpoints.}
#'   \item{col.points} { Color of the scatterpoints.}
#'   \item{bg.points} { Background color of the scatterpoints.}
#'   \item{cex.points} { Numeric; size of the scatterpoints}
#'   \item{density.points} { Logical; whether to add density distribution for
#'      points (univariate ordinations only). Overriden by
#'      \code{density.groups = TRUE}}
#'   }
#' @param _.groups,legend Control how groups of shapes are represented.
#'   \itemize{
#'   \item{groups} { Logical; whether to plot the convex hulls/confidence
#'      ellipses enclosing the groups stored in
#'      \code{mspace$projected$gr_class}.}
#'   \item{mshapes} { Logical; whether to plot the scatter points corresponding
#'      to groups' mean shapes stored in \code{mspace$projected$gr_centroids}.}
#'   \item{col.groups} { Color of the hulls/ellipses and/or scatterpoints
#'      corresponding to groups' mean shapes.}
#'   \item{bg.groups} { Background color of the scatterpoints corresponding to
#'      groups' mean shapes.}
#'   \item{pch.groups} { Numeric; symbol of the scatterpoints corresponding to
#'      groups' mean shapes.}
#'   \item{cex.groups} { Numeric; size of the scatterpoints corresponding to
#'      groups' mean shapes.}
#'   \item{ellipse.groups} { Logical; whether to use confidence ellipses to
#'      delimit groups (if \code{FALSE} convex hulls are used instead).}
#'   \item{conflev.groups} { Numeric; confidence level used for confidence
#'      ellipse(s).}
#'   \item{lwd.groups} { Numeric; width of the lines in groups' ellipses/hulls.}
#'   \item{lty.groups} { Numeric; type of the lines in groups' ellipses/hulls.}
#'   \item{alpha.groups} { Numeric; transparency factor for groups'
#'      ellipses/hulls/density distributions.}
#'   \item{density.groups} { Logical; whether to add density distribution for
#'      groups (univariate ordinations only).}
#'   \item{legend} { Logical; whether to show legend for groups.}
#'   \item{cex.legend} { Numeric; size of legend labels/symbols.}
#'   }
#' @param _.phylo,_.nodes,_.tips Control how phylogenetic relationships are
#'   represented.
#'   \itemize{
#'   \item{phylo} { Logical; whether to plot phylogenetic relationships stored
#'      in \code{mspace$projected$phylo}.}
#'   \item{col.phylo} { Color of the lines depicting phylogenetic
#'      branches.}
#'   \item{lwd.phylo} { Numeric; width of the lines depicting phylogenetic
#'      branches.}
#'   \item{lty.phylo} { Numeric; type of the lines depicting phylogenetic
#'      branches.}
#'   \item{col.nodes} { Color of the scatterpoints representing phylogenetic
#'      nodes.}
#'   \item{bg.nodes} { Background color of the scatterpoints representing
#'      phylogenetic nodes.}
#'   \item{pch.nodes} { Numeric; symbol of the scatterpoints representing
#'      phylogenetic nodes.}
#'   \item{cex.nodes} { Numeric; size of the scatterpoints representing
#'      phylogenetic nodes.}
#'   \item{col.tips} { Color of the scatterpoints representing phylogenetic
#'      tips.}
#'   \item{bg.tips} { Background color of the scatterpoints representing
#'      phylogenetic tips.}
#'   \item{pch.tips} { Numeric; symbol of the scatterpoints representing
#'      phylogenetic tips.}
#'   \item{cex.tips} { Numeric; size of the scatterpoints representing
#'      phylogenetic tips.}
#'   }
#' @param _.axis Control how groups of shapes are represented.
#'   \itemize{
#'   \item{shapeax} { Logical; whether to plot morphometric axes stored in
#'      \code{mspace$projected$shape_axis}.}
#'   \item{type.axis} { Integer; type of arrows ([0] = no arrow; [1] = pointing
#'   towards the maximum; [2] = pointing towards the maximum, [3] = pointing in
#'   both directions).}
#'   \item{col.axis} { Color of the lines depicting a morphometric axis.}
#'   \item{lwd.axis} { Numeric; width of the lines depicting a morphometric
#'      axis.}
#'   \item{lty.axis} { Numeric; type of the lines depicting  a morphometric
#'      axis.}
#'   }
#' @param _.landsc,scalebar Control how landscape surfaces are represented.
#'   \itemize{
#'   \item{landsc} { Logical; whether to plot landscape surface stored in
#'      \code{mspace$projected$landsc}.}
#'   \item{display.landsc} { How to display landscape representation; options
#'      are \code{"contour"} and \code{"filled.contour"}. For bivariate
#'      landscapes only.}
#'   \item{nlevels.landsc} { Number of levels (i.e., contours) to use in
#'      landscape representation.}
#'   \item{palette.landsc} { A function defining a color palette to use for
#'      landscape representation.}
#'   \item{alpha.landsc} { Numeric; transparency factor for filled contours
#'      depicting landscapes.}
#'   \item{lwd.landsc} { Numeric; width of the contour lines depicting
#'      landscapes.}
#'   \item{lty.landsc} { Numeric; type of the contour lines depicting
#'      landscapes.}
#'   \item{drawlabels.landsc} { Logical; should the labels indicating the value
#'      of each surface contour be plotted?}
#'   \item{scalebar} { Logical; whether to show scalebar for landscapes.}
#'   }
#' @param xlim,ylim,xlab,ylab,asp Standard arguments passed to the generic
#'   [graphics::plot()] function.
#'
#' @details This function allows to regenerate/tweak morphospaces contained in
#'   \code{"mspace"} objects already in existence. By default, [plot_mspace]
#'   regenerates the morphospace with all its projected elements, preserving
#'   graphical parameters used originally during the \code{\link{mspace}} +
#'   \code{proj_*} pipeline (and stored in \code{mspace$plotinfo}). However,
#'   all the graphical parameters can be modified to customize graphical
#'   representations. Also, [plot_mspace] can be used to add a legend and/or a
#'   scalebar to aid identification of groups and interpretation of landscapes,
#'   respectively.
#'
#'   In addition, this function expands the range of graphical options available
#'   beyond 'pure' morphospaces. If a numeric non-shape variable (assumed to be
#'   measured for the same specimens in \code{mspace$projected$scores}) is fed
#'   to one of \code{x} or \code{y}, a 'hybrid' morphospace is produced (i.e.
#'   the bivariate plot will be constructed from the combination of \code{x} or
#'   \code{y} and a morphometric axis; background shape models will represent
#'   variation only for the latter). If instead a \code{"phylo"} object (assumed
#'   to describe the phylogenetic relationships among tips scores stored in
#'   \code{mspace$projected$phylo_scores}) is provided for either \code{x} or
#'   \code{y}, a vertical or horizontal phenogram will be deployed (the x/y axis
#'   range will correspond to branch lengths, so caution should be exercised
#'   when interpreting the output).
#'
#'   \emph{Note}: when regenerating landscapes, it's important to keep in mind
#'   that this surface has been calculated for a specific set of shapes and
#'   ordination axes. Hence, its regeneration using [plot_mspace] is only
#'   warranted if the axes being depicted are the same than those used when the
#'   surface landscape was originally computed using \code{\link{mspace}} +
#'   \code{\link{proj_landscape}}. The only exception is when one of the
#'   original axes (i.e., those specified with the \code{axes} argument in
#'   \code{\link{mspace}}) is dropped (i.e., not specified with the \code{axes}
#'   argument of \code{\link{plot_mspace}}). This will result in the collapse of
#'   the 3D landscape projected into a bivariate morphospace into a 2D landscape
#'   projected into a univariate one.
#'
#' @export
#'
#' @examples
#' #load and extract relevant data, packages and information
#' library(magrittr)
#' data("tails")
#' shapes <- tails$shapes
#' species <- tails$data$species
#' sizes <- tails$sizes
#' sp_shapes <- expected_shapes(shapes, species)
#' sp_sizes <- cbind(tapply(sizes, species, mean))
#' tree <- tails$tree
#' links <- tails$links
#'
#' #generate basic morphospace, add sampled shapes, species classification, and
#' #phylogenetic structure
#' msp <- mspace(shapes, axes = c(1,2), plot = FALSE) %>%
#'  proj_shapes(shapes = shapes) %>%
#'  proj_groups(groups = species) %>%
#'  proj_phylogeny(shapes = sp_shapes, tree = tree)
#'
#' ##Plotting 'pure' morphospaces:
#'
#' #plot mspace object as it is
#' plot_mspace(msp)
#'
#' #add colors for points, by species
#' plot_mspace(msp, col.points = species,
#'             col.groups = 1:nlevels(species))
#'
#' #add links for landmark configurations
#' plot_mspace(msp, links = links,
#'             col.points = species, col.groups = 1:nlevels(species))
#'
#' #change number and sizes of shape models in the background
#' plot_mspace(msp, nh = 2, nv = 2, links = links,
#'             size.models = 0.5,
#'             col.points = species, col.groups = 1:nlevels(species))
#'
#' #magnify deformation and highlight landmarks
#' plot_mspace(msp, mag = 1.5, nh = 2, nv = 2, links = links,
#'             size.models = 0.5, col.points = species,
#'             col.groups = 1:nlevels(species), cex.ldm = 5, col.ldm = "red")
#'
#' #change axes 1,2 for 1,3
#' plot_mspace(msp, axes = c(1,3), mag = 1.5, nh = 2, nv = 2, links = links,
#'             size.models = 0.5, col.points = species,
#'             col.groups = 1:nlevels(species), cex.ldm = 5, col.ldm = "red")
#'
#' #change colors for as tips and nodes of the phylogeny
#' plot_mspace(msp, axes = c(1,3), mag = 1.5, nh = 2, nv = 2, links = links,
#'             size.models = 0.5, cex.ldm = 5, col.ldm = "red",
#'             col.tips = "red", col.nodes = "blue")
#'
#' #plot only first PC axis, with general distribution of specimens
#' plot_mspace(msp, axes = 1, mag = 1.5, nh = 2, nv = 2, links = links,
#'             size.models = 0.5, cex.ldm = 5, density.groups = FALSE,
#'             col.ldm = "red", col.tips = "red", col.nodes = "blue")
#'
#' #plot only first PC axis, with distribution of specimens by group
#' plot_mspace(msp, axes = 1, mag = 1.5, nh = 2, nv = 2, links = links,
#'             alpha.groups = 0.5, size.models = 0.5, cex.ldm = 5,
#'             density.groups = TRUE, col.ldm = "red", col.tips = "red",
#'             col.nodes = "blue")
#'
#' #add legend
#' plot_mspace(msp, axes = 1, mag = 1.5, nh = 2, nv = 2, links = links,
#'             alpha.groups = 0.5, size.models = 0.5, cex.ldm = 5,
#'             density.groups = TRUE, col.ldm = "red", col.tips = "red",
#'             col.nodes = "blue", legend = TRUE, cex.legend = 1.3,
#'             pch.groups = 16)
#'
#'
#' #remove all the other elements to see only shape variation captured by the
#' #first two PCs
#' plot_mspace(msp, axes = c(1,2), links = links,
#'             col.points = species, col.groups = 1:nlevels(species),
#'             points = FALSE, groups = FALSE, phylo = FALSE)
#'
#'
#' ##Plotting 'hybrid' morphospaces:
#'
#' #plot size against first PC
#' plot_mspace(msp, x = sizes,  axes = 1, links = links, col.points = species,
#'             col.groups = 1:nlevels(species), pch.points = 16,
#'             xlab = "Centroid size")
#'
#'
#' ##Plotting phenograms:
#'
#' #plot vertical phenogram against PC1
#' plot_mspace(msp, y = tree,  axes = 1, links = links, col.points = species,
#'             col.groups = 1:nlevels(species), pch.points = 16,
#'             xlab = "Branch lengths", pch.tips = 21, bg.tips = "red",
#'             pch.nodes = 21, bg.nodes = "blue")
#'
#' #plot horizontal phenogram against PC2
#' plot_mspace(msp, x = tree,  axes = 2, links = links, col.points = species,
#'             col.groups = 1:nlevels(species), pch.points = 16,
#'             ylab = "Branch lengths", pch.tips = 21, bg.tips = "red",
#'             pch.nodes = 21, bg.nodes = "blue")
plot_mspace <- function(mspace,
                        axes = NULL,
                        links = NULL,
                        template = NULL,
                        x = NULL,
                        y = NULL,
                        nh,
                        nv,
                        mag,
                        rescale,
                        invax = NULL,
                        adj_frame = c(1,1),
                        points = TRUE,
                        models = TRUE,
                        groups = TRUE,
                        phylo = TRUE,
                        shapeax = TRUE,
                        landsc = TRUE,
                        legend = FALSE,
                        scalebar = FALSE,
                        cex.legend = 1,
                        asp = NA,
                        xlab,
                        ylab,
                        xlim = NULL,
                        ylim = NULL,
                        size.models,
                        asp.models,
                        rot.models,
                        col.models,
                        bg.models,
                        lwd.models,
                        alpha.models,
                        cex.ldm,
                        col.ldm,
                        pch.points = 1,
                        col.points = 1,
                        bg.points = 1,
                        cex.points = 1,
                        density.points = TRUE,
                        col.groups = 1,
                        bg.groups = 1,
                        pch.groups = 16,
                        cex.groups = 1,
                        ellipse.groups = mspace$plotinfo$ellipse.groups,
                        conflev.groups = 0.95,
                        lwd.groups = 1,
                        lty.groups = 1,
                        alpha.groups = 0,
                        density.groups = TRUE,
                        col.phylo = 1,
                        lwd.phylo = 1,
                        lty.phylo = 1,
                        col.nodes,
                        pch.nodes,
                        bg.nodes,
                        cex.nodes,
                        col.tips,
                        pch.tips,
                        bg.tips,
                        cex.tips,
                        type.axis,
                        lwd.axis = 1,
                        lty.axis = 1,
                        col.axis = 1,
                        palette.landsc = heat.colors,
                        display.landsc,
                        alpha.landsc,
                        nlevels.landsc = 50,
                        drawlabels.landsc = FALSE,
                        lty.landsc = 1,
                        lwd.landsc = 1) {


  supplied <- names(as.list(match.call()))[-1]
  new_args <- lapply(seq_len(length(supplied)), function(i) {get(supplied[i])})
  names(new_args) <- supplied
  inh_args <- mspace$plotinfo
  merged_args <- utils::modifyList(inh_args, new_args)
  args <- merged_args


  if(!is.null(invax)) {
    mspace$ordination$x[,args$axes[invax]] <- mspace$ordination$x[,args$axes[invax]]  * -1
    mspace$ordination$rotation[,args$axes[invax]] <- mspace$ordination$rotation[,args$axes[invax]] * -1
    if(!is.null(mspace$projected$scores)) {
      mspace$projected$scores[,args$axes[invax]] <- mspace$projected$scores[,args$axes[invax]] * -1
    }
    if(!is.null(mspace$projected$gr_scores)) {
      mspace$projected$gr_scores[,args$axes[invax]] <- mspace$projected$gr_scores[,args$axes[invax]] * -1
    }
    if(!is.null(mspace$projected$phylo_scores)) {
      mspace$projected$phylo_scores[,args$axes[invax]] <- mspace$projected$phylo_scores[,args$axes[invax]] * -1
    }
  }

  ordination <- mspace$ordination

  if(!is.null(x) & !is.null(y)) stop("Only one of x or y can be included")

  if(any(legend, scalebar)) {

    layout_mat <- matrix(c(rep(1, 4), 2), nrow = 1)
    if(all(legend, scalebar)) layout_mat <- rbind(layout_mat, c(rep(1, 4), 3))

    orig_par <- graphics::par(names(graphics::par())[-c(13, 19, 21:23, 54)])
    on.exit(graphics::par(orig_par))

    graphics::layout(layout_mat)
    graphics::par(mar = c(5.1, 4.1, 4.1, 2.1))

  }

  #if neither x nor y were provided, plot pure morphospace #######################################
  if(is.null(x) & is.null(y)) {

    if(ncol(ordination$x) == 1 | length(args$axes) == 1) {
      y <- rep(1, nrow(ordination$x))
      args$ylim <- c(0, 1)
      args$ylab <- "relative density"
    } else {
      y <- NULL
    }

    if(mspace$plotinfo$k == 3 & mspace$ordination$datype == "landm") {

      shapemodels <- morphogrid(ordination = ordination, axes = args$axes, rescale = args$rescale,
                                datype = mspace$ordination$datype, template = NULL, x = NULL,
                                y = y, p = mspace$plotinfo$p, k = mspace$plotinfo$k, nh = args$nh,
                                nv = args$nv, mag = args$mag, asp = args$asp, xlim = args$xlim,
                                ylim = args$ylim, rot.models = args$rot.models,
                                size.models = args$size.models, asp.models = args$asp.models)

      refshape <- matrix(rev_eigen(0,
                                   ordination$rotation[, 1],
                                   ordination$center),
                         nrow = mspace$plotinfo$p, ncol = mspace$plotinfo$k, byrow = TRUE)

      args$xlim <- range(ordination$x[,args$axes[1]])
      if(ncol(ordination$x) > 1 | length(args$axes) > 1) args$ylim <- range(ordination$x[, args$axes[2]])

      rgl::par3d(userMatrix = args$rotation_matrix)

      plot_morphogrid3d(x = NULL, y = y, morphogrid = shapemodels, refshape = refshape,
                        template = args$template, links = args$links,
                        ordtype = mspace$ordination$ordtype, axes = args$axes, xlim = args$xlim,
                        ylim = args$ylim, adj_frame = args$adj_frame, xlab = args$xlab,
                        ylab = args$ylab, cex.ldm = args$cex.ldm, col.ldm = args$col.ldm,
                        col.models = args$col.models, lwd.models = args$lwd.models,
                        bg.models = args$bg.models,  size.models = args$size.models,
                        asp.models = args$asp.models, alpha.models = args$alpha.models,
                        models = args$models, rotmat = args$rotation_matrix)
    } else {

      shapemodels <- morphogrid(ordination = ordination, axes = args$axes, rescale = args$rescale,
                                datype = mspace$ordination$datype, template = args$template,
                                x = NULL, y = y, p = mspace$plotinfo$p, k = mspace$plotinfo$k,
                                nh = args$nh, nv = args$nv, mag = args$mag, asp = args$asp,
                                xlim = args$xlim, ylim = args$ylim, rot.models = args$rot.models,
                                size.models = args$size.models, asp.models = args$asp.models)

      plot_morphogrid2d(x = x, y = y, morphogrid = shapemodels, template = args$template,
                        links = args$links, datype = mspace$ordination$datype,
                        ordtype = mspace$ordination$ordtype, axes = args$axes,
                        adj_frame = args$adj_frame, p = mspace$plotinfo$p,
                        xlab = args$xlab, ylab = args$ylab, cex.ldm = args$cex.ldm,
                        col.ldm = args$col.ldm, col.models = args$col.models,
                        lwd.models = args$lwd.models, bg.models = args$bg.models,
                        models = args$models, )
    }

    #add scatterpoints and hulls/ellipses (intercalated if there are groups present)
    if(points | groups) {

      scores <- mspace$projected$scores

      if(!is.null(mspace$projected$gr_class)){
        gr_class <- mspace$projected$gr_class
        gr_scores <- mspace$projected$gr_scores

        if(!is.null(scores)) {
          index_sc_in_gr <- as.numeric(unlist
                                       (apply(gr_scores, 1, \(x, y) {
                                         which(apply(y, 1, \(z, x) {
                                           all(z == x)}, x))}, scores)))
          if(length(index_sc_in_gr) == 0) index_sc_in_gr <- 1:nrow(scores)
        } else index_sc_in_gr <- 1:nrow(scores <- gr_scores)

      } else {
        index_sc_in_gr <- if(!is.null(scores)) -c(1:nrow(scores)) else NULL
        gr_class <- NULL
      }

      if(length(args$col.groups) == 1) {
        gr_col <- rep(args$col.groups, nlevels(gr_class))
      } else gr_col <- args$col.groups
      if(length(args$lty.groups) == 1) {
        gr_lty <- rep(args$lty.groups, nlevels(gr_class))
      } else gr_lty <- args$lty.groups

      if(length(args$col.points) == 1) {
        gr_col.points <- rep(args$col.points, nrow(mspace$projected$scores))
      } else gr_col.points <- args$col.points
      if(length(args$pch.points) == 1) {
        gr_pch.points <- rep(args$pch.points, nrow(mspace$projected$scores))
      } else gr_pch.points <- args$pch.points
      if(length(args$bg.points) == 1) {
        gr_bg.points <- rep(args$bg.points, nrow(mspace$projected$scores))
      } else gr_bg.points <- args$bg.points
      if(length(args$cex.points) == 1) {
        gr_cex.points <- rep(args$cex.points, nrow(mspace$projected$scores))
      } else gr_cex.points <- args$cex.points


      if(ncol(mspace$ordination$x) > 1) {
        if(length(args$axes) > 1) {
          if(!is.null(scores)) {
            plot_biv_scatter(scores = scores[-index_sc_in_gr,],
                             col = gr_col.points[-index_sc_in_gr],
                             pch = gr_pch.points[-index_sc_in_gr],
                             bg = gr_bg.points[-index_sc_in_gr],
                             cex = gr_cex.points[-index_sc_in_gr])
          }

          if(!is.null(gr_class)) {
            for(i in seq_len(nlevels(gr_class))) {
              plot_biv_scatter(scores = scores[index_sc_in_gr[gr_class == levels(gr_class)[i]], args$axes],
                               col = gr_col.points[index_sc_in_gr][gr_class == levels(gr_class)[i]],
                               pch = gr_pch.points[index_sc_in_gr][gr_class == levels(gr_class)[i]],
                               bg = gr_bg.points[index_sc_in_gr][gr_class == levels(gr_class)[i]],
                               cex = gr_cex.points[index_sc_in_gr][gr_class == levels(gr_class)[i]])
              if(!args$ellipse.groups) {
                hulls_by_group_2D(xy = gr_scores[gr_class == levels(gr_class)[i],],
                                  fac = gr_class[gr_class == levels(gr_class)[i]],
                                  col = gr_col[i], lty = gr_lty[i], lwd = args$lwd.groups,
                                  alpha = args$alpha.groups)
              } else {
                ellipses_by_group_2D(xy = gr_scores[gr_class == levels(gr_class)[i],],
                                     fac = factor(gr_class[gr_class == levels(gr_class)[i]]),
                                     col = gr_col[i], lty = gr_lty[i], lwd = args$lwd.groups,
                                     alpha = args$alpha.groups, conflev = args$conflev.groups)
              }
            }
          }
        } else {
          if(!is.null(gr_class)) {
            if(args$density.groups) {
              density_by_group_2D(gr_scores, gr_class, ax = args$axes[1], alpha = args$alpha.groups,
                                  lwd = args$lwd.groups, lty = args$lty.groups, col = args$col.groups)
            }
            plot_univ_scatter(scores = cbind(scores[index_sc_in_gr, args$axes[1]]), density = FALSE,
                              col = gr_col.points[index_sc_in_gr], pch = gr_pch.points[index_sc_in_gr],
                              cex = gr_cex.points[index_sc_in_gr], bg = gr_bg.points[index_sc_in_gr])
          }
          suppressWarnings(
            plot_univ_scatter(scores = cbind(scores[-index_sc_in_gr, args$axes[1]]), density = FALSE,
                              col = gr_col.points[-index_sc_in_gr], pch = gr_pch.points[-index_sc_in_gr],
                              cex = gr_cex.points[-index_sc_in_gr], bg = gr_bg.points[-index_sc_in_gr])
          )
        }
      } else {
        if(!is.null(gr_class)) {
          if(args$density.groups) {
            density_by_group_2D(gr_scores, gr_class, ax = 1, alpha = args$alpha.groups,
                                lwd = args$lwd.groups, lty = args$lty.groups, col = args$col.groups)
          }
          plot_univ_scatter(scores = cbind(scores[index_sc_in_gr,]), density = FALSE, col = gr_col.points[index_sc_in_gr],
                            pch = gr_pch.points[index_sc_in_gr], cex = gr_cex.points[index_sc_in_gr],
                            bg = gr_bg.points[index_sc_in_gr])
        }
        suppressWarnings(
          plot_univ_scatter(scores = cbind(scores[-index_sc_in_gr,]), density = FALSE,
                            col = gr_col.points[-index_sc_in_gr], pch = gr_pch.points[-index_sc_in_gr],
                            cex = gr_cex.points[-index_sc_in_gr], bg = gr_bg.points[-index_sc_in_gr])
        )
      }
    }

    #add phylogeny
    if(phylo & !is.null(mspace$projected$phylo)) {

      if(length(args$axes) > 1) {
        for(i in seq_len(nrow(mspace$projected$phylo$edge))) {
          graphics::lines(rbind(mspace$projected$phylo_scores[mspace$projected$phylo$edge[i, 1], args$axes],
                                mspace$projected$phylo_scores[mspace$projected$phylo$edge[i, 2], args$axes]),
                          lwd = args$lwd.phylo, lty = args$lty.phylo, col = args$col.phylo)
        }

        tips <- seq_len(length(mspace$projected$phylo$tip.label))
        plot_biv_scatter(scores = mspace$projected$phylo_scores[-tips, args$axes], pch = args$pch.nodes,
                         bg = args$bg.nodes, col = args$col.nodes, cex = args$cex.nodes)
        plot_biv_scatter(scores = mspace$projected$phylo_scores[tips, args$axes], pch = args$pch.tips,
                         bg = args$bg.tips, col = args$col.tips, cex = args$cex.tips)

      } else {
        warning("phylogenetic relationships are not projected into univariate morphospaces")
      }
    }

    #add shape axes
    if(shapeax & !is.null(mspace$projected$shape_axis)) {

      if(length(args$axes) > 1) {

        if(length(args$type.axis) != length(mspace$projected$shape_axis)) {
          args$type.axis <- rep(args$type.axis, length(mspace$projected$shape_axis))
        }
        if(length(args$lwd.axis) != length(mspace$projected$shape_axis)) {
          args$lwd.axis <- rep(args$lwd.axis, length(mspace$projected$shape_axis))
        }
        if(length(args$lty.axis) != length(mspace$projected$shape_axis)) {
          args$lty.axis <- rep(args$lty.axis, length(mspace$projected$shape_axis))
        }
        if(length(args$col.axis) != length(mspace$projected$shape_axis)) {
          args$col.axis <- rep(args$col.axis, length(mspace$projected$shape_axis))
        }

        for(i in seq_len(length(mspace$projected$shape_axis))) {
          graphics::arrows(x0 = mspace$projected$shape_axis[[i]][1, args$axes[1]],
                           y0 = mspace$projected$shape_axis[[i]][1, args$axes[2]],
                           x1 = mspace$projected$shape_axis[[i]][2, args$axes[1]],
                           y1 = mspace$projected$shape_axis[[i]][2, args$axes[2]],
                           code = args$type.axis[i], length = 0.1, col = args$col.axis[i],
                           lwd = args$lwd.axis[i], lty = args$lty.axis[i])
        }
      } else {
        warning("morphometric axes are not projected into univariate morphospaces")
      }
    }

    #add landscape
    if(landsc & !is.null(mspace$projected$landsc)) {

      if(all(args$axes %in% mspace$plotinfo$axes)) { #all specified vectors must be in the original mspace object

        if(length(args$axes) > 1) { # if number of specified axes is 2...
          if(all(args$axes == mspace$plotinfo$axes)) { # ...and axes are in the same order
            cols.landsc <- args$palette.landsc(n = args$nlevels.landsc)

            levs.landsc <- seq(from = min(mspace$projected$landsc$z, na.rm = TRUE),
                               to = max(mspace$projected$landsc$z, na.rm = TRUE),
                               length.out = args$nlevels.landsc)

            plot_biv_landscape(landscape = mspace$projected$landsc, display = args$display.landsc,
                               levels = levs.landsc, col = cols.landsc, lwd = args$lwd.landsc,
                               lty = args$lty.landsc, drawlabels = args$drawlabels.landsc,
                               alpha = args$alpha.landsc, type = mspace$projected$landsc$type)


          } else { # ...but axes are NOT in the same order
            warning("axes used to generate landscape are in the wrong order; won't be regenerated")
          }

        } else { # if the number of specified axes is 1....

          if(length(mspace$plotinfo$axes) > 1) { # ....and the number of original axes is 2
            if(args$axes == mspace$plotinfo$axes[1]) {
              landsc.z <- mspace$projected$landsc$z[, mspace$projected$landsc$which.y0]
              landsc.x <- mspace$projected$landsc$x
            } else {
              landsc.z <- mspace$projected$landsc$z[mspace$projected$landsc$which.x0,]
              landsc.x <- mspace$projected$landsc$y
            }
          } else { # ...and the number of original axes is 1
            landsc.z <- mspace$projected$landsc$z
            landsc.x <- mspace$projected$landsc$x
          }

          cols.landsc <- args$palette.landsc(n = length(landsc.z))[order(order(landsc.z))]
          levs.landsc <- mspace$projected$landsc$z

          plot_univ_landscape(landscape = list(x = landsc.x, z = landsc.z),
                              drawlabels = args$drawlabels.landsc, col = cols.landsc,
                              lwd = args$lwd.landsc)
        }
      } else {
        warning("landscape was originally generated for a different set of axes; won't be regenerated")
      }
    }

  } else {
    #if any of x or y have been provided, deploy hybrid morphospace ################################

    if(is.null(axes)) {
      args$axes <- rep(args$axes[1], 2)
    }

    #if x or y are a phy object, prepare the ground for a phenogram --------------------------------
    if(any(any(class(x) == "phylo"), any(class(y) == "phylo"))) {
      phenogr <- TRUE
    } else {
      phenogr <- FALSE
    }

    if(!is.null(x)) {
      if(any(class(x) == "phylo")) {
        tree <- x
        heights <- phytools::nodeHeights(tree)
        x <- c(rep(max(heights), length(tree$tip.label)),
               unique(heights[, 1]))
        args$xlim <- range(x)
      }
    } else {
      if(any(class(y) == "phylo")) {
        tree <- y
        heights <- phytools::nodeHeights(tree)
        y <- c(rep(max(heights), length(tree$tip.label)),
               unique(heights[, 1]))
        args$ylim <- range(y)
      }
    }

    if(mspace$plotinfo$k == 3 & mspace$ordination$datype == "landm") {

      shapemodels <- morphogrid(ordination = ordination, axes = args$axes, datype = mspace$ordination$datype,
                                template = NULL, x = x, y = y, p = mspace$plotinfo$p,
                                k = mspace$plotinfo$k, nh = args$nh, nv = args$nv, mag = args$mag,
                                asp = args$asp, xlim = args$xlim, ylim = args$ylim, rot.models = args$rot.models,
                                size.models = args$size.models, asp.models = args$asp.models, rescale = args$rescale)

      refshape <- matrix(rev_eigen(0,
                                   ordination$rotation[, 1],
                                   ordination$center),
                         nrow = mspace$plotinfo$p, ncol = mspace$plotinfo$k, byrow = TRUE)

      if(is.null(x)) xlim <- range(ordination$x[,args$axes[1]])
      if(is.null(y)) ylim <- range(ordination$x[,args$axes[1]])

      plot_morphogrid3d(x = x, y = y, morphogrid = shapemodels, refshape = refshape,
                        template = args$template, links = args$links, ordtype = mspace$ordination$ordtype,
                        axes = args$axes, xlim = xlim, ylim = ylim, adj_frame = args$adj_frame,
                        xlab = args$xlab, ylab = args$ylab, cex.ldm = args$cex.ldm, col.ldm = args$col.ldm,
                        col.models = args$col.models, lwd.models = args$lwd.models, bg.models = args$bg.models,
                        size.models = args$size.models, asp.models = args$asp.models, alpha.models = args$alpha.models,
                        models = args$models, rotmat = args$rotation_matrix)
    } else {

      shapemodels <- morphogrid(ordination = ordination, axes = args$axes, datype = mspace$ordination$datype,
                                template = args$template, x = x, y = y, p = mspace$plotinfo$p,
                                k = mspace$plotinfo$k, nh = args$nh, nv = args$nv, mag = args$mag,
                                asp = args$asp, xlim = args$xlim, ylim = args$ylim, rot.models = args$rot.models,
                                size.models = args$size.models, asp.models = args$asp.models, rescale = args$rescale)

      plot_morphogrid2d(x = x, y = y, morphogrid = shapemodels, template = args$template,
                        links = args$links, datype = mspace$ordination$datype, ordtype = mspace$ordination$ordtype,
                        axes = args$axes, adj_frame = args$adj_frame, p = mspace$plotinfo$p,
                        xlab = args$xlab, ylab = args$ylab, cex.ldm = args$cex.ldm, col.ldm = args$col.ldm,
                        col.models = args$col.models, lwd.models = args$lwd.models,
                        bg.models = args$bg.models, models = args$models)
    }


    if(phenogr) { #if x/y is a phy object, plot a phenogram

      tips <- seq_len(length(tree$tip.label))
      if(length(args$col.tips) == 1) args$col.tips <- rep(args$col.tips, length(tips))
      if(length(args$bg.tips) == 1) args$bg.tips <- rep(args$bg.tips, length(tips))
      if(length(args$pch.tips) == 1) args$pch.tips <- rep(args$pch.tips, length(tips))
      if(length(args$cex.tips) == 1) args$cex.tips <- rep(args$cex.tips, length(tips))

      if(length(args$col.nodes) == 1) args$col.nodes <- rep(args$col.nodes, length(tips) - 1)
      if(length(args$bg.nodes) == 1) args$bg.nodes <- rep(args$bg.nodes, length(tips) - 1)
      if(length(args$pch.nodes) == 1) args$pch.nodes <- rep(args$pch.nodes, length(tips) - 1)
      if(length(args$cex.nodes) == 1) args$cex.nodes <- rep(args$cex.nodes, length(tips) - 1)

      maptips <- order(match(rownames(mspace$projected$phylo_scores[tips,]), tree$tip.label))

      plot_phenogram(x = x, y = y, tree = tree, phylo_scores = mspace$projected$phylo_scores,
                     axis = args$axes, points = points, pch.tips = args$pch.tips[maptips],
                     cex.tips = args$cex.tips[maptips], col.tips = args$col.tips[maptips],
                     bg.tips = args$bg.tips[maptips], pch.nodes = args$pch.tips[maptips],
                     cex.nodes = args$cex.nodes[maptips], col.nodes = args$col.nodes[maptips],
                     bg.nodes = args$bg.nodes[maptips], lwd.phylo = args$lwd.phylo,
                     lty.phylo = args$lty.phylo, col.phylo = args$col.phylo)

    } else {
      #if x or y are regular variables, go for a generic hybrid morphospace -----------------------

      if(points | groups) {

        scores <- mspace$projected$scores

        if(!is.null(mspace$projected$gr_class)){
          gr_class <- mspace$projected$gr_class
          gr_scores <- mspace$projected$gr_scores

          if(!is.null(scores)) {
            index_sc_in_gr <- as.numeric(unlist
                                         (apply(gr_scores, 1, \(x, y) {
                                           which(apply(y, 1, \(z, x) {
                                             all(z == x)}, x))}, scores)))
            if(length(index_sc_in_gr) == 0) index_sc_in_gr <- 1:nrow(scores)
          } else index_sc_in_gr <- 1:nrow(scores <- gr_scores)

        } else {
          index_sc_in_gr <- if(!is.null(scores)) -c(1:nrow(scores)) else NULL
          gr_class <- NULL
        }

        xy <- cbind(x, scores[,args$axes[1]], y)

        if(length(args$col.groups) == 1) {
          gr_col <- rep(args$col.groups, nlevels(gr_class))
        } else gr_col <- args$col.groups
        if(length(args$lty.groups) == 1) {
          gr_lty <- rep(args$lty.groups, nlevels(gr_class))
        } else gr_lty <- args$lty.groups

        if(length(args$col.points) == 1) {
          gr_col.points <- rep(args$col.points, nrow(mspace$projected$scores))
        } else gr_col.points <- args$col.points
        if(length(args$pch.points) == 1) {
          gr_pch.points <- rep(args$pch.points, nrow(mspace$projected$scores))
        } else gr_pch.points <- args$pch.points
        if(length(args$bg.points) == 1) {
          gr_bg.points <- rep(args$bg.points, nrow(mspace$projected$scores))
        } else gr_bg.points <- args$bg.points
        if(length(args$cex.points) == 1) {
          gr_cex.points <- rep(args$cex.points, nrow(mspace$projected$scores))
        } else gr_cex.points <- args$cex.points


        if(!is.null(scores)) {
          plot_biv_scatter(scores = xy[-index_sc_in_gr,],
                           col = gr_col.points[-index_sc_in_gr],
                           pch = gr_pch.points[-index_sc_in_gr],
                           bg = gr_bg.points[-index_sc_in_gr],
                           cex = gr_cex.points[-index_sc_in_gr])
        }

        if(!is.null(gr_class)) {
          for(i in seq_len(nlevels(gr_class))) {
            plot_biv_scatter(scores = xy[index_sc_in_gr[gr_class == levels(gr_class)[i]], ],
                             col = gr_col.points[index_sc_in_gr][gr_class == levels(gr_class)[i]],
                             pch = gr_pch.points[index_sc_in_gr][gr_class == levels(gr_class)[i]],
                             bg = gr_bg.points[index_sc_in_gr][gr_class == levels(gr_class)[i]],
                             cex = gr_cex.points[index_sc_in_gr][gr_class == levels(gr_class)[i]])
            if(!args$ellipse.groups) {
              hulls_by_group_2D(xy = xy[index_sc_in_gr[gr_class == levels(gr_class)[i]], ],
                                fac = gr_class[gr_class == levels(gr_class)[i]],
                                col = gr_col[i], lty = gr_lty[i], lwd = args$lwd.groups,
                                alpha = args$alpha.groups)
            } else {
              ellipses_by_group_2D(xy = xy[index_sc_in_gr[gr_class == levels(gr_class)[i]], ],
                                   fac = factor(gr_class[gr_class == levels(gr_class)[i]]),
                                   col = gr_col[i], lty = gr_lty[i], lwd = args$lwd.groups,
                                   alpha = args$alpha.groups, conflev = args$conflev.groups)
            }
          }
        }
      }

      #add phylogeny
      if(phylo & !is.null(mspace$projected$phylo)) {

        tree <- mspace$projected$phylo
        phylo_scores <- mspace$projected$phylo_scores

        if(nrow(cbind(x, y)) == nrow(phylo_scores[tree$tip.label, ]))  {

          if(!is.null(x)) {
            x <- cbind(x)
            if(!is.null(rownames(x))) x <- x[tree$tip.label,] else {
              warning("names for x have not been provided; it will be assumed that x is in the same order as labels in phylogenetic tree")
              x <- c(x)
            }
          }

          if(!is.null(y)) {
            y <- cbind(y)
            if(!is.null(rownames(y))) y <- y[tree$tip.label,] else {
              warning("names for y have not been provided; it will be assumed that y is in the same order as labels in phylogenetic tree")
              y <- c(y)
            }
          }

        } else {

          gr_class <- factor(gsub(x = as.character(mspace$projected$gr_class),
                                  pattern = "_bis.", replacement = ""))

          if(!is.null(gr_class)) {
            if(nrow(cbind(x, y)) == length(gr_class)) {
              if(nlevels(gr_class) == length(tree$tip.label) &
                 all(tree$tip.label %in% levels(gr_class))) {

                if(!is.null(x)) x <- tapply(x, gr_class, mean)[tree$tip.label]
                if(!is.null(y)) y <- tapply(y, gr_class, mean)[tree$tip.label]

              }
            }
          }
        }

        if(nrow(cbind(x, y)) == nrow(phylo_scores[tree$tip.label, ]))  {
          xy.tips <- cbind(x, phylo_scores[tree$tip.label, mspace$plotinfo$axes[1]], y)[tree$tip.label,]
          xy.nodes <- apply(xy.tips[tree$tip.label,], 2, phytools::fastAnc, tree = tree)
          phyloxy <- rbind(xy.tips, xy.nodes)

          for(i in seq_len(nrow(tree$edge))) {
            graphics::lines(rbind(phyloxy[tree$edge[i, 1],], phyloxy[tree$edge[i, 2],]),
                            lwd = args$lwd.phylo, lty = args$lty.phylo, col = args$col.phylo)
          }
          tips <- seq_len(length(tree$tip.label))
          plot_biv_scatter(scores = phyloxy[-tips,], pch = args$pch.nodes,
                           bg = args$bg.nodes, col = args$col.nodes, cex = args$cex.nodes)
          plot_biv_scatter(scores = phyloxy[tips,], pch = args$pch.tips,
                           bg = args$bg.tips, col = args$col.tips, cex = args$cex.tips)

        }
      }
    }
  }

  #add references ###########################################################################
  if(any(legend, scalebar)) {

    if(legend) {
      if(is.null(mspace$projected$gr_class)) {
        stop("Groups levels ($gr_class) are necessary to generate legend labels")
      } else {


        #graphics::par(mar = c(5.1, 1, 4.1, 1))
        graphics::par(mar = c(5.1, 0, 4.1, 0))
        plot(0, type = "n", axes = FALSE, xlab = "", ylab = "",
             ylim = c(0,1), xlim = c(0,1), xaxs = "i", yaxs = "i")


        gr_class <- mspace$projected$gr_class
        gr_scores <- mspace$projected$gr_scores
        scores <- mspace$projected$scores


        index_sc_in_gr <- apply(scores, 1, \(x, y) {any(
          apply(y, 1, \(z, x) {all(z == x)}, x))}, gr_scores)
        if(length(index_sc_in_gr) == 0) index_sc_in_gr <- 1:nrow(scores)

        index_gr_in_sc <- as.numeric(unlist
                                     (apply(scores, 1, \(x, y) {
                                       which(apply(y, 1, \(z, x){
                                         all(z == x)}, x))}, gr_scores)))
        if(length(index_gr_in_sc) == 0) index_gr_in_sc <- 1:nrow(gr_scores)

        pt_params <- NULL
        for(i in seq_len(nrow(cbind(scores[index_sc_in_gr,])))) {
          pt_params <- rbind(pt_params, c(
            as.character(gr_class[index_gr_in_sc][i]),
            args$col.points[index_sc_in_gr][i],
            args$bg.points[index_sc_in_gr][i],
            args$pch.points[index_sc_in_gr][i],
            args$cex.points[index_sc_in_gr][i]))
        }
        gr_params <- data.frame(unique(pt_params))
        if(!is.null(pt_params)) colnames(gr_params) <- c("class", "col", "bg", "pch", "cex")
        gr_params <- gr_params[!gr_params$class %in% names(which(table(gr_params$class) > 1)),]
        gr_params <- gr_params[!is.na(gr_params$class),]

        if(any(!levels(gr_class) %in% gr_params$class)){
          extra_params <- data.frame(cbind(levels(gr_class)[!levels(gr_class) %in% gr_params$class],
                                           NA, NA, NA, NA))
          if(!is.null(extra_params)) colnames(extra_params) <- c("class", "col", "bg", "pch", "cex")
        } else extra_params <- NULL

        gr_params <- rbind(gr_params, extra_params)
        rownames(gr_params) <- gr_params$class
        gr_params <- gr_params[levels(gr_class),]

        gr_params$bg <- args$bg.groups
        gr_params$col2 <- args$col.groups
        gr_params$lty <- args$lty.groups
        gr_params$lwd <- args$lwd.groups

        if(any(grepl(gr_params$class, pattern = "_bis."))) {
          index_bis. <- which(grepl(gr_params$class, pattern = "_bis."))

          gr_params[-index_bis.,][which(is.na(gr_params[-index_bis.,]$col)),] <-
            gr_params[index_bis.,][which(is.na(gr_params[-index_bis.,]$col)),]

          gr_params <- gr_params[-index_bis.,]
          gr_params$class <- gsub(x = gr_params$class, pattern = "_bis.", replacement = "")
          gr_class <- factor(gsub(x = as.character(gr_class), pattern = "_bis.", replacement = ""))

        }

        index_keep <- gr_params$class %in% levels(factor(gr_class))
        if(length(index_keep) != 0) gr_params <- gr_params[index_keep,]
        gr_params$lty[gr_params$lty == 0] <- 1


        graphics::legend("topleft", legend = gr_params$class, pt.cex = as.numeric(gr_params$cex),
                         pch = NA, col = gr_params$col2, cex = cex.legend, pt.bg = gr_params$bg,
                         lty = as.numeric(gr_params$lty), lwd = gr_params$lwd)

        graphics::legend("topleft", legend = gr_params$class, pt.cex = as.numeric(gr_params$cex),
                         pch = as.numeric(gr_params$pch), col = gr_params$col, bty = "n", cex = cex.legend,
                         pt.bg = gr_params$bg, lty = NA, lwd = gr_params$lwd, text.font = NA)


      }
    }

    if(scalebar) {
      if(is.null(mspace$projected$landsc)) {
        stop("Landscape values ($landsc) is necessary to generate a scalebar")
      } else {

        #graphics::par(mar = c(5.1, 1, 4.1, 10))
        graphics::par(mar = c(5.1, 4, 4.1, 6))
        plot(0, type = "n", axes = FALSE, xlab = "", ylab = "",
             ylim = range(levs.landsc), xlim = c(0.2,0.8), xaxs = "i", yaxs = "i")

        if(args$display.landsc == "filled.contour") {
          alpha.landsc <- args$alpha.landsc
        } else {
          alpha.landsc <- 1
        }

        rect(0.2, levs.landsc[-length(levs.landsc)], 0.8, levs.landsc[-1L],
             col = adjustcolor(cols.landsc, alpha = alpha.landsc), border = "#00000000")
        ticks <- quantile(levs.landsc, probs = seq(0, 0.95, length.out = 4))
        axis(side = 4, at = round(ticks, 1))
        box()
        mtext("Landscape", side = 4, line = 3)

      }
    }

    graphics::layout(matrix(1, nrow=1))

  }
}
