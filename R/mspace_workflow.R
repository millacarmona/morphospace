###########################################################################

#' Generate morphospace
#'
#' @description Create an empirical morphospace as an ordination comprising a set
#'   of axes synthesizing shape variation. Allows a variety of multivariate methods
#'   for building the ordination.
#'
#' @param shapes Shapes data.
#' @param axes Numeric of length 2, indicating the axes to be plotted.
#' @param links A list with the indices of the coordinates defining the
#'   wireframe (following the format used in \code{Morpho}).
#' @param template Either a 2-column matrix containing landmarks/semilandmarks
#'   followed by coordinates defining a curve or set of curves describing additional
#'   aspects of morphology (for 2D shape data) or a \code{"mesh3d"} object containing
#'   geometry of the structure the landmarks were placed on (for 3D shape data),
#'   corresponding to the mean shape of the sample. These will be warped using TPS
#'   interpolation to produce the set of background shell models (see
#'   \code{\link{build_template2d}}).
#' @param p Numeric, indicating the number of landmarks/semilandmarks used (for
#'   landmark data only).
#' @param k Numeric, indicating the number of cartesian dimensions of
#'   landmarks/semilandmarks (for landmark data only).
#' @param FUN The function to be used for synthesizing geometric morphometric
#'   variation. Usual alternatives include \code{\link[stats]{prcomp}},
#'   \code{\link{bg_prcomp}} and \code{\link{phy_prcomp}}.
#' @param nh Numeric; number of shape models along the x axis.
#' @param nv Numeric; number of shape models along the y axis.
#' @param mag Numeric; magnifying factor for shape models.
#' @param invax Optional numeric indicating which of the axes provided in
#'   \code{axes} needs to be inverted (optionsare \code{1}, \code{2} or
#'   \code{c(1,2)}).
#' @param adj_frame Numeric of length 2, providing \emph{a posteriori} scaling
#'   factors for the width and height of the frame, respectively.
#' @param rot.models  Numeric; angle (in degrees) to rotate shape models.
#' @param size.models Numeric; size factor for shape models.
#' @param asp.models Numeric; the y/x aspect ratio of shape models.
#' @param col.models Color for wireframes/outlines.
#' @param bg.models Background color for outlines/meshes.
#' @param lwd.models Numeric; width of the lines in wireframes/outlines.
#' @param alpha.models Numeric; transparency factor for background models (3D only).
#' @param points Logical; whether to plot the scatter points.
#' @param cex.ldm Numeric; size of landmarks/semilandmarks in the background
#'   models.
#' @param col.ldm Color of landmarks/semilandmarks in the background models.
#' @param plot Logical; whether to plot morphospace.
#' @param models Logical; whether to plot background shape models.
#' @param xlim,ylim,xlab,ylab,asp Standard arguments passed to the generic plot
#'   function.
#' @param ... Further arguments passed to \code{FUN}.
#'
#' @details This function is the central piece of the \code{morphospace} workflow.
#'   It produces a synthetic space from a sample of normalized shapes using
#'   eigenanalysis-based ordination methods (PCA, between groups PCA, two-block
#'   PLS, and their phylogenetic versions) and generates shape models depicting
#'   the range of realized variation, retaining the information necessary to
#'   project and retrieve new compatible shapes. The output of \code{mspace} is
#'   meant to be expanded using the \code{proj_*} family of functions and the
#'   \code{%>%} operator from \code{magrittr}.
#'
#' @return An object of class \code{"mspace"}, which is a list containing:
#'   \itemize{
#'   \item{$ordtype:} { method used for multivariate ordination.}
#'   \item{$datype:} { type of shape data used.}
#'   \item{$x:} { scores of the sample of shapes in the synthetic axes.}
#'   \item{$rotation:} { eigenvector's coefficients.}
#'   \item{$center:} { the mean values of the original shape variables.}
#'   \item{$plotinfo:} { a list with the information used to create the plot.}
#'   }
#'
#' @seealso \code{\link{proj_shapes}}, \code{\link{proj_consensus}},
#'   \code{\link{proj_groups}}, \code{\link{proj_phylogeny}},
#'   \code{\link{proj_axis}}
#'
#' @export
#'
#' @examples
#'
#' ##2D Landmark data
#'
#' #load and extract relevant data and information
#' data("tails")
#' shapes <- tails$shapes
#' species <- tails$data$species
#' links <- tails$links
#' tree <- tails$tree
#'
#' #generate morphospace using the basic sample of shapes, PCA as ordination method
#' #and the links between landmarks provided for backround models
#' mspace(shapes, links = links, mag = 0.7, axes = c(1,2), points = TRUE)
#'
#' #increase magnification factor x2:
#' mspace(shapes, links = links, mag = 1.5, axes = c(1,2), points = TRUE)
#'
#' #plot PCs 1 and 3
#' mspace(shapes, links = links, mag = 1.5, axes = c(1,3), points = TRUE)
#'
#' #generate morphospace using the basic sample of shapes, bgPCA as ordination method
#' #and links between landmarks for backround models
#' mspace(shapes, links = links, FUN = bg_prcomp, groups = species, mag = 0.7,
#'        axes = c(1,2), invax = 1, points = TRUE)
#'
#' #generate morphospace using species' consensus shapes and phylogenetic tree,
#' #phylogenetic PCA as ordination method, and links between landmarks for backround
#' #models
#' sp_shapes <- expected_shapes(shapes, species)
#' mspace(sp_shapes, links = links, FUN = phy_prcomp, tree = tree, mag = 0.7,
#'        axes = c(1,2), points = TRUE)
#'
#' #just create a morphospace without plotting, save into an object, and inspect
#' morphosp <- mspace(shapes, links = links, mag = 0.7, axes = c(1,2), plot = FALSE)
#' names(morphosp)
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
#' data("shells3D")
#' shapes <- shells3D$shapes
#' mesh_meanspec <- shells3D$mesh_meanspec
#'
#' #generate morphospace. This is interactive, you need to rotate the shape by yourself and
#' #then press enter into the console.
#' mspace(shapes, mag = 1, axes = c(1,2), col.ldm = "black", cex.ldm = 2, points = TRUE)
#'
#' #generate morphospace using a mesh template that improves visualization:
#' #first, get shape corresponding to shells3D$mesh_meanspec using geomorph::findMeanSpec
#' meanspec_id<- findMeanSpec(shapes)
#' meanspec_shape <- shapes[,,meanspec_id]
#'
#' #then get the consensus shape and warp the sample mesh to get the mesh corresponding to the
#' #consensus using Morpho::tps3d
#' meanshape <- expected_shapes(shapes)
#' meanmesh <- tps3d(x = mesh_meanspec , refmat = meanspec_shape, tarmat = meanshape)
#'
#' #finally, generate morphospace providing template (this function used the mesh warped to
#' #the mean shape of the entire sample, hence the previous lines)
#' mspace(shapes, mag = 1, axes = c(1,2), template = meanmesh, bg.models = "gray",
#'        nh = 4, nv = 4, cex.ldm = 0, points = TRUE)
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
#'        asp.models = 1, bg.model = "light gray")
mspace <- function(shapes,
                   axes = c(1,2),
                   links = NULL,
                   template = NULL,
                   p = NULL,
                   k = NULL,
                   FUN = stats::prcomp,
                   nh = 5,
                   nv = 4,
                   mag = 1,
                   invax = NULL,
                   asp = NA,
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

  FUN <- match.fun(FUN)
  ordination <- FUN(data2d, ...)
  ordtype <- class(ordination)

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

  if(k == 3 & datype == "landm") {

    shapemodels <- morphogrid(ordination = ordination, axes = axes, datype = datype, template = NULL,
                              x = NULL, y = y, p = p, k = k, nh = nh, nv = nv, mag = mag,
                              asp = asp, xlim = xlim, ylim = ylim, rot.models = rot.models,
                              size.models = size.models, asp.models = 1)

    refshape <- expected_shapes(shapes)
    xlim <- range(ordination$x[,axes[1]])
    if(ncol(ordination$x) > 1 | length(axes) > 1) ylim <- range(ordination$x[,axes[2]])

    plot_morphogrid3d(x = NULL, y = y, morphogrid = shapemodels, refshape = refshape,
                      template = template, links = links, ordtype = ordtype, adj_frame = adj_frame,
                      axes = axes, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab,
                      cex.ldm = cex.ldm, col.ldm = col.ldm, col.models = col.models,
                      lwd.models = lwd.models, bg.models = bg.models, size.models = size.models,
                      asp.models = asp.models, alpha.models = alpha.models, plot = plot, models = models)
  } else {


    shapemodels <- morphogrid(ordination = ordination, axes = axes, datype = datype, template = template,
                              x = NULL, y = y, p = p, k = k, nh = nh, nv = nv, mag = mag,
                              asp = asp, xlim = xlim, ylim = ylim, rot.models = rot.models,
                              size.models = size.models, asp.models = asp.models)

    plot_morphogrid2d(x = NULL, y = y, morphogrid = shapemodels, template = template,
                      links = links, datype = datype, ordtype = ordtype, axes = axes,
                      adj_frame = adj_frame, p = p, xlab = xlab, ylab = ylab, cex.ldm = cex.ldm,
                      col.ldm = col.ldm, col.models = col.models, lwd.models = lwd.models,
                      bg.models = bg.models, plot = plot, models = models)

  }

  if(points == TRUE) graphics::points(ordination$x[,axes])

  plotinfo <- list(p = p, k = k, links = links, template = template, axes = axes, nh = nh, nv = nv, mag = mag,
                   asp = asp, adj_frame = adj_frame, asp.models = asp.models, rot.models = rot.models,
                   size.models = size.models, lwd.models = lwd.models, bg.models = bg.models, col.models = col.models,
                   alpha.models = alpha.models, cex.ldm = cex.ldm, col.ldm = col.ldm, models = models)

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
#' @param mspace An \code{"mspace"} object.
#' @param density Logical; whether to add density distribution for points
#'   (univariate ordinations only).
#' @param pipe Logical; is the function being included in a pipe?
#' @param ... Further arguments passed to [graphics::points()].
#'
#' @details The purpose of this function is to maintain morphospace
#'   generation and sample representation as independent affairs, and
#'   to add flexibility to graphical representation of scatter points.
#'
#' @return If a plot device with a morphospace is open, shapes feeded to
#'   \code{shapes} are projected into morphospace. If \code{pipe = FALSE}
#'   those scores are returned invisibly. If \code{pipe = TRUE} the supplied
#'   \code{mspace} object will be modified by replacing the original
#'   \code{$x} slot as well as adding some graphical parameters, and returned
#'   invisibly.
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

  if(mspace$datype != datype) stop("shapes and mspace types are not compatible")

  scores <- proj_eigen(x = data2d, vectors = mspace$rotation,
                       center = mspace$center)

  if(.Device != "null device") {

    if(ncol(mspace$x) > 1) {
      if(length(mspace$plotinfo$axes) > 1) {
        graphics::points(scores[, mspace$plotinfo$axes], ...)
      } else {
        if(density == TRUE) {
          dens <- stats::density(scores[, mspace$plotinfo$axes[1]])
          graphics::polygon(dens$x, dens$y / max(dens$y), lwd = 2,
                            col = grDevices::adjustcolor(1, alpha.f = 0.5))
        }
        graphics::abline(h = 0)
        graphics::points(cbind(scores[, mspace$plotinfo$axes[1]], 0), ...)
      }
    } else {
      if(density == TRUE) {
        dens <- stats::density(scores[, 1])
        graphics::polygon(dens$x, dens$y / max(dens$y), lwd = 2,
                          col = grDevices::adjustcolor(1, alpha.f = 0.5))
      }
      graphics::abline(h = 0)
      graphics::points(cbind(scores[, 1], 0), ...)
    }

  }

  mspace$x <- scores
  mspace$plotinfo$pch.points <- args$pch
  mspace$plotinfo$col.points <- args$col
  mspace$plotinfo$bg.points <- args$bg
  mspace$plotinfo$cex.points <- args$cex
  mspace$plotinfo$density.points <- args$density


  if(pipe == FALSE) return(invisible(scores))
  if(pipe == TRUE) return(invisible(mspace))

}


###############################################################################################

#' Project consensus shape(s) into morphospace
#'
#' @description Project one or more mean shapes into an existing morphospace.
#'
#' @param shapes Shapes data.
#' @param mspace An \code{"mspace"} object.
#' @param pipe Logical; is the function being included in a pipe?
#' @param ... Further arguments passed to [graphics::points()].
#'
#' @details The purpose of this function is to add the scores corresponding
#'   to groups' mean shapes to \code{mspace} objects during pipeline. Otherwise,
#'   it does the same than \code{proj_shapes}.
#'
#' @return If a plot device with a morphospace is open, shapes feeded to
#'   \code{shapes} are projected into morphospace. If \code{pipe = FALSE} the
#'   corresponding scores are returned invisibly. If \code{pipe = TRUE} the
#'   supplied \code{mspace} object will be modified by adding a new
#'   \code{$gr_centroids} slot as well as a number of graphical parameters,
#'   and returned invisibly.
#'
#' @seealso \code{\link{expected_shapes}}
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
#' #generate basic morphospace, add sampled and consensus shapes
#' mspace(shapes, mag = 0.7, axes = c(1,2), bg.model = "light gray") %>%
#'   proj_shapes(shapes = shapes, col = c(1:4)[species], pch = 1, cex = 0.7) %>%
#'   proj_consensus(shapes = sp_shapes, pch = 21, bg = 1:4, cex = 2)
proj_consensus <- function(mspace, shapes, pipe = TRUE, ...) {

  args <- c(as.list(environment()), list(...))

  dat <- shapes_mat(shapes)
  data2d <- dat$data2d
  datype <- dat$datype

  if(mspace$datype != datype) stop("shapes and mspace types are not compatible")

  gr_centroids <- proj_eigen(x = data2d, vectors = mspace$rotation,
                             center = mspace$center)

  if(.Device != "null device") {

    if(ncol(mspace$x) > 1) {
      if(length(mspace$plotinfo$axes) > 1) {
        graphics::points(gr_centroids[, mspace$plotinfo$axes], ...)
      } else {
        graphics::abline(h = 0)
        graphics::points(cbind(gr_centroids[, mspace$plotinfo$axes[1]], 0), ...)
      }
    } else {
      graphics::abline(h = 0)
      graphics::points(cbind(gr_centroids[, 1], 0), ...)
    }
  }

  if(is.null(args$col)) args$col <- seq_len(nrow(gr_centroids))
  if(is.null(args$pch)) args$pch <- 1


  mspace$gr_centroids <- gr_centroids
  mspace$plotinfo$pch.groups <- args$pch
  mspace$plotinfo$cex.groups <- args$cex

  if(pipe == FALSE) return(invisible(gr_centroids))
  if(pipe == TRUE) return(invisible(mspace))

}



###############################################################################################

#' Delimit groups in morphospace
#'
#' @description Project convex hulls or confidence ellipses enclosing \emph{a priori}
#'   groups into an existing morphospace.
#'
#' @param mspace An \code{"mspace"} object.
#' @param shapes Optional shapes data.
#' @param groups Factor; classification of observations into groups.
#' @param ellipse Logical; whether to plot confidence ellipses (if \code{FALSE},
#'   convex hulls will be used instead).
#' @param conflev Numeric; confidence level used for confidence ellipse(s).
#' @param density Logical; whether to add density distribution for groups
#'   (univariate ordinations only).
#' @param pipe Logical; is the function being included in a pipe?
#' @param ... Further arguments passed to [hulls_by_group_2D()] or
#'   [ellipses_by_group_2D()].
#'
#' @details The purpose of this function is to add a classification
#'   for shapes populating the morphospace to \code{mspace} objects
#'   during pipeline, as well as to facilitate group visualization.
#'   Otherwise, it is just a wrapper for \code{hulls_by_group_2D} and
#'   \code{ellipses_by_group_2D}.
#'
#' @return If a plot device with a morphospace is open, convex hulls or
#'   confidence ellipses enclosing the scores corresponding to \code{groups}
#'   are projected into morphospace. If \code{pipe = TRUE} the supplied
#'   \code{mspace} object will be modified by adding a new \code{$gr_class}
#'   slot as well as a number of graphical parameters, and returned invisibly.
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
#' mspace(shapes, links = links, mag = 0.7, axes = c(1,2), bg.model = "light gray") %>%
#'   proj_shapes(shapes = shapes, col = c(1:4)[species], pch = 1) %>%
#'   proj_groups(groups = species, col = 1:4, lwd = 2)
#'
#' #generate basic morphospace, add sampled shapes and 95% confidence for species
#' mspace(shapes, links = links, mag = 0.7, axes = c(1,2), bg.model = "light gray") %>%
#'   proj_shapes(shapes = shapes, col = c(1:4)[species], pch = 1) %>%
#'   proj_groups(groups = species, col = 1:4, lwd = 2, ellipse = TRUE, conflev = 0.95)
proj_groups <- function(mspace, shapes = NULL, groups, ellipse = FALSE,
                        conflev = 0.95, density = TRUE, pipe = TRUE, ...) {

  args <- c(as.list(environment()), list(...))

  if(is.null(shapes)) {
    data2d <- mspace$x
  } else {
    dat <- shapes_mat(shapes)
    datype <- dat$datype
    data2d <- proj_eigen(x = dat$data2d, vectors = mspace$rotation, center = mspace$center)

    if(mspace$datype != datype) stop("shapes and mspace types are not compatible")

  }

  if(.Device != "null device") {

    if(ncol(mspace$x) > 1) {
      if(length(mspace$plotinfo$axes) > 1) {
        if(ellipse == FALSE) {
          hulls_by_group_2D(data2d[, mspace$plotinfo$axes], fac = groups, ...)
        } else {
          ellipses_by_group_2D(data2d[, mspace$plotinfo$axes], fac = groups, conflev = conflev, ...)
        }
      } else {
        if(density == TRUE) {
          dens <- lapply(seq_len(nlevels(groups)), function(i) {
            subdens <- stats::density(mspace$x[groups == levels(groups)[i], mspace$plotinfo$axes[1]])
            list(x = subdens$x, y = subdens$y)
          })
          ymax <- max(unlist(lapply(dens, function(x) {x$y})))

          graphics::abline(h = 0)
          for(i in seq_len(nlevels(groups))) {
            graphics::polygon(dens[[i]]$x, dens[[i]]$y/ymax, lwd = 2,
                              col = grDevices::adjustcolor(i, alpha.f = 0.5))
          }
        }
      }

    } else {
      if(density == TRUE) {
        dens <- lapply(seq_len(nlevels(groups)), function(i) {
          subdens <- stats::density(mspace$x[groups == levels(groups)[i], 1])
          list(x = subdens$x, y = subdens$y)
        })
        ymax <- max(unlist(lapply(dens, function(x) {x$y})))

        graphics::abline(h = 0)
        for(i in seq_len(nlevels(groups))) {
          graphics::polygon(dens[[i]]$x, dens[[i]]$y / ymax, lwd = 2,
                            col = grDevices::adjustcolor(i, alpha.f = 0.5))
        }
      }
    }
  }

  if(is.null(args$col)) args$col <- seq_len(nlevels(groups))
  if(is.null(args$alpha)) args$alpha <- 0

  mspace$gr_class <- groups
  mspace$plotinfo$ellipse.groups <- args$ellipse
  mspace$plotinfo$conflev.groups <- args$conflev
  mspace$plotinfo$col.groups <- args$col
  mspace$plotinfo$bg.groups <- args$col
  mspace$plotinfo$lty.groups <- args$lty
  mspace$plotinfo$lwd.groups <- args$lwd
  mspace$plotinfo$alpha.groups <- args$alpha
  mspace$plotinfo$density.groups <- args$density


  if(pipe == TRUE) return(invisible(mspace))
  #if(pipe == FALSE) return(  volume_of_hulls  ) #one for the future

}


###############################################################################################

#' Project a morphometric axis into morphospace
#'
#' @description Project one or more morphometric axes (i.e., linear combinations
#'   of shape variables) into an existing bivariate morphospace.
#'
#' @param obj An object containing either a multivariate ordination of class
#'   \code{"prcomp", "bg_prcomp", "phy_prcomp"} or \code{"pls_shape"} or a
#'   \code{"mlm"} object fitted using [stats::lm()].
#' @param mspace An \code{"mspace"} object.
#' @param axis An optional vector indicating the axis from \code{obj} to be
#'   projected.
#' @param mag Numeric; magnifying factor for representing shape transformation.
#' @param pipe Logical; is the function being included in a pipe?
#' @param ... Further arguments passed to [graphics::lines()].
#'
#' @details This function is primarily aimed at the graphical representation
#'   of morphometric axes (estimated using either linear models or multivariate
#'   ordination methods) into an existing morphospace for heuristic exploration
#'   of patterns in the data. It can also be used to extract theoretical shapes
#'   at the extremes of those axes, although [ax_transformation()] does the
#'   same thing in a more flexible and straightforward way.
#'
#'   Axes computed from linear models and from supervised multivariate ordinations
#'   using the same supervising variable will usually differ in length (i.e. shape
#'   transformations will be magnified or attenuated) due to the assumption of the
#'   former that the explanatory variable has been measured without error. Also,
#'   the former will not be necessarily centered.
#'
#'   For statistical analysis of axes (e.g., trajectory analysis) their vector
#'   coefficients can be extracted directly from slope coefficients stored in
#'   \code{"mlm"} objects or eigenvector coefficients stored in the \code{$rotation}
#'   slot returned by multivariate ordination methods.
#'
#' @return If a plot device with a morphospace is open, a straight line marking the
#'   scores representing shapes at the extremes of the morphometric axis is projected
#'   into morphospace. If \code{pipe = FALSE} those scores are returned invisibly.
#'    If \code{pipe = TRUE} the supplied \code{mspace} object will be modified by adding
#'    a new \code{$shape_axis} slot as well as a number of graphical parameters, and
#'    returned invisibly.
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
#' links <- tails$links
#'
#' #compute intraspecific allometric axis detrend_shapes, using lm and pls_shapes
#' detr_shapes <- arrayspecs(
#'   detrend_shapes(lm(two.d.array(shapes) ~ species)),
#'                           p = 9, k = 2)
#' intrasp_allo_mod <- lm(two.d.array(detr_shapes) ~ logsizes)
#' intrasp_allo_pls <- pls_shapes(shapes = two.d.array(detr_shapes), X = logsizes)
#'
#' #compute intraspecific allometric axis using expected_shapes, tapply, lm and
#' #pls_shapes
#' sp_shapes <- expected_shapes(shapes, species)
#' sp_logsizes <- tapply(logsizes, species, max)
#' intersp_allo_mod <- lm(two.d.array(sp_shapes) ~ sp_logsizes)
#' intersp_allo_pls <- pls_shapes(shapes = two.d.array(sp_shapes), X = sp_logsizes)
#'
#'
#' #generate basic morphospace, add intraspecific (red) and interspecific (blue) axes
#' #for lm models
#' mspace(shapes, links = links, mag = 0.7, axes = c(1,2)) %>%
#'   proj_axis(obj = intrasp_allo_mod, col = "red", lwd = 2) %>%
#'   proj_axis(obj = intersp_allo_mod, col = "blue", lwd = 2)
#' #for pls ordination
#' mspace(shapes, links = links, mag = 0.7, axes = c(1,2)) %>%
#'   proj_axis(obj = intrasp_allo_pls, col = "red", lwd = 2) %>%
#'   proj_axis(obj = intersp_allo_pls, col = "blue", lwd = 2)
proj_axis <- function(mspace, obj, axis = 1, mag = 1, pipe = TRUE, ...) {

  args <- c(as.list(environment()), list(...))

  ext_shapes2d <- ax_transformation(obj = obj, axis = axis, mag = mag)
  ext_scores <- proj_eigen(x = ext_shapes2d, vectors = mspace$rotation,
                           center = mspace$center)

  if(.Device != "null device") {
    if(length(mspace$plotinfo$axes) > 1) {
      graphics::lines(ext_scores[, mspace$plotinfo$axes], ...)
    } else {
      print("\nphenotypic change vectors are omitted from univariate morphospaces")
    }
  }

  mspace$shape_axis <- c(mspace$shape_axis, list(ext_scores))
  mspace$plotinfo$lwd.axis <- c(mspace$plotinfo$lwd.axis, args$lwd)
  mspace$plotinfo$lty.axis <- c(mspace$plotinfo$lty.axis, args$lty)
  mspace$plotinfo$col.axis <- c(mspace$plotinfo$col.axis, args$col)

  if(pipe == FALSE) return(invisible(ext_scores))
  if(pipe == TRUE) return(invisible(mspace))

}



###############################################################################################

#' Project phylogenetic structure into morphospace
#'
#' @description Project phylogenetic relationships among a set of shapes
#'   (representing the consensuses of phhylogenetic terminals) into an existing
#'   bivariate morphospace.
#'
#' @param tree A \code{"phylo"} object containing a phylogenetic tree.
#' @param mspace An \code{"mspace"} object.
#' @param pipe Logical; is the function being included in a pipe?
#' @param ... Further arguments passed to [graphics::lines()].
#'
#' @details The purpose of this function is twofold. First, it is meant to transform
#'   a morphospace into a phylomorphospace by infusing  phylogenetic structure into
#'   the former. To this end, a \code{$gr_centroids} slot matching the tip labels
#'   from \code{tree} needs to be present (either upstream in the pipeline or already
#'   incorporated into an existing \code{mspace} object). Second, this function can be
#'   used to retrieve the scores corresponding to nodes of the phylogenetic tree, which
#'   can in turn be used to compute the associated shapes using [rev_eigen()]. The
#'   position of these shapes in morphospace is estimated using the squared-changes
#'   parsimony algorithm as performed by [phytools::fastAnc()].
#'
#' @return If a plot device with a morphospace is open, shapes representing the nodes
#'   of the phylogenetic tree and lines connecting them and tips are projected into
#'   morphospace. If \code{pipe = FALSE} scores for nodes and tips of the phylogeny
#'   are returned invisibly. If \code{pipe = TRUE} the supplied \code{mspace}
#'   object will be modified by adding a new \code{$phylo_scores} and \code{$phylo}
#'   slots as well as a number of graphical parameters, and returned invisibly.
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
#'   proj_shapes(shapes = shapes, col = c(1:13)[species], pch = 1, cex = 0.7) %>%
#'   proj_consensus(shapes = sp_shapes, pch = 21, bg = 1:13, cex = 2) %>%
#'   proj_phylogeny(tree = tree)
proj_phylogeny <- function(mspace, tree, pipe = TRUE, ...) {

  args <- c(as.list(environment()), list(...))

  if(is.null(mspace$gr_centroids)) stop("Group centroids have not been provided; add proj_consensus() before")

  nodes_scores <- apply(mspace$gr_centroids[tree$tip.label,], 2, phytools::fastAnc, tree = tree)
  phylo_scores <- rbind(mspace$gr_centroids[tree$tip.label,], nodes_scores)


  if(.Device != "null device") {

    if(length(mspace$plotinfo$axes) > 1) {
      for(i in seq_len(nrow(tree$edge))) {
        graphics::lines(rbind(phylo_scores[tree$edge[i, 1], mspace$plotinfo$axes],
                              phylo_scores[tree$edge[i, 2], mspace$plotinfo$axes]), ...)
      }
    } else {
      print("phylogenetic relationships are omitted from univariate morphospaces")
    }
  }

  if(is.null(args$col)) args$col <- 1

  mspace$phylo_scores <- phylo_scores
  mspace$phylo <- tree
  mspace$plotinfo$lwd.phylo <- args$lwd
  mspace$plotinfo$lty.phylo <- args$lty
  mspace$plotinfo$col.phylo <- args$col


  if(pipe == FALSE) {
    return(invisible(phylo_scores))
  } else {
    return(invisible(mspace))
  }

}


#########################################################################################

#' Plot morphospaces and combine them with other variables
#'
#' @description Flexible representation of morphospaces, including their combination
#'   with other variables or a phylogeny.
#'
#' @param mspace An \code{"mspace"} object created using [mspace()].
#' @param axes Numeric of length 1 or 2, indicating the axes to be plotted.
#' @param links A list with the indices of the coordinates defining the
#'   wireframe (following the format used in \code{Morpho}).
#' @param template Either a 2-column matrix containing landmarks/semilandmarks
#'   followed by coordinates defining a curve or set of curves describing additional
#'   aspects of morphology (for 2D shape data) or a \code{"mesh3d"} object containing
#'   geometry of the structure the landmarks were placed on (for 3D shape data),
#'   corresponding to the mean shape of the sample. These will be warped using TPS
#'   interpolation to produce the set of background shell models (see
#'   \code{\link{build_template2d}}).
#' @param x Optional vector with a non-morphometric variable to be plotted in
#'   the x axis. Alternatively, a \code{"phylo"} object can be provided.
#' @param y Optional vector with a non-morphometric variable to be plotted in
#'   the y axis. Alternatively, a \code{"phylo"} object can be provided.
#' @param nh Numeric; the number of shape models along the x axis.
#' @param nv Numeric; the number of shape models along the y axis.
#' @param mag Numeric; magnifying factor for shape models.
#' @param invax Optional numeric indicating which of the axes provided in
#'   \code{axes} needs to be inverted (optionsare \code{1}, \code{2} or
#'   \code{c(1,2)}).
#' @param adj_frame Numeric of length 2, providing \emph{a posteriori} scaling
#'   factors for the width and height of the frame, respectively.
#' @param points Logical; whether to plot the scatter points corresponding to
#'   the sampled shapes stored in \code{mspace$x}.
#' @param models Logical; whether to plot background shape models.
#' @param mshapes Logical; whether to plot the scatter points corresponding to
#'   groups' mean shapes stored in \code{mspace$gr_centroids}.
#' @param groups Logical; whether to plot the convex hulls/confidence ellipses
#'   enclosing the groups stored in \code{mspace$gr_class}.
#' @param phylo Logical; whether to plot phylogenetic relationships stored
#'   in \code{mspace$phylo}.
#' @param shapeax Logical; whether to plot morphometric axes stored in
#'   \code{mspace$shape_axis}.
#' @param legend Logical; whether to show legend for groups (\code{mspace$gr_class}).
#' @param cex.legend Numeric; size of legend labels/symbols.
#' @param size.models Numeric; size factor for shape models.
#' @param asp.models Numeric; y/x aspect ratio of shape models.
#' @param rot.models Numeric; angle (in degrees) to rotate shape models.
#' @param col.models Color for wireframes/outlines.
#' @param bg.models Background color for outlines/meshes.
#' @param lwd.models Numeric; width of the lines in wireframes/outlines.
#' @param alpha.models Numeric; transparency factor for background models (3D only).
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
#'   (univariate ordinations only).
#' @param col.groups Color of the hulls/ellipses and/or scatter points corresponding
#'   to groups mean shapes.
#' @param bg.groups Background color of the scatter points corresponding
#'   to groups mean shapes.
#' @param pch.groups Numeric; symbol of the scatter points corresponding to
#'   groups mean shapes.
#' @param cex.groups Numeric; size of the scatter points corresponding to
#'   groups mean shapes.
#' @param ellipse.groups Logical; whether to plot confidence ellipses (if \code{FALSE},
#'   convex hulls will be used instead).
#' @param conflev.groups Numeric; confidence level used for confidence ellipse(s).
#' @param lwd.groups Numeric; width of the lines in groups' ellipses/hulls.
#' @param lty.groups Numeric; type of the lines in groups' ellipses/hulls.
#' @param alpha.groups Numeric; transparency factor for groups' ellipses/hulls/density
#'   distributions.
#' @param density.groups Logical; whether to add density distribution for groups
#'   (univariate ordinations only).
#' @param lwd.phylo Numeric; width of the lines depicting phylogenetic
#'   branches.
#' @param lty.phylo Numeric; type of the lines depicting phylogenetic branches.
#' @param col.phylo Numeric; color of the lines depicting phylogenetic branches.
#' @param lwd.axis Numeric; width of the lines depicting morphometric axis.
#' @param lty.axis Numeric; type of the lines depicting  morphometric axis.
#' @param col.axis Numeric; color of the lines depicting morphometric axis.
#' @param xlim,ylim,xlab,ylab,asp Standard arguments passed to the generic [graphics::plot()]
#'   function.
#'
#' @details This function allows to regenerate morphospaces contained in \code{mspace}
#'   objects already in existence, either preserving the graphical attributes used
#'   during the pipeline or modifying one or more of them (so there is no need to
#'   execute the pipeline every time morphospaces need to be visualized and/or changed).
#'
#'   Also, \code{plot_mspace} expands the range of graphical options available
#'   beyond 'pure' morphospaces. If a numeric non-shape variable (assumed to have
#'   been measured for the same specimens in \code{mspace$x}) is feeded to one of
#'   \code{x} or \code{y}, a 'hybrid' morphospace is produced (i.e. the bivariate
#'   plot will be constructed from the combination of \code{x} or \code{y} and a
#'   morphometric axis; shape models in the background will represent variation
#'   only for the latter). If instead a \code{"phylo"} object (assumed to describe
#'   the phylogenetic relationships among tips scores stored in
#'   \code{mspace$phylo_scores}) is feeded to one of \code{x} or \code{y}, a vertical
#'   or horizontal phenogram will be deployed (the x/y axis range will correspond to
#'   branch lengths, so caution should be exercised when interpreting the output).
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
plot_mspace <- function(mspace,
                        axes = NULL,
                        links = NULL,
                        template = NULL,
                        x = NULL,
                        y = NULL,
                        nh,
                        nv,
                        mag,
                        invax = NULL,
                        adj_frame = c(1, 1),
                        points = TRUE,
                        models = TRUE,
                        mshapes = TRUE,
                        groups = TRUE,
                        phylo = TRUE,
                        shapeax = TRUE,
                        legend = FALSE,
                        cex.legend = 1,
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
                        lwd.phylo = 1,
                        lty.phylo = 1,
                        col.phylo = 1,
                        lwd.axis = 1,
                        lty.axis = 1,
                        col.axis = 1) {


  supplied <- names(as.list(match.call()))[-1]
  new_args <- lapply(seq_len(length(supplied)), function(i) {get(supplied[i])})
  names(new_args) <- supplied
  inh_args <- mspace$plotinfo
  merged_args <- utils::modifyList(inh_args, new_args)
  args <- merged_args


  if(!is.null(invax)) {
    mspace$x[,axes[invax]] <- mspace$x[,axes[invax]]  * -1
    mspace$rotation[,axes[invax]] <- mspace$rotation[,axes[invax]] * -1
    if(!is.null(mspace$gr_centroids)) {
      mspace$gr_centroids[,axes[invax]] <- mspace$gr_centroids[,axes[invax]] * -1
    }
    if(!is.null(mspace$phylo_scores)) {
      mspace$phylo_scores[,axes[invax]] <- mspace$phylo_scores[,axes[invax]] * -1
    }
  }

  ordination <- list(x = mspace$x, rotation = mspace$rotation, center = mspace$center)

  if(!is.null(x) & !is.null(y)) stop("Only one of x or y can be specified")

  if(legend == TRUE) {

    if(is.null(mspace$gr_class)) {
      stop("Groups ($gr_class) are needed to plot legend")

    } else {

      orig_par <- graphics::par(names(graphics::par())[-c(13, 19, 21:23, 54)])
      on.exit(graphics::par(orig_par))

      graphics::layout(matrix(c(2, 2, 2, 2, 1), nrow = 1))
      graphics::par(mar = c(1, 0, 1, 1))
      plot(0, type = "n", axes = FALSE, xlab = "", ylab = "")

      if(any(args$pch.groups %in% c(21:25))) {
        graphics::legend("left", legend = levels(mspace$gr_class),
                         cex = cex.legend, pch = args$pch.groups,
                         pt.bg = args$bg.groups)
      } else {
        graphics::legend("left", legend = levels(mspace$gr_class),
                         cex = cex.legend, pch = args$pch.groups,
                         col = args$col.groups)
      }

      graphics::par(mar = c(5.1, 4.1, 4.1, 2.1))

    }
  }


  if(is.null(x) & is.null(y)) { #if neither x nor y have been provided, plot pure morphospace

    if(ncol(ordination$x) == 1 | length(args$axes) == 1) {
      y <- rep(1, nrow(ordination$x))
      args$ylim <- c(0, 1)
      args$ylab <- "relative density"
    } else {
      y <- NULL
    }

    if(mspace$plotinfo$k == 3 & mspace$datype == "landm") {

      shapemodels <- morphogrid(ordination = ordination, axes = args$axes, datype = mspace$datype,
                                template = NULL, x = NULL, y = y, p = mspace$plotinfo$p,
                                k = mspace$plotinfo$k, nh = args$nh, nv = args$nv, mag = args$mag,
                                asp = args$asp, xlim = args$xlim, ylim = args$ylim, rot.models = args$rot.models,
                                size.models = args$size.models, asp.models = args$asp.models)

      refshape <- matrix(rev_eigen(0,
                                   ordination$rotation[, 1],
                                   ordination$center),
                         nrow = mspace$plotinfo$p, ncol = mspace$plotinfo$k, byrow = TRUE)

      args$xlim <- range(ordination$x[,args$axes[1]])
      if(ncol(ordination$x) > 1 | length(args$axes) > 1) args$ylim <- range(ordination$x[, args$axes[2]])

      plot_morphogrid3d(x = NULL, y = y, morphogrid = shapemodels, refshape = refshape,
                        template = args$template, links = args$links, ordtype = mspace$ordtype,
                        axes = args$axes, xlim = args$xlim, ylim = args$ylim, adj_frame = args$adj_frame,
                        xlab = args$xlab, ylab = args$ylab, cex.ldm = args$cex.ldm,
                        col.ldm = args$col.ldm, col.models = args$col.models, lwd.models = args$lwd.models,
                        bg.models = args$bg.models,  size.models = args$size.models,
                        asp.models = args$asp.models, alpha.models = args$alpha.models, models = args$models)
    } else {

      shapemodels <- morphogrid(ordination = ordination, axes = args$axes, datype = mspace$datype,
                                template = args$template, x = NULL, y = y, p = mspace$plotinfo$p,
                                k = mspace$plotinfo$k, nh = args$nh, nv = args$nv, mag = args$mag,
                                asp = args$asp, xlim = args$xlim, ylim = args$ylim, rot.models = args$rot.models,
                                size.models = args$size.models, asp.models = args$asp.models)

      plot_morphogrid2d(x = x, y = y, morphogrid = shapemodels, template = args$template,
                        links = args$links, datype = mspace$datype, ordtype = mspace$ordtype,
                        axes = args$axes, adj_frame = args$adj_frame, p = mspace$plotinfo$p,
                        xlab = args$xlab, ylab = args$ylab, cex.ldm = args$cex.ldm, col.ldm = args$col.ldm,
                        col.models = args$col.models, lwd.models = args$lwd.models,
                        bg.models = args$bg.models, models = args$models)
    }

    #add points, hulls, phylogeny, and/or consensus
    if(points == TRUE) {
      if(ncol(mspace$x) > 1) {
        if(length(args$axes) > 1) {
          if(any(args$pch.points %in% c(21:25))) {
            graphics::points(mspace$x[, args$axes], pch = args$pch.points,
                             bg = args$bg.points, cex = args$cex.points)
          } else {
            graphics::points(mspace$x[, args$axes], pch = args$pch.points,
                             col = args$col.points, cex = args$cex.points)
          }
        } else {
          if(args$density.points == TRUE) {
            dens <- stats::density(mspace$x[,args$axes[1]])
            graphics::polygon(dens$x, dens$y / max(dens$y), lwd = 2,
                              col = grDevices::adjustcolor(1, alpha.f = 0.5))
          }
          graphics::abline(h = 0)
          if(any(args$pch.points %in% c(21:25))) {
            graphics::points(cbind(mspace$x[, args$axes[1]], 0), pch = args$pch.points,
                             bg = args$bg.points, cex = args$cex.points)
          } else {
            graphics::points(cbind(mspace$x[, args$axes[1]], 0), pch = args$pch.points,
                             col = args$col.points, cex = args$cex.points)
          }
        }
      } else {
        if(args$density.points == TRUE) {
          dens <- stats::density(mspace$x[, 1])
          graphics::polygon(dens$x, dens$y/max(dens$y), lwd = 2,
                            col = grDevices::adjustcolor(1, alpha.f = 0.5))
        }
        graphics::abline(h = 0)
        if(any(args$pch.points %in% c(21:25))) {
          graphics::points(cbind(mspace$x[, 1], 0), pch = args$pch.points,
                           bg = args$bg.points, cex = args$cex.points)
        } else {
          graphics::points(cbind(mspace$x[, 1], 0), pch = args$pch.points,
                           col = args$col.points, cex = args$cex.points)
        }
      }
    }

    if(groups == TRUE & !is.null(mspace$gr_class)) {
      if(is.null(mspace$gr_class)) {
        stop("groups classification has not been added to mspace object")
      } else {
        if(ncol(mspace$x) > 1) {
          if(length(args$axes) > 1) {
            if(args$ellipse.groups == FALSE) {
              hulls_by_group_2D(mspace$x[, args$axes], fac = mspace$gr_class,
                                col = args$col.groups, lty = args$lty.groups,
                                lwd = args$lwd.groups, alpha = args$alpha.groups)
            } else {
              ellipses_by_group_2D(mspace$x[, args$axes], fac = mspace$gr_class,
                                   conflev = args$conflev.groups, col = args$col.groups,
                                   lty = args$lty.groups, lwd = args$lwd.groups,
                                   alpha = args$alpha.groups)
            }
          } else {
            if(args$density.groups == TRUE) {
              dens <- lapply(seq_len(nlevels(mspace$gr_class)), function(i) {
                subdens <- stats::density(mspace$x[mspace$gr_class == levels(mspace$gr_class)[i],
                                                   args$axes[1]])
                list(x = subdens$x, y = subdens$y)
              })
              ymax <- max(unlist(lapply(dens, function(x) {x$y})))

              graphics::abline(h = 0)
              for(i in seq_len(nlevels(mspace$gr_class))) {
                graphics::polygon(dens[[i]]$x, dens[[i]]$y / ymax, lwd = 2,
                                  col = grDevices::adjustcolor(i, alpha.f = alpha.groups))
              }
            }
          }

        } else {
          if(args$density.groups == TRUE) {
            dens <- lapply(seq_len(nlevels(mspace$gr_class)), function(i) {
              subdens <- stats::density(mspace$x[mspace$gr_class == levels(mspace$gr_class)[i], 1])
              list(x = subdens$x, y = subdens$y)
            })
            ymax <- max(unlist(lapply(dens, function(x) {x$y})))

            graphics::abline(h=0)
            for(i in seq_len(nlevels(mspace$gr_class))) {
              graphics::polygon(dens[[i]]$x, dens[[i]]$y/ymax, lwd = 2,
                                col = grDevices::adjustcolor(i, alpha.f = alpha.groups))
            }
          }
        }
      }
    }

    if(phylo == TRUE & !is.null(mspace$phylo)) {
      if(is.null(mspace$phylo)) {
        stop("phylogenetic relationships have not been added to mspace object")
      } else {
        if(length(args$axes) > 1) {
          for(i in seq_len(nrow(mspace$phylo$edge))) {
            graphics::lines(rbind(mspace$phylo_scores[mspace$phylo$edge[i, 1], args$axes],
                                  mspace$phylo_scores[mspace$phylo$edge[i, 2], args$axes]),
                            lwd = args$lwd.phylo, lty = args$lty.phylo, col = args$col.phylo)
          }
        } else {
          print("phylogenetic relationships are omitted from univariate morphospaces")
        }
      }
    }

    if(mshapes == TRUE & !is.null(mspace$gr_centroids)) {
      if(is.null(mspace$gr_centroids)) {
        stop("groups centroids have not been added to mspace object")
      } else {
        if(ncol(mspace$x) > 1) {
          if(length(args$axes) > 1) {
            if(any(args$pch.groups %in% c(21:25))) {
              graphics::points(mspace$gr_centroids[, args$axes], bg = args$bg.groups,
                               pch = args$pch.groups, cex = args$cex.groups)
            } else {
              graphics::points(mspace$gr_centroids[, args$axes], col = args$col.groups,
                               pch = args$pch.groups, cex = args$cex.groups)
            }
          } else {
            graphics::abline(h = 0)
            if(any(args$pch.groups %in% c(21:25))) {
              graphics::points(cbind(mspace$gr_centroids[, args$axes[1]], 0), bg = args$bg.groups,
                               pch = args$pch.groups, cex = args$cex.groups)
            } else {
              graphics::points(cbind(mspace$gr_centroids[, args$axes[1]], 0), col = args$col.groups,
                               pch = args$pch.groups, cex = args$cex.groups)
            }
          }
        } else {
          graphics::abline(h = 0)
          if(any(args$pch.groups %in% c(21:25))) {
            graphics::points(cbind(mspace$gr_centroids[, 1], 0), bg = args$bg.groups,
                             pch = args$pch.groups, cex = args$cex.groups)
          } else {
            graphics::points(cbind(mspace$gr_centroids[, 1], 0), col = args$col.groups,
                             pch = args$pch.groups, cex = args$cex.groups)
          }
        }
      }
    }

    if(shapeax == TRUE & !is.null(mspace$shape_axis)) {
      if(length(args$axes) > 1) {

        if(length(args$lwd.axis) != length(mspace$shape_axis)) {
          args$lwd.axis <- rep(1, length(mspace$shape_axis))
        }
        if(length(args$lty.axis) != length(mspace$shape_axis)) {
          args$lty.axis <- rep(1, length(mspace$shape_axis))
        }
        if(length(args$col.axis) != length(mspace$shape_axis)) {
          args$col.axis <- rep(1, length(mspace$shape_axis))
        }

        for(i in seq_len(length(mspace$shape_axis))) {
          graphics::lines(mspace$shape_axis[[i]][, args$axes],
                          col = args$col.axis[i], lwd = args$lwd.axis[i], lty = args$lty.axis[i])
        }
      } else {
        print("\nmorphometric axes are omitted from univariate morphospaces")
      }
    }

    ##############################################################################################

  } else { #if x or y have been provided, show hybrid morphospace

    if(is.null(axes)) {
      args$axes <- rep(args$axes[1], 2)
    }

    #if x/y is a phy object, prepare the ground for a phenogram
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



    if(mspace$plotinfo$k == 3 & mspace$datype == "landm") {

      shapemodels <- morphogrid(ordination = ordination, axes = args$axes, datype = mspace$datype,
                                template = NULL, x = x, y = y, p = mspace$plotinfo$p,
                                k = mspace$plotinfo$k, nh = args$nh, nv = args$nv, mag = args$mag,
                                asp = args$asp, xlim = args$xlim, ylim = args$ylim, rot.models = args$rot.models,
                                size.models = args$size.models, asp.models = args$asp.models)

      refshape <- matrix(rev_eigen(0,
                                   ordination$rotation[, 1],
                                   ordination$center),
                         nrow = mspace$plotinfo$p, ncol = mspace$plotinfo$k, byrow = TRUE)

      if(is.null(x)) xlim <- range(ordination$x[,args$axes[1]])
      if(is.null(y)) ylim <- range(ordination$x[,args$axes[1]])

      plot_morphogrid3d(x = x, y = y, morphogrid = shapemodels, refshape = refshape,
                        template = args$template, links = args$links, ordtype = mspace$ordtype,
                        axes = args$axes, xlim = xlim, ylim = ylim, adj_frame = args$adj_frame,
                        xlab = args$xlab, ylab = args$ylab, cex.ldm = args$cex.ldm,
                        col.ldm = args$col.ldm, col.models = args$col.models, lwd.models = args$lwd.models,
                        bg.models = args$bg.models,  size.models = args$size.models,
                        asp.models = args$asp.models, alpha.models = args$alpha.models, models = args$models)
    } else {

      shapemodels <- morphogrid(ordination = ordination, axes = args$axes, datype = mspace$datype,
                                template = args$template, x = x, y = y, p = mspace$plotinfo$p,
                                k = mspace$plotinfo$k, nh = args$nh, nv = args$nv, mag = args$mag,
                                asp = args$asp, xlim = args$xlim, ylim = args$ylim, rot.models = args$rot.models,
                                size.models = args$size.models, asp.models = args$asp.models)

      plot_morphogrid2d(x = x, y = y, morphogrid = shapemodels, template = args$template,
                        links = args$links, datype = mspace$datype, ordtype = mspace$ordtype,
                        axes = args$axes, adj_frame = args$adj_frame, p = mspace$plotinfo$p,
                        xlab = args$xlab, ylab = args$ylab, cex.ldm = args$cex.ldm, col.ldm = args$col.ldm,
                        col.models = args$col.models, lwd.models = args$lwd.models,
                        bg.models = args$bg.models, models = args$models)
    }


    if(phenogr == TRUE) { #if x/y is a phy object, plot a phenogram

      if(length(args$cex.groups) == 1) args$cex.groups <- rep(args$cex.groups, nlevels(mspace$gr_class))
      if(length(args$pch.groups) == 1) args$pch.groups <- rep(args$pch.groups, nlevels(mspace$gr_class))

      maptips <-  order(match(rownames(mspace$gr_centroids), tree$tip.label))
      plot_phenogram(x = x, y = y, tree = tree, axis = args$axes, points = points,
                     pch.groups = args$pch.groups[maptips], phylo_scores = mspace$phylo_scores,
                     cex.groups = args$cex.groups[maptips], col.groups = args$col.groups[maptips],
                     bg.groups = args$bg.groups[maptips], lwd.phylo = args$lwd.phylo,
                     lty.phylo = args$lty.phylo, col.phylo = args$col.phylo)

    } else { #else, go for a generic hybrid morphospace

      xy <- cbind(x, mspace$x[,args$axes[1]], y)

      #add points, hulls, phylogeny, and/or consensus
      if(points == TRUE) {
        if(any(args$pch.points %in% c(21:25))) {
          graphics::points(xy, pch = args$pch.points, bg = args$bg.points,
                           cex = args$cex.points)
        } else {
          graphics::points(xy, pch = args$pch.points, col = args$col.points,
                           cex = args$cex.points)
        }
      }

      if(groups == TRUE & !is.null(mspace$gr_class)) {
        if(is.null(mspace$gr_class)) {
          stop("groups classification has not been added to mspace object")
        } else {

          if(ellipse.groups == FALSE) {
            hulls_by_group_2D(xy, fac = mspace$gr_class, col = args$col.groups,
                              lty = args$lty.groups, lwd = args$lwd.groups, alpha = args$alpha.groups)
          } else {
            ellipses_by_group_2D(xy, fac = mspace$gr_class, conflev = args$conflev.groups,
                                 col = args$col.groups, lty = args$lty.groups,
                                 lwd = args$lwd.groups, alpha = args$alpha.groups)
          }
        }
      }

      if(phylo == TRUE & !is.null(mspace$phylo)) {
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
            nodesxy <- apply(meanxy[mspace$phylo$tip.label,], 2, phytools::fastAnc, tree = mspace$phylo)
            phyloxy <- rbind(meanxy, nodesxy)

          }

          for(i in seq_len(nrow(mspace$phylo$edge))) {
            graphics::lines(rbind(phyloxy[mspace$phylo$edge[i, 1],],
                                  phyloxy[mspace$phylo$edge[i, 2],]),
                            lwd = args$lwd.phylo, lty = args$lty.phylo, col = args$col.phylo)
          }
        }
      }

      if(mshapes == TRUE & !is.null(mspace$gr_centroids)) {
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

          meanxy <- cbind(meanx, mspace$gr_centroids[,args$axes[1]], meany)[mspace$phylo$tip.label,]

          if(any(args$pch.groups %in% c(21:25))) {
            graphics::points(meanxy, bg = args$bg.groups, pch = args$pch.groups,
                             cex = args$cex.groups)
          } else {
            graphics::points(meanxy, col = args$col.groups, pch = args$pch.groups,
                             cex = args$cex.groups)
          }
        }
      }
    }
  }

  if(legend == TRUE) {
    graphics::layout(matrix(1, nrow=1))
  }
}

