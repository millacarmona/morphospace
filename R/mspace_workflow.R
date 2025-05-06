################################################################################

#' Generate morphospace
#'
#' @description Build morphospaces using a variety of multivariate methods,
#'    and depict shape variation represented by the resulting ordination axes.
#'
#' @param shapes Shape data.
#' @param ord Optional object of class \code{"prcomp"}, \code{"bg_prcomp"},
#'    \code{"pls_shapes"}, \code{"phy_pls_shapes"}, \code{"burnaby"},
#'    \code{"phy_burnaby"}, \code{"gm.prcomp"}, \code{"pls"}, \code{"bgPCA"},
#'    \code{"pls2B"}, \code{"phyl.pca"}, \code{"mvgls.pca"} or \code{"PCA"},
#'    containing the results of multivariate ordination of shape data. To be
#'    used instead of the \code{shapes} argument.
#' @param axes Numeric of length 1 (univariate morphospace) or 2 (bivariate
#'    morphospace), indicating the axes to be plotted.
#' @param links A list with the indices of the coordinates defining a wireframe
#'   connecting landmarks (following the format used in \code{Morpho}), or a
#'   2-columns matrix indicating the pairs of landmarks that should be linked
#'    (following the format used in \code{geomorph}).
#' @param template Either a 2-column matrix with landmarks/semilandmarks and
#'    template curves coordinates (for 2D shape data) or a \code{"mesh3d"}
#'    object \strong{representing the mean shape of the sample} (for 3D shape
#'    data). See details below.
#' @param FUN The function (method) to be used for ordination of shape
#'    variation. Supported alternatives include \code{\link[stats]{prcomp}},
#'    \code{\link{bg_prcomp}}, \code{\link{pls_shapes}}, \code{\link{burnaby}},
#'    \code{\link[geomorph]{gm.prcomp}}, \code{\link[geomorph]{two.b.pls}},
#'    \code{\link[Morpho]{groupPCA}}, \code{\link[phytools]{phyl.pca}} and
#'    \code{\link[Momocs]{PCA}}. Ignored if \code{ord} is provided.
#' @param datype Character; the type of shape data used (either \code{"fcoef"}
#'    or \code{"landm"}). Only required if \code{ord} is provided instead of
#'    \code{shapes}.
#' @param p Numeric, indicating the number of landmarks/semilandmarks used (only
#'   for landmark data in 2-margin matrices format).
#' @param k Numeric, indicating the number of Cartesian dimensions of
#'   landmarks/semilandmarks (only for landmark data in 2-margin matrices
#'   format).
#' @param nh Positive integer; number of shape models along the x axis.
#' @param nv Positive integer; number of shape models along the y axis.
#' @param mag Numeric; magnifying factor for shape models.
#' @param invax Optional integer indicating which of the axes provided in
#'    \code{axes} needs to be inverted (options are \code{1}, \code{2} or
#'    \code{c(1,2)}).
#' @param rescale Logical; whether to re-scale background shape models so shape
#'    variation is shown more clearly.
#' @param adj_frame Numeric of length 2, providing \emph{a posteriori} scaling
#'    factors for the width and height of the frame, respectively.
#' @param rot.models Numeric; angle (in degrees) to rotate shape models.
#'    Alternatively, a pre-defined rotation matrix for 3D landmark data
#'    can be provided.
#' @param size.models Numeric; size factor for shape models.
#' @param asp.models Numeric; the y/x aspect ratio of shape models.
#' @param col.models Color for wireframes/outlines.
#' @param bg.models Background color for outlines/meshes.
#' @param lwd.models Integer; width of the lines in wireframes/outlines.
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
#'    a variety of eigenanalysis-based multivariate methods. Alternatively,
#'    the results of multivariate ordinations achieved using functions from
#'    \code{geomorph}, \code{Morpho}, \code{Momocs}, \code{phytools} or
#'    \code{mvMORPH} can be provided directly using the \code{ord} argument
#'    (this requires specification of the \code{datype} argument, and also
#'    of \code{p} and \code{k} for landmark data).
#'
#'    In total, the user can choose from a range of ordination methods to
#'    capture the signal of interest including PCA, Between-groups PCA,
#'    Phylogenetic PCA, 2-blocks PLS, Phylogenetically corrected 2-block PLS,
#'    Phylogenetically-aligned component analysis, and the Burnaby approach for
#'    orthogonalization against nuisance variables (also including a
#'    phylogenetically corrected variant).
#'
#'    Have in mind that when using the \code{geomorph::two.b.pls} or
#'    \code{Morpho::pls2B} in the context of \code{mspace}, it is assumed that
#'    the shape variables the user wants to use to generate the morphospace are
#'    provided in the second (aka response) block (i.e., arguments \code{A1} and
#'    \code{X}, respectively).
#'
#'    \code{mspace} will generate a series of shape models depicting the range
#'    of variation. The resulting \code{"mspace"} object stores all the
#'    information necessary to project and/or retrieve new shapes into the
#'    ordination space. The latter can be populated using the \code{proj_*}
#'    family of functions and the \code{%>%} operator from \code{magrittr}.
#'    Background shapes, shape axes, or other elements projected into
#'    morphospace can be extracted for other uses through
#'    \code{\link{extract_shapes}}.
#'
#'    For landmark data, representation of shape variation can be further
#'    enhanced by providing a template, which contains additional geometric
#'    features from the structure the landmark/semilandmarks were placed upon.
#'    Templates will be warped using TPS interpolation to produce the set of
#'    background shape models. For 2D landmark data, templates must be provided
#'    as a 2-column matrix containing the actual landmarks/semilandmarks,
#'    followed by the coordinates defining a curve or set of curves, separated
#'    from the former and from each other by a row of \code{NA}s (see
#'    \code{\link{build_template2d}}). For 3D landmark data, the template
#'    must be a \code{"mesh3d"} object corresponding to the \strong{actual mean
#'    shape of the sample} (which can be obtained using
#'    \code{\link{expected_shapes}} + \code{\link[Morpho]{tps3d}}; see examples
#'    below).
#'
#'
#' @return An object of class \code{"mspace"}, which is a list containing:
#'   \itemize{
#'   \item \code{$ordination:} a list with the output from the ordination method
#'      used, styled in the \code{\link{prcomp}} format (typically containing at
#'      least an \code{$x}, \code{$rotation}, and \code{$center} slots. Also
#'      contains the type of data (\code{$datype}) and ordination method
#'      (\code{$ordtype}) used.
#'   \item \code{$projected:} a list containing the elements that have been
#'      projected into the morphospace stored in the \code{"mspace"} object.
#'      Initially includes only the background shape models
#'      (\code{$shapemodels}), but more elements can be added using the
#'      \code{proj_*} family of functions.
#'   \item \code{$plotinfo:} a list containing the graphical parameters used to
#'      create the plot. Passed to \code{\link{plot_mspace}} to regenerate
#'      morphospaces.
#'   }
#'
#' @seealso \code{\link{proj_shapes}}, \code{\link{proj_groups}},
#'   \code{\link{proj_phylogeny}}, \code{\link{proj_axis}},
#'   \code{\link{proj_landscape}}, \code{\link{extract_shapes}},
#'   \code{\link{prcomp}}, \code{\link{bg_prcomp}}, \code{\link{pls_shapes}},
#'   \code{\link[geomorph]{gm.prcomp}}, \code{\link[geomorph]{two.b.pls}},
#'   \code{\link[Morpho]{groupPCA}}, \code{\link[Morpho]{pls2B}},
#'   \code{\link[phytools]{phyl.pca}}, \code{\link[mvMORPH]{mvgls.pca}},
#'   \code{\link[geomorph]{gm.prcomp}}, \code{\link[Morpho]{groupPCA}},
#'   \code{\link[Morpho]{pls2B}}, \code{\link[phytools]{phyl.pca}},
#'   \code{\link[mvMORPH]{mvgls.pca}}, \code{\link[Momocs]{PCA}}
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
#' sp_shapes <- expected_shapes(shapes, species)
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
#' #just create a morphospace without plotting, save into an object, and inspect
#' morphosp <- mspace(shapes, links = links, mag = 0.7, axes = c(1,2),
#'                    plot = FALSE)
#' morphosp #general info about the object
#' names(morphosp) #slots
#'
#'
#' # Generate morphospaces using functions from other packages (note that
#' # aside from replacing the shape argument with the ord one, additional
#' #information on number of dimensions and coordinates is required for
#' #landmark data):
#'
#' #phylogenetically-aligned component analysis via geomorph::gm.prcomp
#' library(geomorph)
#' paca <- gm.prcomp(A = two.d.array(sp_shapes), phy = tree,
#'                   align.to.phy = TRUE)
#' mspace(ord = paca, datype = "landm", p = 9, k = 2, links = links,
#'        points = TRUE)
#'
#' #phylogenetic PCA via phytools::phyl.pca
#' library(phytools)
#' phypca <- phyl.pca(Y = two.d.array(sp_shapes), tree = tree)
#' mspace(ord = phypca, datype = "landm", p = 9, k = 2, links = links,
#'        points = TRUE)
#'
#' #between-groups PCA via Morpho::groupPCA
#' library(Morpho)
#' bgpca <- groupPCA(two.d.array(shapes), groups = species, mc.cores = 1,
#'                   rounds = 0, cv = FALSE)
#' mspace(ord = bgpca, datype = "landm", p = 9, k = 2, links = links,
#'        points = TRUE)
#'
#' #phylogenetic PCA via mvMORPH::mvgls.pca **NOTICE a $center slot containing
#' #the center of the original variables need to be appended to the mvgls.pca
#' #object; in this case, this is the phylogenetic mean
#' library(mvMORPH)
#' phypca2 <- mvgls.pca(object = mvgls(two.d.array(sp_shapes) ~ 1, tree = tree,
#'                      model = "BM"), plot = FALSE)
#' phypca2$center <- apply(two.d.array(sp_shapes), 2, phytools::fastAnc, tree = tree)[1,]
#' mspace(ord = phypca2, datype = "landm", p = 9, k = 2, links = links,
#'        points = TRUE)
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
#' if (interactive()) {
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
#'        rescale = FALSE, asp.models = 1, bg.model = "light gray")
#'
#' #shapes in the background can be rescaled to avhieve a (slightly) better
#' #visualization. Also, save the ordination into an object.
#' morphosp_F <- mspace(shapes, mag = 1, axes = c(1,2), nh = 5, nv = 4,
#'                      size.models = 1, asp.models = 1,
#'                      bg.model = "light gray")
#'
#' #inspect the contents of the object
#' morphosp_F
mspace <- function(shapes = NULL,
                   ord = NULL,
                   axes = c(1,2),
                   links = NULL,
                   template = NULL,
                   FUN = stats::prcomp,
                   datype = NULL,
                   p = NULL,
                   k = NULL,
                   nh = 5,
                   nv = 4,
                   mag = 1,
                   invax = NULL,
                   asp = 1,
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

  # ensure only one of the two approaches to entry the data is being used
  if(is.null(shapes) & is.null(ord)) stop("Either (1) shapes (and optionally FUN), or (2) ord and datype must be provided")

  # if shapes are provided, get them into '2D' formar and perform ordination
  if(!is.null(shapes)) {

    if(!is.null(ord)) stop("Drop one of shapes or ord; cannot be used simultaneously")

    dat <- shapes_mat(shapes)
    data2d <- dat$data2d
    datype <- dat$datype

    FUN <- match.fun(FUN)
    ordination <- FUN(data2d, ...)
    ordination <- adapt_ordination(ordination)
    ordination$ordtype <- class(ordination)[1]
    ordination$datype <- datype
  }

  # if an ordination is provided, get data2d and shapes
  if(!is.null(ord)) {

    if(is.null(datype)) stop("provide a value for datype ('landm' or 'fcoef')")

    ordination <- adapt_ordination(ord)

    ordination$ordtype <- class(ordination)
    ordination$datype <- datype

    data2d <- rev_eigen(ordination$x,
                        ordination$rotation,
                        ordination$center)
    shapes <- if(datype == "landm") geomorph::arrayspecs(data2d, k = k, p = p) else data2d

  }

  # if shape data is landmark data, register their number and dimensions; else,
  # prepare grounds for using inv_efourier further down
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

  # invert axes if requested
  if(!is.null(invax)) {
    ordination$x[,axes[invax]] <- ordination$x[,axes[invax]] * -1
    ordination$rotation[,axes[invax]] <- ordination$rotation[,axes[invax]] * -1
  }

  # if only one axis was requested and/or resulted from ordination, prepare
  # ground for univariate morphospace
  if(ncol(ordination$x) == 1 | length(axes) == 1) {
    y <- rep(1, nrow(ordination$x))
    ylim <- c(0, 1)
    ylab <- "relative density"
  } else {
    y <- NULL
  }

  # if only one axis is requested, override default aspect ratio for univariate
  # morphospace
  if(length(axes) == 1) asp <- NA

  # axes is set to == 1 when there is only one ordination axis
  if(ncol(ordination$x) == 1) axes <- 1

  # prepare rotation matrix introduction for 3D data
  if(!is.null(dim(rot.models))) {
    rotation_matrix <- rot.models
    rot.models <- 0
  } else {
    rotation_matrix <- NULL
  }

  # if links were provided in geomorph format, adapt them to Morpho's
  if(!is.null(links)) {
    if(!is.null(dim(links))) {
      links0 <- links
      links <- split(links0, row(links0))
      links <- lapply(1:nrow(links0), function(i) {links0[i, ]})
    }
  }

  # plot 2D / 3D morphospaces
  if(k == 3 & datype == "landm") {

    if(is.null(rotation_matrix)) {
      if(rot.models != 0) stop("Angles are not allowed for 3D shapes; please provide a rotation matrix (see 'set_rotation3D')")
    }

    shapemodels <- morphogrid(ordination = ordination, axes = axes, datype = datype, template = NULL,
                              x = NULL, y = y, p = p, k = k, nh = nh, nv = nv, mag = mag,
                              asp = asp, xlim = xlim, ylim = ylim, rot.models = rot.models,
                              size.models = size.models, asp.models = 1, rescale = rescale)

    refshape <- expected_shapes(shapes)

    rgl::par3d(userMatrix = rotation_matrix)

    plot_morphogrid3d(x = NULL, y = y, morphogrid = shapemodels, refshape = refshape,
                      template = template, links = links, ordtype = ordination$ordtype,
                      adj_frame = adj_frame, axes = axes, xlim = xlim, ylim = ylim, xlab = xlab,
                      ylab = ylab, cex.ldm = cex.ldm, col.ldm = col.ldm, col.models = col.models,
                      lwd.models = lwd.models, bg.models = bg.models, size.models = size.models,
                      asp.models = asp.models, alpha.models = alpha.models, plot = plot,
                      rotmat = rotation_matrix, models = models)

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

  # plot points if requested
  if(points) graphics::points(ordination$x[,axes])

  # prepare shape models for appending
  if(length(axes) == 1) {
    if(datype == "landm") shapemodels <- shapemodels$shapemodels[,,1:nh]
    if(datype == "fcoef") shapemodels <- shapemodels$shapemodels[1:nh,]
  } else {
    shapemodels <- shapemodels$shapemodels
  }

  # save plot's graphical parameters
  plotinfo <- list(p = p, k = k, links = links, template = template, axes = axes, nh = nh, nv = nv, mag = mag,
                   xlim = xlim, ylim = ylim, asp = asp, rescale = rescale, adj_frame = adj_frame,
                   asp.models = asp.models, rot.models = rot.models, size.models = size.models,
                   lwd.models = lwd.models, bg.models = bg.models, col.models = col.models,
                   alpha.models = alpha.models, cex.ldm = cex.ldm, col.ldm = col.ldm,
                   models = models, rotation_matrix = rotation_matrix, plot = plot)

  # prepare output
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
#' @param labels Either logical, indicating whether to include all point labels,
#'   or a character string with the exact names of the points whose labels
#'   should be included.
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
#'
#' #add labels for all the points (needs a named 'shapes' object)
#' mspace(shapes, template = template, mag = 0.7, axes = c(1,2)) %>%
#'   proj_shapes(shapes = shapes, labels = TRUE)
#'
#' #add labels for selected points (needs a named 'shapes' object)
#' mspace(shapes, template = template, mag = 0.7, axes = c(1,2)) %>%
#'   proj_shapes(shapes = shapes) %>%
#'   proj_shapes(shapes = shapes[,,c("Dk M Op 8 Ad_02", "Dk H Tr 11 Ad_03")],
#'               pch = 16, col = 2, labels = c("Dk M Op 8 Ad_02", "Dk H Tr 11 Ad_03"))
proj_shapes <- function(mspace, shapes = NULL, density = TRUE,
                        labels = NULL, pipe = TRUE, ...) {


  args <- c(as.list(environment()), list(...))

  if(is.null(shapes)) {
    scores <- mspace$ordination$x
  } else {
    dat <- shapes_mat(shapes)
    data2d <- dat$data2d
    datype <- dat$datype

    if(mspace$ordination$datype != datype) stop("shapes and mspace types are not compatible")

    scores <- proj_eigen(x = data2d, vectors = mspace$ordination$rotation,
                         center = mspace$ordination$center)
  }


  if(mspace$plotinfo$plot) {

    if(ncol(mspace$ordination$x) > 1) {
      if(length(mspace$plotinfo$axes) > 1) {
        plot_biv_scatter(scores = scores[, mspace$plotinfo$axes], ...)
        add_labels(scores[, mspace$plotinfo$axes], labels)
      } else {
        plot_univ_scatter(scores = cbind(scores[, mspace$plotinfo$axes[1]]),
                          density = density, ...)
        add_labels(cbind(scores[, mspace$plotinfo$axes[1]], 0), labels,
                   srt = 90, adj = c(0,0))
      }
    } else {
      plot_univ_scatter(scores = scores, density = density, ...)
      add_labels(cbind(scores[, 1], 0), labels, srt = 90,
                 srt = 90, adj = c(0,0))
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
  mspace$plotinfo$labels.points <- c(mspace$plotinfo$labels.points, labels)

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
#' @param labels Either logical, indicating whether to include all group labels,
#'   or a character string with the exact names of the groups whose labels
#'   should be included.
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
#' mspace(shapes, links = links, mag = 0.7, axes = c(1,2),
#'        bg.model = "light gray") %>%
#'   proj_shapes(shapes = shapes, col = c(1:4)[species], pch = 1) %>%
#'   proj_groups(groups = species, col = 1:4, lwd = 2, ellipse = TRUE,
#'               conflev = 0.95, alpha = 0.1)
#'
#' #add labels for all the groups
#' mspace(shapes, links = links, mag = 0.7, axes = c(1,2),
#'        bg.model = "light gray") %>%
#'   proj_shapes(shapes = shapes, col = c(1:4)[species], pch = 1) %>%
#'   proj_groups(groups = species, col = 1:4, lwd = 2, ellipse = TRUE,
#'               conflev = 0.95, alpha = 0.1, labels = TRUE)
#'
#' #add labels for selected the groups
#' msp <- mspace(shapes, links = links, mag = 0.7, axes = c(1,2),
#'               bg.model = "light gray") %>%
#'   proj_shapes(shapes = shapes, col = c(1:4)[species], pch = 1) %>%
#'   proj_groups(groups = species, col = 1:4, lwd = 2, ellipse = TRUE,
#'               conflev = 0.95, alpha = 0.1, labels = c("coihuicoensis", "esbelta"))
#'
#' #add legend using plot_mspace()
#' plot_mspace(msp, legend = TRUE, cex.legend = 1, pch.groups = 16)
proj_groups <- function(mspace, shapes = NULL, groups = NULL, ellipse = FALSE,
                        conflev = 0.95, density = TRUE, labels = NULL, pipe = TRUE, ...) {

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

  if(is.null(groups)) groups <- factor(rep(1, nrow(data2d)))

  if(is.null(args$col)) args$col <- sort(unique(as.numeric(groups)))
  args$col <- if(length(args$col) == 1) rep(args$col, nlevels(groups)) else {
    if(nlevels(groups) > length(args$col)) {
      c(rep(NA, nlevels(groups) - length(unique(args$col))),
        unique(args$col))[order(c((1:nlevels(groups))[-which(levels(groups) %in% unique(groups))],
                                  which(levels(groups) %in% unique(groups))))]
    } else args$col
  }
  if(is.factor(args$col)) args$col <- as.numeric(args$col)

  if(mspace$plotinfo$plot) {

    mshapes <- apply(X = data2d, MARGIN = 2, FUN = tapply, groups, mean)
    gcols <- stats::setNames(col2hex(args$col), levels(groups))

    if(ncol(mspace$ordination$x) > 1) {
      if(length(mspace$plotinfo$axes) > 1) {

        if(!ellipse) {
          hulls_by_group_2D(data2d[, mspace$plotinfo$axes], fac = groups, ...)
        } else {
          ellipses_by_group_2D(data2d[, mspace$plotinfo$axes], fac = groups,
                               conflev = conflev, ...)
        }
        add_labels(xy = mshapes[, mspace$plotinfo$axes],
                   labels = labels, col = gcols)
      } else {
        if(density) {
          dens <- density_by_group_2D(data2d, groups, ax = mspace$plotinfo$axes[1], ...)

          gymax <- lapply(dens$dens, function(x){
            if(is.null(x$y)) x$y <- 0
            max(x$y) / dens$ymax
          })
          gxmax <- lapply(dens$dens, function(x){
            if(!is.null(x$x)) x$x[which.max(x$y)] else NA
          })
          gxmax[is.na(gxmax)] <- mshapes[is.na(gxmax), mspace$plotinfo$axes[1]]
          xy <- cbind(gxmax, gymax)
          rownames(xy) <- levels(groups)

          add_labels(xy = xy, labels = labels, col = gcols, pos = 3)
        }
      }

    } else {
      if(density) {
        dens <- density_by_group_2D(data2d, groups, ax = 1, ...)

        gymax <- lapply(dens$dens, function(x){
          if(is.null(x$y)) x$y <- 0
          max(x$y) / dens$ymax
        })
        gxmax <- lapply(dens$dens, function(x){
          if(!is.null(x$x)) x$x[which.max(x$y)] else NA
        })
        gxmax[is.na(gxmax)] <- mshapes[is.na(gxmax), 1]
        xy <- cbind(gxmax, gymax)
        rownames(xy) <- levels(groups)

        add_labels(xy = xy, labels = labels, col = gcols, pos = 3)
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
  mspace$plotinfo$labels.groups <- args$labels
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
#' @param obj An object of class \code{"mlm"}, \code{"procD.lm"},
#'   \code{"lm.rrpp"}, \code{"mvgls"} or \code{"mvols"}, containing a linear
#'   model fit to shape data, or an object of class \code{"prcomp"},
#'   \code{"bg_prcomp"}, \code{"pls_shapes"}, \code{"phy_pls_shapes"},
#'   \code{"burnaby"}, \code{"phy_burnaby"}, \code{"gm.prcomp"}, \code{"pls"},
#'   \code{"bgPCA"}, \code{"pls2B"}, \code{"phyl.pca"}, \code{"mvgls.pca"} or
#'   \code{"PCA"} with the results of multivariate ordination of shape data.
#' @param axis Optional; which axis from \code{obj} should to be projected?
#' @param mag Numeric; magnifying factor for representing shape transformation.
#' @param type Integer; type of arrows (\code{0} = no arrow; \code{1} = pointing
#'   towards the maximum; \code{2} = pointing towards the maximum, \code{3} =
#'   pointing in both directions).
#' @param pipe Logical; is the function being included in a pipe?
#' @param ... Further arguments passed to \code{\link[graphics]{arrows}}.
#'
#' @details This function is primarily aimed at graphically representing
#'   morphometric axes (estimated using either linear models or multivariate
#'   ordination methods) into an existing morphospace for heuristic exploration
#'   of patterns in the data. It can also be used to extract theoretical shapes
#'   at the extremes of those axes, although \code{\link{ax_transformation}} does the
#'   same thing in a more flexible and straightforward way.
#'
#'   Axes computed by fitting linear models to shape data can differ in
#'   extension from axes obtained through supervised ordination using the same
#'   supervising variable (i.e., shape transformations will be either stretched
#'   or truncated). This is because the former assumes that the explanatory
#'   variable has been measured without error. Also, the former will not be
#'   necessarily centered.
#'
#'   For statistical analysis of axes (e.g., trajectory analysis) their vector
#'   coefficients can be extracted directly from slope or eigenvector
#'   coefficients stored in the objects returned by linear model fitting and
#'   multivariate ordination methods.
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
#' @seealso \code{\link{ax_transformation}}, \code{\link[stats]{lm}},
#'   \code{\link[geomorph]{procD.lm}}, \code{\link[geomorph]{procD.pgls}},
#'   \code{\link[RRPP]{lm.rrpp}}, \code{\link[mvMORPH]{mvols}},
#'   \code{\link[mvMORPH]{mvgls}}, \code{\link[geomorph]{gm.prcomp}},
#'   \code{\link[geomorph]{two.b.pls}}, \code{\link[Morpho]{groupPCA}},
#'   \code{\link[Morpho]{pls2B}}, \code{\link[phytools]{phyl.pca}},
#'   \code{\link[mvMORPH]{mvgls.pca}}, \code{\link[Momocs]{PCA}}
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
#' phypca <- phytools::phyl.pca(two.d.array(sp_shapes), tree = tree)
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

  if(mspace$plotinfo$plot) {
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

  if(is.null(args$col)) args$col <- 1
  mspace$plotinfo$col.axis <- c(mspace$plotinfo$col.axis, col2hex(args$col))

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
#' @param evmodel Character, specifying an evolutionary model to perform
#'   ancestral character reconstruction; options are "BM" (Brownian motion;
#'   implemented through \code{\link[phytools]{fastAnc}}), "EB" (Early burst),
#'   "OU" (Ornstein-Uhlenbeck) and "lambda" (Pagel's lambda transformation)
#'   (implemented using \code{\link[mvMORPH]{mvgls}}).
#' @param labels.tips Either logical, indicating whether to include tip labels,
#'   or a character string with the exact names of the tips whose labels
#'   should be included.
#' @param pch.tips Symbol of the scatter points corresponding to the tips of the
#'   phylogeny.
#' @param col.tips Color of the hulls/ellipses and/or scatter points
#'   corresponding to the tips of the phylogeny.
#' @param bg.tips Background color of the scatter points corresponding to the
#'   tips of the phylogeny.
#' @param cex.tips Numeric; size of the scatter points corresponding to the
#'   tips of the phylogeny.
#' @param labels.nodes Either logical, indicating whether to include node
#'   labels, or a character string with the exact names of the nodes whose
#'   labels should be included (with the form of "node_n", e.g., "node_14"
#'   corresponds to the root of a tree with 13 tips).
#' @param pch.nodes Symbol of the scatter points corresponding to the
#'   nodes of the phylogeny.
#' @param col.nodes Color of the hulls/ellipses and/or scatter points
#'   corresponding to the nodes of the phylogeny.
#' @param bg.nodes Background color of the scatter points corresponding to the
#'   nodes of the phylogeny.
#' @param cex.nodes Numeric; size of the scatter points corresponding to the
#'   nodes of the phylogeny.
#' @param pipe Logical; is the function being included in a pipe?
#' @param ... Further arguments passed to \code{\link[graphics]{lines}}
#'   (commonly, \code{pch}, \code{col}/\code{bg} and \code{cex}.
#'
#' @details The purpose of this function is twofold. First, it is meant to
#'   transform a morphospace into a phylomorphospace by projecting node shapes
#'   and phylogenetic relationships. To this end, a set of named shapes must be
#'   provided; \code{dim(shapes)[3]} must match \code{tree$tip.labels}. Second,
#'   this function can be used to retrieve the scores corresponding to nodes of
#'   the phylogenetic tree (\code{$projected$phylo_scores}, which can then be
#'   used to compute the node shapes using \code{\link{extract_shapes}}. The
#'   position of these shapes in morphospace is estimated using
#'   \code{\link[mvMORPH]{mvgls}} which allows specification of different
#'   phenotypic evolutionary models, passed through the \code{evmodel} argument
#'   (with the exception of the default model, Brownian motion, estimated
#'   using \code{\link[phytools]{fastAnc}}).
#'
#' @return If a plot device with a morphospace is open, shapes representing the
#'   tips and nodes of the phylogenetic tree, as well as the lines connecting
#'   them, are projected into morphospace. If \code{pipe = FALSE} scores for
#'   nodes and tips of the phylogeny are returned invisibly.
#'   If \code{pipe = TRUE} the supplied \code{"mspace"} object will be modified
#'   by appending \code{$phylo_scores}, \code{$phylo_evmodel}, \code{$phylo}
#'   slots to \code{$projected}, as well as by adding some graphical parameters
#'   (stored into the \code{$plotinfo} slot), and returned invisibly.
#'
#' @seealso \code{\link[mvMORPH]{mvgls}}, \code{\link[phytools]{fastAnc}}
#'
#' @export
#'
#' @references
#' Clavel, J., Escarguel, G., & Merceron, G. (2015). \emph{mvMORPH: an R package
#'    for fitting multivariate evolutionary models to morphometric data}.
#'    Methods in Ecology and Evolution, 6(11), 1311-1319.
#' Revell, L.J. (2012). \emph{phytools: an R package for phylogenetic
#'    comparative biology (and other things)}. Methods in Ecology and Evolution,
#'    3, 217-223.
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
proj_phylogeny <- function(mspace, shapes = NULL, tree, evmodel = "BM", labels.tips = FALSE,
                           labels.nodes = FALSE, pch.nodes = 16, col.nodes = "gray", bg.nodes = 1,
                           cex.nodes = 0.8, pch.tips = 16, bg.tips = 1, col.tips = 1, cex.tips = 1,
                           pipe = TRUE, ...) {

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

  if(evmodel == "BM") {
    nodes_scores <- apply(tips_scores, 2, phytools::fastAnc, tree = tree)
    rownames(nodes_scores) <- paste0("node_", rownames(nodes_scores))
  } else {
    mvmod <- mvMORPH::mvgls(tips_scores ~ 1, tree = tree, model = evmodel)
    nodes_scores <- mvMORPH::ancestral(mvmod)
  }
  phylo_scores <- rbind(tips_scores, nodes_scores)

  if(mspace$plotinfo$plot) {

    if(length(mspace$plotinfo$axes) > 1) {
      for(i in seq_len(nrow(tree$edge))) {
        graphics::lines(rbind(phylo_scores[tree$edge[i, 1], mspace$plotinfo$axes],
                              phylo_scores[tree$edge[i, 2], mspace$plotinfo$axes]), ...)
      }

      tips <- seq_len(length(tree$tip.label))
      plot_biv_scatter(scores = phylo_scores[-tips, mspace$plotinfo$axes], pch = pch.nodes,
                       bg = bg.nodes, col = col.nodes, cex = cex.nodes)
      plot_biv_scatter(scores = phylo_scores[tips, mspace$plotinfo$axes], pch = pch.tips,
                       bg = bg.tips, col = col.tips, cex = cex.tips)

      add_labels(tips_scores, labels.tips)
      if(all(labels.nodes %in% rownames(tips_scores)))
        labels.nodes <- paste0("node_", ape::getMRCA(phy = tree, tip = labels.nodes))
      add_labels(nodes_scores, labels.nodes)

    } else {
      cat("\nPhylogenetic relationships are not projected into univariate morphospaces")
    }
  }

  mspace$projected$phylo_scores <- phylo_scores
  mspace$projected$phylo <- tree
  mspace$projected$phylo_evmodel <- evmodel

  if(is.null(args$col)) args$col <- 1

  mspace$plotinfo$col.phylo <- col2hex(args$col)
  mspace$plotinfo$lwd.phylo <- args$lwd
  mspace$plotinfo$lty.phylo <- args$lty
  mspace$plotinfo$col.nodes <- col2hex(args$col.nodes)
  mspace$plotinfo$bg.nodes <- col2hex(args$bg.nodes)
  mspace$plotinfo$pch.nodes <- args$pch.nodes
  mspace$plotinfo$cex.nodes <- args$cex.nodes
  mspace$plotinfo$labels.tips <- args$labels.tips
  mspace$plotinfo$labels.nodes <- args$labels.nodes
  mspace$plotinfo$col.tips <- col2hex(args$col.tips)
  mspace$plotinfo$bg.tips <- col2hex(args$bg.tips)
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
#' @param obj An optional object containing a landscape created using
#'   \code{Morphoscape}. If provided, arguments \code{shapes}, \code{FUN}
#'   or \code{X} are ignored.
#' @param shapes Optional shape data. If provided, a landscape will be computed
#'   for the region of the morphospace encompassing that sample of shapes
#'   ("empirical landscape"). If \code{NULL}, the landscape will be computed for
#'   the set of background shape models ("theoretical landscape").
#' @param FUN An optional \emph{ad hoc} function, or a list containing two or
#'   more \emph{ad hoc} functions, to be applied to a set of shapes stored in
#'   "two-dimensional" format. The function/s is/are applied along the first
#'   margin (i.e. to each row) of the set of shapes, and must return a single
#'   numeric value from each. If two or more \emph{ad hoc} functions are
#'   provided, the individual landscapes are combined into an optimality
#'   trade-off landscape (see \code{details}).
#' @param X An optional vector or matrix containing the values assigned to each
#'   shape (vector length / number of rows and their order must match those from
#'   the shapes provided in \code{shapes} or from the background shape models,
#'   depending on whether or not \code{shapes} have been provided). If a matrix
#'   with two or more columns (representing functional metrics) is provided, the
#'   individual landscapes are combined into an optimality trade-off landscape
#'   (see \code{details}).
#' @param opt.increases Logical vector of length \code{o} (where \code{o} is the
#'   number of \emph{ad hoc} functions or variables provided in \code{FUN} or
#'   \code{X}, respectively) indicating optimality direction of each performance
#'   variable. Ignored unless projection of an optimality trade-off landscape is
#'   attempted.
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
#' @param lty Integer; type of the lines depicting contours.
#' @param lwd Integer; width of the lines depicting contours.
#' @param spar Numeric; smoothing parameter used to smooth landscape outline in
#'   univariate representations.
#' @param pipe Logical; is the function being included in a pipe?
#' @param ... Further arguments passed to \code{FUN}.
#'
#' @details The purpose of this function is to generate and depict a 2- (for
#'   univariate morphospaces) or 3-dimensional (for bivariate morphospaces)
#'   surface (i.e., a landscape), interpolated from values assigned to the set
#'   of shapes projected into an existing morphospace. These can be a sample of
#'   shapes specified by the user, producing a surface for a specific region of
#'   the morphospace ("empirical landscapes"). Alternatively, the set of
#'   background shape models can be used to generate a surface for the entire
#'   morphospace ("theoretical morphospace"). Generally, the values that are
#'   interpolated will represent a variable measuring functional performance
#'   (although it can be any kind of continuous variable), and can be either
#'   provided directly through the \code{X} argument or computed automatically
#'   using an \emph{ad hoc} function through the \code{FUN} argument.
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
#'   If more than one \emph{ad hoc} function or set of values are provided,
#'   the multiple landscapes generated are combined into a single optimality
#'   trade-off landscape by computing their Pareto rank ratio, following the
#'   procedures described in Deakin et al. (2022).
#'
#'   Alternatively, a landscape generated using \code{Morphoscape} can be
#'   provided through the argument \code{obj}, to be projected into morphospace.
#'
#' @return If a plot device with a morphospace is open, the landscape surface is
#'   projected into it as a contour map using [akima::interp()] (alternatively,
#'   an object containing a landscape created with  \code{Morphoscape} can be
#'   projected). If \code{pipe = FALSE}, a list containing the x, y and z values
#'   used to plot the landscape (x and z for univariate morphospaces) is
#'   returned invisibly. If \code{pipe = TRUE} the supplied \code{"mspace"}
#'   object will be modified by appending a \code{$landsc} slot to
#'   \code{$projected}, as well as by adding some graphical parameters (stored
#'   into the \code{$plotinfo} slot), and returned invisibly. If two or more
#'   variables or functions are provided, additional \code{$pfront} and
#'   \code{$opt.increases} slots will be appended to \code{$projected}.
#'
#' @seealso \code{\link{morphogrid}}, \code{\link{mspace}},
#'   \code{\link{plot_morphogrid2d}}, \code{\link{plot_morphogrid3d}},
#'   \code{\link{extract_shapes}}, \code{\link{pareto_rank_ratio}},
#'   \code{\link{proj_pfront}}, \code{\link[akima]{interp}}
#'
#' @references
#' Deakin, W. J., Anderson, P. S., den Boer, W., Smith, T. J., Hill, J. J.,
#'   Rücklin, M., Donoghue, P. C. J.  & Rayfield, E. J. (2022). \emph{Increasing
#'   morphological disparity and decreasing optimality for jaw speed and
#'   strength during the radiation of jawed vertebrates}. Science Advances,
#'   8(11), eabl3644.
#' Akima, H., Gebhardt, A. (2022). \emph{akima: interpolation of irregularly and
#'   regularly spaced data}. <https://CRAN.R-project.org/package=akima>
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
#' ###"Theoretical" landscapes
#'
#' #plot morphospace with its associated adaptive landscape
#' mspace(shapes, links = links, nh = 8, nv = 8, size.model = 1.5,
#'        cex.ldm = 0) %>%
#'   proj_shapes(shapes, pch = 16) %>%
#'   proj_landscape(nlevels = 60, FUN = morphospace:::computeLD, expand = 1.2,
#'                  lwd = 2, display = "contour")
#'
#' ##Using the X argument
#'
#' #first, create morphospace and extract background shapes without plotting
#' msp <- mspace(shapes, links = links, nh = 8, nv = 8, size.model = 1.5,
#'               cex.ldm = 0, plot = FALSE)
#' shapemodels2d <- two.d.array(extract_shapes(msp, keep.template = FALSE)$shapes)
#'
#' #run computeLD through the "two-dimensional" matrix of shapes models (this is
#' #the same thing the proj_landscape function is doing internally when FUN is
#' #used, but this vector could be replaced with some other variable obtained in
#' #a different way)
#' LDs <- apply(X = shapemodels2d, FUN = morphospace:::computeLD, MARGIN = 1)
#'
#' #second, plot morphospace with its associated adaptive landscape
#' msp <- mspace(shapes, links = links, nh = 8, nv = 8, size.model = 1.5,
#'               cex.ldm = 0) %>%
#'   proj_shapes(shapes, pch = 16) %>%
#'   proj_landscape(X = LDs, nlevels = 20, expand = 1.2,
#'                  palette = terrain.colors, display = "filled.contour")
#'
#' #if the shapes argument is specified and it matches the values provided in X,
#' #the landscape becomes independent of the background shapes being displayed.
#' msp <- mspace(shapes, links = links, nh = 7, nv = 5, size.model = 1,
#'               cex.ldm = 0) %>%
#'   proj_shapes(shapes, pch = 16) %>%
#'   proj_landscape(shapes = shapemodels2d, X = LDs, nlevels = 50, expand = 1.2,
#'                  palette = terrain.colors, display = "filled.contour")
#'
#' #in this case, the edges of the landscape can be expanded if the initial sample
#' #of shapes covers a larger area of the morphospace
#' msp <- mspace(shapes, links = links, nh = 8, nv = 8, plot = FALSE,
#'               xlim = c(-.6, .2), ylim = c(-.3, .2))
#' shapemodels2d <- two.d.array(extract_shapes(msp, keep.template = FALSE)$shapes)
#' LDs <- apply(X = shapemodels2d, FUN = morphospace:::computeLD, MARGIN = 1)
#'
#' msp <- mspace(shapes, links = links, nh = 7, nv = 5, size.model = 1,
#'               cex.ldm = 0) %>%
#'   proj_shapes(shapes, pch = 16) %>%
#'   proj_landscape(shapes = shapemodels2d, X = LDs, nlevels = 50, expand = 1.2,
#'                  palette = terrain.colors, display = "filled.contour")
#'
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
#'
#' ###"Empirical" landscapes
#'
#' #it is essentially the same, but providing a set of shapes with the shapes
#' #argument of proj_landscape. Let's compute it it only for Tyrannus species with
#' #non-deep forked tail shapes:
#' mspace(shapes, links = links, nh = 8, nv = 8, size.model = 1.5,
#'        cex.ldm = 0) %>%
#'   proj_shapes(shapes, pch = 16) %>%
#'   proj_landscape(shapes = shapes[,,type == "NDF"], nlevels = 20,
#'                  linear = TRUE, FUN = morphospace:::computeLD, expand = 1.2,
#'                  display = "filled.contour")
#'
#' #in this case, the resolution of the projected surface can be improved using
#' #the argument resolution (however, be aware this can be computationally
#' #intense! especially if a theoretical landscape is being computed)
#' mspace(shapes, links = links, nh = 8, nv = 8, size.model = 1.5,
#'        cex.ldm = 0) %>%
#'   proj_shapes(shapes, pch = 16) %>%
#'   proj_landscape(shapes = shapes[,,type == "NDF"], nlevels = 20,
#'                  linear = TRUE, resolution = 200,
#'                  FUN = morphospace:::computeLD, expand = 1.2,
#'                  display = "filled.contour")
#'
#'
#' ##Pareto landscapes
#'
#' #Pareto landscapes are created when providing two or more separate surfaces,
#' #either using a list fed to FUN containing multuple functions, or a matrix with
#' #multiple columns.
#' msp <- mspace(shapes, links = links, nh = 8, nv = 8, plot = FALSE)
#' shapemodels2d <- two.d.array(extract_shapes(msp, keep.template = FALSE)$shapes)
#' LDs <- apply(X = shapemodels2d, FUN = morphospace:::computeLD, MARGIN = 1)
#' MDs <- apply(X = shapemodels2d, FUN = morphospace:::computeMD, MARGIN = 1)
#'
#' #visualize lift/drag landscape (aerodynamic efficiency)
#' mspace(shapes, links = links, nh = 8, nv = 8, size.model = 1.5,
#'        cex.ldm = 0) %>%
#'   proj_landscape(shapes = shapemodels2d, X = LDs,
#'                  nlevels = 20, expand = 1.2, lwd = 3,
#'                  palette = terrain.colors, display = "contour")
#'
#' #visualize moment/drag landscape (maneuverability)
#' mspace(shapes, links = links, nh = 8, nv = 8, size.model = 1.5,
#'        cex.ldm = 0) %>%
#'   proj_landscape(shapes = shapemodels2d, X = MDs,
#'                  nlevels = 20, expand = 1.2, lwd = 3,
#'                  palette = terrain.colors, display = "contour")
#'
#' #visualize combined Pareto landscape (trade-off)
#' mspace(shapes, links = links, nh = 8, nv = 8, size.model = 1.5,
#'        cex.ldm = 0) %>%
#'   proj_landscape(shapes = shapemodels2d, X = cbind(LDs, MDs),
#'                  opt.increases = c(TRUE, TRUE),
#'                  nlevels = 20, expand = 1.2, lwd = 3,
#'                  palette = terrain.colors, display = "contour")
#'
#'
#' #Integration with Morphoscape
#'
#' #Morphoscape objects can be fed to proj_landscape, meaning that
#' #functionalities such as kriging or weighted combinations of surfaces can be
#' #used to create landscapes. For example, let's apply kriging to the Pareto
#' #surface plotted above.
#' require(Morphoscape)
#'
#' #first, let's calculate the Pareto surface (Pareto rank ratio) for LD and MD
#' pareto.vals <- pareto_rank_ratio(performance = cbind(LDs, MDs),
#'                                  opt.increases = c(TRUE, TRUE))$PRR
#'
#' #get scores for background shape models
#' scores <- proj_eigen(shapemodels2d,
#'                      vectors = msp$ordination$rotation,
#'                      center = msp$ordination$center)
#'
#' #store scores and Pareto values in Morphoscape format
#' dat <- data.frame(x = scores[,1],
#'                   y = scores[,2],
#'                   flight = c(pareto.vals))
#' dat_fnc <- as_fnc_df(dat, func.names = c("flight"))
#'
#' # Create alpha-hulled grid for kriging
#' grid <- resample_grid(dat_fnc, hull = "concaveman", alpha = 3)
#' kr_surf <- krige_surf(dat_fnc, grid = grid)
#' plot(kr_surf)
#'
#' #visualise landscape in morphospace using the argument 'obj'
#' mspace(shapes, links = links, nh = 8, nv = 8, size.model = 1.5,
#'        cex.ldm = 0) %>%
#'   proj_landscape(obj = kr_surf,
#'                  nlevels = 20, expand = 1.2, lwd = 3,
#'                  palette = terrain.colors, display = "contour")
proj_landscape <- function(mspace, shapes = NULL, obj = NULL, FUN = NULL, X = NULL, linear = FALSE,
                           resolution = 50, expand = 1, display = "contour", nlevels = 50,
                           palette = grDevices::heat.colors, alpha = 0.5, lwd = 1, lty = 1,
                           drawlabels = FALSE, spar = 0.5, pipe = TRUE, opt.increases = NULL, ...) {

  args <- c(as.list(environment()), list(...))

  if(is.null(obj)) {
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

    gridcoords <- proj_eigen(x = data2d,
                             vectors = mspace$ordination$rotation[, mspace$plotinfo$axes],
                             center = mspace$ordination$center)

    pfront <- NULL

    if(is.null(X)) {

      if(is.list(FUN) & length(FUN) > 1) {

        XX <- vector(mode = "list", length = length(FUN))
        for(i in seq_len(length(FUN)))
          XX[[i]] <- apply(X = data2d, FUN = FUN[[i]], MARGIN = 1, ...)

        XX <- abind::abind(XX, along = 2)

        if(is.null(opt.increases)) stop("Provide values for the 'opt.increases' argument")
        Pareto <- pareto_rank_ratio(performance = XX, opt.increases = opt.increases)
        X <- c(Pareto$PRR)
        cat("\nTwo or more functions have been provided; combining individual landscapes into an optimality trade-off landscape")

        type <- c(type, "Pareto optimality")
        pfront <- Pareto$pfront

      } else X <- apply(X = data2d, FUN = FUN, MARGIN = 1, ...)
    } else {

      if(is.null(dim(X))) {

        if(length(X) != nrow(data2d)) stop("The amount of values in X do not match the number of shapes")

      } else {
        if(nrow(X) != nrow(data2d)) stop("The number of rows of X do not match the number of shapes")

        if(ncol(X) > 1) {
          if(is.null(opt.increases)) stop("Provide values for the 'opt.increases' argument")
          Pareto <- pareto_rank_ratio(performance = X, opt.increases = opt.increases)
          X <- c(Pareto$PRR)
          cat("\nTwo or more variables have been provided; combining individual landscapes into an optimality trade-off landscape")

          type <- c(type, "Pareto optimality")
          pfront <- Pareto$pfront
        }
      }
    }



    if("theoretical" %in% type) {
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


    if("empirical" %in% type & length(mspace$plotinfo$axes) == 1) {
      spline <- stats::smooth.spline(x = landscape$x[!is.na(landscape$z)],
                                     y = landscape$z[!is.na(landscape$z)], spar = spar)
      smoothed_landsc <- Morpho::equidistantCurve(x = cbind(spline$x, spline$y), n = resolution)

      landscape$x <- smoothed_landsc[,1]
      landscape$z <- smoothed_landsc[,2]
    }

  } else {

    if(!is.null(shapes) | !is.null(X) | !is.null(FUN)) warning("'Morphoscape' object provided, overriding 'shapes', 'X' and 'FUN' arguments")

    landscape <- adapt_Morphoscape(obj)
    extrap <- TRUE
    type <- "Morphoscape"


    frame <- list(xlim = graphics::par("usr")[1:2], ylim = graphics::par("usr")[3:4])
    xlim <- frame$xlim * expand
    ylim <- frame$ylim * expand
    xo <- seq(from = xlim[1], to = xlim[2], length.out = resolution)
    yo <- seq(from = ylim[1], to = ylim[2], length.out = resolution)

    xo[which.x0 <- which.min(abs(xo))] <- 0
    yo[which.y0 <- which.min(abs(yo))] <- 0
  }

  if(length(mspace$plotinfo$axes) > 1) {
    landscape$which.y0 <- which.y0
    landscape$which.x0 <- which.x0
  }


  if(mspace$plotinfo$plot) {

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
  if("Pareto optimality" %in% type) {
    mspace$projected$opt.increases <- opt.increases
    mspace$projected$pfront <- data2d[pfront,]
  }

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

#' Project Pareto front optimizing multiple functions
#'
#' @description Project the Pareto front, representing the subset of shapes that
#'   achieve the optimal balance between two antagonistic functional variables,
#'   into morphospace.
#'
#' @param mspace An \code{"mspace"} object.
#' @param shapes Optional shape data. If provided, the associated variables
#'   must be specified using \code{X}.
#' @param X An optional matrix with at least 2 columns, containing the values
#'   assigned to each shape (if provided) in the functional variables the Pareto
#'   is calculated from (vector length / number of rows and their order must
#'   match those from the shapes provided in \code{shapes}).
#' @param opt.increases Logical vector of length \code{o} (where \code{o} is the
#'   number of \emph{ad hoc} functions or variables provided in \code{FUN} or
#'   \code{X}, respectively) indicating optimality direction of each performance
#'   variable. Ignored unless projection of an optimality trade-off landscape is
#'   attempted.
#' @param lwd Numerical, specifying the width of the line used to mark
#'   the Pareto front.
#' @param pipe Logical; is the function being included in a pipe?
#' @param ... Further arguments passed to [graphics::lines()].
#'
#' @details This function projects the Pareto front, formed by the subset of
#'   shapes that achieve a better balance between two or more functional metrics
#'   engaged in a trade-off that cannot be maximized independently, following
#'   Deakin et al. (2022)
#'
#' @return If a plot device with a morphospace is open, the Pareto front is
#'   projected into it. If \code{pipe = FALSE}, a matrix containing the shapes
#'   forming the Pareto front in its "two dimensional" format is returned
#'   invisibly. If \code{pipe = TRUE} the supplied \code{"mspace"} object will
#'   be modified by appending \code{$pfront} and \code{$optimality} slots to
#'   \code{$projected}, as well as by adding some graphical parameters (stored
#'   into the \code{$plotinfo} slot), and returned invisibly.
#'
#' @seealso \code{\link{pareto_front}}, \code{\link{proj_landscape}}
#'
#' @references
#' Deakin, W. J., Anderson, P. S., den Boer, W., Smith, T. J., Hill, J. J.,
#'   Rücklin, M., Donoghue, P. C. J.  & Rayfield, E. J. (2022). \emph{Increasing
#'   morphological disparity and decreasing optimality for jaw speed and
#'   strength during the radiation of jawed vertebrates}. Science Advances,
#'   8(11), eabl3644.
#'
#' @export
#'
#' @examples
#' #load data and packages
#' library(geomorph)
#'
#' data("tails")
#' shapes <- tails$shapes
#' links <- tails$links
#'
#'
#' #Pareto landscapes are created when providing two or more separate surfaces,
#' #either using a list fed to FUN containing multuple functions, or a matrix with
#' #multiple columns.
#' msp <- mspace(shapes, links = links, nh = 8, nv = 8, plot = FALSE)
#' shapemodels2d <- two.d.array(extract_shapes(msp, keep.template = FALSE)$shapes)
#' LDs <- apply(X = shapemodels2d, FUN = morphospace:::computeLD, MARGIN = 1)
#'
#' #create fake functional variable with the opposite pattern to LD
#' antiLDs <- jitter(rev(LDs), factor = 10)
#'
#' #visualize lift/drag landscape (aerodynamic efficiency)
#' mspace(shapes, links = links, nh = 8, nv = 8, size.model = 1.5,
#'        cex.ldm = 0) %>%
#'   proj_landscape(shapes = shapemodels2d, X = LDs,
#'                  nlevels = 20, expand = 1.2, lwd = 3,
#'                  palette = terrain.colors, display = "contour")
#'
#' #visualize antiLD surface
#' mspace(shapes, links = links, nh = 8, nv = 8, size.model = 1.5,
#'        cex.ldm = 0) %>%
#'   proj_landscape(shapes = shapemodels2d, X = antiLDs,
#'                  nlevels = 20, expand = 1.2, lwd = 3,
#'                  palette = terrain.colors, display = "contour")
#'
#' #visualize combined Pareto landscape (trade-off)
#' mspace(shapes, links = links, nh = 8, nv = 8, size.model = 1.5,
#'        cex.ldm = 0) %>%
#'   proj_landscape(shapes = shapemodels2d, X = cbind(LDs, antiLDs),
#'                  opt.increases = c(TRUE, TRUE),
#'                  nlevels = 50, expand = 1.2, lwd = 3,
#'                  palette = terrain.colors, display = "contour")
#'
#'
#' #there are two ways in which proj_pfront can be used:
#'
#' #1. downstream to a proj_landscape instance in which a Pareto
#' #surface was computed
#' mspace(shapes, links = links, nh = 8, nv = 8, size.model = 1.5,
#'        cex.ldm = 0) %>%
#'   proj_landscape(shapes = shapemodels2d, X = cbind(LDs, antiLDs),
#'                  opt.increases = c(TRUE, TRUE),
#'                  nlevels = 50, expand = 1.2, lwd = 3,
#'                  palette = terrain.colors, display = "contour") %>%
#'   proj_shapes(shapemodels2d) %>%
#'   proj_pfront(col = "red")
#'
#' #2. without previous computation of a Pareto surface. In this case,
#' #the shapes, functional values and directions of optimality must be
#' #provided
#' mspace(shapes, links = links, nh = 8, nv = 8, size.model = 1.5,
#'        cex.ldm = 0) %>%
#'   proj_pfront(shapes = shapemodels2d, X = cbind(LDs, antiLDs),
#'               opt.increases = c(TRUE, TRUE), bg = "red", type = "b",
#'               pch = 21)
proj_pfront <- function(mspace, shapes = NULL, X = NULL, opt.increases = NULL,
                        lwd = 2, pipe = TRUE, ...) {


  args <- c(as.list(environment()), list(...))

  if(is.null(X) & is.null(opt.increases)) {

    if(is.null(mspace$projected$pfront)) {
      stop("Either provide values for 'X' and 'opt.increases' or project front downstream to a Pareto optimality surface")
    } else {
      pfront.scores <- proj_eigen(x = mspace$projected$pfront,
                                  vectors = mspace$ordination$rotation,
                                  center = mspace$ordination$center)
    }

  } else {

    dat <- shapes_mat(shapes)
    data2d <- dat$data2d
    datype <- dat$datype

    pfront <- pareto_front(performance = X, opt.increases = opt.increases)
    pfront.scores <- proj_eigen(x = data2d[pfront,],
                                vectors = mspace$ordination$rotation,
                                center = mspace$ordination$center)

  }


  if(mspace$plotinfo$plot) {

    if(length(mspace$plotinfo$axes) > 1) {
      graphics::lines(pfront.scores, lwd = lwd, ...)

    } else {
      cat("\nPareto fronts are not projected into univariate morphospaces")
    }
  }

  if(!is.null(X) & !is.null(opt.increases)) {
    mspace$projected$opt.increases <- opt.increases
    mspace$projected$pfront <- data2d[pfront,]
  }

  mspace$plotinfo$lwd.pfront <- args$lwd
  mspace$plotinfo$lty.pfront <- args$lty
  mspace$plotinfo$col.pfront <- args$col

  if(!pipe) {
    return(invisible(pfront))
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
#' @param ... Further arguments passed to or from other methods.
#'
#' @noRd
#' @export
print.mspace <- function(mspace, ...) {

  if(any(c("prcomp", "mvgls.pca", "PCA") %in% mspace$ordination$ordtype)) ordtype <- "Principal Component Analysis"
  if(mspace$ordination$ordtype == "gm.prcomp") {
    ord <- adapt_ordination(mspace$ordination)
    if(is.null(ord$phy)) ordtype <- "Principal Component Analysis" else {
      if(ord$alignment == "principal") ordtype <- "Phylogenetic Principal Component Analysis"
      if(ord$alignment == "phy") ordtype <- "Phylogenetically aligned Component Analysis"
    }
  }
  if(any(c("bg_prcomp", "bgPCA") %in% mspace$ordination$ordtype)) ordtype <- "Between-Groups Principal Component Analysis"
  if(any(c("phy_prcomp", "phyl.pca") %in% mspace$ordination$ordtype)) ordtype <- "Phylogenetic Principal Component Analysis"
  if(any(c("pls_shapes", "pls2b", "pls2B", "pls") %in% mspace$ordination$ordtype)) ordtype <- "Two-Block Partial Least Squares"
  if(any(c("phy_pls_shapes", "phy_pls_shapes") %in% mspace$ordination$ordtype)) ordtype <- "Phylogenetic Two-Block Partial Least Squares"
  if(any(c("burnaby", "phy_burnaby") %in% mspace$ordination$ordtype)) ordtype <- "Burnaby's approach"


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

  cat("An mspace object containing")
  cat(paste0("\n* an ordination space consisting of ", ncol(mspace$ordination$rotation),
             nax, ", built using:"))
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
    cat(paste0("\n   - a set of shapes for ",
               mspace$projected$phylo$Nnode + 1, " tips and ", mspace$projected$phylo$Nnode, " nodes"))
    cat(paste0("\n   - a phylogenetic tree"))
    cat(paste0("\n   - an evolutionary model (", mspace$projected$phylo_evmodel,
               ") used for ancestral character estimation"))

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
#'   \code{\link{mspace}} %>% \code{proj_*} pipeline.
#' @param axes Numeric of length 1 or 2, indicating the axes to be plotted.
#' @param links A list with the indices of the coordinates defining a wireframe
#'   connecting landmarks (following the format used in \code{Morpho}), or a
#'   2-columns matrix indicating the pairs of landmarks that should be linked
#'    (following the format used in \code{geomorph}).
#' @param template Either a 2-column matrix with landmarks/semilandmarks and
#'   template curves coordinates (for 2D shape data) or a \code{"mesh3d"}
#'   object representing the mean shape of the sample (for 3D shape data).
#' @param x,y Optional vector with a non-morphometric variable (numerical or
#'   categorical) to be plotted in the x or y axis against an ordination axis.
#'   Alternatively, a \code{"phylo"} object can be provided.
#' @param nh,nv Positive integers; the number of shape models along the x
#'   (\code{nh}) and the y (\code{nv}) axes.
#' @param mag Numeric; magnifying factor for shape models.
#' @param rescale Logical; whether to re-scale background shape models so shape
#'    variation is shown more clearly.
#' @param invax Optional numeric indicating which of the axes provided in
#'   \code{axes} needs to be inverted (options are \code{1}, \code{2} or
#'   \code{c(1,2)}).
#' @param adj_frame Numeric of length 2, providing \emph{a posteriori} scaling
#'   factors for the width and height of the frame, respectively.
#'
#'
#' @param models Logical; whether to plot background shape models (stored
#'      in \code{mspace$projected$shapemodels}).
#' @param size.models Numeric; size factor for shape models.
#' @param col.models Color for wireframes/outlines of shape models.
#' @param bg.models Background color for outlines/meshes of shape models.
#' @param asp.models Numeric; y/x aspect ratio of shape models.
#' @param rot.models Numeric; angle (in degrees) to rotate shape models.
#' @param lwd.models Integer; width of border lines in wireframes/outlines
#'      of shape models.
#' @param alpha.models Numeric; transparency factor for background models
#'      (3D only).
#' @param cex.ldm Numeric; size of landmarks/semilandmarks in the background
#'      models.
#' @param col.ldm Color of landmarks/semilandmarks in the background
#'      models.
#'
#'
#' @param points Logical; whether to plot the scatter points corresponding
#'      to the sampled shapes stored in \code{mspace$projected$scores}.
#' @param pch.points Symbol of the scatterpoints.
#' @param col.points Color of the scatterpoints.
#' @param bg.points Background color of the scatterpoints.
#' @param cex.points Numeric; size of the scatterpoints
#' @param density.points Logical; whether to add density distribution for
#'      points (univariate ordinations only). Overriden by
#'      \code{density.groups = TRUE}
#' @param labels.points Either logical, indicating whether to include all point
#'      labels, or a character string with the exact names of the points whose
#'      labels should be included.
#'
#'
#' @param groups Logical; whether to plot the convex hulls/confidence
#'      ellipses enclosing the groups stored in
#'      \code{mspace$projected$gr_class}.
#' @param col.groups Color of the hulls/ellipses and/or scatterpoints
#'      corresponding to groups' mean shapes.
#' @param bg.groups Background color of the scatterpoints corresponding to
#'      groups' mean shapes.
#' @param pch.groups Symbol of the scatterpoints corresponding to
#'      groups' mean shapes.
#' @param cex.groups Numeric; size of the scatterpoints corresponding to
#'      groups' mean shapes.
#' @param ellipse.groups Logical; whether to use confidence ellipses to
#'      delimit groups (if \code{FALSE} convex hulls are used instead).
#' @param conflev.groups Numeric; confidence level used for confidence
#'      ellipse(s).
#' @param lwd.groups Integer; width of the lines in groups' ellipses/hulls.
#' @param lty.groups Integer; type of the lines in groups' ellipses/hulls.
#' @param alpha.groups Numeric; transparency factor for groups'
#'      ellipses/hulls/density distributions.
#' @param boxplot.groups Logical; whether to plot represent categorical data
#'      using boxplots. If \code{FALSE}, violin plots are displayed instead.
#'      Used only if \code{x} or \code{y} are provided.
#' @param density.groups Logical; whether to add density distribution for
#'      groups (univariate ordinations only).
#' @param legend Logical; whether to show legend for groups.
#' @param cex.legend Numeric; size of legend labels/symbols.
#' @param labels.groups Either logical, indicating whether to include all group
#'      labels, or a character string with the exact names of the groups whose
#'      labels should be included.
#'
#'
#' @param phylo Logical; whether to plot phylogenetic relationships stored
#'      in \code{mspace$projected$phylo}.
#' @param col.phylo Color of the lines depicting phylogenetic
#'      branches.
#' @param lwd.phylo Integer; width of the lines depicting phylogenetic
#'      branches.
#' @param lty.phylo Integer; type of the lines depicting phylogenetic
#'      branches.
#' @param labels.nodes Either logical, indicating whether to include node
#'      labels, or a character string with the exact names of the nodes whose
#'      labels should be included (with the form of "node_n", e.g., "node_14"
#'      corresponds to the root of a tree with 13 tips).
#' @param col.nodes Color of the scatterpoints representing the nodes of the
#'      phylogeny.
#' @param bg.nodes Background color of the scatterpoints representing the
#'      nodes of the phylogeny.
#' @param pch.nodes Symbol of the scatterpoints representing the nodes of
#'      the phylogeny.
#' @param cex.nodes Numeric; size of the scatterpoints representing the
#'      nodes of the phylogeny.
#' @param labels.tips Either logical, indicating whether to include tip labels,
#'      or a character string with the exact names of the tips whose labels
#'      should be included.
#' @param col.tips Color of the scatterpoints representing the tips of the
#'      phylogeny.
#' @param bg.tips Background color of the scatterpoints representing the
#'      tips of the phylogeny.
#' @param pch.tips Symbol of the scatterpoints representing the tips of the
#'      phylogeny.
#' @param cex.tips Numeric; size of the scatterpoints representing the tips
#'      of the phylogeny.
#'
#'
#' @param shapeax Logical; whether to plot morphometric axes stored in
#'      \code{mspace$projected$shape_axis}.
#' @param type.axis Integer; type of arrows (\code{0} = no arrow;
#'      \code{1} = pointing towards the maximum; \code{2} = pointing towards
#'      the maximum, \code{3} = pointing in both directions).
#' @param col.axis Color of the lines depicting a morphometric axis.
#' @param lwd.axis Integer; width of the lines depicting a morphometric
#'      axis.
#' @param lty.axis Integer; type of the lines depicting a morphometric
#'      axis.
#'
#'
#' @param landsc Logical; whether to plot landscape surface stored in
#'      \code{mspace$projected$landsc}.
#' @param display.landsc How to display landscape representation; options
#'      are \code{"contour"} and \code{"filled.contour"}. For bivariate
#'      landscapes only.
#' @param nlevels.landsc Number of levels (i.e., contours) to use in
#'      landscape representation.
#' @param palette.landsc A function defining a color palette to use for
#'      landscape representation.
#' @param alpha.landsc Numeric; transparency factor for filled contours
#'      depicting landscapes.
#' @param lwd.landsc Integer; width of the contour lines depicting
#'      landscapes.
#' @param lty.landsc Integer; type of the contour lines depicting
#'      landscapes.
#' @param drawlabels.landsc Logical; should the labels indicating the value
#'      of each surface contour be plotted?
#' @param scalebar Logical; whether to show scalebar for landscapes.
#'
#'
#' @param xlim,ylim,xlab,ylab,asp Standard arguments passed to the generic
#'   [graphics::plot()] function.
#'
#' @details This function allows to regenerate/tweak morphospaces contained in
#'   \code{"mspace"} objects already in existence. By default, [plot_mspace]
#'   regenerates the morphospace with all its projected elements, preserving
#'   graphical parameters used originally during the \code{\link{mspace}} +
#'   \code{proj_*} pipeline (and stored in \code{mspace$plotinfo}). However,
#'   all the graphical parameters can be modified to customize  representation.
#'   Also, [plot_mspace] can be used to add a legend and/or a scalebar to aid
#'   identification of groups and interpretation of landscapes, respectively.
#'
#'   In addition, this function expands the range of graphical options available
#'   beyond 'pure' morphospaces. If a non-shape variable (assumed to be measured
#'   for the same specimens in \code{mspace$projected$scores}) is fed to one of
#'   the \code{x} or \code{y} arguments, a 'hybrid' morphospace is produced
#'   (i.e. the bivariate plot will be constructed from the combination of
#'   \code{x} or \code{y} and a morphometric axis; background shape models
#'   will represent variation only for the latter)If a \code{"phylo"} object
#'   (assumed to describe the phylogenetic relationships among tips scores
#'   stored in \code{mspace$projected$phylo_scores}) is provided instead for
#'   either \code{x} or \code{y}, a vertical or horizontal phenogram will be
#'   deployed (the function assumes the phylogenetic tree is time-calibrated).
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
#' @return None
#'
#' @export
#'
#' @examples
#' \dontrun{
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
#' #generate basic morphospace, add sampled shapes, species classification,
#' #phylogenetic structure and performance landscape
#' msp <- mspace(shapes, axes = c(1,2), plot = FALSE) %>%
#'  proj_shapes(shapes = shapes) %>%
#'  proj_groups(groups = species) %>%
#'  proj_phylogeny(shapes = sp_shapes, tree = tree) %>%
#'  proj_landscape(FUN = morphospace:::computeLD)
#'
#'
#' ##Regenerating/modifying 'pure' morphospaces:
#'
#' #plot mspace object as it is
#' plot_mspace(msp)
#'
#' #remove landscape
#' plot_mspace(msp, landsc = FALSE)
#'
#' #add colors for points, by species
#' plot_mspace(msp, col.points = species, landsc = FALSE,
#'             col.groups = 1:nlevels(species))
#'
#' #add links for landmark configurations
#' plot_mspace(msp, links = links, landsc = FALSE,
#'             col.points = species, col.groups = 1:nlevels(species))
#'
#' #change number and sizes of shape models in the background
#' plot_mspace(msp, nh = 2, nv = 2, links = links,
#'             size.models = 0.5, landsc = FALSE,
#'             col.points = species, col.groups = 1:nlevels(species))
#'
#' #magnify deformation and highlight landmarks
#' plot_mspace(msp, mag = 1.5, nh = 2, nv = 2, links = links,
#'             size.models = 0.5, col.points = species, landsc = FALSE,
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
#' #plot species agains first PC (violin plot, horizontal)
#' plot_mspace(msp, y = species,  axes = 1, links = links, col.points = species,
#'             col.groups = 1:nlevels(species), pch.points = 16)
#'
#' #plot species agains first PC (boxplot, vertical)
#' plot_mspace(msp, x = species,  axes = 1, links = links, col.points = species,
#'             col.groups = 1:nlevels(species), pch.points = 16)
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
#' }
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
                        asp = 1,
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
                        labels.points,
                        col.groups = 1,
                        bg.groups = 1,
                        pch.groups = 16,
                        cex.groups = 1,
                        ellipse.groups = mspace$plotinfo$ellipse.groups,
                        conflev.groups = 0.95,
                        lwd.groups = 1,
                        lty.groups = 1,
                        alpha.groups = 0,
                        boxplot.groups = FALSE,
                        density.groups = TRUE,
                        labels.groups,
                        col.phylo = 1,
                        lwd.phylo = 1,
                        lty.phylo = 1,
                        labels.nodes,
                        col.nodes,
                        pch.nodes,
                        bg.nodes,
                        cex.nodes,
                        col.tips,
                        pch.tips,
                        bg.tips,
                        labels.tips,
                        cex.tips,
                        type.axis,
                        lwd.axis = 1,
                        lty.axis = 1,
                        col.axis = 1,
                        palette.landsc = grDevices::heat.colors,
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

  args$xname <- if(!is.null(x)) deparse(substitute(x)) else NULL
  args$yname <- if(!is.null(y)) deparse(substitute(y)) else NULL

  oldpar <- graphics::par(no.readonly = TRUE)
  if(any(legend, scalebar)) on.exit(graphics::par(oldpar))

  #1 - prepare data ##############################################################

  #invert orientation of variables and vectors
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

  layout <- set_layout(legend = legend, scalebar = scalebar)
  graphics::par(layout$mainpar)

  if(length(args$axes) == 1) args$asp <- NA

  if(is.null(args$density.points)) args$density.points <- FALSE
  if(is.null(args$density.groups)) args$density.groups <- FALSE

  #2 - if neither x nor y were provided, plot pure morphospace #######################################
  if(is.null(x) & is.null(y)) {

    #2.1 - if a single ordination axes is present, set ground for univariate morphospace
    if(ncol(ordination$x) == 1 | length(args$axes) == 1) {
      y <- rep(1, nrow(ordination$x))
      args$ylim <- c(0, 1)
      args$ylab <- "relative density"
    } else {
      y <- NULL
    }

    #2.2 - plot axes and map morphospace -------------------------------------------
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
                        models = args$models)
    }


    #2.3 - project elements -----------------------------------------------------

    #2.3.1 - add scatterpoints and hulls/ellipses (intercalated if there are
    #groups present)
    if(any(points | groups) &
       any(!is.null(mspace$projected$scores), !is.null(mspace$projected$gr_scores))) {

      scores <- mspace$projected$scores

      if(!is.null(mspace$projected$gr_class)) {
        gr_class <- mspace$projected$gr_class
        gr_scores <- mspace$projected$gr_scores

        if(!is.null(scores)) {
          index_sc_in_gr <- as.numeric(unlist
                                       (apply(gr_scores, 1, \(x, y) {
                                         which(apply(y, 1, \(z, x) {
                                           all(z == x)}, x))}, scores)))
          if(length(index_sc_in_gr) == 0) index_sc_in_gr <- 1:nrow(scores)
        } else index_sc_in_gr <- 1:nrow(scores <- gr_scores)

        mshapes <- apply(X = gr_scores, MARGIN = 2, FUN = tapply, gr_class, mean)
        gcols <- stats::setNames(col2hex(args$col.groups), levels(gr_class))

      } else {
        index_sc_in_gr <- if(!is.null(scores)) -c(1:nrow(scores)) else NULL
        gr_class <- NULL

        mshapes <- NULL
        gcols <- NULL

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

      dens <- NULL

      if(ncol(mspace$ordination$x) > 1) { #if there are more than 1 dimensions
        if(length(args$axes) > 1) { #...and more than 1 dimension is specified
          if(!is.null(scores) & points) {
            plot_biv_scatter(scores = scores[-index_sc_in_gr, args$axes],
                             col = gr_col.points[-index_sc_in_gr],
                             pch = gr_pch.points[-index_sc_in_gr],
                             bg = gr_bg.points[-index_sc_in_gr],
                             cex = gr_cex.points[-index_sc_in_gr])
          }

          if(!is.null(gr_class)) {
            for(i in seq_len(nlevels(gr_class))) {
              if(points) {
                plot_biv_scatter(scores = scores[index_sc_in_gr[gr_class == levels(gr_class)[i]], args$axes],
                                 col = gr_col.points[index_sc_in_gr][gr_class == levels(gr_class)[i]],
                                 pch = gr_pch.points[index_sc_in_gr][gr_class == levels(gr_class)[i]],
                                 bg = gr_bg.points[index_sc_in_gr][gr_class == levels(gr_class)[i]],
                                 cex = gr_cex.points[index_sc_in_gr][gr_class == levels(gr_class)[i]])
              }
              if(groups) {
                if(!args$ellipse.groups) {
                  hulls_by_group_2D(xy = gr_scores[gr_class == levels(gr_class)[i], args$axes],
                                    fac = gr_class[gr_class == levels(gr_class)[i]],
                                    col = gr_col[i], lty = gr_lty[i], lwd = args$lwd.groups,
                                    alpha = args$alpha.groups)
                } else {
                  ellipses_by_group_2D(xy = gr_scores[gr_class == levels(gr_class)[i], args$axes],
                                       fac = factor(gr_class[gr_class == levels(gr_class)[i]]),
                                       col = gr_col[i], lty = gr_lty[i], lwd = args$lwd.groups,
                                       alpha = args$alpha.groups, conflev = args$conflev.groups)
                }
              }
            }
          }

          add_labels(scores[, args$axes], args$labels.points)
          add_labels(xy = mshapes[, args$axes], labels = args$labels.groups, col = gcols)

        } else { #...but less than one dimension is specified
          if(!is.null(gr_class)) { #...and there are groups
            if(args$density.groups) {
              if(args$density.points) {
                args$density.points <- FALSE
                cat("\nBoth density.points and density.groups are TRUE; only the latter will be displayed")
              }
              dens <- density_by_group_2D(gr_scores, gr_class, ax = args$axes[1], alpha = args$alpha.groups,
                                          lwd = args$lwd.groups, lty = args$lty.groups, col = args$col.groups)
            }
            if(points) {
              plot_univ_scatter(scores = cbind(scores[index_sc_in_gr, args$axes[1]]), density = args$density.points,
                                col = gr_col.points[index_sc_in_gr], pch = gr_pch.points[index_sc_in_gr],
                                cex = gr_cex.points[index_sc_in_gr], bg = gr_bg.points[index_sc_in_gr])
            }
          }
          if(points) {
            suppressWarnings(
              plot_univ_scatter(scores = cbind(scores[-index_sc_in_gr, args$axes[1]]), density = args$density.points,
                                col = gr_col.points[-index_sc_in_gr], pch = gr_pch.points[-index_sc_in_gr],
                                cex = gr_cex.points[-index_sc_in_gr], bg = gr_bg.points[-index_sc_in_gr])
            )
          }

          add_labels(cbind(scores[, args$axes[1]], 0), args$labels.points, srt = 90, adj = c(0,0))

          gymax <- lapply(dens$dens, function(x){
            if(is.null(x$y)) x$y <- 0
            max(x$y) / dens$ymax
          })
          gxmax <- lapply(dens$dens, function(x){
            if(!is.null(x$x)) x$x[which.max(x$y)] else NA
          })
          gxmax[is.na(gxmax)] <- mshapes[is.na(gxmax), args$axes[1]]
          xy <- cbind(gxmax, gymax)
          rownames(xy) <- levels(gr_class)

          add_labels(xy = xy, labels = args$labels.groups, col = gcols, pos = 3)
        }
      } else {
        if(!is.null(gr_class)) {
          if(args$density.groups) {
            if(args$density.points) {
              args$density.points <- FALSE
              cat("\nBoth density.points and density.groups are TRUE; only the latter will be displayed")
            }
            dens <- density_by_group_2D(gr_scores, gr_class, ax = 1, alpha = args$alpha.groups,
                                        lwd = args$lwd.groups, lty = args$lty.groups, col = args$col.groups)
          }
          if(points) {
            plot_univ_scatter(scores = cbind(scores[index_sc_in_gr,]), density = args$density.points,
                              col = gr_col.points[index_sc_in_gr], pch = gr_pch.points[index_sc_in_gr],
                              cex = gr_cex.points[index_sc_in_gr], bg = gr_bg.points[index_sc_in_gr])
          }
        }
        if(points) {
          suppressWarnings(
            plot_univ_scatter(scores = cbind(scores[-index_sc_in_gr,]), density = args$density.points,
                              col = gr_col.points[-index_sc_in_gr], pch = gr_pch.points[-index_sc_in_gr],
                              cex = gr_cex.points[-index_sc_in_gr], bg = gr_bg.points[-index_sc_in_gr])
          )
        }

        add_labels(cbind(scores[, 1], 0), args$labels.points, srt = 90, adj = c(0,0))

        gymax <- lapply(dens$dens, function(x){
          if(is.null(x$y)) x$y <- 0
          max(x$y) / dens$ymax
        })
        gxmax <- lapply(dens$dens, function(x){
          if(!is.null(x$x)) x$x[which.max(x$y)] else NA
        })
        gxmax[is.na(gxmax)] <- mshapes[is.na(gxmax), 1]
        xy <- cbind(gxmax, gymax)
        rownames(xy) <- levels(gr_class)

        add_labels(xy = xy, labels = args$labels.groups, col = gcols, pos = 3)
      }
    }

    #2.3.2 - add phylogeny
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

        add_labels(mspace$projected$phylo_scores[mspace$projected$phylo$tip.label,], args$labels.tips)
        if(all(args$labels.nodes %in%
               rownames(mspace$projected$phylo_scores[mspace$projected$phylo$tip.label,])))
          args$labels.nodes <- paste0("node_", ape::getMRCA(phy = mspace$projected$phylo,
                                                            tip = args$labels.nodes))
        add_labels(mspace$projected$phylo_scores[-tips,], args$labels.nodes)
      } else {
        cat("\nPhylogenetic relationships are not projected into univariate morphospaces")
      }
    }

    #2.3.3 - add shape axes
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
        cat("\nMorphometric axes are not projected into univariate morphospaces")
      }
    }

    #2.3.4 - add landscape
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
            cat("\nAxes used to generate landscape are in the wrong order; won't be regenerated")
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
        cat("\nLandscape was originally generated for a different set of axes; won't be regenerated")
      }
    }

  } else {

    #3 - if either x or y has been provided, deploy hybrid morphospace #################################
    if(is.null(axes)) {
      args$axes <- rep(args$axes[1], 2)
    }

    args$xlim <- xlim
    args$ylim <- ylim

    args$asp <- NA

    #3.1 - prepare ground ----------------------------------------------------------------

    #3.1.1 - if either x or y is a factor, prepare ground for violin/box plot

    if(is.factor(x) | is.factor(y)) {

      fac <- if(is.factor(x)) x else if(is.factor(y)) y else NULL

      violins <- vector(mode = "list", length = nlevels(fac))
      for(i in seq_len(nlevels(fac))) {
        d <- density_by_group_2D(xy = cbind(mspace$projected$scores[,args$axes[1]], 0),
                                 fac = fac, ax = 1, plot = FALSE)
        d.width <- (nlevels(fac) + 1) / (nlevels(fac) * 4)

        if(is.factor(x)) {
          poly.i.half1 <- cbind((d$dens[[i]]$y / d$ymax) * d.width, d$dens[[i]]$x)
          poly.i.half2 <- cbind(poly.i.half1[nrow(poly.i.half1):1,1] * -1,
                                poly.i.half1[nrow(poly.i.half1):1,2])
          poly.i <- rbind(poly.i.half1, poly.i.half2)
          poly.i <- cbind(poly.i[,1] + i, poly.i[,2])
        } else {
          if(is.factor(y)) {
            poly.i.half1 <- cbind(d$dens[[i]]$x, (d$dens[[i]]$y / d$ymax) * d.width)
            poly.i.half2 <- cbind(poly.i.half1[nrow(poly.i.half1):1,1],
                                  poly.i.half1[nrow(poly.i.half1):1,2] * -1)
            poly.i <- rbind(poly.i.half1, poly.i.half2)
            poly.i <- cbind(poly.i[,1], poly.i[,2] + i)
          }
        }

        violins[[i]] <- poly.i
      }

      dxvals <- NULL
      dyvals <- NULL
      for(i in 1:length(violins)) {
        dxvals <- c(dxvals, violins[[i]][,1])
        dyvals <- c(dyvals, violins[[i]][,2])
      }

      args$xlim <- range(dxvals)
      args$ylim <- range(dyvals)

    }



    #3.1.2 - if either x or y is a phy object, prepare the ground for phenogram
    if(any(any(class(x) == "phylo"), any(class(y) == "phylo"))) {
      phenogr <- TRUE
    } else {
      phenogr <- FALSE
    }

    if(!is.null(x)) {
      if(any(class(x) == "phylo")) {
        tree <- x
        heights <- phytools::nodeHeights(tree)
        x <- NULL
        for(i in 1:(length(tree$tip.label) + tree$Nnode)) x[i] <- unique(heights[tree$edge == i])
        x <- x - max(x)

        args$xlim <- range(x)
        if(is.null(args$xlab)) args$xlab <- "Time"
      }
      if(is.null(args$xlab)) args$xlab <- deparse(substitute(x))

    } else {
      if(any(class(y) == "phylo")) {
        tree <- y
        heights <- phytools::nodeHeights(tree)
        y <- NULL
        for(i in 1:(length(tree$tip.label) + tree$Nnode)) y[i] <- unique(heights[tree$edge == i])
        y <- y - max(y)

        args$ylim <- range(y)
        if(is.null(args$ylab)) args$ylab <- "Time"
      }
      if(is.null(args$ylab)) args$ylab <- deparse(substitute(y))
    }


    #3.2 - plot hybrid morphospace -----------------------------------------------
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

      rgl::par3d(userMatrix = args$rotation_matrix)

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


    #3.3 - add elements and variables

    #3.3.1 - if either x or y is a phylo object, add phenogram -----------------------------
    if(phenogr) {

      tips <- seq_len(length(tree$tip.label))
      if(length(args$col.tips) == 1) args$col.tips <- rep(args$col.tips, length(tips))
      if(length(args$bg.tips) == 1) args$bg.tips <- rep(args$bg.tips, length(tips))
      if(length(args$pch.tips) == 1) args$pch.tips <- rep(args$pch.tips, length(tips))
      if(length(args$cex.tips) == 1) args$cex.tips <- rep(args$cex.tips, length(tips))

      if(length(args$col.nodes) == 1) args$col.nodes <- rep(args$col.nodes, length(tips) - 1)
      if(length(args$bg.nodes) == 1) args$bg.nodes <- rep(args$bg.nodes, length(tips) - 1)
      if(length(args$pch.nodes) == 1) args$pch.nodes <- rep(args$pch.nodes, length(tips) - 1)
      if(length(args$cex.nodes) == 1) args$cex.nodes <- rep(args$cex.nodes, length(tips) - 1)

      maptips <- if(ncol(mspace$projected$phylo_scores) == 1) {
        order(match(names(mspace$projected$phylo_scores[tips,]), tree$tip.label))
      } else {
        order(match(rownames(mspace$projected$phylo_scores[tips,]), tree$tip.label))
      }

      plot_phenogram(x = x, y = y, tree = tree, phylo_scores = mspace$projected$phylo_scores,
                     axis = args$axes, points = points, labels.tips = args$labels.tips,
                     pch.tips = args$pch.tips[maptips], cex.tips = args$cex.tips[maptips],
                     col.tips = args$col.tips[maptips], bg.tips = args$bg.tips[maptips],
                     labels.nodes = args$labels.nodes, pch.nodes = args$pch.tips[maptips],
                     cex.nodes = args$cex.nodes[maptips], col.nodes = args$col.nodes[maptips],
                     bg.nodes = args$bg.nodes[maptips], lwd.phylo = args$lwd.phylo,
                     lty.phylo = args$lty.phylo, col.phylo = args$col.phylo)

    } else {

      #3.3.2 - if either x or y is a regular variable, add typical elements to
      #hybrid morphospace -----------------------------------------------------------------

      skip.phylo <- FALSE

      #3.3.2.1 - add scatterpoints and hulls / ellipses
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

          meanxy <- apply(X = cbind(x, gr_scores[,args$axes[1]], y), MARGIN = 2,
                          FUN = tapply, gr_class, mean)
          gcols <- stats::setNames(col2hex(args$col.groups), levels(gr_class))

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


        if(!is.null(scores) & points) {
          if((is.factor(x) | is.factor(y)) & boxplot.groups) {} else {
            plot_biv_scatter(scores = xy[-index_sc_in_gr,],
                             col = gr_col.points[-index_sc_in_gr],
                             pch = gr_pch.points[-index_sc_in_gr],
                             bg = gr_bg.points[-index_sc_in_gr],
                             cex = gr_cex.points[-index_sc_in_gr])
          }
        }

        if(!is.null(gr_class)) {
          for(i in seq_len(nlevels(gr_class))) {
            if(points) {
              if((is.factor(x) | is.factor(y)) & boxplot.groups) {} else {
                plot_biv_scatter(scores = xy[index_sc_in_gr[gr_class == levels(gr_class)[i]], ],
                                 col = gr_col.points[index_sc_in_gr][gr_class == levels(gr_class)[i]],
                                 pch = gr_pch.points[index_sc_in_gr][gr_class == levels(gr_class)[i]],
                                 bg = gr_bg.points[index_sc_in_gr][gr_class == levels(gr_class)[i]],
                                 cex = gr_cex.points[index_sc_in_gr][gr_class == levels(gr_class)[i]])
              }
            }
            if(groups) {
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

        add_labels(xy, args$labels.points)
        add_labels(xy = meanxy, labels = args$labels.groups, col = gcols)

        if(is.factor(x) | is.factor(y))  {

          if(is.null(args$col.groups)) args$col.groups <- 1:nlevels(fac)

          if(!boxplot.groups) {
            for(i in 1:length(violins)) {
              graphics::polygon(violins[[i]], lwd = args$lwd.groups,
                                border = args$col.groups[i], lty = args$lty.groups,
                                col = grDevices::adjustcolor(args$col.groups[i], alpha.f = 0.2))

              cents <- tapply(mspace$projected$scores[,args$axes[1]], INDEX = fac, FUN = mean)
              graphics::points(if(is.factor(x)) cbind(i, cents[i]) else  cbind(cents[i], i),
                               pch = 21, bg = args$col.groups[i], cex = args$cex.points + .5)
            }
          } else {
            if(is.null(x)) horizontal <- TRUE else if(is.null(y)) horizontal <- FALSE
            graphics::boxplot(mspace$projected$scores[,args$axes[1]] ~ fac, axes = FALSE,
                              add = TRUE, col = grDevices::adjustcolor(args$col.groups, alpha.f = .5),
                              horizontal = horizontal)
          }

          if(phylo & !is.null(mspace$projected$phylo))
            cat("\nPhylogenetic relationships are not projected into violin/box plots")
          if(shapeax & !is.null(mspace$projected$shape_axis))
            cat("\nMorphometric axes are not projected into violin/box plots")
          if(landsc & !is.null(mspace$projected$landsc))
            cat("\nLandscape surfaces are not projected into violin/box plots")

          skip.phylo <- TRUE
        }
      }

      #3.3.2.2 - add phylogeny
      if(phylo & !is.null(mspace$projected$phylo) & !skip.phylo) {

        tree <- mspace$projected$phylo
        phylo_scores <- mspace$projected$phylo_scores
        evmodel <- mspace$projected$phylo_evmodel

        if(nrow(cbind(x, y)) == nrow(cbind(phylo_scores[tree$tip.label, ]))) {

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

        if(nrow(cbind(x, y)) == nrow(cbind(phylo_scores[tree$tip.label, ]))) {
          xy.tips <- cbind(x, phylo_scores[tree$tip.label, args$axes[1]], y)[tree$tip.label,]

          xy <- cbind(x,y)
          xy_anc <- matrix(NA, nrow = tree$Nnode, ncol = 1)
          colnames(xy_anc) <- if(!is.null(args$xname)) args$xname else args$yname

          xy.int <- as.integer(xy[,1])
          if(inherits(xy[,1], "factor") | all(xy.int == xy[,1])) {
            anc. <- ape::ace(xy[,1], tree, model = "ER", type = "discrete")$lik.anc
            for(j in seq_len(nrow(anc.))) xy_anc[j,1] <- as.numeric(colnames(anc.)[which.max(anc.[j,])])
            cat(paste0("\nReconstruction for discrete character '",
                       colnames(xy_anc)[1], "' was performed assuming an equal-rates model using ape::ace"))
          }

          if(inherits(xy[,1], "numeric") & all(xy.int != xy[,1])) {
            xyname <- colnames(xy)[1]
            xy_ <- cbind(xy[,1],xy[,1])
            colnames(xy_) <- c(xyname, xyname)

            xymod <- mvMORPH::mvgls(xy_ ~ 1, tree = tree, model = evmodel)
            xy_anc[,1] <- mvMORPH::ancestral(xymod)[,1]
          }

          if(!is.null(x)) xy.nodes <- cbind(xy_anc,
                                            phylo_scores[!rownames(phylo_scores) %in% tree$tip.label,
                                                         args$axes[1]])
          if(!is.null(y)) xy.nodes <- cbind(phylo_scores[!rownames(phylo_scores) %in% tree$tip.label,
                                                         args$axes[1]],
                                            xy_anc)
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

          add_labels(xy.tips[tree$tip.label,], args$labels.tips)
          if(all(args$labels.nodes %in% rownames(xy.tips[tree$tip.label,])))
            args$labels.nodes <- paste0("node_", ape::getMRCA(phy = tree, tip = args$labels.nodes))
          add_labels(xy.nodes, args$labels.nodes)
        }
      }
    }
  }

  #4 - add references ###########################################################################


  if(any(legend, scalebar)) {

    #4.1 - add legend for groups --------------------------------------------------------
    if(legend) {

      if(is.null(mspace$projected$gr_class)) {
        stop("Groups levels ($gr_class) are necessary to generate legend labels")
      } else {

        graphics::par(layout$legpar)
        plot(0, type = "n", axes = FALSE, xlab = "", ylab = "",
             ylim = c(0,1), xlim = c(0,1), xaxs = "i", yaxs = "i")

        gr_class <- mspace$projected$gr_class
        gr_scores <- mspace$projected$gr_scores
        scores <- if(!is.null(mspace$projected$scores)) mspace$projected$scores else mspace$projected$gr_scores

        # index_sc_in_gr <- apply(scores, 1, \(x, y) {any(
        #   apply(y, 1, \(z, x) {all(z == x)}, x))}, gr_scores)
        index_sc_in_gr <- sapply(seq_len(nrow(scores)), \(i) {
          any(apply(gr_scores, 1, \(z) all(z == scores[i,])))
        })

        if(length(index_sc_in_gr) == 0) index_sc_in_gr <- 1:nrow(scores)

        # index_gr_in_sc <- as.numeric(unlist
        #                              (apply(scores, 1, \(x, y) {
        #                                which(apply(y, 1, \(z, x){
        #                                  all(z == x)}, x))}, gr_scores)))
        index_gr_in_sc <- sapply(seq_len(nrow(scores)), \(i) {
          match(TRUE, apply(gr_scores, 1, \(z) all(z == scores[i,])), nomatch = NA)
        })

        if(length(index_gr_in_sc) == 0) index_gr_in_sc <- 1:nrow(gr_scores)

        if(length(unique(args$cex.points)) > 1)
          args$cex.points <- rep(1, nrow(cbind(scores[index_sc_in_gr,])))

        pt_params <- NULL
        for(i in seq_len(nrow(cbind(scores[index_sc_in_gr,])))) {
          pt_params <- rbind(pt_params, c(
            as.character(gr_class[index_gr_in_sc][i]),
            args$col.points[index_sc_in_gr][i],
            args$bg.points[index_sc_in_gr][i],
            args$pch.points[index_sc_in_gr][i],
            args$cex.points[index_sc_in_gr][i]))
        }

        if(ncol(pt_params) != 5) pt_params <- NULL

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

    #4.2 - add scalebar for landscape ---------------------------------------------------------
    if(scalebar) {
      if(is.null(mspace$projected$landsc)) {
        stop("Landscape values ($landsc) are necessary to generate a scalebar")
      } else {

        graphics::par(layout$scbpar)
        plot(0, type = "n", axes = FALSE, xlab = "", ylab = "",
             ylim = range(levs.landsc), xlim = c(0.2,0.8), xaxs = "i", yaxs = "i")

        if(args$display.landsc == "filled.contour") {
          alpha.landsc <- args$alpha.landsc
        } else {
          alpha.landsc <- 1
        }

        graphics::rect(0.2, levs.landsc[-length(levs.landsc)], 0.8, levs.landsc[-1L],
                       col = grDevices::adjustcolor(cols.landsc, alpha = alpha.landsc),
                       border = "#00000000")
        ticks <- stats::quantile(levs.landsc, probs = seq(0, 0.95, length.out = 4))
        graphics::axis(side = 4, at = round(ticks, 1))
        graphics::box()
        graphics::mtext("Landscape", side = 4, line = 3)

      }
    }
  }
}


################################################################################

#' Plot \code{"mspace"} objects
#'
#' @description Regenerate morphospaces, generic S3 method (to replace
#'   \code{plot_mspace})
#'
#' @param mspace An \code{"mspace"} object.
#' @param ... Further arguments passed to or from other methods.
#'
#' @noRd
#' @export
plot.mspace <- plot_mspace
