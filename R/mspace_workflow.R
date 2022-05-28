
###########################################################################

#' Generate morphospace
#'
#' @description Create an empirical morphospace as an ordination comprising set
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
#' @param nh Numeric; the number of shape models along the x axis.
#' @param nv Numeric; the number of shape models along the y axis.
#' @param mag Numeric; magnifying factor for shape models.
#' @param invax Optional numeric indicating which of the axes provided in
#'   \code{axes} needs to be inverted (optionsare \code{1}, \code{2} or
#'   \code{c(1,2)}).
#' @param adj_frame Numeric of length 2, providing \emph{a posteriori} scaling
#'   factors for the width and height of the frame, respectively.
#' @param rot.models  Numeric; angle (in degrees) to rotate shape models.
#' @param size.models Numeric; size factor for shape models.
#' @param asp.models Numeric; the y/x aspect ratio of shape models.
#' @param col.models The color for wireframes/outlines.
#' @param bg.models Background color for outlines/meshes.
#' @param lwd.models Numeric; the width of the lines in wireframes/outlines.
#' @param alpha.models Numeric; transparency factor for background models (3D only).
#' @param points Logical; whether to plot the scatter points.
#' @param cex.ldm Numeric; size of landmarks/semilandmarks in the background
#'   models.
#' @param col.ldm The color of landmarks/semilandmarks in the background models.
#' @param plot Logical; whether to plot morphospace.
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
#'   \code{%>%} operator from \code{maggritr}.
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
#' sp_shapes <- consensus(shapes, species)
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
#' meanshape <- consensus(shapes)
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
                   adj_frame = c(1,1),
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
    ordination$x[,axes[invax]] <- ordination$x[,axes[invax]]*-1
    ordination$rotation[,axes[invax]] <- ordination$rotation[,axes[invax]]*-1
  }

  if(k == 3 & datype == "landm") {

    shapemodels <- morphogrid(ordination = ordination, axes = axes, datype = datype, template = NULL,
                              x = NULL, y = NULL, p = p, k = k, nh = nh, nv = nv, mag = mag,
                              asp = asp, xlim = xlim, ylim = ylim, rot.models = rot.models,
                              size.models = size.models, asp.models = 1)

    refshape <- consensus(shapes)
    xlim <- range(ordination$x[,axes[1]])
    ylim <- range(ordination$x[,axes[2]])

    plot_morphogrid3d(x = NULL, y = NULL, morphogrid = shapemodels, refshape = refshape,
                      template = template, links = links, ordtype = ordtype, adj_frame = adj_frame,
                      axes = axes, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab,
                      cex.ldm = cex.ldm, col.ldm = col.ldm, col.models = col.models,
                      lwd.models = lwd.models, bg.models = bg.models, size.models = size.models,
                      asp.models = asp.models, alpha.models = alpha.models, plot = plot)
  } else {

    shapemodels <- morphogrid(ordination = ordination, axes = axes, datype = datype, template = template,
                              x = NULL, y = NULL, p = p, k = k, nh = nh, nv = nv, mag = mag,
                              asp = asp, xlim = xlim, ylim = ylim, rot.models = rot.models,
                              size.models = size.models, asp.models = asp.models)

    plot_morphogrid2d(x = NULL, y = NULL, morphogrid = shapemodels, template = template,
                      links = links, datype = datype, ordtype = ordtype, axes = axes, adj_frame = adj_frame,
                      p = p, xlab = xlab, ylab = ylab, cex.ldm = cex.ldm, col.ldm = col.ldm,
                      col.models = col.models, lwd.models = lwd.models, bg.models = bg.models, plot = plot)
  }

  if(points == TRUE) graphics::points(ordination$x[,axes])


  plotinfo <- list(p = p, k = k, links = links, template = template, axes = axes, nh = nh, nv = nv, mag = mag,
                   asp = asp, adj_frame = adj_frame, asp.models = asp.models, rot.models = rot.models,
                   size.models = size.models, lwd.models = lwd.models, bg.models = bg.models, col.models = col.models,
                   alpha.models = alpha.models, cex.ldm = cex.ldm, col.ldm = col.ldm)

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
#' @param pipe Logical; is the function being included in a pipe?
#' @param ... Further arguments passed to [points()].
#'
#' @details The purpose of this function is to mantain morphospace
#'   generation and sample representation as independent affairs, and
#'   to add flexibility to graphical representation of scatter points.
#'
#' @return If a plot device with a morphospace is open, shapes feeded to
#'   \code{shapes} are projected into morphospace. If \code{pipe = FALSE}
#'   those scores are returned invisibly. If \code{pipe = TRUE} the feeded
#'   \code{mspace} object is returned unchanged and invisibly.
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
proj_shapes <- function(shapes, mspace, pipe = TRUE, ...) {

  dat <- shapes_mat(shapes)
  data2d <- dat$data2d
  datype <- dat$datype

  if(mspace$datype != datype) stop("shapes and mspace types are not compatible")

  scores <- proj_eigen(x = data2d, vectors = mspace$rotation,
                       center = mspace$center)

  if(.Device != "null device") graphics::points(scores[, mspace$plotinfo$axes], ...)

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
#' @param ... Further arguments passed to [points()].
#'
#' @details The purpose of this function is to add the scores corresponding
#'   to groups' mean shapes to \code{mspace} objects during pipeline. Otherwise,
#'   it does the same than \code{proj_shapes}.
#'
#' @return If a plot device with a morphospace is open, shapes feeded to
#'   \code{shapes} are projected into morphospace. If \code{pipe = FALSE} the
#'   corresponding scores are returned invisibly. If \code{pipe = TRUE} the
#'   supplied \code{mspace} object will be modified by adding a new
#'   \code{$gr_centroids} slot, and returned invisibly.
#'
#' @seealso \code{\link{consensus}}
#'
#' @export
#'
#' @examples
#' #load and extract relevant data, packages and information
#' library(magrittr)
#' data("shells")
#' shapes <- shells$shapes
#' species <- shells$data$species
#' sp_shapes <- consensus(shapes, species)
#'
#' #generate basic morphospace, add sampled and consensus shapes
#' mspace(shapes, mag = 0.7, axes = c(1,2), bg.model = "light gray") %>%
#'   proj_shapes(shapes = shapes, col = c(1:4)[species], pch = 1, cex = 0.7) %>%
#'   proj_consensus(shapes = sp_shapes, pch = 21, bg = 1:4, cex = 2)
proj_consensus <- function(shapes, mspace, pipe = TRUE, ...) {

  dat <- shapes_mat(shapes)
  data2d <- dat$data2d
  datype <- dat$datype

  if(mspace$datype != datype) stop("shapes and mspace types are not compatible")

  gr_centroids <- proj_eigen(x = data2d, vectors = mspace$rotation,
                             center = mspace$center)

  if(.Device != "null device") graphics::points(gr_centroids[, mspace$plotinfo$axes], ...)

  mspace$gr_centroids <- gr_centroids

  if(pipe == FALSE) return(invisible(gr_centroids))
  if(pipe == TRUE) return(invisible(mspace))

}


###############################################################################################

#' Delimit groups in morphospace
#'
#' @description Project convex hulls ennclosing \emph{a priori} groups
#'   into an existing morphospace.
#'
#' @param mspace An \code{"mspace"} object.
#' @param shapes Optional shapes data.
#' @param groups Factor; classification of observations into groups.
#' @param pipe Logical; is the function being included in a pipe?
#' @param ... Further arguments passed to [hulls_by_group_2D()].
#'
#' @details The purpose of this function is to add a classification
#'   for shapes populating the morphospace to \code{mspace} objects
#'   during pipeline, as well as to facilitate group visualization.
#'   Otherwise, it is just a wrapper for \code{hulls_by_group_2D}.
#'
#' @return If a plot device with a morphospace is open, convex hulls
#'   enclosing the scores corresponding to \code{groups} are projected
#'   into morphospace. If \code{pipe = TRUE} the supplied \code{mspace}
#'   object will be modified by adding a new \code{$gr_class} slot, and
#'   returned invisibly.
#'
#' @seealso \code{\link{hulls_by_group_2D}}
#'
#' @export
#'
#' @examples
#' #load and extract relevant data, packages and information
#' library(magrittr)
#' data("shells")
#' shapes <- shells$shapes
#' species <- shells$data$species
#' sp_shapes <- consensus(shapes, species)
#'
#' #generate basic morphospace, add sampled shapes and convex hulls for species
#' mspace(shapes, links = links, mag = 0.7, axes = c(1,2), bg.model = "light gray") %>%
#'   proj_shapes(shapes = shapes, col = c(1:4)[species], pch = 1) %>%
#'   proj_groups(groups = species, col = 1:4, lwd = 2)
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

  if(.Device != "null device") hulls_by_group_2D(data2d[, mspace$plotinfo$axes], fac = groups, ...)

  mspace$gr_class <- groups

  if(pipe == TRUE) return(invisible(mspace))
  #if(pipe == FALSE) return(  volume_of_hulls  ) #one for the future

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
#' @param mspace An \code{"mspace"} object.
#' @param axis An optional vector indicating the axis from \code{obj} to be
#'   projected.
#' @param mag Numeric; magnifying factor for representing shape transformation.
#' @param pipe Logical; is the function being included in a pipe?
#' @param ... Further arguments passed to [lines()].
#'
#' @details This function is primarily aimed at the graphical representation
#'   of morphometric axes (estimated using either linear models or multivariate
#'   ordination methods) into an existing morphospace for heuristic exploration
#'   of patterns in the data. It can also be used to extract theoretical shapes
#'   at the extremes of those axes, although \code{ax_transformation()} does the
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
#'   \code{mlm} objects or eigenvector coefficients stored in the \code{$rotation}
#'   slot returned by multivariate ordination methods.
#'
#' @return If a plot device with a morphospace is open, a straight line marking the
#'   scores representing shapes at the extremes of the morphometric axis is projected
#'   into morphospace. If \code{pipe = FALSE} those scores are returned invisibly.
#'    If \code{pipe = TRUE} the feeded \code{mspace} object is returned unchanged
#'    and invisibly.
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
#' #compute intraspecific allometric axis using consensus, tapply, lm and pls_shapes
#' sp_shapes <- consensus(shapes, species)
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
proj_axis <- function(obj, mspace, axis = 1, mag = 1, pipe = TRUE, ...) {

  ext_shapes2d <- ax_transformation(obj = obj, axis = axis, mag = mag)
  ext_scores <- proj_eigen(x = ext_shapes2d, vectors = mspace$rotation,
                           center = mspace$center)

  if(.Device != "null device") graphics::lines(ext_scores[, mspace$plotinfo$axes], ...)

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
#' @param mspace An \code{"mspace"} object.
#' @param pipe Logical; is the function being included in a pipe?
#' @param ... Further arguments passed to [lines()].
#'
#' @details The purpose of this function is twofold. First, it is meant to transform
#'   a morphospace into a phylomorphospace by infusing  phylogenetic structure into
#'   the former. To this end, a \code{$gr_centroids} slot matching the tip labels
#'   from \code{tree} needs to be present (either upstream in the pipeline or already
#'   incorporated into an existing \code{mspace} object). Second, this function can be
#'   used to retrieve the scores corresponding to nodes of the phylogenetic tree, which
#'   can in turn be used to compute the associated shapes using \code{rev_eigen()}. The
#'   position of these shapes in morphospace is estimated using the squared-changes
#'   parsimony algorithm as performed by \code{fastAnc} from \code{phytools}.
#'
#' @return If a plot device with a morphospace is open, shapes representing the nodes
#'   of the phylogenetic tree and lines connecting them and tips are projected into
#'   morphospace. If \code{pipe = FALSE} scores for nodes and tips of the phylogeny
#'   are returned invisibly. If \code{pipe = TRUE} the supplied \code{mspace}
#'   object will be modified by adding a new \code{$phylo_scores} and \code{$phylo}
#'   slots, and returned invisibly.
#'
#' @export
#'
#' @examples
#' #load and extract relevant data, packages and information
#' library(magrittr)
#' data("tails")
#' shapes <- tails$shapes
#' species <- tails$data$species
#' sp_shapes <- consensus(shapes, species)
#' tree <- tails$tree
#' links <- tails$links
#'
#' #generate basic morphospace, add sampled shapes, species mean shapes, and
#' #phylogenetic structure
#' mspace(shapes, links = links, mag = 0.7, axes = c(1,2), cex.ldm = 0) %>%
#'   proj_shapes(shapes = shapes, col = c(1:13)[species], pch = 1, cex = 0.7) %>%
#'   proj_consensus(shapes = sp_shapes, pch = 21, bg = 1:13, cex = 2) %>%
#'   proj_phylogeny(tree = tree)
proj_phylogeny <- function(tree, mspace, pipe = TRUE, ...) {

  if(is.null(mspace$gr_centroids)) stop("Group centroids have not been provided; add proj_consensus() before")

  nodes_scores <- apply(mspace$gr_centroids[tree$tip.label,], 2, phytools::fastAnc, tree = tree)
  phylo_scores <- rbind(mspace$gr_centroids[tree$tip.label,], nodes_scores)


  if(.Device != "null device") {
    for(i in 1:nrow(tree$edge)) {
      graphics::lines(rbind(phylo_scores[tree$edge[i, 1], mspace$plotinfo$axes],
                            phylo_scores[tree$edge[i, 2], mspace$plotinfo$axes]), ...)
    }
  }

  mspace$phylo_scores <- phylo_scores
  mspace$phylo <- tree

  if(pipe == FALSE) {
    return(invisible(phylo_scores))
  } else {
    return(invisible(mspace))
  }

}


#########################################################################################

#' Plot morphospaces and combine them with other elements
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
#'   the x axis. Alternatively, a \code{"phy"} object can be provided.
#' @param y Optional vector with a non-morphometric variable to be plotted in
#'   the y axis. Alternatively, a \code{"phy"} object can be provided.
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
#' @param bg.models Background color for outlines/meshes.
#' @param lwd.models Numeric; the width of the lines in wireframes/outlines.
#' @param alpha.models Numeric; transparency factor for background models (3D only).
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
#'   function.
#'
#' @details This function allows to plot the morphospace contained in \code{mspace}
#'   objects already in existence, either preserving the graphical attributes stamped
#'   at during the pipeline or modifying one or more of them, so there is no need to
#'   execute the pipeline every time morphospaces need to be visualized and/or changed.
#'
#'   Also, \code{plot_mspace} expands the range of graphical options available
#'   beyond 'pure' morphospaces. If a numeric non-shape variable (assumed to have
#'   been measured for the same specimens in \code{mspace$x}) is feeded to one of
#'   \code{x} or \code{y}, a 'hybrid' morphospace is produced (i.e. the bivariate
#'   plot will be constructed from the combination of \code{x} or \code{y} and a
#'   morphometric axis; shape models in the background will represent variation
#'   only for the latter). If instead a \code{"phy"} object (assumed to describe
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
#' sp_shapes <- consensus(shapes, species)
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
                        axes,
                        links = NULL,
                        template = NULL,
                        x = NULL,
                        y = NULL,
                        nh,
                        nv,
                        mag,
                        invax = NULL,
                        adj_frame = c(1,1),
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
                        alpha.models,
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


  if(!is.null(invax)) {
    mspace$x[,axes[invax]] <- mspace$x[,axes[invax]]*-1
    mspace$rotation[,axes[invax]] <- mspace$rotation[,axes[invax]]*-1
    if(!is.null(mspace$gr_centroids)) {
      mspace$gr_centroids[,axes[invax]] <- mspace$gr_centroids[,axes[invax]]*-1
    }
    if(!is.null(mspace$phylo_scores)) {
      mspace$phylo_scores[,axes[invax]] <- mspace$phylo_scores[,axes[invax]]*-1
    }
  }

  ordination <- list(x = mspace$x, rotation = mspace$rotation, center = mspace$center)

  if(!is.null(x) & !is.null(y)) stop("Only one of x or y can be specified")



  if(is.null(x) & is.null(y)) { #if neither x nor y have been provided, plot pure morphospace

     if(mspace$plotinfo$k == 3 & mspace$datype == "landm") {

      shapemodels <- morphogrid(ordination = ordination, axes = args$axes, datype = mspace$datype,
                                template = NULL, x = NULL, y = NULL, p = mspace$plotinfo$p,
                                k = mspace$plotinfo$k, nh = args$nh, nv = args$nv, mag = args$mag,
                                asp = args$asp, xlim = args$xlim, ylim = args$ylim, rot.models = args$rot.models,
                                size.models = args$size.models, asp.models = args$asp.models)

      refshape <- matrix(rev_eigen(0,
                                   ordination$rotation[,1],
                                   ordination$center),
                         nrow = mspace$plotinfo$p, ncol = mspace$plotinfo$k, byrow = TRUE)

      xlim <- range(ordination$x[,args$axes[1]])
      ylim <- range(ordination$x[,args$axes[2]])

      plot_morphogrid3d(x = NULL, y = NULL, morphogrid = shapemodels, refshape = refshape,
                        template = args$template, links = args$links, ordtype = mspace$ordtype,
                        axes = args$axes, xlim = xlim, ylim = ylim, adj_frame = args$adj_frame,
                        xlab = args$xlab, ylab = args$ylab, cex.ldm = args$cex.ldm,
                        col.ldm = args$col.ldm, col.models = args$col.models, lwd.models = args$lwd.models,
                        bg.models = args$bg.models,  size.models = args$size.models,
                        asp.models = args$asp.models, alpha.models = args$alpha.models)
    } else {

      shapemodels <- morphogrid(ordination = ordination, axes = args$axes, datype = mspace$datype,
                                template = args$template, x = NULL, y = NULL, p = mspace$plotinfo$p,
                                k = mspace$plotinfo$k, nh = args$nh, nv = args$nv, mag = args$mag,
                                asp = args$asp, xlim = args$xlim, ylim = args$ylim, rot.models = args$rot.models,
                                size.models = args$size.models, asp.models = args$asp.models)

      plot_morphogrid2d(x = x, y = y, morphogrid = shapemodels, template = args$template,
                        links = args$links, datype = mspace$datype, ordtype = mspace$ordtype,
                        axes = args$axes, adj_frame = args$adj_frame, p = mspace$plotinfo$p,
                        xlab = args$xlab, ylab = args$ylab, cex.ldm = args$cex.ldm, col.ldm = args$col.ldm,
                        col.models = args$col.models, lwd.models = args$lwd.models, bg.models = args$bg.models)
    }

    #add points, hulls, phylogeny, and/or consensus
    if(points == TRUE) graphics::points(mspace$x[, args$axes],
                                        pch = pch.points, col = col.points, cex = cex.points)
    if(groups == TRUE & !is.null(mspace$gr_class)) {
      if(is.null(mspace$gr_class)) {
        stop("groups classification has not been added to mspace object")
      } else {
        hulls_by_group_2D(mspace$x[, args$axes], fac = mspace$gr_class,
                          col = col.groups)
      }
    }
    if(phylo == TRUE & !is.null(mspace$phylo)) {
      if(is.null(mspace$phylo)) {
        stop("phylogenetic relationships have not been added to mspace object")
      } else {
        for(i in 1:nrow(mspace$phylo$edge)) {
          graphics::lines(rbind(mspace$phylo_scores[mspace$phylo$edge[i, 1], args$axes],
                                mspace$phylo_scores[mspace$phylo$edge[i, 2], args$axes]),
                          lwd = lwd.branches)
        }
      }
    }
    if(mshapes == TRUE & !is.null(mspace$gr_centroids)) {
      if(is.null(mspace$gr_centroids)) {
        stop("groups centroids have not been added to mspace object")
      } else {
        graphics::points(mspace$gr_centroids[, args$axes],
                         col = col.groups, pch = pch.groups, cex = cex.groups)
      }
    }



  } else { #if x or y have been provided, show hybrid morphospace

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
               unique(heights[,1]))
        args$xlim <- range(x)
      }
    } else {
      if(any(class(y) == "phylo")) {
        tree <- y
        heights <- phytools::nodeHeights(tree)
        y <- c(rep(max(heights), length(tree$tip.label)),
               unique(heights[,1]))
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
                                   ordination$rotation[,1],
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
                        asp.models = args$asp.models, alpha.models = args$alpha.models)
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
                        col.models = args$col.models, lwd.models = args$lwd.models, bg.models = args$bg.models)
    }


    if(phenogr == TRUE) { #if x/y is a phy object, plot a phenogram

      maptips <-  order(match(rownames(mspace$gr_centroids), tree$tip.label))
      plot_phenogram(x = x, y = y, tree = tree, axes = args$axes, points = points,
                     phylo_scores = mspace$phylo_scores, lwd.branches = lwd.branches,
                     cex.groups = cex.groups, col.groups = col.groups[maptips])

    } else { #else, go for a generic hybrid morphospace

      xy <- cbind(x, mspace$x[,args$axes[1]], y)

      #add points, hulls, phylogeny, and/or consensus
      if(points == TRUE) graphics::points(xy, pch = pch.points,
                                          col = col.points, cex = cex.points)
      if(groups == TRUE & !is.null(mspace$gr_class)) {
        if(is.null(mspace$gr_class)) {
          stop("groups classification has not been added to mspace object")
        } else {
          hulls_by_group_2D(xy, fac = mspace$gr_class, col = col.groups)
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
            nodesxy <- apply(meanxy[tree$tip.label,], 2, phytools::fastAnc, tree = tree)
            phyloxy <- rbind(meanxy, nodesxy)

          }

          for(i in 1:nrow(mspace$phylo$edge)) {
            graphics::lines(rbind(phyloxy[mspace$phylo$edge[i, 1],],
                                  phyloxy[mspace$phylo$edge[i, 2],]),
                            lwd = lwd.branches)
          }
        }
      }

      if(mshapes == TRUE & !is.null(mspace$gr_consensus)) {
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

          graphics::points(meanxy, col = col.groups,
                           pch = pch.groups, cex = cex.groups)
        }
      }
    }
  }
}


