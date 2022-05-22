##########################################################################

#' Reverse eigenanalysis-based ordination
#'
#' @description Back-calculate variables in their original scale from existing
#'   ordination data. Used internally.
#'
#' @param scores A matrix of scores resulting from an ordination.
#' @param vectors A matrix of eigenvector coefficients used to build the
#'   ordination.
#' @param center A vector with the mean values of the original variables used in
#'   the ordination.
#'
#' @return A 2-margins matrix of variables on their original scale.
#'
#' @seealso \code{\link{proj_eigen}}
#' @export
#'
#' @examples
#' #load data and packages
#' library(Morpho)
#' library(geomorph)
#' data("tails")
#'
#' #perform principal component analysis on tails data
#' pca <- prcomp(two.d.array(tails$shapes))
#'
#' #transform the scores back to shapes
#' backshapes_mat <- rev_eigen(scores = pca$x,
#'                             vectors = pca$rotation,
#'                             center = pca$center)
#' backshapes_arr <- arrayspecs(backshapes_mat, k = 2, p = 9)
#'
#' #compare
#' pile_shapes(tails$shapes)
#' pile_shapes(backshapes_arr)
#'
#' #obtain shapes at the extremes of PC1
#' extshapes_mat <- rev_eigen(scores = range(pca$x[,1]),
#'                            vectors = pca$rotation[,1],
#'                            center = pca$center)
#' extshapes_arr <- arrayspecs(extshapes_mat, k = 2, p = 9)
#'
#' #plot and compare
#' plot(extshapes_arr[,,1])
#' lineplot(extshapes_arr[,,1], tails$links) ; title("negative")
#' plot(extshapes_arr[,,2])
#' lineplot(extshapes_arr[,,2], tails$links) ; title("positive")
rev_eigen <- function(scores, vectors, center) { t(t(scores %*% t(vectors)) + center) }


###########################################################################

#' Project cases into existing eigenanalysis-based ordination
#'
#' @description Project variables in their original scale into an existing
#'   ordination. Used internally.
#'
#' @param x A matrix with the original variables used for the ordination as
#'   columns, and observations as rows.
#' @param vectors A matrix of eigenvector coefficients used to build the
#'   ordination.
#' @param center A vector with the mean values of the original variables used in
#'   the ordination.
#'
#' @return A 2-margins matrix with the scores on the supplied vector(s).
#'
#' @export
#'
#' @seealso \code{\link{rev_eigen}}
#'
#' @examples
#' #' #load data and packages
#' library(geomorph)
#' data("tails")
#'
#' #perform principal component analysis on tails data
#' pca <- prcomp(two.d.array(tails$shapes))
#'
#' #get project shapes in the pca ordination
#' shapes_mat <- two.d.array(tails$shapes)
#' newscores <- proj_eigen(x = shapes_mat,
#'                         vectors = pca$rotation,
#'                         center = pca$center)
#'
#' #plot and compare
#' plot(pca$x, pch = 16, col = "gray")
#' points(newscores, col = "black")
#'
#' #compute and project species' mean shapes
#' meanshapes_arr <- consensus(tails$shapes, index = tails$data$species)
#' meanshapes_mat <- two.d.array(meanshapes_arr)
#' meanscores <- proj_eigen(x = meanshapes_mat,
#'                         vectors = pca$rotation,
#'                         center = pca$center)
#' points(meanscores, pch = 21, bg = "red", cex = 1.5)
proj_eigen <- function(x, vectors, center) { t(t(rbind(x)) - center) %*% vectors }


##########################################################################

#' Singular value decomposition for 2 blocks of variables
#'
#' @description Just a wrapper for [svd()] that returns an adequate output
#'   when used for blocks of variables. Can deal with phylogenetic data too
#'   using \code{ape} and \code{phytools} functions internally.
#'
#' @param x First block of variables
#' @param y Second block of variables
#' @param tree An optional \code{"phy"} object containing a phylogenetic
#'   tree whose tip.labels match rownames of \code{x} and \code{y}.
#'
#' @return Mimics the [svd()] output.
#'
#' @export
#'
#' @examples
#' #laod data and packages
#' library(geomorph)
#'
#' data("tails")
#' shapes <- tails$shapes
#' sizes <- tails$sizes
#' species <- tails$data$species
#' tree <- tails$tree
#' sp_shapes <- consensus(shapes, species)[,,tree$tip.label]
#' sp_sizes <- cbind(tapply(sizes, species, mean))[tree$tip.label,]
#'
#' #perform partial svd
#' svd_block(x = sizes, y = two.d.array(shapes))
#'
#' #perform partial svd on phylogenetic structure
#' svd_block(x = sp_sizes, y = two.d.array(sp_shapes), tree = tree)
svd_block <- function(x, y, tree = NULL) {

  x <- cbind(x)
  y <- cbind(y)

  if(!is.null(tree)) {
    C <- ape::vcv.phylo(tree)[rownames(x), rownames(x)]
    xy <- cbind(x,y)
    Clambda <- phytools::phyl.vcv(xy, C, 1)$C

    R <- phytools::phyl.vcv(xy, Clambda, 1)$R
    part_R <- R[1:ncol(x), (ncol(x) + 1):(ncol(x) + ncol(y))]

    svd <- svd(part_R)

  } else {
    vcv <- stats::cov(cbind(x,y))
    part_vcv <- vcv[1:ncol(x), (ncol(x) + 1):(ncol(x) + ncol(y))]

    svd <- svd(part_vcv)


  }

  ndims <- min(nrow(x), ncol(x), ncol(y))
  sdev <- (svd$d)[1:ndims]
  if(ncol(x) != ncol(y)) {
    rotations <- list(svd$v, svd$u)
    whichy <- which(unlist(lapply(rotations, nrow)) == ncol(y))
    whichx <- which(unlist(lapply(rotations, nrow)) == ncol(x))
    y_rotation <- cbind(rotations[[whichy]][,1:ndims])
    x_rotation <- cbind(rotations[[whichx]][,1:ndims])
  } else {
    x_rotation <- svd$u
    y_rotation <- svd$v
  }

  results <- list(d = sdev, v = y_rotation, u = x_rotation)
  return(results)

}
##########################################################################

#' Iverse Fourier transform
#'
#' @description A wrapper for [efourier_i()] from \code{Momocs} to transform a
#'   set of Fourier coefficients into (x,y) coordinates. Used internally.
#'
#' @param coe A vector with Fourier coefficients.
#' @param nb.pts Numeric, specifying the number of coordinates for sampling the
#'   outlines.
#'
#' @return A \code{nb.pts x 2} matrix of (x,y) cartesian coordinates defining
#'   single outline shape.
#'
#' @export
#'
#' @seealso \code{\link[Momocs]{efourier_i}}
#'
#' @examples
#' #load data and extract the first outline
#' data("shells")
#' shape_coe <- shells$shapes$coe[1,]
#'
#' #get and plot (x,y) coordinates for using different number of coordinares
#' shape_xy <- inv_efourier(shape_coe, nb.pts = 10)
#' plot(shape_xy)
#' shape_xy <- inv_efourier(shape_coe, nb.pts = 30)
#' plot(shape_xy)
#' shape_xy <- inv_efourier(shape_coe, nb.pts = 100)
#' plot(shape_xy)
inv_efourier <- function(coe, nb.pts = 120) {

  nb.h <- length(coe) / 4
  coords <- matrix(0, nrow = nb.pts, ncol = 2)

  a <- c(coe[1:nb.h])
  b <- c(coe[(nb.h + 1):(nb.h * 2)])
  c <- c(coe[((nb.h * 2) + 1):(nb.h * 3)])
  d <- c(coe[((nb.h * 3) + 1):(nb.h * 4)])

  coefs <- list(an = a, bn = b, cn = c, dn = d)

  coords <- Momocs::efourier_i(coefs, nb.pts = nb.pts)
  return(coords)

}


################################################################################################

#' Identify and arrange shape descriptors
#'
#' @description Identify data type (Fourier coefficients or Procrustes shape
#'   coordinates) and arrange them in 2-margin matrix format, suitable to be
#'   used as input for multivariate analysis. Used internally.
#'
#' @param shapes Shape data.
#'
#' @return A list of length 2 containing:
#' \itemize{
#'   \item \code{$datype:} {the type of geometric morphometrics data.}
#'   \item \code{$data2d:} {the shape descriptors arranged in 2-margins matrix
#'   format.}
#'  }
#'
#' @export
#'
#' @examples
#' #apply on shells dataset
#' data("shells")
#' dat1 <- shapes_mat(shells$shapes)
#'
#' #inspect results
#' head(dat1$data2d)
#' dat1$datype
#'
#' #aply on tails dataset
#' data("tails")
#' dat2 <- shapes_mat(tails$shapes)
#'
#' #inspect results
#' head(dat2$data2d)
#' dat2$datype
shapes_mat <- function(shapes) {

  if(length(dim(shapes)) == 3) {
    datype <- "landm"
    data2d <- geomorph::two.d.array(shapes)
  }

  if(length(dim(shapes)) == 2) {
    if(any(class(shapes) == "OutCoe")) {
      datype <- "fcoef"
      data2d <- shapes$coe
    } else {
      if(!is.null(colnames(shapes))) {
        if(any(colnames(shapes)[1] == "A1", colnames(shapes)[1] == "A2")) {
          datype <- "fcoef"
          data2d <- shapes
        } else {
          datype <- "landm"
          data2d <- shapes
        }
      } else {
        datype <- "landm"
        data2d <- shapes
      }
    }
  }
  return(list(data2d = data2d, datype = datype))
}


################################################################################################

#' Generate background shape models
#'
#' @description Calculate and arrange background shape models for morphospaces.
#'   For internal use.
#'
#' @param ordination An ordination (i.e. a \code{"prcomp"}, \code{"bg_prcomp"},
#'   \code{"phy_prcomp"} or \code{"pls_shape"} object).
#' @param axes Numeric of length 1 or 2, indicating the morphometric axes to be
#'   plotted. If values for either \code{x} or \code{y} are provided, only the
#'   first value of this argument is considered.
#' @param datype Character; the type of shape data (either \code{"fcoef"} or
#'   \code{"landm"}).
#' @param template A 2-column matrix containing 1) the actual
#'   landmarks/semilandmarks being analized, followed by 2) the (x,y) cartesian
#'   coordinates defining a curve or set of curves that will be warped using the
#'   deformation interpolated from the changes between landmarks/semilandmarks
#'   (the actual positions of which must be marked with a row of NA, see the
#'   tails dataset).
#' @param x Optional vector with a non-morphometric variable to be plotted in
#'   the x axis.
#' @param y Optional vector with a non-morphometric variable to be plotted in
#'   the y axis.
#' @param p Numeric, indicating the number of landmarks/semilandmarks used (for
#'   landmark data only).
#' @param k Numeric, indicating the number of cartesian dimensions of
#'   landmarks/semilandmarks (for landmark data only).
#' @param nh Numeric; the number of shape models along the x axis.
#' @param nv Numeric; the number of shape models along the y axis.
#' @param mag Numeric; magnifying factor for shape models.
#' @param rot.models  Numeric; angle (in degrees) to rotate shape models.
#' @param size.models Numeric; size factor for shape models.
#' @param asp.models Numeric; the y/x aspect ratio of shape models.
#' @param xlim,ylim,asp Standard arguments passed to the generic plot function
#'
#' @return A list of length 2 containing:
#' \itemize{
#'   \item \code{$models_mat:} {a 2-column matrix with the (x,y) coordinates of
#'   all the shape models in the background.}
#'   \item \code{$models_arr:} {same as \code{models_mat} but in 3-margins array
#'   format.}
#' }
#'
#' @export
#'
#' @references MacLeod, N. (2009). \emph{Form & shape models}. Palaeontological
#'   Association Newsletter, 72(620), 14-27.
#'
#' @examples
#' #' #load data and packages
#' library(geomorph)
#' data("tails")
#'
#' #perform pca on tails shapes
#' pca <- prcomp(two.d.array(tails$shapes))
#'
#' #generate grid of shape models sampling the range of variation
#' #at 4 locations (the 4 corners of the scatterplot)
#' shapes_grid <- morphogrid(ordination = pca, axes = c(1,2), datype = "landm",
#'                           k = ncol(tails$shapes), p = nrow(tails$shapes),
#'                           nh = 2, nv = 2)
#'
#' #plot grid from $models_mat and project each shape in models_arr
#' plot(shapes_grid$models_mat)
#' points(shapes_grid$models_arr[,,1], pch=16, col = 1)
#' points(shapes_grid$models_arr[,,2], pch=16, col = 2)
#' points(shapes_grid$models_arr[,,3], pch=16, col = 3)
#' points(shapes_grid$models_arr[,,4], pch=16, col = 4)
morphogrid <- function(ordination,
                       axes,
                       datype,
                       template = NULL,
                       x = NULL,
                       y = NULL,
                       p = NULL,
                       k = NULL,
                       nh = 5,
                       nv = 4,
                       mag = 1,
                       asp = NA,
                       xlim = NULL,
                       ylim = NULL,
                       rot.models = 0,
                       size.models = 1,
                       asp.models = 1) {

  if(is.null(x) & is.null(y)) {
    if(is.null(xlim)){
      plotframe_x <- range(ordination$x[,axes[1]])
    } else {
      plotframe_x <- sort(xlim)
    }
    if(is.null(ylim)){
      plotframe_y <- range(ordination$x[,axes[2]])
    } else {
      plotframe_y <- sort(ylim)
    }
  } else {

    if(length(axes) > 1) warning("x or y has been specified, axes[2] will be ignored")
    axes <- axes[1]

    if(!is.null(x)) {
      if(is.null(xlim)){
        plotframe_x <- range(x)
      } else {
        plotframe_x <- sort(xlim)
      }
    } else {
      if(is.null(xlim)){
        plotframe_x <- range(ordination$x[,axes])
      } else {
        plotframe_x <- sort(xlim)
      }
    }

    if(!is.null(y)) {
      if(is.null(ylim)){
        plotframe_y <- range(y)
      } else {
        plotframe_y <- sort(ylim)
      }
    } else {
      if(is.null(ylim)){
        plotframe_y <- range(ordination$x[,axes])
      } else {
        plotframe_y <- sort(ylim)
      }
    }
  }


  scores_x  <-  seq(from = plotframe_x[1],
                    to   = plotframe_x[2],
                    length.out = nh)
  scores_y <-  seq(from = plotframe_y[1],
                   to   = plotframe_y[2],
                   length.out = nv)
  vectors <-  ordination$rotation[,axes]
  center  <-  ordination$center


  gridcoords <- data.matrix(expand.grid(scores_x, scores_y))
  gridcoords_mag <- gridcoords * mag

  if(!is.null(x)) gridcoords_mag <- gridcoords_mag[,2]
  if(!is.null(y)) gridcoords_mag <- gridcoords_mag[,1]


  sh_mat <- rev_eigen(gridcoords_mag, vectors, center)
  if(datype == "fcoef") {
    coords_l <- lapply(1:nrow(sh_mat), function(i){
      inv_efourier(coe = sh_mat[i,], nb.pts = p)
    })
    sh_arr <- abind::abind(coords_l, along = 3) * size.models
  } else {
    sh_arr <- geomorph::arrayspecs(sh_mat, p = p, k = k) * size.models
  }

  if(k < 3) sh_arr[,2,] <- sh_arr[,2,] * asp.models ######

  if(rot.models!=0) for(i in 1:dim(sh_arr)[3]) {
    sh_arr[,,i] <- spdep::Rotation(sh_arr[,,i], rot.models*0.0174532925199)
  }


  if(!is.null(template)) {
    centroid <- matrix(rev_eigen(0, ordination$rotation[,1], ordination$center), ncol = 2, byrow = TRUE)
    temp_cent<-rbind(centroid,
                     Momocs::tps2d(template[-(1:p),],
                                   template[1:p,],
                                   centroid))

    temp_warpd_list <- lapply(1:dim(sh_arr)[3],
                              function(i) {Momocs::tps2d(temp_cent[-(1:p),],
                                                         temp_cent[1:p,],
                                                         sh_arr[,,i])})
    temp_warpd_arr <- abind::abind(temp_warpd_list, along = 3)
    sh_arr <- abind::abind(sh_arr, temp_warpd_arr, along = 1)
    p <- nrow(sh_arr)
  }


  models_mat<-c()
  models_arr<-sh_arr * 0
  for(i in 1:nrow(gridcoords)) {
    descentmat <- matrix(rep(gridcoords[i,], p), p, 2, byrow = TRUE)
    if(k > 2) descentmat <- cbind(descentmat, 0)
    models_mat <- rbind(models_mat, (sh_arr[,,i] * 0.07) + descentmat)
    models_arr[,,i] <- (sh_arr[,,i] * 0.07) + descentmat
  }

  return(list(models_mat = models_mat, models_arr = models_arr))

}

################################################################################################

#' Plot background 2D shape models
#'
#' @description Plot the output from [morphogrid()], for 2-dimensional morphometric
#'   data. Used internally.
#'
#' @param morphogrid An object containing the output of \code{morphogrid}.
#' @param x Optional vector with a non-morphometric variable to be plotted in
#'   the x axis.
#' @param y Optional vector with a non-morphometric variable to be plotted in
#'   the y axis.
#' @param links A list with the indices of the coordinates defining the
#'   wireframe (following the format used in \code{Morpho}).
#' @param template A 2-column matrix containing 1) the actual landmarks or
#'   semilandmarks being analized, followed by 2) the (x,y) cartesian
#'   coordinates defining a curve or set of curves that will be warped using the
#'   deformation interpolated from the changes between landmarks/semilandmarks
#'   (the actual positions of which must be marked with a row of NA, see the
#'   tails dataset).
#' @param datype Character; type of shape data used (\code{"landm"} or
#'   \code{"fcoef"}).
#' @param ordtype Character; method used for multivariate ordination
#'   (\code{"prcomp"}, \code{"bg_prcomp"}, \code{"phy_prcomp"}, \code{"pls_shapes"}
#'   or \code{"phy_pls_shapes"}).
#' @param axes Numeric of length 2, indicating the axes to be plotted.
#' @param p Numeric, indicating the number of landmarks/semilandmarks used (for
#'   landmark data only).
#' @param cex.ldm Numeric; size of landmarks/semilandmarks in the background
#'   models.
#' @param col.ldm The color of landmarks/semilandmarks in the background models.
#' @param col.models The color for wireframes/outlines.
#' @param bg.models Background color for outlines.
#' @param lwd.models Numeric; the width of the lines in wireframes/outlines.
#' @param plot Logical; whether to plot morphospace.
#' @param xlab,ylab Standard arguments passed to the generic plot function.
#'
#' @export
#'
#' @examples
#' #load data and packages
#' library(geomorph)
#' data("tails")
#'
#' #perform pca on tails shapes
#' pca <- prcomp(two.d.array(tails$shapes))
#'
#' #generate grid of shape models sampling the range of variation
#' #at 4 locations (the 4 corners of the scatterplot)
#' shapes_grid <- morphogrid(ordination = pca, axes = c(1,2), datype = "landm",
#'                           k = ncol(tails$shapes), p = nrow(tails$shapes),
#'                           nh = 2, nv = 2)
#'
#' #plot grid
#' plot_morphogrid2d(morphogrid = shapes_grid, links = tails$links,
#'                   datype = "landm", ordtype = "prcomp",
#'                   axes = c(1,2), p = 9, col.ldm = 1, cex.ldm = 1,
#'                   col.models = 1, lwd.models = 1)
#'
#' #project each shape in models_arr
#' points(shapes_grid$models_arr[,,1], pch=16, col = 1)
#' points(shapes_grid$models_arr[,,2], pch=16, col = 2)
#' points(shapes_grid$models_arr[,,3], pch=16, col = 3)
#' points(shapes_grid$models_arr[,,4], pch=16, col = 4)
plot_morphogrid2d <- function(x = NULL,
                              y = NULL,
                              morphogrid,
                              template = NULL,
                              links = NULL,
                              datype,
                              ordtype,
                              axes,
                              p,
                              xlab = NULL,
                              ylab = NULL,
                              cex.ldm,
                              col.ldm,
                              col.models,
                              lwd.models,
                              bg.models,
                              plot = TRUE) {


  xlim <- range(c(na.omit(morphogrid$models_mat[,1])))
  ylim <- range(c(na.omit(morphogrid$models_mat[,2])))

  if(length(axes) == 1) axes <- rep(axes, 2)

  if(is.null(xlab)) {
    if(!is.null(x)) {
      xlab <- "x"
    } else {
      if(ordtype == "prcomp") {
        xlab <- paste0("PC", axes[1])
      }
      if(ordtype == "bg_prcomp") {
        xlab <- paste0("bgPC", axes[1])
      }
      if(ordtype == "phy_prcomp") {
        xlab <- paste0("phyPC", axes[1])
      }
      if(ordtype == "pls_shapes") {
        xlab <- paste0("PLS-", axes[1])
      }
      if(ordtype == "phy_pls_shapes") {
        xlab <- paste0("phyPLS-", axes[1])
      }
    }
  }
  if(is.null(ylab)) {
    if(!is.null(y)) {
      ylab <- "y"
    } else {
      if(ordtype == "prcomp") {
        ylab <- paste0("PC", axes[2])
      }
      if(ordtype == "bg_prcomp") {
        ylab <- paste0("bgPC", axes[2])
      }
      if(ordtype == "phy_prcomp") {
        ylab <- paste0("phyPC", axes[2])
      }
      if(ordtype == "pls_shapes") {
        ylab <- paste0("PLS-", axes[2])
      }
      if(ordtype == "phy_pls_shapes") {
        ylab <- paste0("phyPLS-", axes[2])
      }
    }
  }

  if(plot == TRUE) {
    plot(morphogrid$models_mat, type = "n", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab)

    for(i in 1:dim(morphogrid$models_arr)[3]) {
      if(datype == "landm") {
        graphics::points(morphogrid$models_arr[1:p,,i],
                         pch = 16, cex = cex.ldm * 0.1, col = col.ldm)

        if(!is.null(template)) {
          graphics::lines(morphogrid$models_arr[-c(1:p),,i],
                          col = col.models, lwd = lwd.models)
        } else {
          for(l in 1:length(links)) graphics::lines(morphogrid$models_arr[,,i][links[[l]],],
                                                    col = col.models, lwd = lwd.models)
        }

      } else {
        graphics::polygon(morphogrid$models_arr[,,i], pch = 16,
                          col = bg.models, border = col.models, lwd = lwd.models)
      }
    }
  }
}


##################################################################################

#' Plot background 3D shape models
#'
#' @description Plot the output from [morphogrid()], for 3-dimensional morphometric
#'   (landmark) data. Used internally.
#'
#' @param morphogrid An object containing the output of \code{morphogrid}.
#' @param x Optional vector with a non-morphometric variable to be plotted in
#'   the x axis.
#' @param y Optional vector with a non-morphometric variable to be plotted in
#'   the y axis.
#' @param links A list with the indices of the coordinates defining the
#'   wireframe (following the format used in \code{Morpho}).
#' @param template An optional \code{"mesh3d"} object, which will be warped using
#'   TPS interpolation to produce the set of background shell models.
#' @param refshape reference shape (i.e., the mean landmark configuration)
#'   corresponding to the mesh provided in \code{template}.
#' @param ordtype Character; method used for multivariate ordination
#'   (\code{"prcomp"}, \code{"bg_prcomp"}, \code{"phy_prcomp"}, \code{"pls_shapes"}
#'   or \code{"phy_pls_shapes"}).
#' @param axes Numeric of length 2, indicating the axes to be plotted.
#' @param p Numeric, indicating the number of landmarks/semilandmarks used (for
#'   landmark data only).
#' @param cex.ldm Numeric; size of landmarks/semilandmarks in the background
#'   models.
#' @param col.ldm The color of landmarks/semilandmarks in the background models.
#' @param size.models Numeric; size factor for shape models.
#' @param col.models The color for wireframes/outlines.
#' @param bg.models Background color for meshes.
#' @param lwd.models Numeric; the width of the lines in wireframes/outlines.
#' @param alpha.models Numeric; transparency factor for background models.
#' @param asp.models Numeric; the y/x aspect ratio of shape models.
#' @param plot Logical; whether to plot morphospace.
#' @param xlim,ylim,xlab,ylab Standard arguments passed to the generic plot
#'   function.
#'
#' @details This function allows the user to choose the orientation of the 3D
#'   models by interactively rotating a shape model. Do not close the \code{rgl}
#'   window, or minimize it actively (just bring back Rstudio to the front and
#'   let the device get minimized pasively). The process of morphospace generation
#'    is rather slow, specially if a mesh is provided for \code{template}, a large
#'    number of shape models is asked, and/or \code{alpha.models} value is lower
#'    than \code{1}.
#'
#' @export
#'
#' @examples
#' #load data and packages
#' library(geomorph)
#' data("shells3D")
#' shapes <- shells3D$shapes
#'
#' #perform pca on tails shapes
#' pca <- prcomp(two.d.array(shapes))
#'
#' #generate grid of shape models sampling the range of variation
#' #at 4 locations (the 4 corners of the scatterplot)
#' shapes_grid <- morphogrid(ordination = pca, axes = c(1,2), datype = "landm",
#'                           k = ncol(shapes), p = nrow(shapes),
#'                           nh = 2, nv = 2)
#'
#' #get meanshape
#' meanshape <- consensus(shapes)
#'
#' \dontrun{
#' #plot grid (shape coordinates only)
#' plot_morphogrid3d(morphogrid = shapes_grid, refshape = meanshape,
#'                   ordtype = "prcomp", axes = c(1,2), p = 9, col.ldm = 1, cex.ldm = 1,
#'                   col.models = 1, lwd.models = 1, bg.models = "gray", size.models = 2,
#'                   asp.models = 3)
#'
#' #get shape corresponding to shells3D$mesh_meanspec using geomorph::findMeanSpec,
#' #then get mesh corresponding to mean shape using Morpho::tps3d
#' meanspec_id<- findMeanSpec(shapes)
#' meanspec_shape <- shapes[,,meanspec_id]
#' meanmesh <- tps3d(x = shells3D$mesh_meanspec , refmat = meanspec_shape, tarmat = meanshape)
#'
#' #plot grid (includinh mesh template)
#' plot_morphogrid3d(morphogrid = shapes_grid, template = meanmesh, refshape = meanshape,
#'                   ordtype = "prcomp", axes = c(1,2), p = 9, col.ldm = 1, cex.ldm = 1,
#'                   col.models = 1, lwd.models = 1, bg.models = "gray", size.models = 2,
#'                   asp.models = 3)
#' }
plot_morphogrid3d <- function(x = NULL,
                              y = NULL,
                              morphogrid,
                              refshape,
                              template = NULL,
                              links = NULL,
                              ordtype,
                              axes,
                              p,
                              xlim = NULL,
                              ylim = NULL,
                              xlab = NULL,
                              ylab = NULL,
                              cex.ldm,
                              col.ldm,
                              col.models,
                              lwd.models,
                              size.models,
                              alpha.models = 1,
                              bg.models,
                              asp.models,
                              plot = TRUE) {


  if(is.null(xlim)) xlim <- range(c(morphogrid$models_mat[,1]))
  if(is.null(ylim)) ylim <- range(c(morphogrid$models_mat[,2]))

  if(length(axes) == 1) axes <- rep(axes, 2)

  if(is.null(xlab)) {
    if(!is.null(x)) {
      xlab <- "x"
    } else {
      if(ordtype == "prcomp") {
        xlab <- paste0("PC", axes[1])
      }
      if(ordtype == "bg_prcomp") {
        xlab <- paste0("bgPC", axes[1])
      }
      if(ordtype == "phy_prcomp") {
        xlab <- paste0("phyPC", axes[1])
      }
      if(ordtype == "pls_shapes") {
        xlab <- paste0("PLS-", axes[1])
      }
      if(ordtype == "phy_pls_shapes") {
        xlab <- paste0("phyPLS-", axes[1])
      }
    }
  }
  if(is.null(ylab)) {
    if(!is.null(y)) {
      ylab <- "y"
    } else {
      if(ordtype == "prcomp") {
        ylab <- paste0("PC", axes[2])
      }
      if(ordtype == "bg_prcomp") {
        ylab <- paste0("bgPC", axes[2])
      }
      if(ordtype == "phy_prcomp") {
        ylab <- paste0("phyPC", axes[2])
      }
      if(ordtype == "pls_shapes") {
        ylab <- paste0("PLS-", axes[2])
      }
      if(ordtype == "phy_pls_shapes") {
        ylab <- paste0("phyPLS-", axes[2])
      }
    }
  }


  wd <- tempdir()
  enter <- NULL
  if(!is.null(template)) {

    while(is.null(enter)) {

      refmesh <- template

      rgl::plot3d(refmesh, col = bg.models, specular = "black", axes = FALSE, aspect = FALSE,
                  xlab = "", ylab = "", zlab = "", alpha = alpha.models)

      rgl::plot3d(refshape, col = col.ldm, specular = "black", axes = FALSE, aspect = FALSE,
                  xlab = "", ylab = "", zlab = "", type = "s", size = cex.ldm, add = TRUE)

      for(l in 1:length(links)) rgl::lines3d(refshape[links[[l]],],
                                             col = col.models, lwd = lwd.models)

      cat("Preparing for snapshot: rotate mean shape to the desired orientation\n (don't close or minimize the rgl device).")

      enter <- readline("Press <Enter> in the console to continue:")

      cat("This will take a minute")

    }

    for(i in 1:dim(morphogrid$models_arr)[3]) {

      modelmesh <- Morpho::tps3d(x = refmesh , refmat = refshape, tarmat = morphogrid$models_arr[,,i])

      rgl::plot3d(modelmesh, col = bg.models, specular = "black", axes = FALSE, aspect = FALSE,
                  xlab = "", ylab = "", zlab = "", alpha = alpha.models)

      rgl::plot3d(morphogrid$models_arr[,,i], col = col.ldm, specular = "black", axes = FALSE, aspect = FALSE,
                  xlab = "", ylab = "", zlab = "", type = "s", size = cex.ldm, add = TRUE)

      for(l in 1:length(links)) rgl::lines3d(morphogrid$models_arr[,,i][links[[l]],],
                                             col = col.models, lwd = lwd.models)

      rgl::rgl.snapshot(paste0(wd, "model", i, ".png"))
    }

  } else {
    while(is.null(enter)) {

      rgl::plot3d(refshape, col = col.ldm, specular = "black", axes = FALSE, aspect = FALSE,
                  xlab = "", ylab = "", zlab = "", type = "s", size = cex.ldm)

      for(l in 1:length(links)) rgl::lines3d(refshape[links[[l]],],
                                             col = col.models, lwd = lwd.models)

      cat("Preparing for snapshot: rotate mean shape to the desired orientation\n (don't close or minimize the rgl device).")

      enter <- readline("Press <Enter> in the console to continue:")

    }

    for(i in 1:dim(morphogrid$models_arr)[3]) {

      rgl::plot3d(morphogrid$models_arr[,,i], col = col.ldm, specular = "black",
                  axes = FALSE, aspect = FALSE, xlab = "", ylab = "", zlab = "",
                  type = "s", size = cex.ldm)

      for(l in 1:length(links)) rgl::lines3d(morphogrid$models_arr[,,i][links[[l]],],
                                             col = col.models, lwd = lwd.models)

      rgl::rgl.snapshot(paste0(wd, "model", i, ".png"))
    }
  }



  for(i in 1:dim(morphogrid$models_arr)[3]) {
    model_i <- magick::image_read(paste0(wd, "model", i, ".png"))
    model_i_clean <- magick::image_fill(model_i, color = "transparent", refcolor = "white", fuzz = 4, point = "+1+1")
    model_i_clean_trimmed <- magick::image_trim(model_i_clean)
    magick::image_write(model_i_clean_trimmed, path = paste0(wd,"model",i,".png"), format = "png")
  }


  models <- lapply(1:dim(morphogrid$models_arr)[3], function (i) {
    model_i <- png::readPNG(paste0(wd, "model", i, ".png"), native = TRUE)
  })


  model_centers <- t(apply(morphogrid$models_arr, 3, colMeans))
  model_ranges <- lapply(1:length(models), function(i) {

    halfsize <- size.models * 0.015

    matrix(c(c((model_centers[i,1] - halfsize), (model_centers[i,1] + halfsize)),
             c((model_centers[i,2] - halfsize * (asp.models/2)),
               (model_centers[i,2] + halfsize * (asp.models/2)))),
           nrow = 2, byrow = FALSE)
  })

  xlim <- range(c(lapply(model_ranges, function(x) {x[,1]}), xlim))
  ylim <- range(c(lapply(model_ranges, function(x) {x[,2]}), ylim))

  if(plot == TRUE) {
    plot(0, type = "n", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab)

    for(i in 1:length(models)) {
      rasterImage(models[[i]], model_ranges[[i]][1,1], model_ranges[[i]][1,2],
                  model_ranges[[i]][2,1], model_ranges[[i]][2,2])
    }
  }
}



################################################################################################


#' Rotate Fourier shape 180 degrees
#'
#' @description Correct 180-degrees spurious rotation in closed outline shapes.
#'   Used internally.
#'
#' @param fcoef The set of Fourier coefficients measuring the shape(s) to be
#'   rotated.
#'
#' @return A set of Fourier coefficients describing the rotated outline(s).
#'
#' @export
#'
#' @references Iwata, H., & Ukai, Y. (2002). \emph{SHAPE: a computer program
#'   package for quantitative evaluation of biological shapes based on elliptic
#'   Fourier descriptors}. Journal of Heredity, 93(5), 384-385.
#'
#' @examples
#' #load shells data, plot the first shape
#' data("shells")
#' shape1 <- shells$shapes$coe[1,]
#' plot(inv_efourier(shape1, 300))
#'
#' #rotate and plot again
#' shape1_rot <- rotate_fcoef(shape1)
#' plot(inv_efourier(shape1_rot, 300))
rotate_fcoef <- function(fcoef) {

  fcoef <- rbind(fcoef)
  nb.h <- ncol(fcoef) / 4

  a <- rbind(fcoef[,1:nb.h])
  b <- rbind(fcoef[,(nb.h + 1):(nb.h * 2)])
  c <- rbind(fcoef[,((nb.h * 2) + 1):(nb.h * 3)])
  d <- rbind(fcoef[,((nb.h * 3) + 1):(nb.h * 4)])

  a[,1:nb.h %% 2 == 0] <- a[,1:nb.h %% 2 == 0] * -1
  b[,1:nb.h %% 2 == 0] <- b[,1:nb.h %% 2 == 0] * -1
  c[,1:nb.h %% 2 == 0] <- c[,1:nb.h %% 2 == 0] * -1
  d[,1:nb.h %% 2 == 0] <- d[,1:nb.h %% 2 == 0] * -1

  rot_fcoef <- rbind(cbind(a, b, c, d))
  colnames(rot_fcoef) <- paste0(rep(c("A", "B", "C", "D"), each = nb.h), 1:nb.h)
  return(rot_fcoef)
}


##########################################################################################


#' Plot phenogram
#'
#' @description Plot phenogram using phylogenetic and morphometric information. Used internally.
#'
#' @param x Optional vector with a non-morphometric variable to be plotted in
#'   the x axis.
#' @param y Optional vector with a non-morphometric variable to be plotted in
#'   the y axis.
#' @param tree A \code{"phy"} object containing a phylogenetic tree. Tip labels
#'   should match the row names from \code{x} or \code{y}.
#' @param phylo_scores A matrix containing the scores from tips and nodes of
#'   the phylogeny provided in \code{tree}.
#' @param axes Numeric of length 2, indicating the axes to be plotted.
#' @param col.groups The color of the scatter points corresponding to groups
#'   mean shapes.
#' @param cex.groups Numeric; the size of the scatter points corresponding to
#'   groups mean shapes.
#' @param lwd.branches Numeric; the width of the lines depicting phylogenetic
#'   branches.
#' @param points Logical; whether to plot the scatter points.
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
#' msp <- mspace(shapes, mag = 0.7, axes = c(1,2), cex.ldm = 0, plot = FALSE) %>%
#'   proj_shapes(shapes = shapes, col = c(1:13)[species], pch = 1, cex = 0.7) %>%
#'   proj_consensus(shapes = sp_shapes, pch = 21, bg = 1:13, cex = 2) %>%
#'   proj_phylogeny(tree = tree)
#'
#' #get node heights
#' heights <- phytools::nodeHeights(tree)
#' node_heights <- c(rep(max(heights), length(tree$tip.label)), unique(heights[,1]))
#'
#' #plot simple phengram
#' plot(node_heights, msp$phylo_scores[,1])
#' plot_phenogram(tree = tree, x = node_heights, phylo_scores = msp$phylo_scores,
#'                axes = 1, lwd.branches = 1, cex.groups = 1, col.groups = 1, points = TRUE)
plot_phenogram <- function(x = NULL,
                           y = NULL,
                           tree,
                           phylo_scores,
                           axes,
                           lwd.branches,
                           cex.groups,
                           col.groups,
                           points) {

  for(i in 1:nrow(tree$edge)) {
    phyloxy <- cbind(x, phylo_scores[,axes[1]], y)
    graphics::lines(rbind(phyloxy[tree$edge[i, 1],],
                          phyloxy[tree$edge[i, 2],]), lwd = lwd.branches)
  }
  if(points == TRUE) {
    ntips <- length(tree$tip.label)
    graphics::points(phyloxy[-c(1:ntips),], pch = 16)
    graphics::points(phyloxy[c(1:ntips),][rownames(phylo_scores)[1:ntips],],
                     bg = col.groups, pch = 21, cex = cex.groups)

  }
}
