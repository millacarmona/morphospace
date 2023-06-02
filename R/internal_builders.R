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
#' meanshapes_arr <- expected_shapes(tails$shapes, x = tails$data$species)
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
#'   when used on blocks of variables. Can deal with phylogenetic data too
#'   using \code{ape} and \code{phytools} functions. Used internally.
#'
#' @param x First block of variables
#' @param y Second block of variables
#' @param tree An optional \code{"phylo"} object containing a phylogenetic
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
#' sp_shapes <- expected_shapes(shapes, species)[,,tree$tip.label]
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
    part_R <- R[seq_len(ncol(x)), (ncol(x) + 1):(ncol(x) + ncol(y))]

    svd <- svd(part_R)

  } else {
    vcv <- stats::cov(cbind(x,y))
    part_vcv <- vcv[seq_len(ncol(x)), (ncol(x) + 1):(ncol(x) + ncol(y))]

    svd <- svd(part_vcv)


  }

  ndims <- min(nrow(x), ncol(x), ncol(y))
  sdev <- (svd$d)[seq_len(ndims)]
  if(ncol(x) != ncol(y)) {
    rotations <- list(svd$v, svd$u)
    whichy <- which(unlist(lapply(rotations, nrow)) == ncol(y))
    whichx <- which(unlist(lapply(rotations, nrow)) == ncol(x))
    y_rotation <- cbind(rotations[[whichy]][,seq_len(ndims)])
    x_rotation <- cbind(rotations[[whichx]][,seq_len(ndims)])
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
#' @description A wrapper for [Momocs::efourier_i()] to transform a
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

  a <- c(coe[seq_len(nb.h)])
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

#' Adjust aspect and scale of background shape models for 2D data
#'
#' @description Avoid background shape models distortion caused by differences in ranges
#'   of x and y axes. Used internally.
#'
#' @param models An array containing the background shape models.
#' @param frame The frame in which shape models are to be plotted.
#' @param model_width Numeric; the width of a reference shape model (usually the consensus).
#' @param model_height Numeric; the height of a reference shape model (usually the consensus).
#'
#' @return An array containing the adjusted shape models.
#'
#' @export
#'
#' @examples
#' #load package and data
#' library(geomorph)
#' library(Morpho)
#' data("wings")
#' shapes <- wings$shapes
#'
#' #perform pca, extract ranges for PC1 and 2 and plot
#' pca <- prcomp(two.d.array(shapes))
#' xlim <- range(pca$x[,1])
#' ylim <- range(pca$x[,2])
#' plot(pca$x, col = "gray", xlim = xlim, ylim = ylim)
#'
#' #calculate and plot shape grid
#' shapes_grid0 <- morphogrid(ordination = pca, axes = c(1,2), datype = "landm",
#'                           k = ncol(shapes), p = nrow(shapes),
#'                           nh = 4, nv = 4, xlim = xlim, ylim = ylim)
#' for(i in 1:dim(shapes_grid0$models_arr)[3]) {
#'   points(shapes_grid0$models_arr[,,i], type = "l", col = "red")
#' }
#'
#' #amplify range of x axis and replot grid
#' newframe <- cbind(pca$x[,1] * 20, pca$x[,2])
#' plot(newframe, col = "gray")
#' for(i in 1:dim(shapes_grid0$models_arr)[3]) {
#'   points(shapes_grid0$models_arr[,,i], type = "l", col = "red")
#' }
#'
#' #calculate width and height of consensus shape
#' wh <- abs(apply(apply(expected_shapes(shapes_grid0$models_arr), 2, range), 2, diff))
#'
#' #ajust grid to new frame and plot
#' adj_grid <- adjust_models2d(models = shapes_grid0$models_arr, frame = newframe,
#'                             model_width = wh[1], model_height = wh[2])
#' for(i in 1:dim(adj_grid)[3]) points(adj_grid[,,i], type = "l", col = "blue", lwd = 2)
adjust_models2d <- function(models, frame, model_width, model_height) {

  frame_xydiffrange <- abs(apply(apply(frame, 2, range), 2, diff))
  frame_width <- frame_xydiffrange[1]
  frame_height <- frame_xydiffrange[2]

  models_adj_l <- lapply(seq_len(dim(models)[3]), function(i) {

    model <- models[,,i]
    if(frame_width >= frame_height) {
      model_ratio <- model_width / model_height
      frame_ratio <- frame_width / frame_height

      model[,1] <- model[,1] * (frame_ratio / model_ratio)

    } else {
      model_ratio <- model_height / model_width
      frame_ratio <- frame_height / frame_width

      model[,2] <- model[,2] * (frame_ratio / model_ratio)
    }

    model
  })
  models_adj <- abind::abind(models_adj_l, along = 3)


  model_width2  <- abs(apply(apply(expected_shapes(models_adj), 2, range), 2, diff))[1]
  model_height2 <- abs(apply(apply(expected_shapes(models_adj), 2, range), 2, diff))[2]

  models_adj_l2 <- lapply(seq_len(dim(models)[3]), function(i) {

    model <- models_adj[,,i]

    model[,1] <- model[,1] * ((frame_width / dim(models)[3]) / model_width2)
    model[,2] <- model[,2] * ((frame_height / dim(models)[3]) / model_height2)

    model * 3
  })
  models_adj <- abind::abind(models_adj_l2, along = 3)

  return(models_adj)
}


################################################################################################

#' Adjust aspect and scale of background shape models for 3D data
#'
#' @description Avoid background shape models distortion caused by differences in ranges
#'   of x and y axes. Used internally.
#'
#' @param models A list containing PNG images depicting background shape models.
#' @param frame The frame in which shape models are to be plotted.
#' @param size.models Numeric; size factor for shape models.
#' @param asp.models Numeric; the y/x aspect ratio of shape models.
#'
#' @return A list containing the adjusted frame of each element in \code{models}, and
#'   new, adjusted limits for the x and y axes from \code{frame}.
#'
#' @export
#'
#' @seealso \code{\link{morphogrid}}, \code{\link{plot_morphogrid3d}}
#'
#' @examples
#' #load package and data
#' library(geomorph)
#' library(Morpho)
#' library(rgl)
#' library(png)
#' data("shells3D")
#' shapes <- shells3D$shapes
#'
#' \dontrun{
#' #perform pca, extract ranges for PC1 and 2
#' pca <- prcomp(two.d.array(shapes))
#' xlim <- range(pca$x[,1])
#' ylim <- range(pca$x[,2])
#'
#' #calculate shape grid and compute model centroid
#' shapes_grid0 <- morphogrid(ordination = pca, axes = c(1,2), datype = "landm",
#'                           k = ncol(shapes), p = nrow(shapes),
#'                           nh = 2, nv = 2, xlim = xlim, ylim = ylim)
#' model_centers <- t(apply(shapes_grid0$models_arr, 3, colMeans))[,1:2]
#'
#' #create a temporal directory, plot 3D models with rgl, and take a snapshot
#' wd <- tempdir()
#' for(i in 1:dim(shapes_grid0$models_arr)[3]) {
#'   plot3d(shapes_grid0$models_arr[,,i], aspect = FALSE, axes = FALSE,
#'               xlab = "", ylab = "", zlab = "")
#'   rgl.snapshot(paste0(wd,"model",i,".png"))
#' }
#'
#' #re-import snapshots
#' models <- lapply(1:dim(shapes_grid0$models_arr)[3], function (i) {
#'   model_i <- readPNG(paste0(wd, "model", i, ".png"), native = TRUE)
#' })
#'
#' #adjust snapshots
#' adj <- adjust_models3d(models = models, frame = model_centers, size.models = 1, asp.models = 1)
#'
#' #plot everything
#' plot(pca$x, xlim = adj$xlim, ylim = adj$ylim)
#' for(i in 1:length(models)) {
#'   rasterImage(models[[i]], adj$model_frames[[i]][1,1], adj$model_frames[[i]][1,2],
#'               adj$model_frames[[i]][2,1], adj$model_frames[[i]][2,2])
#' }
#'
#' #kind of. Anyway, use plot_morphogrid3d that do the full process.
#' }
adjust_models3d <- function(models, frame, size.models, asp.models) {

  mw <- sapply(models, function(x) {dim(x)})[1,]
  mh <- sapply(models, function(x) {dim(x)})[2,]

  fw <- abs(diff(range(frame[,1])))
  fh <- abs(diff(range(frame[,2])))
  N <- nrow(frame)

  mw_std <- (mw/max(mw)) * (fw / N) * (size.models * 7)
  mh_std <- (mh/max(mh)) * (fh / N) * (size.models * 7) * asp.models

  xlim_min <- xlim_max <- ylim_min <- ylim_max <- NULL
  model_frames <- list()
  for(i in seq_len(nrow(frame))) {

    pw_min <- frame[i,1] - (mw_std[i] / 2)
    pw_max <- frame[i,1] + (mw_std[i] / 2)
    xlim_min <- min(xlim_min, pw_min)
    xlim_max <- max(xlim_max, pw_max)


    ph_min <- frame[i,2] - (mh_std[i] / 2)
    ph_max <- frame[i,2] + (mh_std[i] / 2)
    ylim_min <- min(ylim_min, ph_min)
    ylim_max <- max(ylim_max, ph_max)


    w <- rbind(pw_min, pw_max)
    h <- rbind(ph_min, ph_max)
    model_frames[[i]] <- cbind(w, h)
  }

  new_xlim <- c(xlim_min, xlim_max)
  new_ylim <- c(ylim_min, ylim_max)

  results <- list(model_frames = model_frames,
                  xlim = new_xlim, ylim = new_ylim)
  return(results)

}

####################################################################################

#' Generate background shape models
#'
#' @description Calculate and arrange background shape models for morphospaces.
#'   Used internally.
#'
#' @param ordination An ordination (i.e. a \code{"prcomp"}, \code{"bg_prcomp"},
#'   \code{"phy_prcomp"} or \code{"pls_shape"} object).
#' @param axes Numeric of length 1 or 2, indicating the morphometric axes to be
#'   plotted. If values for either \code{x} or \code{y} are provided, only the
#'   first value of this argument is considered.
#' @param datype Character; the type of shape data (either \code{"fcoef"} or
#'   \code{"landm"}).
#' @param template A 2-column matrix containing landmarks/semilandmarks followed
#'   by coordinates defining a curve or set of curves describing additional aspects
#'   of morphology, which will be warped using TPS interpolation to produce the set
#'   of background shell models(see \code{\link{build_template2d}}).
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
#' @return A list of length 3 containing:
#' \itemize{
#'   \item \code{$models_mat:} {a 2-column matrix with the (x,y) coordinates of
#'   all the shape models in the background.}
#'   \item \code{$models_arr:} {same as \code{models_mat} but in 3-margins array
#'   format.}
#'   \item \code{$grid:} {coordinates marking the centroid of each shape model.
#'   Intended for internal use.}
#' }
#'
#' @export
#'
#' @seealso \code{\link{plot_morphogrid2d}}, \code{\link{plot_morphogrid3d}}
#'
#' @references MacLeod, N. (2009). \emph{Form & shape models}. Palaeontological
#'   Association Newsletter, 72(620), 14-27.
#'
#' @examples
#'  #load data and packages
#'  library(geomorph)
#'  data("tails")
#'
#'  #perform pca on tails shapes
#'  pca <- prcomp(two.d.array(tails$shapes))
#'
#'  #generate grid of shape models sampling the range of variation
#'  #at 4 locations (the 4 corners of the scatterplot)
#'  shapes_grid <- morphogrid(ordination = pca, axes = c(1,2), datype = "landm",
#'                            k = ncol(tails$shapes), p = nrow(tails$shapes),
#'                            nh = 2, nv = 2)
#'
#'  #plot grid from $models_mat and project each shape in models_arr
#'  plot(shapes_grid$models_mat)
#'  points(shapes_grid$models_arr[,,1], pch=16, col = 1)
#'  points(shapes_grid$models_arr[,,2], pch=16, col = 2)
#'  points(shapes_grid$models_arr[,,3], pch=16, col = 3)
#'  points(shapes_grid$models_arr[,,4], pch=16, col = 4)
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
      if(is.null(xlim)) {
        plotframe_x <- range(x)
      } else {
        plotframe_x <- sort(xlim)
      }
    } else {
      if(is.null(xlim)) {
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
      if(is.null(ylim)) {
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
    coords_l <- lapply(seq_len(nrow(sh_mat)), function(i) {
      sh <- inv_efourier(coe = sh_mat[i,], nb.pts = p)
      sh / Morpho::cSize(sh)
    })
    sh_arr <- abind::abind(coords_l, along = 3)
  } else {
    sh_arr <- geomorph::arrayspecs(sh_mat, p = p, k = k)
    for(i in seq_len(dim(sh_arr)[3])) sh_arr[,,i] <- sh_arr[,,i] / Morpho::cSize(sh_arr[,,i])
  }


  if(rot.models != 0) for(i in seq_len(dim(sh_arr)[3])) {
    sh_arr[,,i] <- spdep::Rotation(sh_arr[,,i], rot.models * 0.0174532925199)
  }


  if(k < 3) {
    wh <- abs(apply(apply(expected_shapes(sh_arr), 2, range), 2, diff))
    sh_arr <- adjust_models2d(sh_arr, gridcoords, wh[1], wh[2])

    sh_arr[,2,] <- sh_arr[,2,] * asp.models
  }
  sh_arr <- sh_arr * size.models


  if(!is.null(template)) {
    centroid <- matrix(rev_eigen(0, ordination$rotation[,1], ordination$center), ncol = 2, byrow = TRUE)
    temp_cent<-rbind(centroid,
                     Momocs::tps2d(template[-(seq_len(p)),],
                                   template[(seq_len(p)),],
                                   centroid))

    temp_warpd_list <- lapply(1:dim(sh_arr)[3],
                              function(i) {Momocs::tps2d(temp_cent[-(seq_len(p)),],
                                                         temp_cent[(seq_len(p)),],
                                                         sh_arr[,,i])})
    temp_warpd_arr <- abind::abind(temp_warpd_list, along = 3)
    sh_arr <- abind::abind(sh_arr, temp_warpd_arr, along = 1)
    p <- nrow(sh_arr)
  }


  models_mat <- NULL
  models_arr <- sh_arr * 0
  for(i in seq_len(nrow(gridcoords))) {
    descentmat <- matrix(rep(gridcoords[i,], p), p, 2, byrow = TRUE)
    if(k > 2) descentmat <- cbind(descentmat, 0)
    models_arr[,,i] <- (sh_arr[,,i]) + descentmat
    models_mat <- rbind(models_mat, (sh_arr[,,i]) + descentmat)
  }


  return(list(models_mat = models_mat, models_arr = models_arr, grid = gridcoords))

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
#' @param template A 2-column matrix containing landmarks/semilandmarks followed
#'   by coordinates defining a curve or set of curves describing additional aspects
#'   of morphology, which will be warped using TPS interpolation to produce the set
#'   of background shell models(see \code{\link{build_template2d}}).
#' @param datype Character; type of shape data used (\code{"landm"} or
#'   \code{"fcoef"}).
#' @param ordtype Character; method used for multivariate ordination
#'   (\code{"prcomp"}, \code{"bg_prcomp"}, \code{"phy_prcomp"}, \code{"pls_shapes"}
#'   or \code{"phy_pls_shapes"}).
#' @param axes Numeric of length 2, indicating the axes to be plotted.
#' @param adj_frame Numeric of length 2, providing \emph{a posteriori} scaling
#'   factors for the width and height of the frame, respectively.
#' @param p Numeric, indicating the number of landmarks/semilandmarks used (for
#'   landmark data only).
#' @param cex.ldm Numeric; size of landmarks/semilandmarks in the background
#'   models.
#' @param col.ldm The color of landmarks/semilandmarks in the background models.
#' @param col.models The color for wireframes/outlines.
#' @param bg.models Background color for outlines.
#' @param lwd.models Numeric; the width of the lines in wireframes/outlines.
#' @param plot Logical; whether to plot morphospace.
#' @param models Logical; whether to plot background shape models.
#' @param xlab,ylab Standard arguments passed to the generic plot function.
#'
#' @export
#'
#' @seealso \code{\link{morphogrid}}, \code{\link{adjust_models2d}}
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
                              adj_frame = c(1,1),
                              p,
                              xlab = NULL,
                              ylab = NULL,
                              cex.ldm,
                              col.ldm,
                              col.models,
                              lwd.models,
                              bg.models,
                              plot = TRUE,
                              models = TRUE) {


  xlim <- range(c(stats::na.omit(morphogrid$models_mat[,1])))
  ylim <- range(c(stats::na.omit(morphogrid$models_mat[,2])))

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

    xlim <- c(xlim[1] + (diff(xlim) * (1 - adj_frame[1]) / 2),
              xlim[2] - (diff(xlim) * (1 - adj_frame[1]) / 2))
    ylim <- c(ylim[1] + (diff(ylim) * (1 - adj_frame[2]) / 2),
              ylim[2] - (diff(ylim) * (1 - adj_frame[2]) / 2))

    plot(morphogrid$models_mat, type = "n", xlim = xlim, ylim = ylim , xlab = xlab, ylab = ylab)

    if(models == TRUE) {

      for(i in seq_len(dim(morphogrid$models_arr)[3])) {
        if(datype == "landm") {
          graphics::points(morphogrid$models_arr[seq_len(p),,i],
                           pch = 16, cex = cex.ldm * 0.1, col = col.ldm)

          if(!is.null(template)) {
            graphics::lines(morphogrid$models_arr[-c(seq_len(p)),,i],
                            col = col.models, lwd = lwd.models)
          } else {
            for(l in seq_len(length(links))) graphics::lines(morphogrid$models_arr[,,i][links[[l]],],
                                                      col = col.models, lwd = lwd.models)
          }

        } else {
          graphics::polygon(morphogrid$models_arr[,,i], pch = 16,
                            col = bg.models, border = col.models, lwd = lwd.models)
        }
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
#' @param template An optional \code{"mesh3d"} object containing
#'   geometry of the structure the landmarks were placed on (for 3D shape data),
#'   corresponding to the mean shape of the sample, which will be warped using
#'   TPS interpolation to produce the set of background shell models.
#' @param refshape reference shape (i.e., the mean landmark configuration)
#'   corresponding to the mesh provided in \code{template}.
#' @param ordtype Character; method used for multivariate ordination
#'   (\code{"prcomp"}, \code{"bg_prcomp"}, \code{"phy_prcomp"}, \code{"pls_shapes"}
#'   or \code{"phy_pls_shapes"}).
#' @param axes Numeric of length 2, indicating the axes to be plotted.
#' @param adj_frame Numeric of length 2, providing \emph{a posteriori} scaling
#'   factors for the width and height of the frame, respectively.
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
#' @param models Logical; whether to plot background shape models.
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
#' @seealso \code{\link{morphogrid}}, \code{\link{adjust_models3d}}
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
#' meanshape <- expected_shapes(shapes)
#'
#' \dontrun{
#' #plot grid (shape coordinates only)
#' plot_morphogrid3d(morphogrid = shapes_grid, refshape = meanshape,
#'                   ordtype = "prcomp", axes = c(1,2), col.ldm = 1, cex.ldm = 1,
#'                   col.models = 1, lwd.models = 1, bg.models = "gray", size.models = 2,
#'                   asp.models = 1)
#'
#' #get shape corresponding to shells3D$mesh_meanspec using geomorph::findMeanSpec,
#' #then get mesh corresponding to mean shape using Morpho::tps3d
#' meanspec_id<- findMeanSpec(shapes)
#' meanspec_shape <- shapes[,,meanspec_id]
#' meanmesh <- tps3d(x = shells3D$mesh_meanspec , refmat = meanspec_shape, tarmat = meanshape)
#'
#' #plot grid (includinh mesh template)
#' plot_morphogrid3d(morphogrid = shapes_grid, template = meanmesh, refshape = meanshape,
#'                   ordtype = "prcomp", axes = c(1,2), col.ldm = 1, cex.ldm = 1,
#'                   col.models = 1, lwd.models = 1, bg.models = "gray", size.models = 2,
#'                   asp.models = 1)
#' }
plot_morphogrid3d <- function(x = NULL,
                              y = NULL,
                              morphogrid,
                              refshape,
                              template = NULL,
                              links = NULL,
                              ordtype,
                              axes,
                              xlim = NULL,
                              ylim = NULL,
                              xlab = NULL,
                              ylab = NULL,
                              adj_frame = 1,
                              cex.ldm,
                              col.ldm,
                              col.models,
                              lwd.models,
                              size.models,
                              alpha.models = 1,
                              bg.models,
                              asp.models,
                              plot = TRUE,
                              models = TRUE) {

  plot.models <- models

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

      for(l in seq_len(length(links))) rgl::lines3d(refshape[links[[l]],],
                                             col = col.models, lwd = lwd.models)

      cat("Preparing for snapshot: rotate mean shape to the desired orientation\n (don't close or minimize the rgl device).")

      enter <- readline("Press <Enter> in the console to continue:")

      rgl::rgl.snapshot(paste0(wd, "model", dim(morphogrid$models_arr)[3] + 1, ".png"))

      cat("This will take a minute...")
    }

    for(i in seq_len(dim(morphogrid$models_arr)[3])) {
      modelmesh <- Morpho::tps3d(x = refmesh , refmat = refshape, tarmat = morphogrid$models_arr[,,i])

      rgl::plot3d(modelmesh, col = bg.models, specular = "black", axes = FALSE, aspect = FALSE,
                  xlab = "", ylab = "", zlab = "", alpha = alpha.models)

      rgl::plot3d(morphogrid$models_arr[,,i], col = col.ldm, specular = "black", axes = FALSE, aspect = FALSE,
                  xlab = "", ylab = "", zlab = "", type = "s", size = cex.ldm, add = TRUE)

      for(l in seq_len(length(links))) rgl::lines3d(morphogrid$models_arr[,,i][links[[l]],],
                                             col = col.models, lwd = lwd.models)

      rgl::rgl.snapshot(paste0(wd, "model", i, ".png"))
    }
    cat("\nDONE.")

  } else {
    while(is.null(enter)) {
      rgl::plot3d(refshape, col = col.ldm, specular = "black", axes = FALSE, aspect = FALSE,
                  xlab = "", ylab = "", zlab = "", type = "s", size = cex.ldm)

      for(l in seq_len(length(links))) rgl::lines3d(refshape[links[[l]],],
                                             col = col.models, lwd = lwd.models)

      cat("Preparing for snapshot: rotate mean shape to the desired orientation\n (don't close or minimize the rgl device).")

      enter <- readline("Press <Enter> in the console to continue:")

      rgl::rgl.snapshot(paste0(wd, "model", dim(morphogrid$models_arr)[3] + 1, ".png"))

    }

    for(i in seq_len(dim(morphogrid$models_arr)[3])) {
      rgl::plot3d(morphogrid$models_arr[,,i], col = col.ldm, specular = "black",
                  axes = FALSE, aspect = FALSE, xlab = "", ylab = "", zlab = "",
                  type = "s", size = cex.ldm)

      for(l in seq_len(length(links))) rgl::lines3d(morphogrid$models_arr[,,i][links[[l]],],
                                             col = col.models, lwd = lwd.models)

      rgl::rgl.snapshot(paste0(wd, "model", i, ".png"))
    }
  }


  for(i in seq_len((dim(morphogrid$models_arr)[3] + 1))) {
    model_i <- magick::image_read(paste0(wd, "model", i, ".png"))
    model_i_clean <- magick::image_fill(model_i, color = "transparent", refcolor = "white", fuzz = 4, point = "+1+1")
    magick::image_write(model_i_clean, path = paste0(wd,"model",i,".png"), format = "png")
  }


  models <- lapply(seq_len((dim(morphogrid$models_arr)[3] + 1)), function (i) {
    model_i <- png::readPNG(paste0(wd, "model", i, ".png"), native = TRUE)
  })

  model_centers <- t(apply(morphogrid$models_arr, 3, colMeans))[,1:2]

  adj <- adjust_models3d(models = models, frame = model_centers,
                         size.models = size.models, asp.models = asp.models)

  model_frames <- adj$model_frames
  new_xlim <- c(adj$xlim[1] + (diff(adj$xlim) * (1 - adj_frame[1]) / 2),
                adj$xlim[2] - (diff(adj$xlim) * (1 - adj_frame[1]) / 2))
  new_ylim <- c(adj$ylim[1] + (diff(adj$ylim) * (1 - adj_frame[2]) / 2),
                adj$ylim[2] - (diff(adj$ylim) * (1 - adj_frame[2]) / 2))


  if(plot == TRUE) {
    plot(0, type = "n", xlim = new_xlim, ylim = new_ylim, xlab = xlab, ylab = ylab)

    if(plot.models == TRUE) {

      for(i in seq_len((length(models) - 1))) {
      graphics::rasterImage(models[[i]], model_frames[[i]][1,1], model_frames[[i]][1,2],
                  model_frames[[i]][2,1], model_frames[[i]][2,2])
    }
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

  a <- rbind(fcoef[,seq_len(nb.h)])
  b <- rbind(fcoef[,(nb.h + 1):(nb.h * 2)])
  c <- rbind(fcoef[,((nb.h * 2) + 1):(nb.h * 3)])
  d <- rbind(fcoef[,((nb.h * 3) + 1):(nb.h * 4)])

  a[,seq_len(nb.h) %% 2 == 0] <- a[,seq_len(nb.h) %% 2 == 0] * -1
  b[,seq_len(nb.h) %% 2 == 0] <- b[,seq_len(nb.h) %% 2 == 0] * -1
  c[,seq_len(nb.h) %% 2 == 0] <- c[,seq_len(nb.h) %% 2 == 0] * -1
  d[,seq_len(nb.h) %% 2 == 0] <- d[,seq_len(nb.h) %% 2 == 0] * -1

  rot_fcoef <- rbind(cbind(a, b, c, d))
  colnames(rot_fcoef) <- paste0(rep(c("A", "B", "C", "D"), each = nb.h), 1:nb.h)
  return(rot_fcoef)
}


##########################################################################################

#' Plot phenogram
#'
#' @description Plot phenogram using phylogenetic and morphometric information.
#'   Used internally.
#'
#' @param x Optional vector with a non-morphometric variable to be plotted in
#'   the x axis.
#' @param y Optional vector with a non-morphometric variable to be plotted in
#'   the y axis.
#' @param tree A \code{"phylo"} object containing a phylogenetic tree. Tip labels
#'   should match the row names from \code{x} or \code{y}.
#' @param phylo_scores A matrix containing the scores from tips and nodes of
#'   the phylogeny provided in \code{tree}.
#' @param axis Numeric; the axis to be plotted.
#' @param pch.groups Numeric; the symbol of the scatter points corresponding to
#'   groups mean shapes.
#' @param col.groups The color of the scatter points corresponding to groups
#'   mean shapes.
#' @param bg.groups The background color of the scatter points corresponding
#'   to groups mean shapes.
#' @param cex.groups Numeric; the size of the scatter points corresponding to
#'   groups mean shapes.
#' @param lwd.phylo Numeric; the width of the lines depicting phylogenetic
#'   branches.
#' @param lty.phylo Numeric; the type of the lines depicting phylogenetic
#'   branches.
#' @param col.phylo Numeric; the color of the lines depicting phylogenetic
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
#' sp_shapes <- expected_shapes(shapes, species)
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
#'                axis = 1, lwd.phylo = 1, lty.phylo = 1, col.phylo = 1, cex.groups = 1,
#'                col.groups = 1, pch.groups = 1, points = TRUE)
plot_phenogram <- function(x = NULL,
                           y = NULL,
                           tree,
                           phylo_scores,
                           axis,
                           lwd.phylo,
                           lty.phylo,
                           col.phylo,
                           cex.groups,
                           pch.groups,
                           col.groups,
                           bg.groups,
                           points) {

  for(i in seq_len(nrow(tree$edge))) {
    phyloxy <- cbind(x, phylo_scores[,axis[1]], y)
    graphics::lines(rbind(phyloxy[tree$edge[i, 1],],
                          phyloxy[tree$edge[i, 2],]),
                    lwd = lwd.phylo, lty = lty.phylo, col = col.phylo)
  }
  if(points == TRUE) {
    ntips <- length(tree$tip.label)
    graphics::points(phyloxy[-c(seq_len(ntips)),], pch = 16)

    if(any(pch.groups %in% c(21:25))) {
      graphics::points(phyloxy[seq_len(ntips),][rownames(phylo_scores)[seq_len(ntips)],],
                       bg = bg.groups, pch = pch.groups, cex = cex.groups)
    } else {
      graphics::points(phyloxy[seq_len(ntips),][rownames(phylo_scores)[seq_len(ntips)],],
                       col = col.groups, pch = pch.groups, cex = cex.groups)
    }

  }
}
