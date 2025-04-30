################################################################################

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
#'
#' @export
#' @keywords internal
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


################################################################################

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
#' @keywords internal
#'
#' @seealso \code{\link{rev_eigen}}
#'
#' @examples
#' #load data and packages
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


################################################################################

#' Singular value decomposition for 2 blocks of variables
#'
#' @description Just a wrapper for [svd()] that returns an adequate output
#'   when used on blocks of variables. Can deal with phylogenetic data too
#'   using \code{ape}, \code{phytools} and \code{mvMORPH} functions.
#'   Used internally.
#'
#' @param x First block of variables.
#' @param y Second block of variables.
#' @param tree An optional \code{"phylo"} object containing a phylogenetic
#'   tree whose tip.labels match rownames of \code{x} and \code{y}.
#' @param evmodel Character, specifying an evolutionary model to perform
#'   ancestral character reconstruction; options are "BM" (Brownian motion),
#'   "EB" (Early burst) and "lambda" (Pagel's lambda transformation) (see
#'   \code{\link[mvMORPH]{mvgls}} for more details).
#'
#' @export
#' @keywords internal
#'
#' @return Mimics the \code{\link{svd}} output.
#'
#' @examples
#' #load data and packages
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
#' svd_block(x = sizes, y = two.d.array(shapes), evmodel = "BM")
#'
#' #perform partial svd on phylogenetic structure
#' svd_block(x = sp_sizes, y = two.d.array(sp_shapes),
#'           tree = tree, evmodel = "BM")
svd_block <- function(x, y, tree = NULL, evmodel) {

  x <- cbind(x)
  y <- cbind(y)

  if(!is.null(tree)) {
    adj_tree <- mvMORPH::mvgls(y ~ x, tree = tree, model = evmodel)$corrSt$phy
    C <- ape::vcv.phylo(adj_tree)[rownames(x), rownames(x)]
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


################################################################################

#' Inverse Fourier transform
#'
#' @description A wrapper for [Momocs::efourier_i()] to transform a
#'   set of Fourier coefficients into (x,y) coordinates. Used internally.
#'
#' @param coe A vector with Fourier coefficients.
#' @param nb.pts Numeric, specifying the number of coordinates for sampling the
#'   outlines.
#'
#' @return A \code{nb.pts x 2} matrix of (x,y) Cartesian coordinates defining
#'   single outline shape.
#'
#' @export
#' @keywords internal
#'
#' @seealso \code{\link[Momocs]{efourier_i}}
#'
#' @references
#' Bonhomme, V., Picq, S., Gaucherel, C., & Claude, J. (2014). \emph{Momocs:
#'   outline analysis using R}. 56(13). <https://www.jstatsoft.org/v56/i13/>.
#'
#' @examples
#' #load data and extract the first outline
#' data("shells")
#' shape_coe <- shells$shapes$coe[1,]
#'
#' #get and plot (x,y) coordinates for using different number of coordinates
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


################################################################################

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
#'   \item \code{$datype:} the type of geometric morphometrics data.
#'   \item \code{$data2d:} the shape descriptors arranged in 2-margins matrix
#'   format.
#'  }
#'
#' @export
#' @keywords internal
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
  #if is just a vector, arrange as a one-row matrix (will be treated as "2D" data)
  if(is.null(dim(shapes))) {
    shapes <- rbind(shapes)
    rownames(shapes) <- NULL
  }

  if(inherits(shapes, "OutCoe")) { #if it's a Momocs object
    datype <- "fcoef" #then it's Fourier data
    data2d <- shapes$coe #extract coefficients
  } else { #otherwise....

    if(inherits(shapes, "matrix")) { #if it is a matrix,
      if(any(colnames(shapes)[1] == "A1", colnames(shapes)[1] == "A2")) { #and columns have coefficients names
        datype <- "fcoef" #then its Fourier data (already in "2D" format)
        data2d <- shapes
      } else { #if there are no names,
        if(ncol(shapes) <= 3) { #and the matrix has 3 columns or less
          datype <- "landm" #then it is a single landmark shape in "3D" format.
          data2d <- t(matrix(t(shapes))) #stretch into a single row
        } else { #if there are more than 3 columns
          datype <- "landm" #then is landmark data (already in "2D" format)
          data2d <- shapes
        }
      }
    } else { #if it's not a matrix
      if(inherits(shapes, "array")) { #but it is an array
        if(dim(shapes)[3] > 1) { #with several slides
          datype <- "landm" #then is landmark data in "3D" format.
          data2d <- geomorph::two.d.array(shapes) #transform into "2D" data
        } else { #but if it only has one slide,
          datype <- "landm" #then is a single landmark shape in "3D" format.
          data2d <- t(matrix(t(shapes[,,1]))) #stretch into a single row
        }
      }
    }
  }
  return(list(data2d = data2d, datype = datype))
}


################################################################################

#' Adjust aspect and scale of background shape models for 2D data
#'
#' @description Avoid background shape models distortion caused by differences
#'   in ranges of x and y axes. Used internally.
#'
#' @param models An array containing the background shape models.
#' @param frame The frame in which shape models are to be plotted.
#' @param model_width Numeric; the width of a reference shape model (usually the
#'   consensus).
#' @param model_height Numeric; the height of a reference shape model (usually
#'   the consensus).
#'
#' @export
#' @keywords internal
#'
#' @return An array containing the adjusted shape models.
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
#' wh <- abs(apply(apply(expected_shapes(shapes_grid0$models_arr), 2, range), 2,
#'           diff))
#'
#' #ajust grid to new frame and plot
#' adj_grid <- adjust_models2d(models = shapes_grid0$models_arr,
#'                             frame = newframe, model_width = wh[1],
#'                             model_height = wh[2])
#' for(i in 1:dim(adj_grid)[3]) points(adj_grid[,,i], type = "l", col = "blue",
#'                                     lwd = 2)
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


################################################################################

#' Adjust aspect and scale of background shape models for 3D data
#'
#' @description Avoid background shape models distortion caused by differences
#'   in ranges of x and y axes. Used internally.
#'
#' @param models A list containing PNG images depicting background shape models.
#' @param frame The frame in which shape models are to be plotted.
#' @param size.models Numeric; size factor for shape models.
#' @param asp.models Numeric; the y/x aspect ratio of shape models.
#'
#' @return A list containing the adjusted frame of each element in
#'   \code{models}, and new, adjusted limits for the x and y axes from
#'   \code{frame}.
#'
#' @export
#' @keywords internal
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
#' if (interactive()) {
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
#' adj <- adjust_models3d(models = models, frame = model_centers,
#'                        size.models = 1, asp.models = 1)
#'
#' #plot everything
#' plot(pca$x, xlim = adj$xlim, ylim = adj$ylim)
#' for(i in 1:length(models)) {
#'   rasterImage(models[[i]], adj$model_frames[[i]][1,1],
#'               adj$model_frames[[i]][1,2], adj$model_frames[[i]][2,1],
#'               adj$model_frames[[i]][2,2])
#' }
#'
#' #kind of. Anyway, use plot_morphogrid3d that go thru the full process.
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


################################################################################

#' Adapt format of foreign multivariate analyses
#'
#' @description Reorganize results from \code{\link[geomorph]{gm.prcomp}},
#'   \code{\link[Morpho]{groupPCA}}, \code{\link[Morpho]{pls2B}},
#'   \code{\link[phytools]{phyl.pca}} and \code{\link[mvMORPH]{mvgls.pca}} for
#'   their use in the mspace workflow. Used internally.
#'
#' @param ordination An object of class \code{"gm.prcomp"}, \code{"bgPCA"},
#'    \code{"PCA"}, \code{"pls2B"}, \code{"pls"}, \code{"phyl.pca"} or \code{"mvgls.pca"}.
#'
#' @return An object of equivalent class with scores, eigenvectors, eigenvalues/
#'    standard deviations, and original centroid arranged so they can be used
#'    in downstream \code{morphospace} operations.
#'
#' @export
#' @keywords internal
#'
#' @examples
#' #load data
#' data("tails")
#'
#' #perform PLS
#' pls <- Morpho::pls2B(x = tails$sizes, y = geomorph::two.d.array(tails$shapes))
#'
#' #adapt format
#' pls_adapted <- adapt_ordination(pls)
#'
#' #compare contents
#' names(pls)
#' names(pls_adapted)
adapt_ordination <- function(ordination) {

  if(class(ordination)[1] == "PCA") {
    ordination_ <- list()

    ordination_$x <- ordination$x
    ordination_$rotation <- ordination$rotation
    ordination_$sdev <- ordination$sdev
    ordination_$center <- ordination$center
    ordination_$raw_pca <- ordination

    ordination <- ordination_
    class(ordination) <- "PCA"
  }


  if(class(ordination)[1] == "gm.prcomp") {
    if(ordination$GLS) ordination$center <- c(ordination$center)
    class(ordination) <- "gm.prcomp"
  }

  if(class(ordination)[1] == "phyl.pca") {
    ordination$a <- c(ordination$a)

    names(ordination)[which(names(ordination) == "S")] <- "x"
    names(ordination)[which(names(ordination) == "Evec")] <- "rotation"
    names(ordination)[which(names(ordination) == "Eval")] <- "values"
    names(ordination)[which(names(ordination) == "a")] <- "center"
  }

  if(class(ordination)[1] == "bgPCA") {
    ordination$eigenvalues <- sqrt(ordination$eigenvalues)
    ordination$Scores <- cbind(ordination$Scores)

    names(ordination)[which(names(ordination) == "Scores")] <- "x"
    names(ordination)[which(names(ordination) == "groupPCs")] <- "rotation"
    names(ordination)[which(names(ordination) == "eigenvalues")] <- "sdev"
    names(ordination)[which(names(ordination) == "Grandmean")] <- "center"
    names(ordination)[which(names(ordination) == "groupmeans")] <- "grcenters"
  }

  if(class(ordination)[1] == "mvgls.pca") {
    if(!("center" %in% names(ordination))) stop("Please assign a $center slot containing the centroid/phylgenetic mean to the 'mvgls.pca' object (see examples in ?adapt_model)")
    ordination$values <- sqrt(ordination$values)

    names(ordination)[which(names(ordination) == "scores")] <- "x"
    names(ordination)[which(names(ordination) == "vectors")] <- "rotation"
    names(ordination)[which(names(ordination) == "values")] <- "sdev"
  }


  if(any(class(ordination)[1] == c("pls2B", "pls"))) {
    if(class(ordination)[1] == "pls2B") {
      ordination$xrotation <- ordination$svd$u
      ordination$yrotation <- ordination$svd$v
    } else {
      ordination$xrotation <- ordination$left.pls.vectors
      ordination$yrotation <- ordination$right.pls.vectors
    }
    ordination$sdev <- sqrt(ordination$svd$d)
    ordination$svd <- NULL

    class <- if(class(ordination)[1] == "pls2B") "pls2B" else "pls"

    if(class(ordination)[1] == "pls2B") {
      names(ordination)[which(names(ordination) == "Xscores")] <- "xscores"
      names(ordination)[which(names(ordination) == "Yscores")] <- "yscores"
    } else {
      names(ordination)[which(names(ordination) == "XScores")] <- "xscores"
      names(ordination)[which(names(ordination) == "YScores")] <- "yscores"
    }

    if(class(ordination)[1] == "pls2B") {
      ycenter <- ordination$ycenter
    } else {
      ycenter <- colMeans(ordination$A2.matrix)
    }


    ordination <- list(sdev = apply(ordination$yscores, 2, stats::sd),
                       rotation = ordination$yrotation,
                       x = ordination$yscores,
                       x2 = ordination$xscores,
                       center = ycenter,
                       totvar = ordination$ytotvar,
                       raw_pls = ordination)

    class(ordination) <- class
  }
  return(ordination)
}


################################################################################

#' Adapt format of different linear model functions
#'
#' @description Extract and format results from \code{\link[stats]{lm}},
#'   \code{\link[geomorph]{procD.lm}}, \code{\link[geomorph]{procD.pgls}},
#'   \code{\link[RRPP]{lm.rrpp}}, \code{\link[mvMORPH]{mvols}} and
#'   \code{\link[mvMORPH]{mvgls}} for their use in the mspace workflow and shape
#'   operations. Used internally.
#'
#' @param model An object of class \code{"mlm"}, \code{"procD.lm"},
#'    \code{"lm.rrpp"}, \code{"mvgls"} or \code{"mvols"}.
#'
#' @return A list containing the response and explanatory variables,
#'    coefficients, original centroid of response variables (or phylogenetic
#'    mean in the case of phylogenetic linear models), fitted values, and method
#'    used for linear model fitting.
#'
#' @export
#' @keywords internal
#'
#' @examples
#' #load data
#' data("tails")
#'
#' #perform PGLS
#' gmdf <- geomorph::geomorph.data.frame(shapes = expected_shapes(tails$shapes, tails$data$species),
#'                                       sizes = cbind(tapply(tails$sizes, tails$data$species, mean)),
#'                                       phy = tails$tree)
#'
#' #adapt PGLS format
#' mod <- geomorph::procD.pgls(shapes ~ sizes, phy = phy, data = gmdf)
#' mod_adapted <- adapt_model(mod)
#'
#' #compare contents
#' names(mod)
#' names(mod_adapted)
adapt_model <- function(model) {

  if(class(model)[1] == "mlm") {
    coefs <- model$coefficients
    y <- model$model[[1]]
    x <- cbind(model$model[-1])
    grandmean <- colMeans(y)
    fitted <- model$fitted.values
    modtype <- "ols"
  }

  if(class(model)[1] == "procD.lm") {
    y <- model$Y
    x <- model$data[-1]
    if("pgls.coefficients" %in% names(model)) {
      coefs <- model$pgls.coefficients
      grandmean <- model$pgls.mean
      fitted <- model$pgls.fitted
      modtype <- "pgls"
    } else {
      coefs <- model$coefficients
      grandmean <- colMeans(y)
      fitted <- model$fitted
      modtype <- "ols"
    }
  }

  if(any(class(model)[1] == c("mvgls", "mvols"))) {
    rnames <- rownames(stats::model.matrix(model))
    coefs <- model$coefficients
    rownames(coefs) <- colnames(model$variables$X)
    y <- model$variables$Y

    resptype <- attr(model$terms, "dataClasses")[-1]
    x <- data.frame(matrix(NA, nrow = nrow(model$variables$Y), ncol = length(resptype)))
    colnames(x) <- names(resptype) ; rownames(x) <- rownames(model$variables$X)
    for(i in seq_len(length(resptype))) {
      if(resptype[i] == "factor") {


        respindex <- which(grepl(x = colnames(model$variables$X)[-1], names(resptype)[i]))
        subX <- cbind(cbind(model$variables$X[,-1])[,respindex])
        if(is.null(colnames(subX))) colnames(subX) <- colnames(model$variables$X)[-1][i]

        response <- matrix("Intercept", ncol = 1, nrow = nrow(model$variables$Y))
        for(j in seq_len(ncol(subX))) {
          respiname_raw <- colnames(subX)[j]
          responsei <- cbind(subX[,j])
          levi <- gsub(x = respiname_raw, pattern = names(resptype)[i], replacement = "")
          response[which(responsei == 1)] <- levi

        }
        x[,i] <- factor(response)
      } else x[,i] <- model$variables$X[,colnames(model$variables$X) == names(resptype)[i]]
    }

    fitted <- model$fitted
    grandmean <- colMeans(model$fitted)

    modtype <- if(class(model)[1] == "mvgls") "pgls" else "ols"
    if(modtype == "pgls") {
      fitted <- fitted[rnames,]
      y <- y[rnames,]
      x <- x[rnames,]
      if(is.null(dim(x))) {
        x <- data.frame(x)
        colnames(x) <- names(resptype)
        rownames(x) <- rnames
      }
    }
  }


  if(class(model)[1] == "lm.rrpp") {
    if(model$LM$ols) {
      coefs <- model$LM$coefficients
      grandmean <- model$LM$mean
      fitted <- model$LM$fitted
      modtype <- "ols"
    }
    if(model$LM$gls) {
      coefs <- model$LM$gls.coefficients
      grandmean <- model$LM$gls.mean
      modtype <- "gls"
    }
    y <- model$LM$Y
    x <- model$LM$data[-1]
  }

  adapted <- list(coefs = coefs, grandmean = grandmean, fitted = fitted,
                  y = y, x = x, modtype = modtype)
  return(adapted)
}

################################################################################

#' Adapt format of Morphoscape objects
#'
#' @description Extract surface information from landscapes created
#'   with \code{Morphoscape} and adapt format to be used by
#'   [graphics::contour()] and [graphics::image()].
#'
#' @param obj An object of class \code{"kriged_surfaces"},
#'   \code{"wtd_lscp"}, \code{"poly_surf"} or \code{"multi_poly"}.
#'
#' @return A list containing the x and y values of the grid, and a matrix
#'   with the z values of each combination of x and y.
#'
#' @keywords internal
adapt_Morphoscape <- function(obj) {

  if(inherits(obj, "kriged_surfaces")) {
    grid <- obj$dataframes$grid[,1:3]

    if(ncol(obj$dataframes$grid) > 3)
      cat(paste("Only", colnames(obj$dataframes$grid)[3]), "landscape is plotted")
  }

  if(inherits(obj, "wtd_lscp")) {
    grid <- obj$Wprime$grid[,c("x", "y", "Z")]
  }

  if(inherits(obj, "poly_surf")) {
    grid <- obj$grid[,1:3]
  }

  if(inherits(obj, "multi_poly")) {
    grid <- obj[[1]]$grid[,1:3]
    if(names(obj$dataframes$grid)[1] > 1)
      cat(paste("Only", names(obj$dataframes$grid)[1], "landscape is plotted"))
  }


  x <- unique(grid[,1])
  y <- unique(grid[,2])
  z <- matrix(grid[,3], nrow = length(y), ncol = length(x))

  return(list(x = x, y= y, z = z))
}


################################################################################

#' Generate background shape models
#'
#' @description Calculate and arrange background shape models for morphospaces.
#'   Used internally.
#'
#' @param ordination An object containing and ordination, formatted in the style
#'   of the \code{"prcomp"} class
#' @param axes Numeric of length 1 or 2, indicating the morphometric axes to be
#'   plotted. If values for either \code{x} or \code{y} are provided, only the
#'   first value of this argument is considered.
#' @param datype Character; the type of shape data (either \code{"fcoef"} or
#'   \code{"landm"}).
#' @param rescale Logical; whether to re-scale background shape models so shape
#'    variation is shown more clearly.
#' @param template A 2-column matrix containing landmarks/semilandmarks followed
#'   by coordinates defining a curve or set of curves describing additional
#'   aspects of morphology, which will be warped using TPS interpolation to
#'   produce the set of background shell models (see
#'   \code{\link{build_template2d}}).
#' @param x Optional vector with a non-morphometric variable to be plotted in
#'   the x axis.
#' @param y Optional vector with a non-morphometric variable to be plotted in
#'   the y axis.
#' @param p Numeric, indicating the number of landmarks/semilandmarks used (for
#'   landmark data only).
#' @param k Numeric, indicating the number of Cartesian dimensions of
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
#'   \item \code{$models_mat:} a 2-column matrix with the (x,y) coordinates of
#'   all the shape models in the background.
#'   \item \code{$models_arr:} same as \code{models_mat} but in 3-margins array
#'   format.
#'   \item \code{$grid:} coordinates marking the centroid of each shape model.
#'   Intended for internal use.
#' }
#'
#' @export
#' @keywords internal
#'
#' @seealso \code{\link{plot_morphogrid2d}}, \code{\link{plot_morphogrid3d}}
#'
#' @references
#' MacLeod, N. (2009). \emph{Form & shape models}. Palaeontological Association
#'   Newsletter, 72(620), 14-27.
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
#' #plot grid from $models_mat and project each shape in models_arr
#' plot(shapes_grid$models_mat)
#' points(shapes_grid$models_arr[,,1], pch=16, col = 1)
#' points(shapes_grid$models_arr[,,2], pch=16, col = 2)
#' points(shapes_grid$models_arr[,,3], pch=16, col = 3)
#' points(shapes_grid$models_arr[,,4], pch=16, col = 4)
morphogrid <- function(ordination,
                       axes,
                       datype,
                       rescale = TRUE,
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

    if(length(axes) > 1) cat(paste0("\n", c("x","y")[which(c(!is.null(x), !is.null(y)))],
                                        " has been specified, axes[2] will be ignored"))
    axes <- axes[1]

    if(!is.null(x)) {
      if(is.null(xlim)) {
        plotframe_x <- range(as.numeric(x))
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
        plotframe_y <- range(as.numeric(y))
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

  if(!is.na(asp)) {
    adjframe <- adjust_asp(xlim = plotframe_x, ylim = plotframe_y, asp = asp)

    plotframe_x <- adjframe$xlim
    plotframe_y <- adjframe$ylim
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
    sh_originals <- sh_mat
    coords_l <- lapply(seq_len(nrow(sh_mat)), function(i) {
      sh <- inv_efourier(coe = sh_mat[i,], nb.pts = p)
    })
    sh_arr <- abind::abind(coords_l, along = 3)
  } else {
    sh_arr <- sh_originals <- geomorph::arrayspecs(sh_mat, p = p, k = k)

  }

  if(rescale) {
    for(i in seq_len(dim(sh_arr)[3])) sh_arr[,,i] <- sh_arr[,,i] / Morpho::cSize(sh_arr[,,i])
  }

  if(rot.models != 0) for(i in seq_len(dim(sh_arr)[3])) {
    sh_arr[,,i] <- rotate_coords(sh_arr[,,i], rot.models)
  }


  if(k < 3) {
    wh <- abs(apply(apply(expected_shapes(sh_arr), 2, range), 2, diff))
    sh_arr <- adjust_models2d(sh_arr, gridcoords, wh[1], wh[2])

    sh_arr[,2,] <- sh_arr[,2,] * asp.models

    for(i in seq_len(dim(sh_arr)[3])) sh_arr[,2,i] <- sh_arr[,2,i] * (wh[2] / wh[1])
  }
  sh_arr <- sh_arr * size.models


  if(!is.null(template)) {
    centroid <- matrix(rev_eigen(0, ordination$rotation[,1], ordination$center), ncol = 2, byrow = TRUE)
    temp_cent <- rbind(centroid,
                       Momocs::tps2d(template[-(seq_len(p)), ],
                                     template[(seq_len(p)), ],
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

  return(list(models_mat = models_mat, models_arr = models_arr, shapemodels = sh_originals))

}


################################################################################

#' Plot background 2D shape models
#'
#' @description Plot the output from [morphogrid()], for 2-dimensional
#'   morphometric data. Used internally.
#'
#' @param morphogrid An object containing the output of \code{morphogrid}.
#' @param x Optional vector with a non-morphometric variable to be plotted in
#'   the x axis.
#' @param y Optional vector with a non-morphometric variable to be plotted in
#'   the y axis.
#' @param links A list with the indices of the coordinates defining the
#'   wireframe (following the format used in \code{Morpho}).
#' @param template A 2-column matrix containing landmarks/semilandmarks followed
#'   by coordinates defining a curve or set of curves describing additional
#'   aspects of morphology, which will be warped using TPS interpolation to
#'   produce the set of background shell models (see
#'   \code{\link{build_template2d}}).
#' @param datype Character; type of shape data used (\code{"landm"} or
#'   \code{"fcoef"}).
#' @param ordtype Character; method used for multivariate ordination (options
#'   available are \code{"prcomp"}, \code{"gm.prcomp"}, \code{"PCA"},
#'   \code{"mvgls.pca"}, \code{"phyl.pca"}, \code{"bg_prcomp"}, \code{"bgPCA"},
#'   \code{"pls_shapes"}, \code{"phy_pls_shapes"}, \code{"pls2B"} and
#'   \code{"pls"}).
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
#' @keywords internal
#'
#' @return None
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
      if(any(c("prcomp", "PCA") %in% ordtype)) {
        xlab <- paste0("PC", axes[1])
      }
      if(any(c("bg_prcomp", "bgPCA") %in% ordtype)) {
        xlab <- paste0("bgPC", axes[1])
      }
      if(any(c("phy_prcomp", "phyl.pca") %in% ordtype)) {
        xlab <- paste0("phyPC", axes[1])
      }
      if(ordtype == "phyalign_comp") {
        xlab <- paste0("PAC", axes[1])
      }
      if(any(c("pls_shapes", "pls2B", "pls") %in% ordtype)) {
        xlab <- paste0("PLS-", axes[1])
      }
      if(ordtype == "phy_pls_shapes") {
        xlab <- paste0("phyPLS-", axes[1])
      }
      if(any(c("burnaby", "phy_burnaby", "gm.prcomp", "mvgls.pca") %in% ordtype)) {
        xlab <- paste0("Axis ", axes[1])
      }
    }
  }
  if(is.null(ylab)) {
    if(!is.null(y)) {
      ylab <- "y"
    } else {
      if(any(c("prcomp", "PCA") %in% ordtype)) {
        ylab <- paste0("PC", axes[2])
      }
      if(any(c("bg_prcomp", "bgPCA") %in% ordtype)) {
        ylab <- paste0("bgPC", axes[2])
      }
      if(any(c("phy_prcomp", "phyl.pca") %in% ordtype)) {
        ylab <- paste0("phyPC", axes[2])
      }
      if(ordtype == "phyalign_comp") {
        ylab <- paste0("PAC", axes[2])
      }
      if(any(c("pls_shapes", "pls2B", "pls") %in% ordtype)) {
        ylab <- paste0("PLS-", axes[2])
      }
      if(ordtype == "phy_pls_shapes") {
        ylab <- paste0("phyPLS-", axes[2])
      }
      if(any(c("burnaby", "phy_burnaby", "gm.prcomp", "mvgls.pca") %in% ordtype)) {
        ylab <- paste0("Axis ", axes[2])
      }
    }
  }

  if(plot) {

    xlim <- c(xlim[1] + (diff(xlim) * (1 - adj_frame[1]) / 2),
              xlim[2] - (diff(xlim) * (1 - adj_frame[1]) / 2))
    ylim <- c(ylim[1] + (diff(ylim) * (1 - adj_frame[2]) / 2),
              ylim[2] - (diff(ylim) * (1 - adj_frame[2]) / 2))

    graphics::plot(morphogrid$models_mat, type = "n", xlim = xlim, ylim = ylim ,
                   xlab = "", ylab = "", axes = FALSE)
    graphics::box()
    if(is.factor(x)) {
      nch <- sapply(levels(x), nchar)[which.max(sapply(levels(x), nchar))]
      cex.axis <- if(nch > 10) 10/nch else 1
      graphics::axis(side = 1, x, at = seq_len(nlevels(x)),
                     labels = levels(x), las = 2, cex.axis = cex.axis)
    } else {
      graphics::axis(side = 1)
      graphics::mtext(side = 1, line = 3, text = xlab)
    }

    if(is.factor(y)) {
      nch <- sapply(levels(y), nchar)[which.max(sapply(levels(y), nchar))]
      cex.axis <- if(nch > 10) 10/nch else 1
      graphics::axis(side = 2, y, at = seq_len(nlevels(y)),
                     labels = levels(y), las = 2, cex.axis = cex.axis)
    } else {
      graphics::axis(side = 2)
      graphics::mtext(side = 2, line = 3, text = ylab)
    }

    if(models) {

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


################################################################################

#' Plot background 3D shape models
#'
#' @description Plot the output from [morphogrid()], for 3-dimensional
#'   morphometric (landmark) data. Used internally.
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
#' @param ordtype Character; method used for multivariate ordination (options
#'   available are \code{"prcomp"}, \code{"gm.prcomp"}, \code{"PCA"},
#'   \code{"mvgls.pca"}, \code{"phyl.pca"}, \code{"bg_prcomp"}, \code{"bgPCA"},
#'   \code{"pls_shapes"}, \code{"phy_pls_shapes"}, \code{"pls2B"} and
#'   \code{"pls"}).
#' @param axes Numeric of length 2, indicating the axes to be plotted.
#' @param rotmat Optional rotation matrix for background shape models. If
#'   \code{NULL}, the user will be asked to define a preferred orientation.
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
#'   let the device get minimized passively). The process of morphospace
#'   generation is rather slow, specially if a mesh is provided for
#'   \code{template}, a large number of shape models is asked, and/or
#'   \code{alpha.models} value is lower than \code{1}.
#'
#' @export
#' @keywords internal
#'
#' @return None
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
#' if (interactive()) {
#' #plot grid (shape coordinates only)
#' plot_morphogrid3d(morphogrid = shapes_grid, refshape = meanshape,
#'                   ordtype = "prcomp", axes = c(1,2), col.ldm = 1,
#'                   cex.ldm = 1, col.models = 1, lwd.models = 1,
#'                   bg.models = "gray", size.models = 2, asp.models = 1)
#'
#' #get shape corresponding to shells3D$mesh_meanspec using
#' #geomorph::findMeanSpec, then get mesh corresponding to mean shape using
#' #Morpho::tps3d
#' meanspec_id<- findMeanSpec(shapes)
#' meanspec_shape <- shapes[,,meanspec_id]
#' meanmesh <- tps3d(x = shells3D$mesh_meanspec , refmat = meanspec_shape,
#'                   tarmat = meanshape)
#'
#' #plot grid (includinh mesh template)
#' plot_morphogrid3d(morphogrid = shapes_grid, template = meanmesh,
#'                   refshape = meanshape, ordtype = "prcomp", axes = c(1,2),
#'                   col.ldm = 1, cex.ldm = 1, col.models = 1, lwd.models = 1,
#'                   bg.models = "gray", size.models = 2, asp.models = 1)
#' }
plot_morphogrid3d <- function(x = NULL,
                              y = NULL,
                              morphogrid,
                              refshape,
                              template = NULL,
                              links = NULL,
                              ordtype,
                              axes,
                              rotmat = NULL,
                              xlim = NULL,
                              ylim = NULL,
                              xlab = NULL,
                              ylab = NULL,
                              adj_frame = c(1,1),
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
      if(any(c("prcomp", "PCA") %in% ordtype)) {
        xlab <- paste0("PC", axes[1])
      }
      if(any(c("bg_prcomp", "bgPCA") %in% ordtype)) {
        xlab <- paste0("bgPC", axes[1])
      }
      if(any(c("phy_prcomp", "phyl.pca") %in% ordtype)) {
        xlab <- paste0("phyPC", axes[1])
      }
      if(ordtype == "phyalign_comp") {
        xlab <- paste0("PAC", axes[1])
      }
      if(any(c("pls_shapes", "pls2B", "pls") %in% ordtype)) {
        xlab <- paste0("PLS-", axes[1])
      }
      if(ordtype == "phy_pls_shapes") {
        xlab <- paste0("phyPLS-", axes[1])
      }
      if(any(c("burnaby", "phy_burnaby", "gm.prcomp", "mvgls.pca") %in% ordtype)) {
        xlab <- paste0("Axis ", axes[1])
      }
    }
  }
  if(is.null(ylab)) {
    if(!is.null(y)) {
      ylab <- "y"
    } else {
      if(any(c("prcomp", "PCA") %in% ordtype)) {
        ylab <- paste0("PC", axes[2])
      }
      if(any(c("bg_prcomp", "bgPCA") %in% ordtype)) {
        ylab <- paste0("bgPC", axes[2])
      }
      if(any(c("phy_prcomp", "phyl.pca") %in% ordtype)) {
        ylab <- paste0("phyPC", axes[2])
      }
      if(ordtype == "phyalign_comp") {
        ylab <- paste0("PAC", axes[2])
      }
      if(any(c("pls_shapes", "pls2B", "pls") %in% ordtype)) {
        ylab <- paste0("PLS-", axes[2])
      }
      if(ordtype == "phy_pls_shapes") {
        ylab <- paste0("phyPLS-", axes[2])
      }
      if(any(c("burnaby", "phy_burnaby", "gm.prcomp", "mvgls.pca") %in% ordtype)) {
        ylab <- paste0("Axis ", axes[2])
      }
    }
  }

  wd <- tempdir()
  enter <- NULL
  while(is.null(enter)) {

    if(!is.null(template)) {
      refmesh <- template
      add <- TRUE
      rgl::plot3d(refmesh, col = bg.models, specular = "black", axes = FALSE, aspect = FALSE,
                  xlab = "", ylab = "", zlab = "", alpha = alpha.models)
    } else add <- FALSE

    rgl::plot3d(refshape, col = col.ldm, specular = "black", axes = FALSE, aspect = FALSE,
                xlab = "", ylab = "", zlab = "", type = "s", size = cex.ldm, add = add)

    for(l in seq_len(length(links))) rgl::lines3d(refshape[links[[l]],],
                                                  col = col.models, lwd = lwd.models)

    if(is.null(rotmat)) {
      cat("Preparing for snapshot: rotate mean shape to the desired orientation\n (don't close or minimize the rgl device).")
      enter <- readline("Press <Enter> in the console to continue:")
    } else {
      enter <- 1
    }

    rgl::rgl.snapshot(paste0(wd, "model", dim(morphogrid$models_arr)[3] + 1, ".png"))

    cat("\nThis can take a few seconds...")
  }

  for(i in seq_len(dim(morphogrid$models_arr)[3])) {

    if(!is.null(template)) {
      modelmesh <- Morpho::tps3d(x = refmesh , refmat = refshape, tarmat = morphogrid$models_arr[,,i])
      rgl::plot3d(modelmesh, col = bg.models, specular = "black", axes = FALSE, aspect = FALSE,
                  xlab = "", ylab = "", zlab = "", alpha = alpha.models)
    }

    rgl::plot3d(morphogrid$models_arr[,,i], col = col.ldm, specular = "black", axes = FALSE, aspect = FALSE,
                xlab = "", ylab = "", zlab = "", type = "s", size = cex.ldm, add = add)

    for(l in seq_len(length(links))) rgl::lines3d(morphogrid$models_arr[,,i][links[[l]],],
                                                  col = col.models, lwd = lwd.models)

    rgl::rgl.snapshot(paste0(wd, "model", i, ".png"))
  }
  cat("\nDONE.")

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


  if(plot) {
    graphics::plot(0, type = "n", xlim = new_xlim, ylim = new_ylim,
                   xlab = "", ylab = "", axes = FALSE)
    graphics::box()
    if(is.factor(x)) {
      nch <- sapply(levels(x), nchar)[which.max(sapply(levels(x), nchar))]
      cex.axis <- if(nch > 10) 10/nch else 1
      graphics::axis(side = 1, x, at = seq_len(nlevels(x)),
                     labels = levels(x), las = 2, cex.axis = cex.axis)
    } else {
      graphics::axis(side = 1)
      graphics::mtext(side = 1, line = 3, text = xlab)
    }

    if(is.factor(y)) {
      nch <- sapply(levels(y), nchar)[which.max(sapply(levels(y), nchar))]
      cex.axis <- if(nch > 10) 10/nch else 1
      graphics::axis(side = 2, y, at = seq_len(nlevels(y)),
                     labels = levels(y), las = 2, cex.axis = cex.axis)
    } else {
      graphics::axis(side = 2)
      graphics::mtext(side = 2, line = 3, text = ylab)
    }
    if(plot.models) {

      for(i in seq_len((length(models) - 1))) {
      graphics::rasterImage(models[[i]], model_frames[[i]][1,1], model_frames[[i]][1,2],
                  model_frames[[i]][2,1], model_frames[[i]][2,2])
      }
    }
  }
}



################################################################################

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
#' @keywords internal
#'
#' @references
#' Iwata, H., & Ukai, Y. (2002). \emph{SHAPE: a computer program package for
#'   quantitative evaluation of biological shapes based on elliptic Fourier
#'   descriptors}. Journal of Heredity, 93(5), 384-385.
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


################################################################################

#' Rotate x,y coordinates
#'
#' @description Rotate x,y coordinates by an arbitrary amount of degrees
#'
#' @param xy (x,y) coordinates
#' @param degrees Numeric; angle (in degrees) to rotate \code{(x,y)}
#'
#'
#' @return Matrix of rotated c(x,y) coordinates.
#'
#' @export
#' @keywords internal
#'
#' @examples
#' #load data
#' data(wings)
#'
#' shapes <- wings$shapes
#' links <- wings$links
#'
#' # rotate first shape from the data set
#' rot_shape <- rotate_coords(shapes[,,1], 180)
#'
#' #plot and compare
#' plot(shapes[,,1], pch = 1)
#' Morpho::lineplot(shapes[,,1], links)
#'
#' points(rot_shape, pch = 16, col = "red")
#' Morpho::lineplot(rot_shape, links, col = "red")
rotate_coords <- function(xy, degrees) {

  radians <- degrees * pi / 180

  x <- xy[,1]
  y <- xy[,2]

  cos_angle <- cos(radians)
  sin_angle <- sin(radians)

  rotated_x <- x * cos_angle - y * sin_angle
  rotated_y <- x * sin_angle + y * cos_angle

  rotated_coordinates <- cbind(x = rotated_x, y = rotated_y)
  return(rotated_coordinates)
}


################################################################################

#' Plot phenogram
#'
#' @description Plot phenogram using phylogenetic and morphometric information.
#'   Used internally.
#'
#' @param x Optional vector with a non-morphometric variable to be plotted in
#'   the x axis.
#' @param y Optional vector with a non-morphometric variable to be plotted in
#'   the y axis.
#' @param tree A \code{"phylo"} object containing a phylogenetic tree. Tip
#'   labels should match the row names from \code{x} or \code{y}.
#' @param phylo_scores A matrix containing the scores from tips and nodes of
#'   the phylogeny provided in \code{tree}.
#' @param axis Numeric; the axis to be plotted.
#' @param labels.tips Either logical, indicating whether to include tip labels,
#'   or a character string with the exact names of the tips whose labels
#'   should be included.
#' @param pch.tips Symbol of the scatterpoints representing the tips of the
#'   phylogeny.
#' @param col.tips Color of the scatterpoints representing the tips of the
#'   phylogeny.
#' @param bg.tips Background color of the scatterpoints representing the
#'   tips of the phylogeny.
#' @param cex.tips Numeric; size of the scatterpoints representing the tips
#'   of the phylogeny.
#' @param labels.nodes Either logical, indicating whether to include node
#'   labels, or a character string with the exact names of the nodes whose
#'   labels should be included (with the form of "node_n", e.g., "node_14"
#'   corresponds to the root of a tree with 13 tips).
#' @param pch.nodes Symbol of the scatterpoints representing the nodes of the
#'   phylogeny.
#' @param col.nodes Color of the scatterpoints representing the nodes of the
#'   phylogeny.
#' @param bg.nodes Background color of the scatterpoints representing the
#'   nodes of the phylogeny.
#' @param cex.nodes Numeric; size of the scatterpoints representing the nodes
#'   of the phylogeny.
#' @param lwd.phylo Integer; the width of the lines depicting phylogenetic
#'   branches.
#' @param lty.phylo Integer; the type of the lines depicting phylogenetic
#'   branches.
#' @param col.phylo Numeric; the color of the lines depicting phylogenetic
#'   branches.
#' @param points Logical; whether to plot the scatter points.
#'
#' @export
#' @keywords internal
#'
#' @return None
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
#' msp <- mspace(shapes, mag = 0.7, axes = c(1,2), cex.ldm = 0,
#'               plot = FALSE) %>%
#'   proj_shapes(shapes = shapes, col = c(1:13)[species], pch = 1,
#'               cex = 0.7) %>%
#'   proj_shapes(shapes = sp_shapes, pch = 21, bg = 1:13, cex = 2) %>%
#'   proj_phylogeny(shapes = sp_shapes, tree = tree)
#'
#' #get node heights
#' heights <- phytools::nodeHeights(tree)
#' node_heights <- NULL
#' for(i in 1:(length(tree$tip.label) + tree$Nnode)) node_heights[i] <- unique(heights[tree$edge == i])
#'
#' #plot simple phengram
#' plot(node_heights, msp$projected$phylo_scores[,1])
#' plot_phenogram(tree = tree, x = node_heights,
#'                phylo_scores = msp$projected$phylo_scores, axis = 1, lwd.phylo = 1,
#'                lty.phylo = 1, col.phylo = 1, cex.tips = 1, col.tips = 1,
#'                pch.tips = 1, cex.nodes = 1, col.nodes = 1, pch.nodes = 1,
#'                points = TRUE, labels.tips = NULL, labels.nodes = NULL)
plot_phenogram <- function(x = NULL,
                           y = NULL,
                           tree,
                           phylo_scores,
                           axis,
                           lwd.phylo,
                           lty.phylo,
                           col.phylo,
                           labels.nodes,
                           cex.nodes,
                           pch.nodes,
                           col.nodes,
                           bg.nodes,
                           labels.tips,
                           cex.tips,
                           pch.tips,
                           col.tips,
                           bg.tips,
                           points) {

  for(i in seq_len(nrow(tree$edge))) {
    phyloxy <- cbind(x, phylo_scores[,axis[1]], y)
    graphics::lines(rbind(phyloxy[tree$edge[i, 1],],
                          phyloxy[tree$edge[i, 2],]),
                    lwd = lwd.phylo, lty = lty.phylo, col = col.phylo)
  }
  if(points) {
    ntips <- length(tree$tip.label)
    plot_biv_scatter(phyloxy[-c(seq_len(ntips)),], bg = bg.nodes, pch = pch.nodes,
                     cex = cex.nodes, col = col.nodes)
    plot_biv_scatter(phyloxy[seq_len(ntips),][rownames(phylo_scores)[seq_len(ntips)],],
                     bg = bg.tips, pch = pch.tips, cex = cex.tips, col = col.tips)

  }
  add_labels(phyloxy[tree$tip.label,], labels.tips)
  add_labels(phyloxy[-seq_len(ntips),], labels.nodes)
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
#' @keywords internal
#'
#' @return None
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
    if(sum(fac == levels(fac)[i]) > 1) {
      x <- xy[fac == levels(fac)[i], 1]
      y <- xy[fac == levels(fac)[i], 2]
      hullp <- grDevices::chull(x = x, y = y)
      graphics::polygon(x[hullp], y[hullp], border = col[i],
                        col = grDevices::adjustcolor(col[i], alpha.f = alpha), lty = lty[i], ...)
    }
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
#' @keywords internal
#'
#' @return None
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
                          radius = sqrt(stats::qchisq(conflev, df = 2)),
                          draw = FALSE)
      graphics::polygon(ell, border = col[i],
                        col = grDevices::adjustcolor(col[i], alpha.f = alpha), lty = lty[i], ...)
    }
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
#' @keywords internal
#'
#' @return None
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
density_by_group_2D <- function(xy, fac, ax, alpha = 0.2, lwd = 1, lty = 1, plot = TRUE,
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

  if(plot) {
    for(i in seq_len(nlevels(fac))) {
      graphics::polygon(dens[[i]]$x, dens[[i]]$y / ymax, lwd = lwd, border = col[i],
                        lty = lty[i], col = grDevices::adjustcolor(col[i], alpha.f = alpha))
      graphics::abline(h = 0)
    }
  }
  return(invisible(list(ymax = ymax, dens = dens)))
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
#'
#' @noRd
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
#'
#' @noRd
plot_biv_scatter <- function(scores, col = 1, bg = 1, pch = 1, cex = 1, ...) {

  if(is.null(dim(scores))) scores <- rbind(scores)

  if(any(pch %in% c(21:25))) {
    graphics::points(scores, pch = pch, bg = bg, cex = cex, ...)
  } else {
    graphics::points(scores, pch = pch, col = col, cex = cex, ...)
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
#'
#' @noRd
plot_univ_landscape <- function(landscape, drawlabels, col, lwd) {

  if(drawlabels) {
    w.transp <-round(stats::quantile(x = 1:length(landscape$z), probs = c(0.25, 0.5, 0.75)))
    w.transp <- sort(c(w.transp - 1, w.transp, w.transp + 1))

    label_col <- col[w.transp[c(2,5,8)]]
    col[w.transp] <- NA
  }

  landscape$z <- landscape$z / max(landscape$z)
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

#' Plot bivariate landscape into morphospace
#'
#' @description An \code{ad hoc} wrapper for [graphics::contour()] /
#'   [graphics::.filled.contour()] /  [graphics::image()]. Used internally for
#'   projecting landscapes into bivariate morphospaces.
#'
#' @param landscape A list containing the results of [akima::interp()] or .
#'   [adapt_Morphoscape()]
#' @param display Either \code{"contour"} or \code{"filled.contour"}.
#' @param type Type of landscape to be plotted. Options are \code{"theoretical"}
#'   \code{"empirical"}, or \code{"Morphoscape"}.
#' @param levels Levels to be used to create the contours.
#' @param lty Integer; type of the lines depicting contours.
#' @param lwd Integer; width of the lines depicting contours.
#' @param col Colors used to represent the landscape contours.
#' @param drawlabels Logical; should the labels indicating the value of each
#'   surface contour be plotted?
#' @param alpha Transparency factor for filled contours.
#'
#' @noRd
plot_biv_landscape <- function(landscape, display, type, levels, lwd, lty, col, drawlabels, alpha) {
  if(display == "contour") {
    graphics::contour(landscape$x, landscape$y, landscape$z, levels = levels,
                      lwd = lwd, lty = lty, col = col, labels = round(levels, digits = 3),
                      drawlabels = drawlabels, add = TRUE)
    graphics::box()
  }

  if(display == "filled.contour") {
    if(any(c("theoretical", "Morphoscape") %in% type)) {
      graphics::.filled.contour(landscape$x, landscape$y, landscape$z, levels = levels,
                                col = grDevices::adjustcolor(col = col, alpha = alpha))
    }
    if("empirical" %in% type) {
      graphics::image(landscape$x, landscape$y, landscape$z, add = TRUE,
                      col = grDevices::adjustcolor(col = col, alpha = alpha))
    }
  }
}


################################################################################

#' Add labels to scatterplots
#'
#' @description A wrapper for [graphics::text()]. Used internally.
#'
#' @param xy A matrix with scatter point coordinates.
#' @param labels Either logical, indicating whether to add labels to all the
#'   points in the scatter plot (taken from row names), or character string
#'   containing specific names to be added.
#' @param col Either numeric or character, specifying the color(s) to be used.
#' @param ... Further arguments passed to [graphics::text()].
#'
#' @noRd
add_labels <- function (xy, labels = NULL, col = 1, ...) {
  if(!is.null(labels)) {
    if(is.character(labels)) {
      label.scores <- xy[labels, ]
      label.text <- labels
      graphics::text(rbind(label.scores), label.text, col = col[labels], ...)
    } else {
      if(labels) {
        label.scores <- xy
        label.text <- rownames(xy)
        graphics::text(rbind(label.scores), label.text, col = col, ...)
      }
    }
  }
}


################################################################################

#' Adjust aspect ratio of morphospaces
#'
#' @description Adjusts xlim and ylim so they reflect the desired aspect ratio.
#'   Used internally.
#'
#' @param xlim A vector of length 2 indicating the negative and positive
#'   extremes of values along the x axis.
#' @param ylim A vector of length 2 indicating the negative and positive
#'   extremes of values along the y axis.
#' @param asp Desired y/x proportions.
#'
#' @noRd
adjust_asp <- function(xlim, ylim, asp) {

  plotw <- graphics::par("pin")[1]
  ploth <- graphics::par("pin")[2]

  asp0 <- asp
  asp <- if(plotw > ploth) asp else (1 / asp)

  xperin <- diff(xlim) / plotw #x units per inch
  yperin <- diff(ylim) / ploth #y units per inch

  if(yperin > xperin) { #if more units per inch vertically, adjust xlim

    aspfac <- asp * (yperin * plotw) / (xperin * plotw)
    halfnewxlim <- (diff(xlim) * aspfac) / 2

    newxlim <- c(mean(xlim) - halfnewxlim,
                 mean(xlim) + halfnewxlim)
    newylim <- ylim

    #if correction misframes xlim, reset and adjust ylim instead
    if(min(xlim) <= min(newxlim) | max(xlim) >= max(newxlim)) {
      newxlim <- xlim

      aspfac <- (1 / asp) * (xperin * ploth) / (yperin * ploth)
      halfnewylim <- (diff(ylim) * aspfac) / 2

      newylim <- c(mean(ylim) - halfnewylim,
                   mean(ylim) + halfnewylim)
    }
  }

  if(yperin < xperin) { #if more units per inch horizontally, adjust ylim

    aspfac <- asp * (xperin * ploth) / (yperin * ploth)
    halfnewylim <- (diff(ylim) * aspfac) / 2

    newylim <- c(mean(ylim) - halfnewylim,
                 mean(ylim) + halfnewylim)
    newxlim <- xlim

    #if correction misframes ylim, reset and adjust xlim instead
    if(min(ylim) <= min(newylim) | max(ylim) >= max(newylim)) {
      newylim <- ylim

      aspfac <- (1 / asp) * (yperin * plotw) / (xperin * plotw)
      halfnewxlim <- (diff(xlim) * aspfac) / 2

      newxlim <- c(mean(xlim) - halfnewxlim,
                   mean(xlim) + halfnewxlim)
    }
  }

  return(list(xlim = newxlim, ylim = newylim))

}


################################################################################

#' Set layout for morphospaces
#'
#' @description Sets distribution of morphospaces, legends and scale bars when
#'   using \code{\link{plot_mspace}}. Used internally.
#'
#' @param legend Logical; whether to include legend for groups.
#' @param scalebar Logical; whether to include scale bars for landscapes.
#'
#' @noRd
set_layout <- function(legend = TRUE, scalebar = TRUE) {

  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))

  mainpar <- legpar <- scbpar <- NULL


  if(any(legend, scalebar)) {
    plot(0, type = "n", axes = F, xlab = "", ylab = "")

    mainpar <- list(fig = c(0, 0.7, 0, 1), mar = c(5, 4, 4, 1), new = TRUE)

    if(legend & !scalebar) legpar <- list(fig = c(0.7, 1, 0, 1), mar = c(5, 1, 4, 2), new = TRUE)
    if(!legend & scalebar) scbpar <- list(fig = c(0.7, 0.8, 0, 1), mar = c(5, 1, 4, 2), new = TRUE)

    if(all(legend, scalebar)) {
      legpar <- list(fig = c(0.7, 1, 0.5, 1), mar = c(1.5, 1, 4, 2), new = TRUE)
      scbpar <- list(fig = c(0.7, 0.8, 0, 0.5), mar = c(5, 1, 1.5, 2), new = TRUE)
    }

  } else {
    mainpar <- oldpar
  }
  return(invisible(list(mainpar = mainpar, legpar = legpar, scbpar = scbpar)))
}

