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
#' @seealso \code{\link{project_eigen}}
#' @export
#'
#' @examples
#' #perform principal component analysis of tails data
#' data("tails")
#' pca <- prcomp(geomorph::two.d.array(tails$shapes))
#'
#' #transform the scores back to shapes
#' backshapes_mat <- rev_eigen(scores = pca$x, vectors = pca$rotation, center = pca$center)
#' backshapes_arr <- geomorph::arrayspecs(backshapes_mat, k = 2, p = 9)
#'
#' #compare
#' pile_shapes(tails$shapes)
#' pile_shapes(backshapes_arr)
#'
#' #obtain shapes at the extremes of PC1
#' extshapes_mat <- rev_eigen(scores = range(pca$x[,1]), vectors = pca$rotation[,1], center = pca$center)
#' extshapes_arr <- geomorph::arrayspecs(extshapes_mat, k = 2, p = 9)
#'
#' #plot and compare
#' plot(extshapes_arr[,,1])
#' Morpho::lineplot(extshapes_arr[,,1], tails$links) ; title("negative")
#' plot(extshapes_arr[,,2])
#' Morpho::lineplot(extshapes_arr[,,2], tails$links) ; title("positive")
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
#' #perform principal component analysis of tails data
#' data("tails")
#' pca <- prcomp(geomorph::two.d.array(tails$shapes))
#'
#' #get project shapes in the pca ordination
#' shapes_mat <- geomorph::two.d.array(tails$shapes)
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
#' meanshapes_mat <- geomorph::two.d.array(meanshapes)
#' meanscores <- proj_eigen(x = meanshapes_mat,
#'                         vectors = pca$rotation,
#'                         center = pca$center)
#' points(meanscores, pch = 21, bg = "red", cex = 1.5)
proj_eigen <- function(x, vectors, center) { t(t(rbind(x)) - center) %*% vectors }


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
#'   \code{"phy_prcomp"} or \code{"pls"} object).
#' @param axes Numeric of length 1 or 2, indicating the morphometric axes to be
#'   plotted. If values for either \code{x} or \code{y} are provided, only the
#'   first value of this argument is considered.
#' @param datype type of shape data (either "fcoef" or "landm").
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
#' @examples
#' #perform pca on tails shapes
#' data("tails")
#' pca <- prcomp(geomorph::two.d.array(tails$shapes))
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

  sh_arr[,2,] <- sh_arr[,2,] * asp.models

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
    sh_arr <- abind::abind(temp_warpd_list, along = 3)
    p <- nrow(sh_arr)
  }


  models_mat<-c()
  models_arr<-sh_arr * 0
  for(i in 1:nrow(gridcoords)) {
    descentmat <- matrix(rep(gridcoords[i,], p), p, k, byrow = TRUE)
    models_mat <- rbind(models_mat, (sh_arr[,,i] * 0.07) + descentmat)
    models_arr[,,i] <- (sh_arr[,,i] * 0.07) + descentmat
  }

  return(list(models_mat = models_mat, models_arr = models_arr))

}




################################################################################################


#' Rotate Fourier shape 180 degrees
#'
#' @description Correct 180-degrees spurious rotation in closed outline shapes.
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

