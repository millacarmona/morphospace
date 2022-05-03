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
#' @return
#' @export
#'
#' @examples
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
#' @return
#' @export
#'
#' @examples
proj_eigen <- function(x, vectors, center) { t(t(rbind(x)) - center) %*% vectors }


###########################################################################

#' Sample shapes along a morphometric axis
#'
#' @param scores
#' @param vector
#' @param center
#' @param p
#' @param k
#'
#' @return
#' @export
#'
#' @examples
# sample_ax <- function(scores, vector, center, p, k) {
#
#   sh_mat <- rev_eigen(scores, vector, center)
#   sh_arr <- geomorph::arrayspecs(sh_mat, p = p, k = k)
#   return(sh_arr)
#
# }

##########################################################################

#' Iverse Fourier transform
#'
#' @description A wrapper for [efourier_i()] from \code{Momocs} to transform a
#'   set of Fourier coefficients into (x,y) coordinates. Used internally.
#'
#' @param coe A vector with Fourier coefficients.
#' @param nb.p Numeric, specifying the number of coordinates for sampling the
#'   outlines.
#'
#' @return
#' @export
#'
#' @examples
inv_efourier<-function(coe, nb.p = 120){

  nb.h <- length(coe) / 4
  coords <- matrix(0, nrow = nb.p, ncol = 2)

  a <- c(coe[1:nb.h])
  b <- c(coe[(nb.h + 1):(nb.h * 2)])
  c <- c(coe[((nb.h * 2) + 1):(nb.h * 3)])
  d <- c(coe[((nb.h * 3) + 1):(nb.h * 4)])

  coefs <- list(an = a, bn = b, cn = c, dn = d)

  coords <- Momocs::efourier_i(coefs, nb.pts = nb.p)
  return(coords)

}


################################################################################################

#' Identify and arrange shape descriptors
#'
#' @description Identify data type (Fourier coefficients or Procrustes shape
#'   coordinates) and arrange them in matrix format, suitable to be used as
#'   input for multivariate analysis. Used internally.
#'
#' @param shapes Shape data.
#'
#' @return
#' @export
#'
#' @examples
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
        if(colnames(shapes)[1] == "A1") {
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

#' Background shape models
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
#' @return
#' @export
#'
#' @examples
morphogrid <- function(ordination,
                       axes,
                       datype,
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

    if(length(axes) > 1) warning("x or y has been specified, axes[2] is ignored")
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
      inv_efourier(coe = sh_mat[i,], nb.p = p)
    })
    sh_arr <- abind::abind(coords_l, along = 3) * size.models
  } else {
    sh_arr <- geomorph::arrayspecs(sh_mat, p = p, k = k) * size.models
  }

  sh_arr[,2,] <- sh_arr[,2,] * asp.models

  if(rot.models!=0) for(i in 1:dim(sh_arr)[3]) {
    sh_arr[,,i] <- spdep::Rotation(sh_arr[,,i], rot.models*0.0174532925199)
  }

  models_mat<-c()
  models_arr<-sh_arr * 0
  for(i in 1:nrow(gridcoords)) {
    descentmat <- matrix(rep(gridcoords[i,], p), p, k, byrow = TRUE)
    models_mat <- rbind(models_mat, (sh_arr[,,i] * 0.07) + descentmat)
    models_arr[,,i] <- (sh_arr[,,i] * 0.07) + descentmat
  }

  return(list(models_mat = models_mat, models_arr = models_arr,
              plotframe = cbind(plotframe_x, plotframe_y)))

}


