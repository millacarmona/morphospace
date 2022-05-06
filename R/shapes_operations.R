
###########################################################################

#' Compute consensus shape(s)
#'
#' @description Compute the mean shape from the entire sample, the mean shape of
#'   a subset of samples, or the mean shape of the levels of a factor (for
#'   landmark configurations).
#'
#' @param shapes Shape data.
#' @param index Either a numeric vector indicating the configurations to be
#'   averaged or a factor whose levels are used to average groups of
#'   configurations.
#'
#' @return For landmark data, either a \code{p x k} matrix defining a single
#'   mean shape or a \code{p x k x n} array containing \code{n} mean shapes, one
#'   for each level of \code{index}. For Fourier data, a \code{n x (4 x nb.h)}
#'   matrix of Fourier coefficients (with \code{n} being either 1 or the number
#'   of levels of \code{index} and \code{nb.h} being the number of harmonics
#'   used during elliptic Fourier analysis).
#'
#' @export
#'
#' @examples
#' #load tails data
#' data("tails")
#'
#' #compute and plot mean shape of the entire sample
#' mshape <- consensus(shapes)
#' mshape
#' plot(mshape)
#' Morpho::lineplot(mshape, tails$links)
#'
#' #getting mean shape for levels of a factor: compute and plot mean shape of
#' #each of the 13 species
#' (index <- tails$data$species)
#' sp_mshapes <- consensus(shapes, index = index)
#' sp_mshapes
#' pile_shapes(sp_mshapes, links = tails$links, mshape = FALSE)
#' for(i in 1:13) Morpho::lineplot(sp_mshapes[,,i], tails$links)
#'
#' #getting mean shape for a subset of specimens: compute and plot mean shape of
#' #deep-forked species
#' (index <- which(tails$data$type == "DF"))
#' df_mshape <- consensus(shapes, index = index)
#' df_mmsape
#' plot(df_mshape)
#' Morpho::lineplot(df_mshape, tails$links)

#' #quick demo for Fourier data:
#' data("shells")
#'
#' #mean shape of the entire sample
#' mshape <- consensus(shells$shapes)
#' mshape
#' plot(inv_efourier(mshape, nb.pts = 200), type = "l")
#'
#' #mean shape of each of the four species
#' sp_mshapes <- consensus(shells$shapes, index = shells$data$species)
#' sp_mshapes
#' pile_shapes(sp_mshapes, mshape = FALSE)
consensus <- function(shapes, index = NULL) {

  dat <- shapes_mat(shapes)
  data2d <- dat$data2d
  datype <- dat$datype

  if(is.null(index)) index <- 1:nrow(data2d)

  if(is.numeric(index)) {
    cons <- rbind(colMeans(data2d[index,]))
  } else {
    if(is.character(index)) index <- factor(index)
    if(is.factor(index)) {
      cons <- apply(X = data2d, MARGIN = 2, FUN = tapply, index, mean)
      rownames(cons) <- levels(index)
    }
  }

  if(datype == "landm") {
    p <- nrow(shapes)
    k <- ncol(shapes)
    if(nrow(cons) == 1) {
      cons <- matrix(cons, nrow = p, byrow = TRUE)
    } else {
      cons <- geomorph::arrayspecs(cons, p = p, k = k)
    }
  }

  return(cons)
}


###########################################################################

#' Remove shape variation associated to external variables
#'
#' @description Detrend shape data using the functional relationship between
#'   shape data and some external explanatory variable(s) (works for both
#'   factors and numerics), estimated using a linear model.
#'
#' @param model A \code{mlm} object created using [lm()].
#' @param xvalue A value (numeric) or level (character) at which shape data is
#'   to be standardized (i.e. centered); If NULL, the mean of the complete
#'   sample is used.
#' @param newx,newy New data to be standardized instead of the used in
#'   \code{model}. Coefficients are taken from the linear model and applied to
#'   the new data to predict shapes at the desired \code{xvalue}.
#'
#' @details This function detrends (or standardizes, or corrects) shapes
#'   (landmarks or Fourier coefficients) from variation associated with
#'   non-shape variables, using a lm model fitting the former to the latter. It
#'   returns a 2-margins matrix of shapes (which means an extra step has to be
#'   taken for larndmark data in order to retrieve configurations in 3-margins
#'   array format) corresponding to the detrended versions of the shapes
#'   specified in the left size of the formula in \code{model}.
#'
#'   However, if \code{newx} and \code{newy} are provided, the function will
#'   instead return the shapes provided in \code{newy}, detrended using the
#'   relationship described in \code{model}.
#'
#'   The grand mean of the sample of shapes is used by default to center shape
#'   variation, although a \code{xvalue} specifying a level or numeric value of
#'   the explanatory variable in \code{model} to center shapes at can be
#'   provided. This 'displacement' can only be applied for one explanatory
#'   variable at the time.
#'
#' @return A 2-margins matrix, of dimensions \code{n x (k x p)} for the case of
#'   landmark data and \code{n x (4 x nb.h)} for the case of Fourer data (where
#'   nb. h is the number of harmonics used in elliptic Fourier analysis).
#'
#' @export
#'
#' @references Klingenberg, C. P. (2016). \emph{Size, shape, and form: concepts
#'   of allometry in geometric morphometrics}. Development Genes and Evolution,
#'   226(3), 113-137.
#'
#' @examples
#' #### Landmark data
#'
#' #load tails data
#' data(tails)
#' shapes <- tails$shapes
#' species <- tails$data$species
#' sex <- tails$data$sex
#' logsizes <- log(tails$sizes)
#' msp <- mspace(shapes, links = tails$links, points = TRUE)
#' hulls_by_group_2D(msp$x, fac = species)
#'
#'
#' ### For numeric variables
#'
#' #fit linear model between shapes and sizes, then center at the grand mean of the sample
#' model <- lm(two.d.array(shapes) ~ logsizes)
#' detr_shapes_mat <- detrend_shapes(model)
#'
#' detr_shapes_nosize <- geomorph::arrayspecs(detr_shapes_mat, k = 2, p = 9)
#'
#' msp_nosize <- mspace(detr_shapes_nosize, links = tails$links, points = TRUE)
#' hulls_by_group_2D(msp_nosize$x, fac = species)
#'
#' ## Using xvalue
#'
#' #fit linear model between shapes and sizes, then center at the shape corresponding to the
#' #maximum size of the sample
#' model <- lm(two.d.array(shapes) ~ logsizes)
#' detr_shapes_mat2 <- detrend_shapes(model,
#'                                    xvalue = max(logsizes))
#'
#' detr_shapes_nosize2 <- geomorph::arrayspecs(detr_shapes_mat2, k = 2, p = 9)
#'
#' msp_nosize2 <- mspace(detr_shapes_nosize2, links = tails$links, points = TRUE)
#' hulls_by_group_2D(msp_nosize2$x, fac = species)
#'
#' ## Using newx, newy
#'
#' #fit linear model between shapes and sizes for NDF species, then use the NDF allometry to
#' #detrend DF shapes from alometric variation
#' #maximum size of the sample
#' index <- tails$data$type == "NDF"
#' shapes_ndf <- shapes[,,index]
#' logsizes_ndf <- logsizes[index]
#' shapes_df <- shapes[,,!index]
#' logsizes_df <- logsizes[!index]
#'
#' model <- lm(geomorph::two.d.array(shapes_ndf) ~ logsizes_ndf)
#' detr_shapes_mat3 <- detrend_shapes(model,
#'                                    newx = logsizes_df,
#'                                    newy = geomorph::two.d.array(shapes_df))
#'
#' detr_shapes_nosize3 <- geomorph::arrayspecs(detr_shapes_mat3, k = 2, p = 9)
#'
#' msp_nosize3 <- mspace(detr_shapes_nosize3, links = tails$links, points = TRUE)
#' hulls_by_group_2D(msp_nosize3$x, fac = factor(species[!index]))
#'
#'
#' ### For factors
#'
#' #fit linear model between shapes and species, then center at the grand mean of the sample
#' model <- lm(two.d.array(shapes) ~ species)
#' detr_shapes_mat <- detrend_shapes(model)
#'
#' detr_shapes_nospp <- geomorph::arrayspecs(detr_shapes_mat, k = 2, p = 9)
#'
#' msp_nospp <- mspace(detr_shapes_nospp, links = tails$links, points = TRUE)
#' hulls_by_group_2D(msp_nospp$x, fac = species)
#'
#' ## Using xvalue
#'
#' #fit linear model between shapes and species, then center at the shape corresponding to the
#' #mean shape of T. savana
#' model <- lm(two.d.array(shapes) ~ species)
#' detr_shapes_mat2 <- detrend_shapes(model,
#'                                    xvalue = "T. savana")
#'
#' detr_shapes_nospp2 <- geomorph::arrayspecs(detr_shapes_mat2, k = 2, p = 9)
#'
#' msp_nospp2 <- mspace(detr_shapes_nospp2, links = tails$links, points = TRUE)
#' hulls_by_group_2D(msp_nospp2$x, fac = species)
#'
#' ## Using newx, newy
#'
#' #fit linear model between shapes and sex for NDF species, then use the NDF sexual dimorphism to
#' #detrend DF shapes from variation between sexes
#' #maximum size of the sample
#' index <- tails$data$type == "NDF"
#' shapes_ndf <- shapes[,,index]
#' sex_ndf <- sex[index]
#' shapes_df <- shapes[,,!index]
#' sex_df <- sex[!index]
#'
#' model_ndf <- lm(geomorph::two.d.array(shapes_ndf) ~ sex_ndf)
#' detr_shapes_mat3 <- detrend_shapes(model_ndf,
#'                                    newx = sex_df,
#'                                    newy = geomorph::two.d.array(shapes_df))
#'
#' detr_shapes_nosize3 <- geomorph::arrayspecs(detr_shapes_mat3, k = 2, p = 9)
#'
#' msp_nosize3 <- mspace(detr_shapes_nosize3, links = tails$links, points = TRUE)
#' hulls_by_group_2D(msp_nosize3$x, fac = factor(species[!index]))
#'
#' #### Fourier data (quick demo)
#'
#' #load shells data
#' data(shells)
#' shapes <- shells$shapes$coe
#' species <- shells$data$species
#' logsizes <- log(shells$sizes)
#' msp <- mspace(shapes, mag = 0.5, size.models = 0.4, asp.models = 0.5, points = TRUE)
#' hulls_by_group_2D(msp$x, fac = species)
#'
#' #fit linear model between shapes and sizes, then center at the shape corresponding to the grand mean
#' model <- lm(shapes ~ logsizes)
#' detr_shapes_nosize <- detrend_shapes(model)
#'
#' msp_nosize <- mspace(detr_shapes_nosize, mag = 0.5, size.models = 0.4, asp.models = 0.5, points = TRUE)
#' hulls_by_group_2D(msp_nosize$x, fac = species)
#'
#' #fit linear model between shapes and sizes, then center at the shape corresponding to the maximum size
#' model <- lm(shapes ~ logsizes)
#' detr_shapes_nosize2 <- detrend_shapes(model,
#'                                       xvalue = max(logsizes))
#'
#' msp_nosize2 <- mspace(detr_shapes_nosize2, mag = 0.5, size.models = 0.4, asp.models = 0.5, points = TRUE)
#' hulls_by_group_2D(msp_nosize2$x, fac = species)
#'
#' #fit linear model between shapes and sizes, then center at the maximum size
#' shapes_koeneni <- shapes[species == "koeneni",]
#' logsizes_koeneni <- logsizes[species == "koeneni"]
#' shapes_esbelta <- shapes[species == "esbelta",]
#' logsizes_esbelta <- logsizes[species == "esbelta"]
#'
#' model_koeneni <- lm(shapes_koeneni ~ logsizes_koeneni)
#' detr_shapes_nosize3 <- detrend_shapes(model_koeneni,
#'                                       newx = logsizes_esbelta,
#'                                       newy = shapes_esbelta)
#'
#' msp_nosize3 <- mspace(shapes_esbelta, mag = 0.5, size.models = 0.15, asp.models = 0.7, points = TRUE)
#' title("raw P. esbelta morphospace")
#' msp_nosize4 <- mspace(detr_shapes_nosize3, mag = 0.5, size.models = 0.15, asp.models = 0.7, points = TRUE)
#' title("P. esbelta morphospace, refined using \n allometric variation from P. koeneni")


detrend_shapes <- function(model, xvalue = NULL, newx = NULL, newy = NULL) {

  coefs <- model$coefficients
  resids <- model$resid

  x <- model$model[,ncol(model$model)]

  if(!is.null(newx) & is.null(newy)) stop("Provide values for newy")
  if(is.null(newx) & !is.null(newy)) stop("Provide values for newx")

  if(!is.null(newx) & !is.null(newy)) {
    designmat0 <- cbind(1, newx)
    resids <- newy - designmat0 %*% coefs
  }

  n <- nrow(resids)

  if(is.null(xvalue)) {
    grandmean <- colMeans(model$resid + model$fitted.values)

    grandmean_vec <- rep(1, n) %*% t(grandmean)
    predicted_mat <- resids + grandmean_vec

  } else {

    if(is.numeric(x) == TRUE) {
      designmat <- cbind(1, xvalue)
    }

    if(is.factor(x) == TRUE) {
      designmat <- rep(0, nlevels(x))
      designmat[which(levels(x) == xvalue)] <- 1

      if(isFALSE(designmat[1] == 1)) designmat[1] <- 1
    }

    fitted <- as.numeric(designmat %*% coefs)

    fitted_vec <- rep(1, n) %*% t(fitted)
    predicted_mat <- resids + fitted_vec
  }

  return(predicted_mat)

}


##################################################################################

#' Correct artificially rotated set of Fourier shapes interactively
#'
#' @description Choose and correct 180-degrees spurious rotation in a sample of
#'   closed outline shapes.
#'
#' @param ef An \code{"OutCoe"} object.
#' @param index An optional numeric vector providing the indices identifying the
#'   shapes to be rotated. If declared, overrides interative selection of
#'   shapes.
#'
#' @details The usage of this function is inspired in the SHAPE program for
#'   elliptic Fourier analysis (Iwata & Ukai 2002). The \code{index} argument is
#'   intended to facilitate its inclusion in scripts without having to manually
#'   select the outlines every time.
#'
#' @return An \code{"OutCoe"} object containing the corrected outlines.
#'
#' @export
#'
#' @references Iwata, H., & Ukai, Y. (2002). \emph{SHAPE: a computer program
#'   package for quantitative evaluation of biological shapes based on elliptic
#'   Fourier descriptors}. Journal of Heredity, 93(5), 384-385.
#'
#' @examples
#' #load shells data
#' data("shells")
#' ef <- shells$shapes
#'
#' #correct first 3 shapes automatically with index
#' ef_corr1 <- correct_efourier(ef, index = 1:3)
#' pile_shapes(ef_corr1, mshape = FALSE)
#'
#' #correct interactively (the process will remain open until
#' #you click 'finish' or press the Esc key)
#' \dontrun{
#'
#' ef_corr2 <- correct_efourier(ef_corr1)
#' }
correct_efourier<-function(ef, index = NULL) {

  orig_frame <- par("mar", "oma")
  on.exit(par(orig_frame))
  par(mar = c(1,1,1,1), oma = c(0,0,0,0))

  if(is.null(index)) print("Click Finish (top-right corner of plot pane) or enter <Esc> to finish selection")

  options(warn = -1)

  outs_coords<-array(0, c(300, 2, nrow(ef$coe)))
  for(i in 1:nrow(ef$coe)) outs_coords[,,i] <- inv_efourier(ef$coe[i,], 300)
  outs <- Momocs::Out(outs_coords)

  pos <- Momocs::coo_listpanel(
    Momocs::coo_template(outs$coo))


  if(is.null(index)) {

    n <- length(ef)
    sel <- rep(FALSE, length(ef))
    while(sum(sel) < n) {
      index <- identify(pos[!sel,], n = 1, plot = FALSE)
      if(!length(index)) break
      index <- which(!sel)[index]
      sel[index] <- TRUE

      ef$coe[index,] <- rotate_fcoef(ef$coe[index,])

      outs_coords[,,index] <- inv_efourier(ef$coe[index,], 300)
      new_outs <- Momocs::Out(outs_coords)

      Momocs::coo_listpanel(
        Momocs::coo_template(new_outs$coo))

    }
  } else {

    ef$coe[index,] <- rotate_fcoef(ef$coe[index,])

  }

  options(warn=0)
  return(ef)

}




