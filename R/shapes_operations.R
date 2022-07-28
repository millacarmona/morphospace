###########################################################################

#' Compute mean/expected shape(s)
#'
#' @description Compute the mean shape from the entire sample, or the shape(s)
#'   expected for one or more levels (factors) or values (numerics) of an external
#'   explanatory variable as fitted by a linear model.
#'
#' @param shapes Shape data.
#' @param x A vector or column vector containing a single explanatory variable
#'    (can be either a factor or a numeric). If \code{NULL}, the grand mean of
#'    the entire sample is computed.
#' @param xvalue One or more numeric value(s) or factor level(s) of \code{x} at
#'   which calculate expected shape(s). If \code{NULL}, all the value(s) or
#'   level(s) are used.
#' @param tree A \code{"phylo"} object containing a phylogenetic tree. Tip labels
#'   should match names in \code{x} and \code{shapes}.
#'
#' @details If a phylogenetic tree is supplied for interspecific shape data, the
#'   procedure is performed using the phylogenetically-corrected regression coefficients
#'   (see Revell, 2009) assuming a Brownian motion model of evolution.
#'
#' @return For landmark data, either a \code{p x k} matrix defining a single
#'   mean shape or a \code{p x k x n} array containing \code{n} mean shapes.
#'   For Fourier data, a \code{n x (4 x nb.h)} matrix of Fourier coefficients
#'   (with \code{nb.h} being the number of harmonics used during elliptic Fourier
#'   analysis).
#'
#' @export
#'
#' @references Revell, L. J. (2009). \emph{Size-correction and principal components
#'   for interspecific comparative studies}. Evolution, 63, 3258-3268.
#'
#' @examples
#' #load tails data and packages
#' library(Morpho)
#' library(Momocs)
#' data("tails")
#' shapes <- tails$shapes
#' sizes <- log(tails$sizes)
#' species <- tails$data$species
#' type <- tails$data$type
#' tree <- tails$tree
#'
#' #compute and plot mean shape of the entire sample
#' mshape <- expected_shapes(shapes)
#' plot(mshape)
#' lineplot(mshape, tails$links)
#'
#' #getting mean shapes for the levels of a factor: compute and plot the mean shape
#' #of each of the 13 species
#' sp_shapes <- expected_shapes(shapes, x = species)
#' pile_shapes(sp_shapes, links = tails$links, mshape = FALSE)
#'
#' #getting the mean shape for a specific level of a factor: compute and plot the
#' #mean shape of deep-forked specimens
#' df_shape <- expected_shapes(shapes, x = type, xvalue = "DF")
#' plot(df_shape)
#' lineplot(df_shape, tails$links)
#'
#' #getting the mean shape for a specific level of a factor, correcting for phylogeny:
#' #compute and plot mean the shape of deep-forked species
#' sp_type <- factor(c(tapply(as.character(type), species, unique)))
#' df_sp_shape <- expected_shapes(sp_shapes, x = sp_type, xvalue = "DF", tree = tree)
#' plot(df_sp_shape)
#' lineplot(df_sp_shape, tails$links)
#'
#' #getting the shapes expected for a covariate: compute and plot the shapes expected
#' #under the linear regression size on of shape
#' exp_shapes <- expected_shapes(shapes, x = sizes)
#' pile_shapes(exp_shapes, links = tails$links, mshape = FALSE)
#'
#' #getting the shape expected for specific values of a covariate: compute and plot the
#' #shapes expected at the maximum size
#' large_shape <- expected_shapes(shapes, x = sizes, xvalue = max(sizes))
#' plot(large_shape)
#' lineplot(large_shape, tails$links)
#'
#' #getting the shape expected for specific values of a covariate, correcting for phylogeny:
#' #compute and plot the shapes expected at the maximum size
#' sp_sizes <- c(tapply(sizes, species, mean))
#' large_sp_shape <- expected_shapes(sp_shapes, x = sp_sizes, xvalue = max(sp_sizes), tree = tree)
#' plot(large_sp_shape)
#' lineplot(large_sp_shape, tails$links)
#'
#'
#' #quick demo for Fourier data:
#' data("shells")
#' shapes <- shells$shapes
#'
#' #mean shape of the entire sample
#' mshape <- expected_shapes(shapes)
#' plot(inv_efourier(mshape, nb.pts = 200), type = "l")
#'
#' #mean shape of each of the four species
#' sp_shapes <- expected_shapes(shapes, x = shells$data$species)
#' pile_shapes(sp_shapes, mshape = FALSE)
#'
#' #mean shape of P. esbelta
#' esbelta_shape <- expected_shapes(shapes, x = shells$data$species, xvalue = "esbelta")
#' plot(inv_efourier(esbelta_shape, nb.pts = 200), type = "l")
#'
#' #shapes expected by the linear regression of size on shape
#' exp_shapes <- expected_shapes(shapes, x = shells$sizes)
#' pile_shapes(exp_shapes, mshape = FALSE)
#'
#' #shapes expected at the minimum size
#' large_shape <- expected_shapes(shapes, x = shells$sizes, xvalue = min(shells$sizes))
#' plot(inv_efourier(large_shape, nb.pts = 200), type = "l")
expected_shapes <- function(shapes, x = NULL, xvalue = NULL, tree = NULL) {

  dat <- shapes_mat(shapes)
  data2d <- dat$data2d
  datype <- dat$datype

  if(is.null(rownames(data2d))) rownames(data2d) <- seq_len(nrow(data2d))
  nams <- rownames(data2d)

  if(is.null(x)) {
    designmat <- stats::model.matrix(~ 1, data = data.frame(data2d))
  } else {
    if(!is.null(dim(x)) ) {if(ncol(x) > 1) stop("Multiple explanatory variables are not allowed")}
    designmat <- stats::model.matrix(~ x, data = data.frame(data2d))
  }


  if(is.null(tree)) {
    coefs <- solve(t(designmat) %*% designmat) %*% t(designmat) %*% data2d
  } else {
    designmat <- cbind(designmat[tree$tip.label,])
    data2d <- cbind(data2d[tree$tip.label,])
    C <- ape::vcv.phylo(tree)
    coefs <- solve(t(designmat) %*% solve(C) %*% designmat) %*%
      t(designmat) %*% solve(C) %*% data2d
  }


  if(!is.null(x) & !is.null(xvalue)) {
    if(is.numeric(x) == TRUE) {
      designmat <- cbind(1, xvalue)
    }
    if(is.factor(x) == TRUE) {
      designmat <- matrix(rep(0, nlevels(x) * length(xvalue)),
                          nrow = length(xvalue), byrow = TRUE)
      for(i in seq_len(length(xvalue))) designmat[i, which(levels(x) == xvalue[i])] <- 1
      designmat[,1] <- 1
      rownames(designmat) <- xvalue
    }
  }

  predicted_mat <- designmat %*% coefs
  if(!is.null(tree) & is.null(xvalue)) predicted_mat <- predicted_mat[nams,]


  if(is.null(x)) {
    predicted_mat <- unique(predicted_mat)
  } else {
    if(is.factor(x)) {
      if(is.null(xvalue)) {
        predicted_mat <- unique(predicted_mat[order(x),])
        rownames(predicted_mat) <- levels(x)
      } else {
        rownames(predicted_mat) <- xvalue
      }
    }
  }

  if(datype == "landm") {
    p <- nrow(shapes)
    k <- ncol(shapes)
    predicted_mat <- geomorph::arrayspecs(predicted_mat, p = p, k = k)
    if(dim(predicted_mat)[3] == 1) predicted_mat <- predicted_mat[,,1]
  }

  return(predicted_mat)
}


###########################################################################

#' Remove shape variation associated to external variables
#'
#' @description Detrend shape data using the functional relationship between
#'   shape data and some external explanatory variable(s) (works for both
#'   factors and numerics), estimated using a linear model (and potentially
#'   corrected using phylogenetic relationships).
#'
#' @param model A \code{"mlm"} object created using [stats::lm()].
#' @param xvalue A value (numeric) or level (character) at which shape data is
#'   to be standardized (i.e. centered); If NULL, the mean of the complete
#'   sample is used.
#' @param tree A \code{"phylo"} object containing a phylogenetic tree. Tip labels
#'   should match row names from \code{x}.
#' @param newdata New data to be standardized instead of the used in \code{model},
#'   provided as a \code{"mlm"} object created using [stats::lm()]. Explanatory
#'   variables must be the same as those from \code{model}. Coefficients are taken
#'   from the linear model from \code{model} and applied to the new data to predict
#'   shapes at the desired \code{xvalue}.
#'
#' @details This function detrends (or standardizes, or corrects) shapes
#'   (landmarks or Fourier coefficients) from variation associated with
#'   non-shape variables, using a lm model fitting the former to the latter. It
#'   returns a 2-margins matrix of shapes (which means an extra step has to be
#'   taken for larndmark data in order to retrieve configurations in 3-margins
#'   array format) corresponding to the detrended versions of the shapes
#'   specified in the left size of the formula in \code{model}.
#'
#'   However, if \code{newdata} is provided, the function will instead return
#'   the shapes provided in the new model, detrended using the relationship described
#'   in \code{model}.
#'
#'   The grand mean of the sample of shapes is used by default to center shape
#'   variation, although a \code{xvalue} specifying a level or numeric value of
#'   the explanatory variable in \code{model} to center shapes at can be
#'   provided. This shift can only be applied for one explanatory variable
#'   at a time.
#'
#'   If a phylogenetic tree is supplied for interspecific shape data, the procedure
#'   is performed using the phylogenetically-corrected regression coefficients (and
#'   the phylogenetic mean is used instead of the grand mean for re-centering data;
#'   see Revell, 2009) assuming a Brownian motion model of evolution.
#'
#' @return A 2-margins matrix, of dimensions \code{n x (k x p)} for the case of
#'   landmark data and \code{n x (4 x nb.h)} for the case of Fourer data (where
#'   \code{nb.h} is the number of harmonics used in elliptic Fourier analysis).
#'
#' @export
#'
#' @references Revell, L. J. (2009). \emph{Size-correction and principal components
#'   for interspecific comparative studies}. Evolution, 63, 3258-3268.
#'
#' Klingenberg, C. P. (2016). \emph{Size, shape, and form: concepts
#'   of allometry in geometric morphometrics}. Development Genes and Evolution,
#'   226(3), 113-137.
#'
#' @examples
#'  #### Landmark data
#'
#'  #load tails data and packages
#'  library(geomorph)
#'  data(tails)
#'  shapes <- tails$shapes
#'  species <- tails$data$species
#'  sex <- tails$data$sex
#'  logsizes <- log(tails$sizes)
#'  msp <- mspace(shapes, links = tails$links, points = TRUE)
#'  hulls_by_group_2D(msp$x, fac = species)
#'
#'
#'  ### For numeric variables
#'
#'  #fit linear model between shapes and sizes, then center at the grand mean of the sample
#'  model <- lm(two.d.array(shapes) ~ logsizes)
#'  detr_shapes_mat <- detrend_shapes(model)
#'
#'  detr_shapes_nosize <- arrayspecs(detr_shapes_mat, k = 2, p = 9)
#'
#'  msp_nosize <- mspace(detr_shapes_nosize, links = tails$links, points = TRUE)
#'  hulls_by_group_2D(msp_nosize$x, fac = species)
#'
#'  ## using (phylogenetic) tree
#'
#'  #fit linear model between shapes and sizes, then center at the shape corresponding to the
#'  #maximum size of the sample
#'  sp_shapes <- expected_shapes(shapes, species)
#'  sp_logsizes <- c(tapply(logsizes, species, mean))
#'  model <- lm(two.d.array(sp_shapes) ~ sp_logsizes)
#'  detr_shapes_mat1 <- detrend_shapes(model)
#'
#'  detr_shapes_nosize1 <- arrayspecs(detr_shapes_mat1, k = 2, p = 9)
#'
#'  msp_nosize1 <- mspace(detr_shapes_nosize1, links = tails$links, points = TRUE)
#'  points(msp_nosize1$x, pch = 21, bg = c(1:13)[species])
#'
#'  ## Using xvalue
#'
#'  #fit linear model between shapes and sizes, then center at the shape corresponding to the
#'  #maximum size of the sample
#'  model <- lm(two.d.array(shapes) ~ logsizes)
#'  detr_shapes_mat2 <- detrend_shapes(model,
#'                                     xvalue = max(logsizes))
#'
#'  detr_shapes_nosize2 <- arrayspecs(detr_shapes_mat2, k = 2, p = 9)
#'
#'  msp_nosize2 <- mspace(detr_shapes_nosize2, links = tails$links, points = TRUE)
#'  hulls_by_group_2D(msp_nosize2$x, fac = species)
#'
#'  ## Using newdata
#'
#'  #fit linear model between shapes and sizes for NDF species, then use the NDF allometry to
#'  #detrend DF shapes from alometric variation
#'  #maximum size of the sample
#'  index <- tails$data$type == "NDF"
#'  shapes_ndf <- shapes[,,index]
#'  logsizes_ndf <- logsizes[index]
#'  shapes_df <- shapes[,,!index]
#'  logsizes_df <- logsizes[!index]
#'
#'  model_ndf <- lm(two.d.array(shapes_ndf) ~ logsizes_ndf)
#'  model_df <- lm(two.d.array(shapes_df) ~ logsizes_df)
#'  detr_shapes_mat3 <- detrend_shapes(model_ndf, newdata = model_df)
#'
#'  detr_shapes_nosize3 <- arrayspecs(detr_shapes_mat3, k = 2, p = 9)
#'
#'  msp_nosize3 <- mspace(detr_shapes_nosize3, links = tails$links, points = TRUE)
#'  hulls_by_group_2D(msp_nosize3$x, fac = factor(species[!index]))
#'
#'  ### For factors
#'
#'  #fit linear model between shapes and species, then center at the grand mean of the sample
#'  model <- lm(two.d.array(shapes) ~ species)
#'  detr_shapes_mat <- detrend_shapes(model)
#'
#'  detr_shapes_nospp <- arrayspecs(detr_shapes_mat, k = 2, p = 9)
#'
#'  msp_nospp <- mspace(detr_shapes_nospp, links = tails$links, points = TRUE)
#'  hulls_by_group_2D(msp_nospp$x, fac = species)
#'
#'  ## Using xvalue
#'
#'  #fit linear model between shapes and species, then center at the shape corresponding to
#'  #the mean shape of T. savana
#'  model <- lm(two.d.array(shapes) ~ species)
#'  detr_shapes_mat2 <- detrend_shapes(model,
#'                                     xvalue = "T. savana")
#'
#'  detr_shapes_nospp2 <- arrayspecs(detr_shapes_mat2, k = 2, p = 9)
#'
#'  msp_nospp2 <- mspace(detr_shapes_nospp2, links = tails$links, points = TRUE)
#'  hulls_by_group_2D(msp_nospp2$x, fac = species)
#'
#'  ## Using newdata
#'
#'  #fit linear model between shapes and sex for NDF species, then use the NDF sexual
#'  #dimorphism to detrend DF shapes from variation between sexes
#'  #maximum size of the sample
#'  index <- tails$data$type == "NDF"
#'  shapes_ndf <- shapes[,,index]
#'  sex_ndf <- sex[index]
#'  shapes_df <- shapes[,,!index]
#'  sex_df <- sex[!index]
#'
#'  model_ndf <- lm(two.d.array(shapes_ndf) ~ sex_ndf)
#'  detr_shapes_mat3 <- detrend_shapes(model_ndf, newdata = model_df)
#'
#'  detr_shapes_nosize3 <- arrayspecs(detr_shapes_mat3, k = 2, p = 9)
#'
#'  msp_nosize3 <- mspace(detr_shapes_nosize3, links = tails$links, points = TRUE)
#'  hulls_by_group_2D(msp_nosize3$x, fac = factor(species[!index]))
#'
#'  #### Fourier data (quick demo)
#'
#'  #load shells data
#'  data(shells)
#'  shapes <- shells$shapes$coe
#'  species <- shells$data$species
#'  logsizes <- log(shells$sizes)
#'  msp <- mspace(shapes, mag = 0.5, points = TRUE)
#'  hulls_by_group_2D(msp$x, fac = species)
#'
#'  #fit linear model between shapes and sizes, then center at the shape corresponding to
#'  #the grand mean
#'  model <- lm(shapes ~ logsizes)
#'  detr_shapes_nosize <- detrend_shapes(model)
#'
#'  msp_nosize <- mspace(detr_shapes_nosize, mag = 0.5, points = TRUE)
#'  hulls_by_group_2D(msp_nosize$x, fac = species)
#'
#'  #fit linear model between shapes and sizes, then center at the shape corresponding
#'  #to the maximum size
#'  model <- lm(shapes ~ logsizes)
#'  detr_shapes_nosize2 <- detrend_shapes(model,
#'                                        xvalue = max(logsizes))
#'
#'  msp_nosize2 <- mspace(detr_shapes_nosize2, mag = 0.5, points = TRUE)
#'  hulls_by_group_2D(msp_nosize2$x, fac = species)
#'
#'  #fit linear model between shapes and sizes, then center at the maximum size
#'  shapes_koeneni <- shapes[species == "koeneni",]
#'  logsizes_koeneni <- logsizes[species == "koeneni"]
#'  shapes_esbelta <- shapes[species == "esbelta",]
#'  logsizes_esbelta <- logsizes[species == "esbelta"]
#'
#'  model_koeneni <- lm(shapes_koeneni ~ logsizes_koeneni)
#'  model_esbelta <- lm(shapes_esbelta ~ logsizes_esbelta)
#'  detr_shapes_nosize3 <- detrend_shapes(model_koeneni, newdata = model_esbelta)
#'
#'  msp_nosize3 <- mspace(shapes_esbelta, mag = 0.5, points = TRUE)
#'  title("raw P. esbelta morphospace")
#'  msp_nosize4 <- mspace(detr_shapes_nosize3, mag = 0.5, points = TRUE)
#'  title("P. esbelta morphospace, refined using \n allometric variation from P. koeneni")
detrend_shapes <- function(model, xvalue = NULL, tree = NULL, newdata = NULL) {

  x <- model$model[-1]
  y <- model$model[[1]]

  if(is.null(rownames(x))) rownames(x) <- seq_len(nrow(x))
  if(is.null(rownames(y))) rownames(y) <- seq_len(nrow(y))
  namesx <- rownames(x)
  namesy <- rownames(y)

  designmat <- stats::model.matrix(stats::as.formula(paste0("~ ",
                                              strsplit(as.character(
                                                model$call[2]), split = "~")[[1]][2])), data = x)

  if(!is.null(tree)){
    if(!all(length(tree$tip.label) == nrow(x), length(tree$tip.label) == nrow(y))) {
      stop("Number of tips in the tree does not match the number of observations in x and/or y data sets")
    }
    if(!all(tree$tip.label %in% namesy, tree$tip.label %in% namesx)) {
      stop("Names in phylogenetic tree does not match names in x and/or y data sets")
    } else {
      x <- cbind(x[tree$tip.label,])
      y <- cbind(y[tree$tip.label,])
      designmat <- cbind(designmat[tree$tip.label,])
    }

    grandmean <- apply(y, 2, phytools::fastAnc, tree = tree)[1,]
    C <- ape::vcv.phylo(tree)
    coefs <- solve(t(designmat) %*% solve(C) %*% designmat) %*%
      t(designmat) %*% solve(C) %*% y

  } else {
    grandmean <- colMeans(y)
    coefs <- model$coefficients
  }

  which_na <- is.na(apply(coefs,1,sum))
  resids <- y - designmat[,!which_na] %*% coefs[!which_na,]

  if(!is.null(newdata)) {
    newx <- newdata$model[-1]
    newy <- newdata$model[[1]]
    if(is.null(rownames(newy))) rownames(newy) <- seq_len(nrow(newy))
    namesy <- rownames(newy)

    newdesignmat <- stats::model.matrix(stats::as.formula(paste0("~ ",
                                                   strsplit(as.character(newdata$call[2]),
                                                            split = "~")[[1]][2])),
                                 data = newx)

    which_na <- is.na(apply(coefs,1,sum))
    resids <- newy - newdesignmat[,!which_na] %*% coefs[!which_na,]

  }

  n <- nrow(resids)

  if(is.null(xvalue)) {

    grandmean_vec <- rep(1, n) %*% t(grandmean)
    predicted_mat <- resids + grandmean_vec

  } else {

    if(is.numeric(x[[1]]) == TRUE) {
      designmat <- cbind(1, xvalue)
    }

    if(is.factor(x[[1]]) == TRUE) {
      designmat <- rep(0, nlevels(x[[1]]))
      designmat[which(levels(x[[1]]) == xvalue)] <- 1

      if(isFALSE(designmat[1] == 1)) designmat[1] <- 1
    }

    fitted <- as.numeric(designmat %*% coefs)

    fitted_vec <- rep(1, n) %*% t(fitted)
    predicted_mat <- resids + fitted_vec
  }

  predicted_mat <- predicted_mat[namesy,]
  return(rbind(predicted_mat))

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

  orig_frame <- graphics::par("mar", "oma")
  on.exit(graphics::par(orig_frame))
  graphics::par(mar = c(1,1,1,1), oma = c(0,0,0,0))

  if(is.null(index)) print("Click Finish (top-right corner of Plots pane) or enter <Esc> in the console to finish selection")

  options(warn = -1)

  outs_coords<-array(0, c(300, 2, nrow(ef$coe)))
  for(i in seq_len(nrow(ef$coe))) outs_coords[,,i] <- inv_efourier(ef$coe[i,], 300)
  outs <- Momocs::Out(outs_coords)

  pos <- Momocs::coo_listpanel(
    Momocs::coo_template(outs$coo))


  if(is.null(index)) {

    n <- length(ef)
    sel <- rep(FALSE, length(ef))
    while(sum(sel) < n) {
      index <- graphics::identify(pos[!sel,], n = 1, plot = FALSE)
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


##################################################################################

#' Compute shapes at the extremes of a morphometric axis
#'
#' @description This function computes the theoretical shapes corresponding to
#'   the extremes of a morphometric axis, which can be supplied as an object
#'   containing a linear model or a multivariate ordination analysis.
#'
#' @param obj An object containing either a multivariate ordination of class
#'   \code{"prcomp", "bg_prcomp", "phy_prcomp", "pls_shapes"} or
#'   \code{"phy_pls_shapes"} or a \code{"mlm"} object fitted using [stats::lm()].
#' @param axis An optional numeric value specifying the axis of the multivariate
#'   ordination which is to be represented.
#' @param mag Numeric; magnifying factor for representing shape transformation.
#'
#' @return A 2-margins matrix, of dimensions \code{2 x (k x p)} for the case of
#'   landmark data and \code{2 x (4 x nb.h)} for the case of Fourer data (where
#'   \code{nb.h} is the number of harmonics used in elliptic Fourier analysis).
#'
#' @details If an object of class \code{"mlm"} fitting shape to a factor (only two
#'   levels allowed) is supplied in \code{obj}, magnification of the axis range
#'   is attained through bgPCA.
#'
#' @export
#'
#' @seealso \code{\link{expected_shapes}}, \code{\link{rev_eigen}}
#'
#' @references MacLeod, N. (2009). \emph{Form & shape models}. Palaeontological
#'   Association Newsletter, 72(620), 14-27.
#'
#' @examples
#' #load tail data and packages
#' library(geomorph)
#' data("tails")
#' shapes <- tails$shapes
#' links <- tails$links
#'
#' #perform PCA, compute and plot extreme shapes of PC1 at its natural range
#' pca <- prcomp(two.d.array(shapes))
#' extshapes2d <- ax_transformation(obj = pca, axis = 1, mag = 1)
#' extshapes <- arrayspecs(extshapes2d, k = ncol(shapes), p = nrow(shapes))
#' pile_shapes((extshapes), links = links, mshape = FALSE)
#'
#' #compute and plot extreme shapes of PC2 at its natural range
#' extshapes2d <- ax_transformation(obj = pca, axis = 2, mag = 1)
#' extshapes <- arrayspecs(extshapes2d, k = ncol(shapes), p = nrow(shapes))
#' pile_shapes(extshapes, links = links, mshape = FALSE)
#'
#' #compute and plot extreme shapes of PC2 magnified x2
#' extshapes2d <- ax_transformation(obj = pca, axis = 2, mag = 2)
#' extshapes <- arrayspecs(extshapes2d, k = ncol(shapes), p = nrow(shapes))
#' pile_shapes(extshapes, links = links, mshape = FALSE)
#'
#' #perform lm of shape on size, compute and plot extreme shapes at its natural range
#' model <- lm(two.d.array(tails$shapes) ~ log(tails$sizes))
#' extshapes2d <- ax_transformation(obj = model, mag = 1)
#' extshapes <- arrayspecs(extshapes2d, k = ncol(shapes), p = nrow(shapes))
#' pile_shapes(extshapes, links = links, mshape = FALSE)
#'
#' #perform lm of shape on size, compute and plot extreme shapes at its natural range
#' model <- lm(two.d.array(tails$shapes) ~ tails$data$sex)
#' extshapes2d <- ax_transformation(obj = model, mag = 1)
#' extshapes <- arrayspecs(extshapes2d, k = ncol(shapes), p = nrow(shapes))
#' pile_shapes(extshapes, links = links, mshape = FALSE)
#'
#' #perform lm of shape on size, compute and plot extreme shapes magnified x2
#' model <- lm(two.d.array(tails$shapes) ~ tails$data$sex)
#' extshapes2d <- ax_transformation(obj = model, mag = 2)
#' extshapes <- arrayspecs(extshapes2d, k = ncol(shapes), p = nrow(shapes))
#' pile_shapes(extshapes, links = links, mshape = FALSE)
ax_transformation <- function(obj, axis = 1, mag = 1) {

  if(any(class(obj) %in% c("prcomp", "bg_prcomp", "phy_prcomp",
                           "pls_shapes", "phy_pls_shapes"))) {
    extshapes_mat <- rev_eigen(range(obj$x[,axis]) * mag,
                               obj$rotation[,axis],
                               obj$center)
  }

  if(any(class(obj) %in% "mlm")) {

    x <- obj$model[,ncol(obj$model)]

    if(is.numeric(x)) {
      cent <- mean(range(x))
      halfdif <- diff(range(x)) / 2
      newrange <- c(cent - (halfdif * mag),
                    cent + (halfdif * mag))

      coefs <- obj$coefficients
      designmat <- cbind(1, newrange)
      extshapes_mat <- rbind(designmat %*% coefs)
    }

    if(is.factor(x)) {

      if(nlevels(x) > 2) stop("Only two levels are allowed for extracting axes from mlm objects; try with a bg_prcomp object")

      Y <- obj$model[,1]
      mshapes <- apply(X = Y, MARGIN = 2, FUN = tapply, x, mean)
      bgpca <- stats::prcomp(mshapes)
      bgax1 <- bgpca$x[,1]

      cent <- mean(range(bgax1))
      halfdif <- diff(range(bgax1)) / 2
      newrange <- c(cent - (halfdif * mag),
                    cent + (halfdif * mag))

      extshapes_mat <- rev_eigen(newrange, bgpca$rotation[,1], bgpca$center)

    }

  }

  return(rbind(extshapes_mat))

}


