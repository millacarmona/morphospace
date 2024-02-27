################################################################################

#' Compute mean/expected shape(s)
#'
#' @description Compute the mean shape from the entire sample, or the shape(s)
#'   expected for one or more levels (factors) or values (numerics) of an
#'   external explanatory variable as fitted by a linear model.
#'
#' @param shapes Shape data.
#' @param x A vector or column vector containing a single explanatory variable
#'    (can be either a factor or a numeric). If \code{NULL}, the grand mean of
#'    the entire sample is computed.
#' @param xvalue One or more numeric value(s) or factor level(s) of \code{x} at
#'   which calculate expected shape(s). If \code{NULL}, all the value(s) or
#'   level(s) are used.
#' @param tree A \code{"phylo"} object containing a phylogenetic tree. Tip
#'   labels should match names in \code{x} and \code{shapes}.
#' @param evmodel Character, specifying an evolutionary model to perform
#'   ancestral character reconstruction; options are "BM" (Brownian motion),
#'   "EB" (Early burst) and "lambda" (Pagel's lambda transformation) (see
#'   \code{\link[mvMORPH]{mvgls}} for more details).
#' @param returnarray Logical, indicating whether shapes should be returned
#'   in "3D" array format (landmark shapes only). Mostly intended for internal
#'   use.
#'
#' @details If a phylogenetic tree is supplied for interspecific shape data, the
#'   procedure is performed using the phylogenetically-corrected regression
#'   coefficients (see Revell, 2009) under different possible evolutionary
#'   models using \code{\link[mvMORPH]{mvgls}}.
#'
#' @return For landmark data, either a \code{p x k} matrix defining a single
#'   mean shape or a \code{p x k x n} array containing \code{n} mean shapes,
#'   unless \code{returnarray = TRUE} (in which case a \code{n x (p x k)} matrix
#'   will be returned. For Fourier data, a \code{n x (4 x nb.h)} matrix of
#'   Fourier coefficients (with \code{nb.h} being the number of harmonics used
#'   during elliptic Fourier analysis).
#'
#' @export
#'
#' @references
#' Revell, L. J. (2009). \emph{Size-correction and principal components for
#'   interspecific comparative studies}. Evolution, 63, 3258-3268.
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
#' #getting mean shapes for the levels of a factor: compute and plot the mean
#' #shape of each of the 13 species
#' sp_shapes <- expected_shapes(shapes, x = species)
#' pile_shapes(sp_shapes, links = tails$links, mshape = FALSE)
#'
#' #getting the mean shape for a specific level of a factor: compute and plot
#' #the mean shape of deep-forked specimens
#' df_shape <- expected_shapes(shapes, x = type, xvalue = "DF")
#' plot(df_shape)
#' lineplot(df_shape, tails$links)
#'
#' #getting the mean shape for a specific level of a factor, correcting for
#' #phylogeny: compute and plot mean the shape of deep-forked species
#' sp_type <- factor(c(tapply(as.character(type), species, unique)))
#' df_sp_shape <- expected_shapes(sp_shapes, x = sp_type, xvalue = "DF",
#'                                tree = tree)
#' plot(df_sp_shape)
#' lineplot(df_sp_shape, tails$links)
#'
#' #getting the shapes expected for a covariate: compute and plot the shapes
#' #expected under the linear regression size on of shape
#' exp_shapes <- expected_shapes(shapes, x = sizes)
#' pile_shapes(exp_shapes, links = tails$links, mshape = FALSE)
#'
#' #getting the shape expected for specific values of a covariate: compute and
#' #plot the shapes expected at the maximum size
#' large_shape <- expected_shapes(shapes, x = sizes, xvalue = max(sizes))
#' plot(large_shape)
#' lineplot(large_shape, tails$links)
#'
#' #getting the shape expected for specific values of a covariate, correcting
#' #for phylogeny: compute and plot the shapes expected at the maximum size
#' sp_sizes <- c(tapply(sizes, species, mean))
#' large_sp_shape <- expected_shapes(sp_shapes, x = sp_sizes,
#'                                   xvalue = max(sp_sizes), tree = tree)
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
#' esbelta_shape <- expected_shapes(shapes, x = shells$data$species,
#'                                  xvalue = "esbelta")
#' plot(inv_efourier(esbelta_shape, nb.pts = 200), type = "l")
#'
#' #shapes expected by the linear regression of size on shape
#' exp_shapes <- expected_shapes(shapes, x = shells$sizes)
#' pile_shapes(exp_shapes, mshape = FALSE)
#'
#' #shapes expected at the minimum size
#' large_shape <- expected_shapes(shapes, x = shells$sizes,
#'                                xvalue = min(shells$sizes))
#' plot(inv_efourier(large_shape, nb.pts = 200), type = "l")
expected_shapes <- function(shapes, x = NULL, xvalue = NULL,
                            tree = NULL, evmodel = "BM", returnarray = TRUE) {


  # shapes = shapes
  # x =logsizes#species
  # xvalue = NULL
  # tree = NULL
  # evmodel = "BM"
  # returnarray = TRUE



  dat <- shapes_mat(shapes)
  data2d <- dat$data2d
  datype <- dat$datype

  if(length(dim(x)) == 1) x <- cbind(x)
  #if(length(dim(x)) == 1 | is.null(dim(x))) x <- cbind(x) #!

  if(is.null(rownames(data2d))) rownames(data2d) <- seq_len(nrow(data2d))
  nams <- rownames(data2d)

  ################ old ########################
  if(is.null(x)) {
    designmat <- stats::model.matrix(~ 1, data = data.frame(data2d))
  } else {
    if(!is.null(dim(x))) if(ncol(x) > 1) stop("Multiple explanatory variables are not allowed")
    designmat <- stats::model.matrix(~ x, data = data.frame(data2d))
  }

  if(is.null(tree)) {
    coefs <- solve(t(designmat) %*% designmat) %*% t(designmat) %*% data2d
  } else {
    #!
    # designmat <- cbind(designmat[tree$tip.label,])
    # data2d <- cbind(data2d[tree$tip.label,])
    # C <- ape::vcv.phylo(tree)
    # coefs <- solve(t(designmat) %*% solve(C) %*% designmat) %*%
    #   t(designmat) %*% solve(C) %*% data2d
    if(is.null(x)) mvmod <- mvMORPH::mvgls(data2d ~ 1, tree = tree, model = evmodel) else {
      if(!is.null(dim(x))) if(ncol(x) > 1) stop("Multiple explanatory variables are not allowed")
      mvmod <- mvMORPH::mvgls(data2d ~ x, tree = tree, model = evmodel)
    }

    designmat <- cbind(stats::model.matrix(mvmod)[tree$tip.label,])
    data2d <- cbind(data2d[tree$tip.label,])
    C <- ape::vcv.phylo(mvmod$corrSt$phy)[rownames(designmat), rownames(designmat)]
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
    rownames(predicted_mat) <- NULL
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

  # ################# new #################
  # if(is.null(x))  mvmod <- stats::lm(data2d ~ 1) else {
  #   if(!is.null(dim(x))) if(ncol(x) > 1) stop("Multiple explanatory variables are not allowed")
  #   mvmod <- stats::lm(data2d ~ x)
  #   designmat <- stats::model.matrix(mvmod)
  #   coefs <- mvmod$coefficients
  # }
  #
  # if(!is.null(tree)) {
  #   if(is.null(x)) mvmod <- mvMORPH::mvgls(data2d ~ 1, tree = tree, model = evmodel) else {
  #     if(!is.null(dim(x))) if(ncol(x) > 1) stop("Multiple explanatory variables are not allowed")
  #     mvmod <- mvMORPH::mvgls(data2d ~ x, tree = tree, model = evmodel)
  #   }
  #
  #   designmat <- cbind(stats::model.matrix(mvmod)[tree$tip.label,])
  #   data2d <- cbind(data2d[tree$tip.label,])
  #   C <- ape::vcv.phylo(mvmod$corrSt$phy)[rownames(designmat), rownames(designmat)]
  #   coefs <- solve(t(designmat) %*% solve(C) %*% designmat) %*%
  #     t(designmat) %*% solve(C) %*% data2d
  #
  # }
  # predicted_mat <- mvmod$fitted
  #
  # if(!is.null(x) & !is.null(xvalue)) {
  #   if(is.numeric(x) == TRUE) {
  #     designmat <- cbind(1, xvalue)
  #   }
  #   if(is.factor(x) == TRUE) {
  #     designmat <- matrix(rep(0, nlevels(x) * length(xvalue)),
  #                         nrow = length(xvalue), byrow = TRUE)
  #     for(i in seq_len(length(xvalue))) designmat[i, which(levels(x) == xvalue[i])] <- 1
  #     designmat[,1] <- 1
  #     rownames(designmat) <- xvalue
  #   }
  #
  #   predicted_mat <- designmat %*% coefs
  # }
  # if(!is.null(tree) & is.null(xvalue)) predicted_mat <- predicted_mat[nams,]
  #
  # predicted_mat <- if(is.null(x)) colMeans(predicted_mat) else {
  #   if(is.factor(x)) apply(predicted_mat,2,tapply,INDEX=x,mean) else predicted_mat
  # }
  # ##########################################


  if(datype == "landm" & returnarray) {
    if(length(dim(shapes)) == 3) {
      p <- nrow(shapes)
      k <- ncol(shapes)
      predicted_mat <- geomorph::arrayspecs(predicted_mat, p = p, k = k)
      if(dim(predicted_mat)[3] == 1) predicted_mat <- predicted_mat[,,1]
    } else warning("landmark data has been provided as a 2-margin matrix; expected shapes cannot be returned as an array")
  }

  return(predicted_mat)
}


################################################################################

#' Remove shape variation associated to external variables
#'
#' @description Detrend shape variation using the relationship between shape
#'   data and some external explanatory variable(s) (works for both factors and
#'   numerics).
#'
#' @param model An object of class \code{"mlm"}, \code{"procD.lm"},
#'   \code{"lm.rrpp"}, \code{"mvgls"} or \code{"mvols"}, containing a linear
#'   model fit to shape data.
#' @param xvalue A value (numeric) or level (character) at which shape data is
#'   to be standardized (i.e., centered); If NULL, the mean of the complete
#'   sample is used (only available if there is a single explanatory variable).
#' @param method Method used for detrending; options are \code{"orthogonal"}
#'   (the default) and \code{"residuals"} (see details).
#' @param newdata New data to be standardized. It should be provided as either
#'   a linear model object fitting the same variables used in \code{model}
#'   measured in a new sample, or a 2-margin matrix of shape descriptors (only
#'   for \code{method = "orthogonal"}). See details.
#' @param tree,evmodel Further arguments needed only for orthogonal detrending
#'   with phylogenetic correction for a specific \code{xvalue}. Ignored
#'   otherwise.
#'
#' @details This function detrends (or standardizes, or corrects) shapes from
#'   variation associated with non-shape variables, using either the residuals
#'   computed from a linear model fitting the former to the latter
#'   (\code{method = "residuals"}), or the projection of shapes into a subspace
#'   that is orthogonal to variation associated to non-shape variables, computed
#'   using the Burnaby approach (\code{method = "orthogonal"}).
#'
#'   If \code{newdata} is provided as an object containing a linear model fit,
#'   shapes from \code{newdata} are detrended by "correcting" their relationship
#'   to the explanatory variables using the one estimated for the data provided
#'   in \code{model} ("partial detrending"). Specifically, if
#'   \code{method = "residuals"}, shape residuals will be computed by
#'   subtracting values fitted to the coefficients from \code{model} to the
#'   shapes provided in \code{newdata}; whereas if \code{method = "orthogonal"},
#'   shapes are projected into the subspace resulting from the subtraction of
#'   the orthogonal subspaces computed for \code{newdata} and \code{model}.
#'
#'   If \code{method = "orthogonal"} and \code{newdata} is provided as a
#'   2-margin matrix of shape descriptors, the function will instead return the
#'   new set of shapes, detrended using the relationship estimated for the data
#'   provided in \code{model} (i.e., the shapes provided in \code{newdata} are
#'   projected directly into the orthogonal subspace computed for the data
#'   provided in \code{model}).
#'
#'   The grand mean of the sample of shapes is used by default to center shape
#'   variation, although an \code{xvalue} specifying a level or numeric value of
#'   the explanatory variable in \code{model} to center shapes at can be
#'   provided. This shift can only be applied for one explanatory variable.
#'
#'   If a phylogenetic linear model is supplied to \code{model} (e.g., a linear
#'   model fitted using \code{\link[geomorph]{procD.pgls}} or
#'   \code{\link[mvMORPH]{mvgls}}, the procedure will use
#'   phylogenetically-corrected coefficients and phylogenetic mean (see Revell,
#'   2009) as calculated by the source function.
#'
#' @return A 2-margins matrix, of dimensions \code{n x (k x p)} for the case of
#'   landmark data and \code{n x (4 x nb.h)} for the case of Fourier data (where
#'   \code{nb.h} is the number of harmonics used in elliptic Fourier analysis).
#'
#' @seealso \code{\link[stats]{lm}}, \code{\link{burnaby}},
#'   \code{\link{expected_shapes}}, \code{\link[stats]{lm}},
#'   \code{\link[geomorph]{procD.lm}}, \code{\link[geomorph]{procD.pgls}},
#'   \code{\link[RRPP]{lm.rrpp}}, \code{\link[mvMORPH]{mvols}} and
#'   \code{\link[mvMORPH]{mvgls}}
#'
#' @export
#'
#' @references
#' Burnaby, T. P. (1966) \emph{Growth-invariant discriminant functions and
#'   generalized distances. Biometrics, 22, 96–110.}
#'
#' Revell, L. J. (2009). \emph{Size-correction and principal components for
#'   interspecific comparative studies}. Evolution, 63, 3258-3268.
#'
#' Klingenberg, C. P. (2016). \emph{Size, shape, and form: concepts of allometry
#'   in geometric morphometrics}. Development Genes and Evolution, 226(3),
#'   113-137.
#'
#' @examples
#'  #### Landmark data
#'
#'  #load tails data and packages
#'  library(geomorph)
#'  library(Morpho)
#'
#'  data(tails)
#'  shapes <- tails$shapes
#'  species <- tails$data$species
#'  sex <- tails$data$sex
#'  logsizes <- log(tails$sizes)
#'  msp <- mspace(shapes, links = tails$links, points = TRUE)
#'  hulls_by_group_2D(msp$ordination$x, fac = species)
#'
#'
#'  ### For numeric variables
#'
#'  #fit linear model between shapes and sizes, then center at the grand mean of
#'  #the sample
#'  model <- lm(two.d.array(shapes) ~ logsizes)
#'  detr_shapes_mat <- detrend_shapes(model)
#'
#'  detr_shapes_nosize <- arrayspecs(detr_shapes_mat, k = 2, p = 9)
#'
#'  msp_nosize <- mspace(detr_shapes_nosize, links = tails$links, points = TRUE)
#'  hulls_by_group_2D(msp_nosize$ordination$x, fac = species)
#'
#'  ## Using xvalue
#'
#'  #fit linear model between shapes and sizes, then center at the shape
#'  #corresponding to the maximum size of the sample
#'  model <- lm(two.d.array(shapes) ~ logsizes)
#'  detr_shapes_mat2 <- detrend_shapes(model,
#'                                     xvalue = max(logsizes))
#'
#'  detr_shapes_nosize2 <- arrayspecs(detr_shapes_mat2, k = 2, p = 9)
#'
#'  msp_nosize2 <- mspace(detr_shapes_nosize2, links = tails$links,
#'                        points = TRUE)
#'  hulls_by_group_2D(msp_nosize2$ordination$x, fac = species)
#'
#'  ## Using a phylogenetic tree
#'
#'  #fit linear model between shapes and sizes, then center at the shape
#'  #corresponding to the grand size of the sample
#'  sp_shapes <- expected_shapes(shapes, species)
#'  sp_logsizes <- c(tapply(logsizes, species, mean))
#'  model <- lm(two.d.array(sp_shapes) ~ sp_logsizes)
#'  detr_shapes_mat1 <- detrend_shapes(model)
#'
#'  detr_shapes_nosize1 <- arrayspecs(detr_shapes_mat1, k = 2, p = 9)
#'
#'  msp_nosize1 <- mspace(detr_shapes_nosize1, links = tails$links,
#'                        points = TRUE)
#'  points(msp_nosize1$ordination$x, pch = 21, bg = c(1:13))
#'
#'  ## Using newdata
#'
#'  #fit linear model between shapes and sizes for NDF species, then use the NDF
#'  #allometry to detrend DF shapes from allometric variation (just for
#'  #illustration, not implying this makes any sense)
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
#'  msp_nosize3 <- mspace(detr_shapes_nosize3, links = tails$links,
#'                        points = TRUE)
#'  hulls_by_group_2D(msp_nosize3$ordination$x, fac = factor(species[!index]))
#'
#'  ### For factors
#'
#'  #load wings data
#'
#'  data(wings)
#'  shapes <- wings$shapes
#'  cactus <- wings$data$cactus
#'  sex <- wings$data$sex
#'  species <- wings$data$species
#'  msp <- mspace(shapes, template = wings$template, points = TRUE)
#'  hulls_by_group_2D(msp$ordination$x, fac = cactus)
#'
#'
#'  #fit linear model between shapes and sex, then center at the grand mean
#'  #of the sample.
#'  model <- lm(two.d.array(shapes) ~ cactus)
#'  detr_shapes_mat <- detrend_shapes(model)
#'
#'  detr_shapes_nocac <- arrayspecs(detr_shapes_mat, k = 2, p = 9)
#'
#'  msp_nocac <- mspace(detr_shapes_nocac, template = wings$template,
#'                      points = TRUE)
#'  hulls_by_group_2D(msp_nocac$ordination$x, fac = cactus)
#'
#'  ## Using xvalue
#'
#'  #fit linear model between shapes and species, then center at the shape
#'  #corresponding to the mean shape of T. savana
#'  model <- lm(two.d.array(shapes) ~ cactus)
#'  detr_shapes_mat <- detrend_shapes(model, xvalue = "Tr")
#'
#'  detr_shapes_nocac2 <- arrayspecs(detr_shapes_mat, k = 2, p = 9)
#'
#'  msp_nocac2 <- mspace(detr_shapes_nocac2, template = wings$template,
#'                       points = TRUE)
#'  hulls_by_group_2D(msp_nocac2$ordination$x, fac = cactus)
#'
#'  ## Using newdata
#'
#'  #fit linear model between shapes and cactus for Db species, then use this
#'  #relationship to detrend Dk shapes from variation between cactus
#'  index <- species == "Db"
#'  shapes_Db <- shapes[,,index]
#'  cactus_Db <- cactus[index]
#'  shapes_Dk <- shapes[,,!index]
#'  cactus_Dk <- cactus[!index]
#'
#'  model_Db <- lm(two.d.array(shapes_Db) ~ cactus_Db)
#'  model_Dk <- lm(two.d.array(shapes_Dk) ~ cactus_Dk)
#'  detr_shapes_mat3 <- detrend_shapes(model_Db, newdata = model_Dk)
#'
#'  detr_shapes_nocactus <- arrayspecs(detr_shapes_mat3, k = 2, p = 9)
#'
#'  msp_cactus <- mspace(shapes_Dk, template = wings$template, points = TRUE)
#'  hulls_by_group_2D(msp_cactus$ordination$x, fac = factor(cactus[!index]))
#'
#'  msp_nocactus <- mspace(detr_shapes_nocactus, template = wings$template,
#'                      points = TRUE)
#'  hulls_by_group_2D(msp_nocactus$ordination$x, fac = factor(cactus[!index]))
#'
#'
#'  ### Comparing residuals vs orthogonal methods
#'
#'  \dontrun{
#'  #load shells3D data, retain only specimens belonging to S. vacaensis
#'  data("shells3D")
#'  index <- species == levels(species)[7]
#'  shapes <- shells3D$shapes
#'  refmesh <- shells3D$mesh_meanspec
#'  species <- shells3D$data$species
#'  sizes <- log(shells3D$sizes)
#'  template <- Morpho::tps3d(x = refmesh,
#'                            refmat = shapes[,,findMeanSpec(shapes)],
#'                            tarmat = expected_shapes(shapes[,,index]))
#'
#'  #compute allometric axis (i.e., shape variation associated to changes in
#'  #size)
#'  alloax <- lm(two.d.array(shapes[,,index]) ~ sizes[index])
#'
#'  #project vacaensis specimens into the overall morphospace together with
#'  #allometric axis
#'  mspace(shapes, template = template, bg.models = "gray",
#'         cex.ldm = 0, alpha.models = 0.7, adj_frame = c(0.93,0.93)) %>%
#'    proj_shapes(shapes[,,index], pch = 21, bg = 7) %>%
#'    proj_axis(alloax, type = 2, lwd = 2, lty = 1)
#'
#'  #compute non-allometric variation using linear models (center on the shape
#'  #corresponding to the maximum size)
#'  detr_shapes_lm <- lm(two.d.array(shapes[,,index]) ~ sizes[index]) %>%
#'    detrend_shapes(method = "residuals", xvalue = max(sizes[index])) %>%
#'    arrayspecs(k = 3, p = 90)
#'
#'  #visualize
#'  mspace(shapes, template = template, bg.models = "gray",
#'         cex.ldm = 0, alpha.models = 0.7, adj_frame = c(0.93,0.93)) %>%
#'    proj_shapes(shapes[,,index], pch = 1, col = 7) %>%
#'    proj_axis(alloax, type = 2, lwd = 2, lty = 1) %>%
#'    proj_shapes(detr_shapes_lm, pch = 21, bg = 7)
#'
#'  #compute non-allometric variation using a orthogonal subspace (center on the
#'  #shape corresponding to the maximum size)
#'  detr_shapes_os <- lm(two.d.array(shapes[,,index]) ~ sizes[index]) %>%
#'    detrend_shapes(method = "orthogonal", xvalue = max(sizes[index])) %>%
#'    arrayspecs(k = 3, p = 90)
#'
#'  #visualize
#'  mspace(shapes, template = template, bg.models = "gray",
#'         cex.ldm = 0, alpha.models = 0.7, adj_frame = c(0.93,0.93)) %>%
#'    proj_shapes(shapes[,,index], pch = 1, col = 7) %>%
#'    proj_axis(alloax, type = 2, lwd = 2, lty = 1) %>%
#'    proj_shapes(detr_shapes_os, pch = 21, bg = 7)
#'  }
#'
#'
#'  #### Fourier data (quick demo)
#'
#'  #load shells data
#'  data(shells)
#'  shapes <- shells$shapes$coe
#'  species <- shells$data$species
#'  logsizes <- log(shells$sizes)
#'  msp <- mspace(shapes, mag = 0.5, points = TRUE)
#'  hulls_by_group_2D(msp$ordination$x, fac = species)
#'
#'  #fit linear model between shapes and sizes, then center at the shape
#'  #corresponding to the grand mean
#'  model <- lm(shapes ~ logsizes)
#'  detr_shapes_nosize <- detrend_shapes(model)
#'
#'  msp_nosize <- mspace(detr_shapes_nosize, mag = 0.5, points = TRUE)
#'  hulls_by_group_2D(msp_nosize$ordination$x, fac = species)
#'
#'  #fit linear model between shapes and sizes, then center at the shape
#'  #corresponding to the maximum size
#'  model <- lm(shapes ~ logsizes)
#'  detr_shapes_nosize2 <- detrend_shapes(model,
#'                                        xvalue = max(logsizes))
#'
#'  msp_nosize2 <- mspace(detr_shapes_nosize2, mag = 0.5, points = TRUE)
#'  hulls_by_group_2D(msp_nosize2$ordination$x, fac = species)
#'
#'  #fit linear model between shapes and sizes, then center at the maximum size
#'  shapes_koeneni <- shapes[species == "koeneni",]
#'  logsizes_koeneni <- logsizes[species == "koeneni"]
#'  shapes_esbelta <- shapes[species == "esbelta",]
#'  logsizes_esbelta <- logsizes[species == "esbelta"]
#'
#'  model_koeneni <- lm(shapes_koeneni ~ logsizes_koeneni)
#'  model_esbelta <- lm(shapes_esbelta ~ logsizes_esbelta)
#'  detr_shapes_nosize3 <- detrend_shapes(model_koeneni,
#'                                        newdata = model_esbelta)
#'
#'  msp_nosize3 <- mspace(shapes_esbelta, mag = 0.5, points = TRUE)
#'  title("raw P. esbelta morphospace")
#'  msp_nosize4 <- mspace(detr_shapes_nosize3, mag = 0.5, points = TRUE)
#'  title("P. esbelta morphospace, refined using \n allometric variation from P. koeneni")
detrend_shapes <- function(model, method = "residuals", xvalue = NULL, newdata = NULL,
                           tree = NULL, evmodel = NULL) {

  # # # load_all()
  # # #
  # model = model2
  # model$call
  #model = mvMORPH::mvols(sp_shapes ~ sp_logsizes)

  # model <- mvMORPH::mvgls(sp_shapes ~ sp_type, tree = tree, model = "BM")
  #

  # model <- geomorph::procD.pgls(sp_shapes ~ sp_type, phy = phy, data = gdf)
  # #model <- geomorph::procD.lm(sp_shapes ~ sp_type, data = gdf)
  # model$call
  # method = "orthogonal"
  # newdata = NULL
  # tree = NULL
  # xvalue = NULL
  # evmodel = NULL
  # # # # #############################


  #data preparation
  if(!any(method == "residuals", method == "orthogonal")) stop("method should be one of 'residuals' or 'orthogonal'")

  #!
  # x <- model$model[-1]
  # y <- model$model[[1]]
  designmat <- stats::model.matrix(model)

  admod <- adapt_model(model)
  coefs <- admod$coefs
  grandmean <- admod$grandmean
  fitted <- admod$fitted
  y <- admod$y
  x <- admod$x
  modtype <- admod$modtype


  if(ncol(x) > 1 & !is.null(xvalue)) stop("xvalue can only be specified for a single explanatory variable")

  if(is.null(rownames(x))) rownames(x) <- seq_len(nrow(x))
  if(is.null(rownames(y))) rownames(y) <- seq_len(nrow(y))
  namesx <- rownames(x)
  namesy <- rownames(y)

  #!
  # designmat <- stats::model.matrix(stats::as.formula(paste0("~ ",
  #                                                           strsplit(as.character(
  #                                                             model$call[2]), split = "~")[[1]][2])), data = x)

  # if(!is.null(tree)) {
  #   if(!all(length(tree$tip.label) == nrow(x), length(tree$tip.label) == nrow(y))) {
  #     stop("\nNumber of tips in the tree does not match the number of observations in x and/or y data sets")
  #   }
  #   if(!all(tree$tip.label %in% namesy, tree$tip.label %in% namesx)) {
  #     stop("\nNames in phylogenetic tree does not match names in x and/or y data sets")
  #   } else {
  #     x <- cbind(x[tree$tip.label,])
  #     y <- cbind(y[tree$tip.label,])
  #     designmat <- cbind(designmat[tree$tip.label,])
  #   }
  # }

  #detrending
  if(method == "orthogonal") {

    # if(modtype == "pgls" & !is.null(tree))
    #   stop("method = 'orthogonal', together with a phylogenetic linear model, requires the phylogenetic tree to be provided")

    axmat <- if(nrow(coefs) < 3) cbind(coefs[-1,]) else cbind(t(coefs[-1,]))
    center <- colMeans(fitted)

    #ortho_space <- suppressWarnings(burnaby(x = y, axmat = axmat, tree = tree)) #!
    ortho_space <- burnaby(x = y, axmat = axmat, center = center)
    ortho_scores <- ortho_space$x
    ortho_rotation <- ortho_space$rotation
    ax <- sum(ortho_space$sdev^2 > 1e-15)

    if(!is.null(xvalue)) {
      if(is.factor(x[[1]]) & !xvalue %in% levels(x[[1]]) & inherits(model, "mvgls")) {
        xvalue <- "Intercept"
        warning("Be sure that the level indicated in 'xvalue' has been spelled correctly; it is assumed it was, and that it corresponds to the Intercept of the mvgls model")
      }

      # ortho_center <- c(expected_shapes(shapes = y, x = x[[1]], xvalue = xvalue,
      #                                   tree = tree, returnarray = FALSE)) #!

      if(modtype == "pgls" & any(is.null(tree), is.null(evmodel)))
        stop("Please feed the tree and evmodel arguments")

      ortho_center <- c(expected_shapes(shapes = y, x = x[[1]], xvalue = xvalue,
                                        tree = tree, evmodel = evmodel, returnarray = FALSE))
    } else {
      ortho_center <- ortho_space$center
    }

    if(!is.null(newdata)) {

      #if(inherits(newdata, "mlm")) { #!
      if(any(class(newdata)[1] %in% c("mlm", "procD.lm", "lm.rrpp", "mvgls", "mvols"))) {

        # newx <- newdata$model[-1]
        # newy <- newdata$model[[1]] #!
        newadmod <- adapt_model(newdata)
        #newx <- newadmod$x #!
        newy <- newadmod$y

        if(is.null(rownames(newy))) rownames(newy) <- seq_len(nrow(newy))
        namesy <- rownames(newy)

        ortho_scores <- proj_eigen(x = newy, vectors = ortho_space$rotation,
                                   center = ortho_space$center)

        # #newortho_space <- burnaby(x = newy, vars = newx) #!
        # # newortho_space <- burnaby(x = newy, axmat = newadmod$coefs[2,],
        # #                           center = newadmod$grandmean) #!
        # # newortho_space <- burnaby(x = newy, axmat = admod$coefs[2,],
        # #                           center = admod$grandmean)
        # newortho_space <- burnaby(x = newy, axmat = axmat, center = center)
        # ax <- min(sum(newortho_space$sdev^2 > 1e-15), ax)
        # ortho_scores <- newortho_space$x
        # ortho_rotation <- newortho_space$rotation#[,seq_len(ax)]# - ortho_space$rotation[,seq_len(ax)] #!
        #
        #
        # if(!is.null(xvalue)) {
        #
        #   if(is.factor(x[[1]]) & !xvalue %in% levels(x[[1]]) & inherits(model, "mvgls")) {
        #     xvalue <- "Intercept"
        #     warning("Be sure that the level indicated in 'xvalue' has been spelled correctly; it is assumed it was, and that it corresponds to the Intercept of the mvgls model")
        #   }
        #   # ortho_center <- c(expected_shapes(shapes = newy, x = newx[[1]],
        #   #                                   xvalue = xvalue, returnarray = FALSE)) #!
        #   # ortho_center <- c(expected_shapes(shapes = newy, x = newx[[1]], evmodel = evmodel,
        #   #                                   xvalue = xvalue, returnarray = FALSE))
        #   # ortho_center <- c(expected_shapes(shapes = y, x = x[[1]], xvalue = xvalue,
        #   #                                   tree = tree, evmodel = evmodel, returnarray = FALSE))
        #   ortho_center <- ortho_center
        # } else {
        #   ortho_center <- ortho_center#newortho_space$center
        # }

      } else {
        ortho_scores <- proj_eigen(x = newdata, vectors = ortho_space$rotation,
                                   center = ortho_space$center)
      }
    }

    # ortho_shapes2d <- rev_eigen(scores = ortho_scores[, seq_len(ax)],
    #                             vectors = ortho_rotation[, seq_len(ax)],
    #                             center = ortho_center) #!
    ortho_shapes2d00 <- rev_eigen(scores = ortho_scores[, seq_len(ax)],
                                vectors = ortho_rotation[, seq_len(ax)],
                                center = ortho_center)
    ortho_shapes2d0 <- scale(ortho_shapes2d00, scale = FALSE, center = TRUE)
    ortho_shapes2d <- matrix(t(t(ortho_shapes2d0) + ortho_center),
                             nrow = nrow(ortho_shapes2d00),
                             ncol = ncol(ortho_shapes2d00), byrow = FALSE)
    colnames(ortho_shapes2d) <- colnames(ortho_shapes2d00)
    rownames(ortho_shapes2d) <- rownames(ortho_scores)
    return(ortho_shapes2d)
  }

  if(method == "residuals") {
    which_na <- is.na(apply(coefs,1,sum))
    # resids <- y - designmat[,!which_na] %*% coefs[!which_na,] #!
    resids <- y - fitted

    if(!is.null(newdata)) {

      if(any(class(newdata)[1] %in% c("mlm", "procD.lm", "lm.rrpp", "mvgls", "mvols"))) {
        # newx <- newdata$model[-1]
        # newy <- newdata$model[[1]] #!
        newadmod <- adapt_model(newdata)
        newx <- newadmod$x
        newy <- newadmod$y
        #newfitted <- newadmod$fitted

        if(is.null(rownames(newy))) rownames(newy) <- seq_len(nrow(newy))
        namesy <- rownames(newy)

        # formula <- stats::as.formula(paste0("~ ", strsplit(as.character(newdata$call[2]),
        #                                                    split = "~")[[1]][2])) #!
        newdesignmat <- stats::model.matrix(newdata)

        which_na <- is.na(apply(coefs,1,sum))
        resids <- newy - newdesignmat[,!which_na] %*% coefs[!which_na,]
        #resids <- y - newfitted
      } else stop("Only linear models can be supplied as newdata for method = 'residuals'")
    }

    n <- nrow(resids)

    if(is.null(xvalue)) {

      fitted_vec <- rep(1, n) %*% t(grandmean)
      predicted_mat00 <- resids + fitted_vec
      fitted_vec <- colMeans(fitted_vec)
      # fitted_vec <- c(colMeans(admod$fitted))
      # predicted_mat <- t(t(resids) + fitted_vec)

    } else {

      if(is.factor(x[[1]])) {
        if(!xvalue %in% levels(x[[1]]) & inherits(model, "mvgls")) {
          xvalue <- "Intercept"
          warning("Be sure that the level indicated in 'xvalue' has been spelled correctly; it is assumed it was, and that it corresponds to the Intercept of the mvgls model")

        }
      }

      # if(is.numeric(x[[1]]) == TRUE) {
      #   designmat <- cbind(1, xvalue)
      # }
      #
      # if(is.factor(x[[1]]) == TRUE) {
      #   designmat <- rep(0, nlevels(x[[1]]))
      #   designmat[which(levels(x[[1]]) == xvalue)] <- 1
      #
      #   if(isFALSE(designmat[1] == 1)) designmat[1] <- 1
      # }
      #
      # fitted <- as.numeric(designmat %*% coefs)
      #
      # fitted_vec <- rep(1, n) %*% t(fitted)
      # predicted_mat <- resids + fitted_vec

      if(modtype == "pgls" & any(is.null(tree), is.null(evmodel)))
        stop("Please feed  the tree and evmodel arguments")

      fitted_vec <- c(expected_shapes(shapes = y, x = x[[1]], xvalue = xvalue,
                                    tree = tree, evmodel = evmodel, returnarray = FALSE))
      # fitted <- if(is.null(newdata)) fitted[which(x == xvalue),] else newfitted[which(newx == xvalue),]
      # fitted_vec <- if(is.null(dim(fitted))) fitted else colMeans(fitted)
      predicted_mat00 <- t(t(resids) + fitted_vec)
      # fitted_vec <- rep(1, n) %*% t(fitted_vec)
      # predicted_mat <- resids + fitted_vec
    }

    predicted_mat0 <- scale(predicted_mat00, scale = FALSE, center = TRUE)
    predicted_mat <- matrix(t(t(predicted_mat0) + fitted_vec),
                            nrow = nrow(predicted_mat00),
                            ncol = ncol(predicted_mat00), byrow = FALSE)
    colnames(predicted_mat) <- colnames(predicted_mat00)
    rownames(predicted_mat) <- rownames(predicted_mat00)

    predicted_mat <- predicted_mat[namesy,]
    return(rbind(predicted_mat))

  }
}


################################################################################

#' Correct artificially rotated set of Fourier shapes interactively
#'
#' @description Choose and correct 180-degrees spurious rotation in a sample of
#'   closed outline shapes.
#'
#' @param ef An \code{"OutCoe"} object.
#' @param index An optional integer vector providing the indices identifying the
#'   shapes to be rotated. If declared, overrides interactive selection of
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
#' @references
#' Iwata, H., & Ukai, Y. (2002). \emph{SHAPE: a computer program package for
#'   quantitative evaluation of biological shapes based on elliptic Fourier
#'   descriptors}. Journal of Heredity, 93(5), 384-385.
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


################################################################################

#' Compute shapes at the extremes of a morphometric axis
#'
#' @description This function computes the theoretical shapes corresponding to
#'   the extremes of a morphometric axis, which can be supplied as an object
#'   containing a linear model or a multivariate ordination.
#'
#' @param obj An object of class \code{"mlm"}, \code{"procD.lm"},
#'   \code{"lm.rrpp"}, \code{"mvgls"} or \code{"mvols"}, containing a linear
#'   model fitting a single explanatory variable to shape data, or an object of
#'   class \code{"prcomp"}, \code{"bg_prcomp"}, \code{"pls_shapes"},
#'   \code{"phy_pls_shapes"}, \code{"burnaby"}, \code{"phy_burnaby"},
#'   \code{"gm.prcomp"}, \code{"bgPCA"}, \code{"pls2B"}, \code{"phyl.pca"} or
#'   \code{"mvgls.pca"} with the results of multivariate ordination of shape
#'   data.
#' @param axis An optional integer value specifying the axis of the multivariate
#'   ordination which is to be represented.
#' @param mag Numeric; magnifying factor for representing shape transformation.
#'
#' @return A 2-margins matrix, of dimensions \code{2 x (k x p)} for the case of
#'   landmark data and \code{2 x (4 x nb.h)} for the case of Fourier data (where
#'   \code{nb.h} is the number of harmonics used in elliptic Fourier analysis).
#'
#' @details If an object of class \code{"mlm"} fitting shape to a factor (only
#'   two levels allowed) is supplied in \code{obj}, magnification of the axis
#'   range is attained through bgPCA.
#'
#' @export
#'
#' @seealso \code{\link{expected_shapes}}, \code{\link{rev_eigen}},
#'   \code{\link[stats]{lm}}, \code{\link[geomorph]{procD.lm}},
#'   \code{\link[geomorph]{procD.pgls}}, \code{\link[RRPP]{lm.rrpp}},
#'   \code{\link[mvMORPH]{mvols}}, \code{\link[mvMORPH]{mvgls}},
#'   \code{\link[geomorph]{gm.prcomp}}, \code{\link[Morpho]{groupPCA}},
#'   \code{\link[Morpho]{pls2B}}, \code{\link[phytools]{phyl.pca}},
#'   \code{\link[mvMORPH]{mvgls.pca}}
#'
#' @references
#' MacLeod, N. (2009). \emph{Form & shape models}. Palaeontological Association
#'   Newsletter, 72(620), 14-27.
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
#' #perform lm of shape on size, compute and plot extreme shapes at its natural
#' #range
#' model <- lm(two.d.array(tails$shapes) ~ log(tails$sizes))
#' extshapes2d <- ax_transformation(obj = model, mag = 1)
#' extshapes <- arrayspecs(extshapes2d, k = ncol(shapes), p = nrow(shapes))
#' pile_shapes(extshapes, links = links, mshape = FALSE)
#'
#' #perform lm of shape on size, compute and plot extreme shapes at its natural
#' #range
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

  # if(any(class(obj) %in% c("prcomp", "bg_prcomp", "phy_prcomp", "phyalign_comp",
  #                          "pls_shapes", "phy_pls_shapes", "burnaby", "phy_burnaby"))) { #!
  if(any(class(obj) == c("pls2b", "phy_pls2b"))) stop("pls2b objects are not allowed; use pls_shapes instead")
  if(any(class(obj)[1] %in% c("prcomp", "bg_prcomp", "phy_prcomp", "phyalign_comp",
                              "pls_shapes", "phy_pls_shapes", "burnaby", "phy_burnaby",
                              "gm.prcomp", "pls2B", "bgPCA", "phyl.pca", "mvgls.pca", "PCA", "pls"))) {

    obj <- adapt_ordination(obj)
    extshapes_mat <- rev_eigen(range(obj$x[,axis]) * mag,
                               obj$rotation[,axis],
                               obj$center)
  }

  #if(any(class(obj) %in% "mlm")) {  #!
  if(any(class(obj)[1] %in% c("mlm", "procD.lm", "lm.rrpp", "mvgls", "mvols"))) {

    # x <- obj$model[,ncol(obj$model)]  #!
    obj <- adapt_model(obj)
    x <- obj$x[,1]

    if(is.numeric(x)) {
      cent <- mean(range(x))
      halfdif <- diff(range(x)) / 2
      newrange <- c(cent - (halfdif * mag),
                    cent + (halfdif * mag))

      # coefs <- obj$coefficients  #!
      coefs <- obj$coefs
      designmat <- cbind(1, newrange)
      extshapes_mat <- rbind(designmat %*% coefs)
    }

    if(is.factor(x)) {

      if(nlevels(x) > 2) stop("Only two levels are allowed for extracting axes from mlm objects; try providing a bg_prcomp object")

      # Y <- obj$model[,1]  #!
      Y <- obj$y
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


################################################################################

#' Retrieve / compute shapes from an existing morphospace
#'
#' @description Extracts shapes from \code{"mspace"} objects in various ways
#'   (background shape models, specific axes, or particular coordinates - either
#'   provided or chosen interactively by the user).
#'
#' @param mspace An \code{"mspace"} object.
#' @param axis Optional integer indicating an axis along which shapes should
#'   be sampled.
#' @param range Optional integer vector of length 2, indicating the range of
#'   values the axis should be sampled over.
#' @param nshapes Optional integer indicating the number of shapes the user
#'   wishes to extract.
#' @param scores An optional vector of length 2 or 2-column matrix indicating
#'   the (x,y) coordinates in the morphospace that the user wishes to extract
#'   as shapes. If \code{NULL}, a new device will open and the user will be
#'   asked choose the coordinates interactively. Ignored when \code{axis} is
#'   provided.
#' @param keep.template Logical; should warped templates be returned as well?
#' @param mag Optional numeric indicating a magnifying factor for representing
#'   shape transformation.
#'
#' @details This function provides the user with an easy way to extract
#'   theoretical shapes from an existing morphospace. If only an \code{"mspace"}
#'   object is provided, the set of background shape models (optionally
#'   amplified by a factor of \code{mag}) will be returned.
#'
#'   If \code{axis} is provided, a sample of \code{nshapes} shapes computed at
#'   regular intervals along the specified ordination axis (either over its
#'   empirical range -optionally amplified by a factor of \code{mag}- or, if
#'   provided, between the extremes of \code{range}) will be returned.
#'
#'   If \code{axis = NULL}, this function will let the user to select arbitrary
#'   coordinates in the morphospace to be back-transformed into shapes. There
#'   are two alternatives: 1) if \code{scores = NULL} (the default option) the
#'   user will be asked to interactively select the location(s) of
#'   \code{nshapes} points in a new graphical device. 2) Otherwise, the set of
#'   shapes represented by the (x,y) coordinates provided in \code{scores} will
#'   be returned.
#'
#' @return A list containing sampled shapes (\code{$shapes}), as well as their
#'   associated templates (\code{$templates}), when warranted.
#'
#' @export
#'
#' @examples
#' #load all the relevant data and packages
#' library(Morpho)
#' library(geomorph)
#'
#' data("tails")
#' shapes <- tails$shapes
#' sizes <- tails$sizes
#' species <- tails$data$species
#' type <- tails$data$type
#' links <- tails$links
#' sp_shapes <- expected_shapes(shapes, species)
#' tree <- tails$tree
#'
#' #build phylomorphospace
#' phylomsp <- mspace(shapes, links = links) %>%
#'   proj_phylogeny(sp_shapes, tree = tree)
#'
#'
#' ##Extracting background shape models
#'
#' #extract background shape models
#' background_shapes <- extract_shapes(phylomsp)
#'
#' #pile shapes and visualise the corresponding coordinates sampled in the
#' #morphospace
#' pile_shapes(background_shapes$shapes, links = links)
#'
#' plot_mspace(phylomsp, phylo = FALSE)
#' background_scores <- proj_eigen(two.d.array(background_shapes$shapes),
#'                                 phylomsp$ord$rotation, phylomsp$ord$center)
#' points(background_scores, pch = 21, bg = "red")
#'
#'
#' ##Sampling a particular ordination axis
#'
#' #extract shapes along PC2
#' PC2_shapes <- extract_shapes(phylomsp, axis = 2, nshapes = 8)
#'
#' #pile shapes and visualise the corresponding coordinates sampled in the
#' #morphospace
#' pile_shapes(PC2_shapes$shapes, links = links, mshape = FALSE)
#'
#' plot_mspace(phylomsp, phylo = FALSE)
#' PC2_scores <- proj_eigen(two.d.array(PC2_shapes$shapes),
#'                          phylomsp$ord$rotation, phylomsp$ord$center)
#' points(PC2_scores, pch = 21, bg = "blue")
#'
#'
#' ##Sampling particular (x,y) locations
#'
#' #1. Interactively
#'
#' \dontrun{
#'
#' #select 1 shape in the new window
#' arbitrary_shape <- extract_shapes(phylomsp, nshapes = 1)
#'
#' #plot shape and visualise the corresponding coordinates sampled in the
#' #morphospace
#' plot(arbitrary_shape$shapes[,,1], axes = FALSE, xlab = "", ylab = "")
#' lineplot(arbitrary_shape$shapes[,,1], links)
#'
#' plot_mspace(phylomsp, phylo = FALSE)
#' arbitrary_scores <- proj_eigen(two.d.array(arbitrary_shape$shapes),
#'                                 phylomsp$ord$rotation, phylomsp$ord$center)
#' points(arbitrary_scores, pch = 21, bg = "magenta")
#'
#' }
#'
#'
#' #2. Specifying coordinates
#'
#' #get scores of the nodes of the phylogeny for the first two PCs
#' nodes_scores0 <- phylomsp$projected$phylo_scores[14:25,1:2]
#'
#' #extract shapes from morphospace
#' nodes_shapes <- extract_shapes(phylomsp, scores = nodes_scores0)
#'
#' plot_mspace(phylomsp, phylo = TRUE)
#' nodes_scores <- proj_eigen(two.d.array(nodes_shapes$shapes),
#'                          phylomsp$ord$rotation, phylomsp$ord$center)
#' points(nodes_scores, pch = 21, bg = "green")
extract_shapes <- function(mspace,
                           axis = NULL,
                           nshapes = NULL,
                           scores = NULL,
                           range = NULL,
                           keep.template = TRUE,
                           mag = 1) {


  if(all(is.null(nshapes), is.null(axis), is.null(scores))) {
    axis <- mspace$plotinfo$axes
    scores <- proj_eigen(x = shapes_mat(mspace$projected$shapemodels)$data2d,
                         vectors = mspace$ordination$rotation[,axis],
                         center = mspace$ordination$center) * mag
  } else {
    if(!is.null(axis)) {
      if(!is.null(nshapes)) {
        if(!is.null(scores)) warning("an axis has been provided; scores will be ignored")
        if(nshapes > 1) {
          if(!is.null(range)) {
            xrange <- range
          }  else {
            xrange <- range(mspace$ordination$x[,axis]) * mag
          }
          scores <- seq(xrange[1], xrange[2], length.out = nshapes)
        } else stop("At least two shapes are necessary to represent variation along a PC axis")
      } else stop("nshapes must be provided")
    } else {
      axis <- mspace$plotinfo$axes
      if(!is.null(scores)) {
        if(!is.null(nshapes)) {
          if((is.null(dim(scores)) & nshapes != 1) | (nshapes != nrow(scores))) {
            warning("number of scores provided and nshapes differ, the latter will be ignored")
          }
        }
        if(is.null(dim(scores))) scores <- rbind(scores)
      } else {
        if(Sys.info()["sysname"] == "Darwin") grDevices::quartz() else grDevices::X11()
        cat(paste0("Select ", nshapes, " points in the morphospace"))
        plot_mspace(mspace)
        scores <- matrix(unlist(graphics::locator(n = nshapes)),
                         nrow = nshapes, ncol = 2, byrow = TRUE)
        grDevices::dev.off()
      }
    }
  }

  data2d <- rev_eigen(scores = scores,
                      vectors = mspace$ordination$rotation[,axis],
                      center = mspace$ordination$center)

  if(mspace$ordination$datype == "landm") {
    shapes <- geomorph::arrayspecs(data2d,
                                   p = mspace$plotinfo$p,
                                   k = mspace$plotinfo$k)

    if(keep.template) {
      if(is.null(mspace$plotinfo$template)) {
        warning("there are no templates to warp; won't be returned")
      } else {
        centroid <- matrix(rev_eigen(0, mspace$ordination$rotation[,1], mspace$ordination$center),
                           ncol = mspace$plotinfo$k, byrow = TRUE)
        if(mspace$plotinfo$k == 2) {

          temp_cent <- rbind(centroid,
                             Momocs::tps2d(mspace$plotinfo$template[-(seq_len(mspace$plotinfo$p)), ],
                                           mspace$plotinfo$template[(seq_len(mspace$plotinfo$p)), ],
                                           centroid))

          temp_warpd_list <- lapply(seq_len(dim(shapes)[3]),
                                    function(i) {Momocs::tps2d(temp_cent[-(seq_len(mspace$plotinfo$p)),],
                                                               temp_cent[(seq_len(mspace$plotinfo$p)),],
                                                               shapes[,,i])})
          templates <- temp_warpd_arr <- abind::abind(temp_warpd_list, along = 3)

        } else {
          templates <- vector(length = dim(shapes)[3], mode = "list")
          for(i in seq_len(dim(shapes)[3])) {
            templates[[i]] <- Morpho::tps3d(x = mspace$plotinfo$template ,
                                            refmat = centroid, tarmat = shapes[,,i])
          }
        }
      }
    }
  } else {
    shapes <- data2d
  }

  results <- list(shapes = shapes)
  if(keep.template & !is.null(mspace$plotinfo$template)) results$templates <- templates

  return(invisible((results)))
}

