
########################################################################################

#' Between-groups Principal Component Analysis
#'
#' @description Performs between group PCA allowing for leave-one-out cross-validation, which is useful one the number of variables
#'   exceeds the number of observations (i.e., alleviates spurious separation
#'   between groups).
#'
#' @param x A matrix with variables as columns and observations as rows.
#' @param groups Factor; classification of observations of \code{x} into a
#'   priori groups.
#' @param gweights Logical; whether to weight each group by its number of
#'   observations.
#' @param LOOCV Logical; whether to apply leave-one-out cross-validation.
#' @param recompute Logical; whether to re-compute rotation matrix using the
#'   scores resulting from LOOCV.
#' @param corr Logical; whether to use correlation instead of covariance matrix
#'   as input.
#'
#' @return
#'
#' @seealso \code{\link[base]{prcomp}}
#'
#' @references
#'
#' @export
#'
#' @examples
bg_prcomp <- function(x, groups, gweights = TRUE,
                      LOOCV = FALSE, recompute = FALSE,
                      corr = FALSE) {

  grandmean <- colMeans(x)
  totvar <- sum(apply(x, 2, stats::var))

  if(LOOCV == FALSE) {

    x_centered <- scale(x, scale = FALSE, center = TRUE)
    x_gmeans <- apply(X = x_centered, MARGIN = 2, FUN = tapply, groups, mean)

    n_g <- tapply(groups, groups, length)
    if(gweights == TRUE) {
      wts <- as.numeric(n_g / sum(n_g))
    } else {
      wts <- rep(1, nlevels(groups))
    }

    vcv_g <- stats::cov.wt(x_gmeans, wt = wts, cor = corr)$cov
    svd <- svd(x = vcv_g)

    scores <- proj_eigen(x, svd$v, center = colMeans(x))
    scores <- cbind(scores[, 1:(nlevels(groups) - 1)])

    results <- list(sdev = sqrt(svd$d), rotation = svd$v, x = scores,
                    center = grandmean, grcenters = x_gmeans, totvar = totvar)
    class(results) <- "bg_prcomp"
    return(results)

  } else {

    n_g <- tapply(groups, groups, length)
    if(gweights == TRUE) {
      wts <- as.numeric(n_g / sum(n_g))
    } else {
      wts <- rep(1, nlevels(groups))
    }

    x_centered <- scale(x, scale = FALSE, center = TRUE)
    x_gmeans <- apply(X = x_centered, MARGIN = 2, FUN = tapply, groups, mean)

    vcv_g <- cov.wt(x_gmeans, wt = wts, cor = corr)$cov
    svd <- svd(vcv_g)

    rotation <- svd$v
    values <- svd$d
    refaxis <- svd$v[, 1]

    scores_l <- lapply(1:nrow(x), function(i) {
      subx <- x[-i,]
      subgroups <- groups[-i]

      subx_centered <- scale(subx, scale = FALSE, center = TRUE)
      subx_gmeans <- apply(X = subx_centered, MARGIN = 2, FUN = tapply, subgroups, mean)

      n_g <- tapply(subgroups, subgroups, length)
      if(gweights == TRUE) {
        wts <- as.numeric(n_g / sum(n_g))
      } else {
        wts <- rep(1, nlevels(subgroups))
      }

      vcv_g <- stats::cov.wt(subx_gmeans, wt = wts, cor = corr)$cov
      svd <- svd(vcv_g)

      scores <- proj_eigen(x[i,], svd$v, center = colMeans(x))
      scores <- scores * sign(stats::cor(svd$v[, 1], refaxis))

    })
    scores <- do.call("rbind", scores_l)

    scores <- cbind(scale(scores[, 1:(nlevels(groups) - 1)], scale = FALSE, center = TRUE))
    if(recompute == TRUE) rotation <- matrix(t(t(scores) %*% t(t(x))), ncol = (nlevels(groups) - 1))

    results <- list(sdev = sqrt(values), rotation = rotation, x = scores,
                    center = grandmean, grcenters = x_gmeans, totvar = totvar)
    class(results) <- "bg_prcomp"
    return(results)
  }

}


########################################################################################

#' Phylogenetic Principal Component Analysis
#'
#' @description A wrapper for [pyl.pca()] from \code{phytools}.
#'
#' @param x A matrix with one or more variables as columns and observations as
#'   rows. Row must be named and match tip labels from the phylogenetic
#'   \code{tree}.
#' @param tree A \code{"phy"} object containing a phylogenetic tree. Tip labels
#'   should match the row names from \code{x}.
#' @param corr Logical; whether to use correlation instead of covariance matrix
#'   as input.
#' @param ...
#'
#' @return
#'
#' @seealso \code{\link[phytools]{phyl.pca}}, \code{\link[base]{prcomp}}
#'
#' @export
#'
#' @examples
phy_prcomp <- function(x, tree, corr = FALSE, ...) {

  if(corr == FALSE) {
    mode <- "cov"
  } else {
    mode <- "corr"
  }
  center <- colMeans(x)

  phypca <- phytools::phyl.pca(Y = x, tree = tree, mode = mode, ...)

  results <- list(sdev = sqrt(phypca$Eval), rotation = phypca$Evec,  x = phypca$S,
                  center = center)
  if(!is.null(phypca$lambda)) results$lambda <- phypca$lambda
  if(!is.null(phypca$logL)) results$logL <- phypca$logL

  class(results) <- "phy_prcomp"
  return(results)

}


########################################################################################

#' Two-block Partial Least Squares
#'
#' @description Performs Partial Least Squares allowing for leave-one-out
#'   cross-validation, which is useful one the number of variables exceeds the
#'   number of observations (i.e., alleviates spurious covariation between
#'   variables).
#'
#' @param y A matrix with one or more variables as columns and observations as
#'   rows, representing the first block.
#' @param x A matrix with one or more variables as columns and observations as
#'   rows, representing the second block.
#' @param LOOCV Logical; whether to apply leave-one-out cross-validation.
#' @param recompute Logical; whether to re-compute rotation matrix using the
#'   scores resulting from LOOCV.
#'
#' @return
#'
#' @seealso \code{\link{pls_shapes}}
#'
#' @export
#'
#' @examples
pls2b <- function(y, x, LOOCV = FALSE, recompute = FALSE) {

  if(is.vector(x)) x <- matrix(x)
  if(is.vector(y)) y <- matrix(y)

  totvar_x <- sum(apply(x, 2, stats::var))
  totvar_y <- sum(apply(y, 2, stats::var))

  y_center <- colMeans(y)
  x_center <- colMeans(x)

  y_centered <- scale(y, scale = FALSE, center = TRUE)
  x_centered <- scale(x, scale = FALSE, center = TRUE)

  svd <- Morpho:::svd2B(x_centered, y_centered)
  values <- svd$d
  rotation_y <- svd$v
  rotation_x <- svd$u

  if(LOOCV == FALSE) {

    yscores <- y_centered %*% svd$v
    xscores <- x_centered %*% svd$u

    results <- list(yscores = yscores, yrotation = svd$v, ycenter = y_center, ytotvar = totvar_y,
                    xscores = xscores, xrotation = svd$u, xcenter = x_center, xtotvar = totvar_x,
                    values = svd$d)
    class(results) <- "pls"
    return(results)

  } else {

    ndim <- min(ncol(y), ncol(x))
    yscores <- matrix(NA, nrow = nrow(y), ncol = ndim)
    xscores <- matrix(NA, nrow = nrow(x), ncol = ndim)

    for(i in 1:nrow(y)) {

      suby_centered <- scale(y[-i,], scale = FALSE, center = TRUE)
      subx_centered <- scale(x[-i,], scale = FALSE, center = TRUE)

      svd <- Morpho:::svd2B(subx_centered, suby_centered)

      if(i == 1) refyaxis <- svd$v[, 1]
      if(i == 1) refxaxis <- svd$u[, 1]

      if(length(refyaxis) > 1) {
        yscores[i,] <- (y[i,] %*% svd$v) * sign(stats::cor(svd$v[,1], refyaxis))
      } else {
        yscores[i,] <- (y[i,] %*% svd$v) * sign(svd$v[,1] * refyaxis)
      }

      if(length(refxaxis) > 1) {
        xscores[i,] <- (x[i,] %*% svd$u) * sign(stats::cor(svd$u[,1], refxaxis))
      } else {
        xscores[i,] <- (x[i,] %*% svd$u) * sign(svd$u[,1] * refxaxis)
      }

    }

    yscores_centered <- scale(yscores, scale = FALSE, center = TRUE)
    xscores_centered <- scale(xscores, scale = FALSE, center = TRUE)

    if(recompute == TRUE) {

      rotation_y <- t(t(yscores) %*% t(t(y)))
      rotation_x <- t(t(xscores) %*% t(t(x)))

    }

    results <- list(yscores = yscores, yrotation = rotation_y, ycenter = y_center, ytotvar = totvar_y,
                    xscores = xscores, xrotation = rotation_x, xcenter = x_center, xtotvar = totvar_x,
                    values = values)
    class(results) <- "pls"
    return(results)

  }

}


########################################################################################

#' 2B Partial Least Squares for shape data
#'
#' @description A wrapper for [pls2b()] aimed specifically at synthesizing
#'   covariation between shape data and other external, non-shape variable(s).
#'
#' @param shapes Shape data.
#' @param x A matrix with variables as columns and observations as rows,
#'   representing the external variables that supervise ordination.
#' @param LOOCV Logical; whether to apply leave-one-out cross-validation.
#' @param recompute Logical; whether to re-compute rotation matrix using the
#'   scores resulting from LOOCV.
#'
#' @return
#'
#' @seealso \code{\link{pls2b}}
#'
#' @export
#'
#' @examples
pls_shapes <- function(shapes, x, LOOCV = FALSE, recompute = FALSE) {

  y <- shapes_mat(shapes)$data2d
  parlesqu <- pls2b(y = y, x = x, LOOCV = LOOCV, recompute = recompute)

  results <- list(sdev = stats::sd(parlesqu$yscores),
                  rotation = parlesqu$yrotation,
                  x = parlesqu$yscores,
                  center = parlesqu$ycenter,
                  totvar = parlesqu$ytotvar)
  class(results) <- "pls_shape"
  return(results)

}


##########################################################################

#' Calculate percentages of variation accounted by synthetic axes
#'
#' @description Calculate the percentage of total original variation accounted
#'   by syntethtic axes generated by different multivariate ordination methods.
#'
#' @param ord An ordination (i.e. a \code{"prcomp"}, \code{"bg_prcomp"},
#'   \code{"phy_prcomp"} or \code{"pls_shape"} object).
#'
#' @return
#' @export
#'
#' @examples
exp_var <- function(ord) {

  ax_var <- apply(ord$x, 2, stats::var)

  if(class(ord) == "prcomp") {
    totvar <- sum(ax_var)
  } else {
    totvar <- ord$totvar
  }

  acc_var <- 100 * (ax_var / totvar)
  tab <- as.data.frame(round(x = cbind(variance=acc_var,
                                       cummulative=cumsum(acc_var)),
                             digits = 5))

  if(class(ord) == "prcomp") axname <- "PC"
  if(class(ord) == "bg_prcomp") axname <- "bgPC"
  if(class(ord) == "phy_prcomp") axname <- "phyPC"
  if(class(ord) == "pls") axname <- "PLS-"

  rownames(tab) <- paste0(axname, 1:nrow(tab))
  return(tab)

}
