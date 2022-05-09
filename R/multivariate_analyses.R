
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
#' @details Phylogenetic PCA finds the linear combination of variables (in the
#'   the context of \code{morphospace} will generally be a series of shapes
#'   arranged as 2-margin matrix) maximizing the residual variation left after
#'   removing covariation explained by phylogenetic history (i.e. they reflect
#'   the covariance that would correspond to a star phylogeny).
#'
#'   Phylogenetic PCA has some important differences to regular PCA. First, the
#'   resulting ordination is centered around the phylogenetic mean (i.e. the
#'   values estimated for the root of the tree) instead of the overall centroid
#'   of the original variables. More importantly, both phylogenetic PCA's
#'   eigenvectors and eigenvalues have been constructed using covariation that
#'   is independent of phylogenetic structure. However only the orientation of
#'   the scores on those axes, and not their variances, reflect this adjustment.
#'   In other words, orientation of the phylogenetic PC axes are devoid of
#'   phylogenetic structure, but magnitudes measured in the resulting morphospace
#'   (e.g. distances, variances) still retain phylogenetic information. Stemming
#'   from this same fact, the variances computed using phylogenetic scores differ
#'   from the ones calculated from the phylogenetic eigenvalues (which represent
#'   the total amount of variance among the taxa after removing covariance
#'   accounted by phylogeny, although their magnitude depends on the units used
#'   to measure branch length), and do not necessarily decrease for subordinate
#'   axes. Also, the set of phylogenetic scores are not uncorrelated, meaning
#'   they can contain redundant information. For more details, see Polly et al.
#'   2013.
#'
#'
#' @return A \code{"phy_prcomp"} object formatted following the \code{"prcomp"}
#' class:
#' \itemize{
#'   \item \code{$sdev:} {the standard deviations of the principal components
#'   (i.e., the square roots of the eigenvalues of the covariance/correlation
#'   matrix).}
#'   \item \code{$rotation:} {a matrix of eigenvector coefficients.}
#'   \item \code{$center:} {the phylogenetic mean (i.e. the shape estimated
#'   for the root of the tree).}
#'   \item \code{$totvar:} {the sum of the variances from all the original
#'   variables.}
#'   \item \code{$lambda, $logL:} {fitted value of lambda and log-likelihood
#'   of the model; see \code{\link[phytools]{phyl.pca}}.}
#'
#' @seealso \code{\link[phytools]{phyl.pca}}, \code{\link[base]{prcomp}},
#'   \code{\link{exp_var}}
#'
#' @export
#'
#' @references Revell, L. J. (2009). \emph{Size-correction and principal
#' components for interspecific comparative studies}. Evolution, 63, 3258-3268.
#'
#' Polly, P. D., Lawing, A. M., Fabre, A. C., & Goswami, A. (2013).
#' \emph{Phylogenetic principal components analysis and geometric
#' morphometrics}. Hystrix, 24(1), 33.
#'
#' Monteiro, L. R. (2013). \emph{Morphometrics and the comparative method:
#' studying the evolution of biological shape}. Hystrix, the Italian Journal
#' of Mammalogy, 24(1), 25-32.
#'
#' @examples
phy_prcomp <- function(x, tree, corr = FALSE, ...) {

  if(corr == FALSE) {
    mode <- "cov"
  } else {
    mode <- "corr"
  }
  center <- colMeans(x)
  totvar <- sum(apply(x, 2, stats::var))

  phypca <- phytools::phyl.pca(Y = x, tree = tree, mode = mode, ...)

  if(is.null(phypca$lambda)) {
    lambda <- 1
  } else {
    lambda <- phypca$lambda
  }

  C <- ape::vcv.phylo(tree)[rownames(x), rownames(x)]
  anc <- c(phytools::phyl.vcv(x, C, lambda)$alpha)

  results <- list(sdev = sqrt(phypca$Eval), rotation = phypca$Evec,  x = phypca$S,
                  center = anc, totvar = totvar)
  if(!is.null(phypca$lambda)) results$lambda <- phypca$lambda
  if(!is.null(phypca$logL)) results$logL <- phypca$logL

  class(results) <- "phy_prcomp"
  return(results)

}


########################################################################################

#' Between-groups Principal Component Analysis
#'
#' @description Performs between group PCA allowing for leave-one-out
#'   cross-validation, which is useful one the number of variables exceeds the
#'   number of observations (i.e., alleviates spurious separation between
#'   groups).
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
#' @details bgPCA finds the liner combination of variables (which in the
#'   context of \code{morphospace} will generally be a series of shapes
#'   arranged as 2-margin matrix) maximizing variation between groups'
#'   centroids, and then project the actual observation into the resulting
#'   synthetic axes. This method is preferred here to LDA/CVA as a way to
#'   produce ordinations maximizing separation between groups because it avoids
#'   spherization of shape variation carried out for the former methods.
#'
#'   Recently, it has been discovered that bgPCA produces spurious separation
#'   between groups when the number of variables exceeds the number of
#'   observations (which is a common situation in geometric morphometrics
#'   analyses). This problem can be alleviated by carrying out a
#'   leave-one-out cross-validation (LOOCV; i.e. each observation is
#'   excluded from the calculation of bgPCA prior to its projection in the
#'   resulting ordination as a way to calculate its score).
#'
#' @return A \code{"bg_prcomp"} object formatted following the \code{"prcomp"}
#' class:
#' \itemize{
#'   \item \code{$sdev:} {the standard deviations of the principal components
#'   (i.e., the square roots of the eigenvalues of the covariance/correlation
#'   matrix).}
#'   \item \code{$rotation:} {a \code{n x (g - 1)} matrix of eigenvector
#'   coefficients (with \code{g} being the number of groups.}
#'   \item \code{$center:} {the mean values of the original variables for the
#'   entire sample (i.e. the grand mean).}
#'   \item \code{$totvar:} {the mean values of the original variables for
#'   each group.}
#'   \item \code{$grcenters:} {the sum of the variances from all the original
#'   variables.}
#'
#' @seealso \code{\link[base]{prcomp}}, \code{\link{exp_var}}
#'
#' @references Mitteroecker, P., & Bookstein, F. (2011). \emph{Linear
#' discrimination, ordination, and the visualization of selection gradients in
#' modern morphometrics}. Evolutionary Biology, 38(1), 100-114.
#'
#' Bookstein, F. L. (2019). \emph{Pathologies of between-groups principal
#' components analysis in geometric morphometrics}. Evolutionary Biology, 46(4),
#' 271-302.
#'
#' Cardini, A., O’Higgins, P., & Rohlf, F. J. (2019). \emph{Seeing distinct
#' groups where there are none: spurious patterns from between-group PCA}.
#' Evolutionary Biology, 46(4), 303-316.
#'
#' Cardini, A., & Polly, P. D. (2020). \emph{Cross-validated between group PCA
#' scatterplots: A solution to spurious group separation?}. Evolutionary Biology,
#' 47(1), 85-95.
#'
#' Rohlf, F. J. (2021). \emph{Why clusters and other patterns can seem to be
#' found in analyses of high-dimensional data}. Evolutionary Biology, 48(1), 1-16.
#'
#' Thioulouse, J., Renaud, S., Dufour, A. B., & Dray, S. (2021). \emph{Overcoming
#' the spurious groups problem in between-group PCA}. Evolutionary Biology, 48(4),
#' 458-471.
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

#' Two-blocks Partial Least Squares
#'
#' @description Performs 2B Partial Least Squares allowing for leave-one-out
#'   cross-validation, which is useful one the number of variables exceeds the
#'   number of observations (i.e., alleviates spurious covariation between
#'   variables).
#'
#' @param x A matrix with one or more variables as columns and observations as
#'   rows, representing the first block.
#' @param y A matrix with one or more variables as columns and observations as
#'   rows, representing the second block.
#' @param LOOCV Logical; whether to apply leave-one-out cross-validation.
#' @param recompute Logical; whether to re-compute rotation matrix using the
#'   scores resulting from LOOCV.
#'
#' @details Starting with two blocks of variables measured for the same
#'   cases, two-blocks PLS finds the linear combination of variables on each
#'   block maximizing covariation with the variables on the other. In the
#'   context of \code{morphospace}, one or both of these blocks will usually
#'   be a series of shapes arranged as a 2-margins matrix.
#'
#'   It has been reported that PLS (as an algebraic equivalent of bgPCA)
#'   produces spurious covariation between blocks when the number of variables
#'   exceeds the number of observations (which is a common situation in
#'   geometric morphometrics analyses). This problem can be alleviated by
#'   carrying out a leave-one-out cross-validation (LOOCV; i.e. each
#'   observation is excluded from the calculation of PLS axes before its
#'   projection in the resulting ordination as a way to calculate its score).
#'
#' @return A \code{"pls2b"} object, containing:
#' \itemize{
#'   \item \code{$values:} {vector of singular values accounting for the
#'   covariation among blocks explained by each par of axes.}
#'   \item \code{$xrotation:} {matrix of vector coefficients for the first
#'   block.}
#'   \item \code{$yrotation:} {matrix of vector coefficients for the second
#'   block.}
#'   \item \code{$xcenter:} {the mean values of the original variables from
#'   the first block.}
#'   \item \code{$ycenter:} {the mean values of the original variables from
#'   the second block.}
#'   \item \code{$xtotvar:} {the sum of the variances from the original
#'   variables from the second block.}
#'   \item \code{$ytotvar:} {the sum of the variances from the original
#'   variables from the second block.}
#'
#' @seealso \code{\link{pls_shapes}}, \code{\link{phy_pls2b}},
#'   \code{\link{exp_var}}
#'
#' @export
#'
#' @references Rohlf, F. J., & Corti, M. (2000). \emph{Use of two-block partial
#' least-squares to study covariation in shape}. Systematic Biology, 49(4),
#' 740-753.
#'
#' Zelditch, M. L., Swiderski, D. L., & Sheets, H. D. (2012).
#' \emph{Geometric morphometrics for biologists: A primer}. 2nd ed. Academic
#' Press.
#'
#' Bookstein, F. L. (2019). \emph{Pathologies of between-groups principal
#' components analysis in geometric morphometrics}. Evolutionary Biology, 46(4),
#' 271-302.
#'
#' @examples
pls2b <- function(x, y, LOOCV = FALSE, recompute = FALSE) {

  if(is.vector(x)) x <- matrix(x)
  if(is.vector(y)) y <- matrix(y)

  totvar_x <- sum(apply(x, 2, stats::var))
  totvar_y <- sum(apply(y, 2, stats::var))

  y_center <- colMeans(y)
  x_center <- colMeans(x)

  y_centered <- scale(y, scale = FALSE, center = TRUE)
  x_centered <- scale(x, scale = FALSE, center = TRUE)

  svd <- Morpho:::svd2B(x_centered, y_centered)
  values <- svd$d^2
  rotation_y <- svd$v
  rotation_x <- svd$u

  if(LOOCV == FALSE) {

    yscores <- y_centered %*% svd$v
    xscores <- x_centered %*% svd$u

  } else {

    ndim <- min(ncol(y), ncol(x))
    yscores <- matrix(NA, nrow = nrow(y), ncol = ndim)
    xscores <- matrix(NA, nrow = nrow(x), ncol = ndim)

    for(i in 1:nrow(y)) {

      suby_centered <- scale(y[-i,], scale = FALSE, center = TRUE)
      subx_centered <- scale(x[-i,], scale = FALSE, center = TRUE)

      subsvd <- Morpho:::svd2B(subx_centered, suby_centered)

      if(i == 1) refyaxis <- subsvd$v[, 1]
      if(i == 1) refxaxis <- subsvd$u[, 1]

      if(length(refyaxis) > 1) {
        yscores[i,] <- (y[i,] %*% subsvd$v) * sign(stats::cor(subsvd$v[,1], refyaxis))
      } else {
        yscores[i,] <- (y[i,] %*% subsvd$v) * sign(subsvd$v[,1] * refyaxis)
      }

      if(length(refxaxis) > 1) {
        xscores[i,] <- (x[i,] %*% subsvd$u) * sign(stats::cor(subsvd$u[,1], refxaxis))
      } else {
        xscores[i,] <- (x[i,] %*% subsvd$u) * sign(subsvd$u[,1] * refxaxis)
      }

    }

    yscores <- scale(yscores, scale = FALSE, center = TRUE)
    xscores <- scale(xscores, scale = FALSE, center = TRUE)

    if(recompute == TRUE) {
      rotation_y <- t(t(yscores) %*% t(t(y)))
      rotation_x <- t(t(xscores) %*% t(t(x)))
    }
  }

  results <- list(yscores = yscores, yrotation = rotation_y, ycenter = y_center, ytotvar = totvar_y,
                  xscores = xscores, xrotation = rotation_x, xcenter = x_center, xtotvar = totvar_x,
                  values = values)
  class(results) <- "pls2b"
  return(results)

}


##########################################################################

#' Phylogenetic Two-blocks Partial Least Squares
#'
#' @description Performs phylogenetic 2B Partial Least Squares allowing for
#'   leave-one-out cross-validation. Experimental.
#'
#' @param x A matrix with one or more variables as columns and observations as
#'   rows, representing the first block.
#' @param y A matrix with one or more variables as columns and observations as
#'   rows, representing the second block.
#' @param tree A \code{"phy"} object containing a phylogenetic tree. Tip labels
#'   should match the row number and names from \code{x} and \code{y}.
#' @param LOOCV Logical; whether to apply leave-one-out cross-validation.
#' @param recompute Logical; whether to re-compute rotation matrix using the
#'   scores resulting from LOOCV.
#'
#' @details Similar to \code{\link{pls2b}} but cases are linked by phylogenetic
#'   relationships. Like in the case of \code{\link{phy_prcomp}}, the resulting
#'   axes maximize the residual covariation among blocks left after removing
#'   covariation among blocks accounted by phylogenetic history (assuming a
#'   Brownian model of evolution an 100% phylogenetic signal, which is
#'   equivalent to setting \code{method = "BM"} in [phytools::phyl.pca()]).
#'   This method display the same variational properties than phylogenetic PCA,
#'   (i.e. centering on the phylogenetic mean; orientation of scores reflect
#'   non-phylogenetic covariation but their variance is not scaled and thus
#'   contain phylogenetic information; for the latter reason, variance of
#'   scores and eigenvalues differ; scores are correlated; etc; see Polly
#'   et al. 2013).
#'
#'   As in other supervised multivariate ordination functions from
#'   \code{morphospace}, leave-one-out cross-validation can be implemented to
#'   alleviate the spurious covariation among blocks that is produced when
#'   the number of variables exceeds the number of cases.
#'
#' @return A \code{"phy_pls2b"} object, containing:
#' \itemize{
#'   \item \code{$values:} {vector of singular values accounting for the
#'   covariation among blocks explained by each par of axes.}
#'   \item \code{$xrotation:} {matrix of vector coefficients for the first
#'   block.}
#'   \item \code{$yrotation:} {matrix of vector coefficients for the second
#'   block.}
#'   \item \code{$xcenter:} {the mean values of the original variables from
#'   the first block.}
#'   \item \code{$ycenter:} {the mean values of the original variables from
#'   the second block.}
#'   \item \code{$xtotvar:} {the sum of the variances from the original
#'   variables from the second block.}
#'   \item \code{$ytotvar:} {the sum of the variances from the original
#'   variables from the second block.}
#'
#' @seealso \code{\link{pls_shapes}}, \code{\link{pls2b}},
#'   \code{\link{phy_prcomp}}, \code{\link{exp_var}}
#'
#' @export
#'
#' @references
#' Rohlf, F. J., & Corti, M. (2000). \emph{Use of two-block partial
#' least-squares to study covariation in shape}. Systematic Biology, 49(4),
#' 740-753.
#'
#' Bookstein, F. L. (2019). \emph{Pathologies of between-groups principal
#' components analysis in geometric morphometrics}. Evolutionary Biology, 46(4),
#' 271-302.
#'
#' Polly, P. D., Lawing, A. M., Fabre, A. C., & Goswami, A. (2013).
#' \emph{Phylogenetic principal components analysis and geometric
#' morphometrics}. Hystrix, 24(1), 33.
#'
#' Monteiro, L. R. (2013). \emph{Morphometrics and the comparative method:
#' studying the evolution of biological shape}. Hystrix, the Italian Journal
#' of Mammalogy, 24(1), 25-32.
#'
#' @examples
phy_pls2b <- function(x, y, tree, LOOCV = FALSE, recompute = FALSE) {

  if(is.vector(x)) x <- cbind(x)
  if(is.vector(y)) y <- cbind(y)

  namesx <- rownames(x)
  namesy <- rownames(y)

  if(!all(length(tree$tip.label) == nrow(x), length(tree$tip.label) == nrow(y))) {
    stop("Number of tips in the tree does not match the number of observations in x and/or y data sets")
  }

  if(!all(tree$tip.label %in% namesy, tree$tip.label %in% namesx)) {
    stop("Names in phylogenetic tree does not match names in x and/or y data sets")
  } else {
    x <- cbind(x[tree$tip.label,])
    y <- cbind(y[tree$tip.label,])
  }

  totvar_x <- sum(apply(x, 2, stats::var))
  totvar_y <- sum(apply(y, 2, stats::var))

  C <- ape::vcv.phylo(tree)[rownames(x), rownames(x)]

  xy <- cbind(x,y)
  C <- phytools::phyl.vcv(xy, C, 1)$C
  anc <- phytools::phyl.vcv(xy, C, 1)$alpha
  R <- phytools::phyl.vcv(xy, C, 1)$R

  part_R <- R[1:ncol(x), (ncol(x) + 1):(ncol(x) + ncol(y))]
  svd <- svd(part_R)

  x_anc <- anc[1:ncol(x)]
  y_anc <- anc[(ncol(x) + 1):(ncol(x) + ncol(y))]

  x_centered <- t(t(x) - x_anc)
  y_centered <- t(t(y) - y_anc)

  rotations <- list(svd$v, svd$u)
  whichy <- which(unlist(lapply(rotations, nrow)) == ncol(y))
  whichx <- which(unlist(lapply(rotations, nrow)) == ncol(x))
  y_rotation <- rotations[[whichy]]
  x_rotation <- rotations[[whichx]]

  values <- svd$d^2

  if(LOOCV == FALSE) {

    yscores <- (y_centered %*% y_rotation)[namesy,]
    xscores <- (x_centered %*% x_rotation)[namesx,]

  } else {

    ndim <- min(ncol(y), ncol(x))
    yscores <- matrix(NA, nrow = nrow(y), ncol = ndim)
    xscores <- matrix(NA, nrow = nrow(x), ncol = ndim)

    for(i in 1:nrow(y)) {

      subx <- cbind(x[-i,])
      suby <- cbind(y[-i,])
      subtree <- ape::drop.tip(phy = tree, tip = i)

      subC <- ape::vcv.phylo(subtree)[rownames(subx), rownames(subx)]

      subxy <- cbind(subx, suby)
      subC <- phytools::phyl.vcv(subxy, subC, 1)$C
      subanc <- phytools::phyl.vcv(subxy, subC, 1)$alpha
      subR <- phytools::phyl.vcv(subxy, subC, 1)$R

      subpart_R <- subR[1:ncol(subx), (ncol(subx) + 1):(ncol(subx) + ncol(suby))]
      subsvd <- svd(subpart_R)

      subx_anc <- anc[1:ncol(subx)]
      suby_anc <- anc[(ncol(subx) + 1):(ncol(subx) + ncol(suby))]

      subrotations <- list(subsvd$v, subsvd$u)
      whichy <- which(unlist(lapply(subrotations, nrow)) == ncol(suby))
      whichx <- which(unlist(lapply(subrotations, nrow)) == ncol(subx))
      suby_rotation <- subrotations[[whichy]]
      subx_rotation <- subrotations[[whichx]]

      if(i == 1) refyaxis <- suby_rotation[, 1]
      if(i == 1) refxaxis <- subx_rotation[, 1]


      iy_centered <- (t(cbind(y[i,])) - y_anc)
      reorienty <- sign(stats::cor(suby_rotation[,1], refyaxis))
      if(length(refyaxis) > 1) {
        yscores[i,] <- (iy_centered %*% suby_rotation) * reorienty
      } else {
        yscores[i,] <- (iy_centered %*% suby_rotation) * sign(suby_rotation[,1] * refyaxis)
      }

      ix_centered <- (t(cbind(x[i,])) - x_anc)
      reorientx <- sign(stats::cor(subx_rotation[,1], refxaxis))
      if(length(refxaxis) > 1) {
        xscores[i,] <- (ix_centered %*% subx_rotation) * reorientx
      } else {
        xscores[i,] <- (ix_centered %*% subx_rotation) * sign(subx_rotation[,1] * refxaxis)
      }

    }

    if(recompute == TRUE) {
      yrotation <- t(t(yscores) %*% t(t(y)))
      xrotation <- t(t(xscores) %*% t(t(x)))
    }

    yscores <- yscores[namesy,]
    xscores <- xscores[namesx,]

  }

  results <- list(yscores = yscores, yrotation = y_rotation, ycenter = y_anc, ytotvar = totvar_y,
                  xscores = xscores, xrotation = y_rotation, xcenter = x_anc, xtotvar = totvar_x,
                  values = values)
  class(results) <- "phy_pls2b"
  return(results)
}



########################################################################################

#' 2B Partial Least Squares for shape data
#'
#' @description A wrapper for [pls2b()] or [phy_pls2b()] aimed specifically
#'   at synthesizing covariation between shape data and other external,
#'   non-shape variable(s).
#'
#' @param X A matrix with variables as columns and observations as rows,
#'   representing the external variables that supervise ordination
#'   (corresponding to the first block).
#' @param shapes Shape data (corresponding to the second block of
#' \code{pls2b}.
#' @param tree A \code{"phy"} object containing a phylogenetic tree. Tip labels
#'   should match the row number and names from \code{x} and \code{y}.
#' @param LOOCV Logical; whether to apply leave-one-out cross-validation.
#' @param recompute Logical; whether to re-compute rotation matrix using the
#'   scores resulting from LOOCV.
#'
#' @details This function finds the linear combination maximizing covariation
#'   between a block of variables assumed to be shape variables formatted as
#'   a 2-margins matrix and another block which can be either another set of
#'   shape variables or one or more non-shape variables for supervizing the
#'   analysis. If a phylogenetic tree is supplied, observations from \code{x}
#'   and \code{shapes} blocks are treated as data measured for its tips.
#'
#'   It has been reported that PLS (as an algebraic equivalent of bgPCA)
#'   produces spurious covariation between blocks when the number of variables
#'   exceeds the number of observations (which is a common situation in
#'   geometric morphometrics analyses). This problem can be alleviated by
#'   carrying out a leave-one-out cross-validation (LOOCV; i.e. each
#'   observation is excluded from the calculation of PLS axes before its
#'   projection in the resulting ordination as a way to calculate its score).
#'
#' @return A \code{"pls_shapes"} or \code{"phy_pls_shapes"} object formatted
#' following the \code{"prcomp"} class:
#' \itemize{
#'   \item \code{$sdev:} {the standard deviation of the PLS axis/axis of the
#'    shape block.}
#'   \item \code{$rotation:} {a matrix of vector coefficients for the shape
#'   block.}
#'   \item \code{$center:} {the mean values of the original variables from
#'   the shape block.}
#'   \item \code{$totvar:} {the sum of the variances from the original
#'   variables in the shape block.}
#'   \item \code{$x2:} {the scores from the supervizing block (i.e. the
#'   \code{x} scores.}
#'
#' @seealso \code{\link{pls2b}}, \code{\link{phy_pls2b}}
#'
#' @export
#'
#' @references Rohlf, F. J., & Corti, M. (2000). \emph{Use of two-block partial
#' least-squares to study covariation in shape}. Systematic Biology, 49(4),
#' 740-753.
#'
#' Zelditch, M. L., Swiderski, D. L., & Sheets, H. D. (2012).
#' \emph{Geometric morphometrics for biologists: A primer}. 2nd ed. Academic
#' Press.
#'
#' Bookstein, F. L. (2019). \emph{Pathologies of between-groups principal
#' components analysis in geometric morphometrics}. Evolutionary Biology, 46(4),
#' 271-302.
#'
#' @examples
pls_shapes <- function(X, shapes, tree = NULL, LOOCV = FALSE, recompute = FALSE) {

  y <- shapes_mat(shapes)$data2d

  if(is.null(tree)) {
    pls <- pls2b(y = y, x = X, LOOCV = LOOCV, recompute = recompute)
  } else {
    pls <- phy_pls2b(y = y, x = X, tree = tree, LOOCV = LOOCV, recompute = recompute)
  }

  yscores <- cbind(pls$yscores)
  xscores <- cbind(pls$xscores)

  results <- list(sdev = apply(yscores, 2, stats::sd),
                  rotation = pls$yrotation,
                  x = yscores,
                  x2 = xscores,
                  center = pls$ycenter,
                  totvar = pls$ytotvar)

  if(class(pls) == "pls2b") {
    class(results) <- "pls_shapes"
  } else {
    class(results) <- "phy_pls_shapes"
  }
  return(results)

}


##########################################################################

#' Calculate percentages of variation accounted by synthetic axes
#'
#' @description Calculate the percentage of total original variation accounted
#'   by syntethtic axes generated by different multivariate ordination methods.
#'
#' @param ordination An ordination (i.e. a \code{"prcomp"}, \code{"bg_prcomp"},
#'   \code{"phy_prcomp"}, \code{"pls_shapes"} or \code{"phy_pls_shapes"}
#'   object).
#'
#' @return A table informing the percentages and cumulative percentages of
#'   original variation accounted by each synthetic axis of the multivariate
#'   ordination.
#'
#' @export
#'
#' @seealso \code{\link[base]{prcomp}}, \code{\link{bg_prcomp}},
#' \code{\link{phy_prcomp}}, \code{\link{pls_shapes}}
#'
#' @examples
exp_var <- function(ordination) {

  ax_var <- apply(ordination$x, 2, stats::var)

  if(class(ordination) == "prcomp") {
    totvar <- sum(ax_var)
  } else {
    totvar <- ordination$totvar
  }

  acc_var <- 100 * (ax_var / totvar)
  tab <- as.data.frame(round(x = cbind(variance=acc_var,
                                       cummulative=cumsum(acc_var)),
                             digits = 5))

  if(class(ordination) == "prcomp") axname <- "PC"
  if(class(ordination) == "bg_prcomp") axname <- "bgPC"
  if(class(ordination) == "phy_prcomp") axname <- "phyPC"
  if(class(ordination) == "pls_shapes") axname <- "PLS-"
  if(class(ordination) == "phy_pls_shapes") axname <- "phyPLS-"

  rownames(tab) <- paste0(axname, 1:nrow(tab))
  return(tab)

}
