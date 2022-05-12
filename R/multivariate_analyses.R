
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
#'   is independent of phylogenetic structure. However, only the orientation of
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
#'   (i.e. the square roots of the eigenvalues of the covariance/correlation
#'   matrix).}
#'   \item \code{$rotation:} {a matrix of eigenvector coefficients.}
#'   \item \code{$center:} {the phylogenetic mean (i.e. the shape estimated
#'   for the root of the tree).}
#'   \item \code{$totvar:} {the sum of the variances from all the original
#'   variables.}
#'   \item \code{$lambda, $logL:} {fitted value of lambda and log-likelihood
#'   of the model; see \code{\link[phytools]{phyl.pca}}.}
#'   }
#'
#' @seealso \code{\link[phytools]{phyl.pca}}, \code{\link[base]{prcomp}},
#'   \code{\link{exp_var}}
#'
#' @export
#'
#' @references
#' Revell, L. J. (2009). \emph{Size-correction and principal components
#'   for interspecific comparative studies}. Evolution, 63, 3258-3268.
#'
#' Polly, P. D., Lawing, A. M., Fabre, A. C., & Goswami, A. (2013).
#'   \emph{Phylogenetic principal components analysis and geometric
#'   morphometrics}. Hystrix, 24(1), 33-41.
#'
#' Monteiro, L. R. (2013). \emph{Morphometrics and the comparative method:
#'   studying the evolution of biological shape}. Hystrix, the Italian Journal
#'   of Mammalogy, 24(1), 25-32.
#'
#' @examples
#' #load data
#' data("tails")
#'
#' #compute mean shapes for all species and extract the phylogenetic tree
#' sp_shapes <- consensus(shapes = tails$shapes, index = tails$data$species)
#' tree <- tails$tree
#'
#' #perform phylogenetic PCA
#' ppca <- phy_prcomp(x = geomorph::two.d.array(sp_shapes), tree = tree)
#'
#' #look at the results
#' names(ppca) #the contents of the resulting object
#' exp_var(ppca) #variance explained by each axis
#' plot(ppca$x) #ordination
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
#'   number of observations (i.e. alleviates spurious separation between
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
#'   (i.e. the square roots of the eigenvalues of the covariance/correlation
#'   matrix).}
#'   \item \code{$rotation:} {a \code{n x (g - 1)} matrix of eigenvector
#'   coefficients (with \code{g} being the number of groups.}
#'   \item \code{$center:} {the mean values of the original variables for the
#'   entire sample (i.e. the grand mean).}
#'   \item \code{$totvar:} {the mean values of the original variables for
#'   each group.}
#'   \item \code{$grcenters:} {the sum of the variances from all the original
#'   variables.}
#'   }
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
#' #load data
#' data("shells")
#'
#' #extract species classification and shapes
#' species <- shells$data$species
#' shapes <- shells$shapes$coe
#'
#' #perform between-groups PCA
#' bgpca <- bg_prcomp(x = shapes, groups = species)
#'
#' #look at the results
#' names(bgpca) #the contents of the resulting object
#' exp_var(bgpca) #variance explained by each axis
#' plot(bgpca$x) #ordination
#' hulls_by_group_2D(bgpca$x, species) #add convex hulls for species
bg_prcomp <- function(x, groups, gweights = TRUE,
                      LOOCV = FALSE, recompute = FALSE,
                      corr = FALSE) {

  grandmean <- colMeans(x)
  groupmeans <- apply(X = x, MARGIN = 2, FUN = tapply, groups, mean)
  totvar <- sum(apply(x, 2, stats::var))

  n_g <- tapply(groups, groups, length)
  if(gweights == TRUE) {
    wts <- as.numeric(n_g / sum(n_g))
  } else {
    wts <- rep(1, nlevels(groups))
  }

  vcv_g <- stats::cov.wt(groupmeans, wt = wts, cor = corr)$cov
  svd <- svd(x = vcv_g)

  ndims <- min(nrow(x), nlevels(groups) - 1)
  rotation <- svd$v[,1:ndims]
  values <- svd$d[1:ndims]

  if(LOOCV == FALSE) {

    scores <- proj_eigen(x, rotation, colMeans(x))

  } else {

    refaxis <- rotation[, 1]

    scores_l <- lapply(1:nrow(x), function(i) {
      subx <- x[-i,]
      subgroups <- groups[-i]

      subgroupmeans <- apply(X = subx, MARGIN = 2, FUN = tapply, subgroups, mean)

      subn_g <- tapply(subgroups, subgroups, length)
      if(gweights == TRUE) {
        subwts <- as.numeric(subn_g / sum(subn_g))
      } else {
        subwts <- rep(1, nlevels(subgroups))
      }

      subvcv_g <- stats::cov.wt(subgroupmeans, wt = subwts, cor = corr)$cov
      subsvd <- svd(subvcv_g)

      subrotation <- subsvd$v[,1:ndims]

      iscore <- proj_eigen(x[i,], subrotation, center = colMeans(x)) *
        sign(stats::cor(subrotation[,1], refaxis))

    })

    scores <- do.call("rbind", scores_l)
    scores <- scale(scores, center = TRUE, scale = FALSE)

    if(recompute == TRUE) rotation <- matrix(t(t(scores) %*% t(t(x))), ncol = (nlevels(groups) - 1))

  }

  results <- list(sdev = sqrt(values), rotation = rotation, x = cbind(scores),
                  center = grandmean, grcenters = groupmeans, totvar = totvar)
  class(results) <- "bg_prcomp"
  return(results)
}


########################################################################################

#' Two-blocks Partial Least Squares
#'
#' @description Performs 2B Partial Least Squares allowing for leave-one-out
#'   cross-validation, which is useful one the number of variables exceeds the
#'   number of observations (i.e. alleviates spurious covariation between
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
#'   }
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
#' #load data
#' data("tails")
#'
#' #extract shapes and log sizes
#' shapes <- tails$shapes
#' sizes <- log(tails$sizes)
#'
#' #perform PLS between shape and size
#' pls <- pls2b(y = geomorph::two.d.array(shapes), x = sizes)
#'
#' #look at the results
#' names(pls) #the contents of the resulting object
#' plot(pls$xscores, pls$yscores) #ordination
#'
#' #pls_shapes achieves identical results, but it is formatted to be used more
#' #easily when analyzing shape variables and its output is compatible with other
#' #functions from morphospace
#' pls2 <- pls_shapes(shapes = shapes, X = sizes)
#'
#' names(pls2) #the contents of the resulting object
#' exp_var(pls2) #variance explained by each axis
#' plot(pls2$x2, pls2$x) #ordination
pls2b <- function(x, y, LOOCV = FALSE, recompute = FALSE) {

  x <- cbind(x)
  y <- cbind(y)

  totvar_x <- sum(apply(x, 2, stats::var))
  totvar_y <- sum(apply(y, 2, stats::var))

  y_center <- colMeans(y)
  x_center <- colMeans(x)

  ndims <- min(nrow(x), ncol(x), ncol(y))
  svd <- Morpho:::svd2B(x, y)
  values <- (svd$d^2)[1:ndims]
  y_rotation <- cbind((svd$v)[,1:ndims])
  x_rotation <- cbind((svd$u)[,1:ndims])

  if(LOOCV == FALSE) {

    yscores <- proj_eigen(y, y_rotation, y_center)
    xscores <- proj_eigen(x, x_rotation, x_center)

  } else {

    yscores <- c()
    xscores <- c()

    refyaxis <- y_rotation[,1]
    refxaxis <- x_rotation[,1]

    for(i in 1:nrow(y)) {

      suby <- y[-i,]
      subx <- x[-i,]

      subsvd <- Morpho:::svd2B(subx, suby)

      if(length(refyaxis) > 1) {
        yscores <- rbind(yscores,
                         proj_eigen(y[i,], subsvd$v[,1:ndims], y_center) *
                           sign(stats::cor(subsvd$v[,1], refyaxis)) )
      } else {
        yscores <- rbind(yscores,
                         proj_eigen(y[i,], subsvd$v[,1:ndims], y_center) *
                           sign(subsvd$v[,1] * refyaxis) )
      }

      if(length(refxaxis) > 1) {
        xscores <- rbind(xscores,
                         proj_eigen(x[i,], subsvd$u[,1:ndims], x_center) *
                           sign(stats::cor(subsvd$u[,1], refxaxis)) )
      } else {
        xscores <- rbind(xscores,
                         proj_eigen(x[i,], subsvd$u[,1:ndims], x_center) *
                           sign(subsvd$u[,1] * refxaxis) )
      }

    }

    yscores <- scale(yscores, center = TRUE, scale = FALSE)
    xscores <- scale(xscores, center = TRUE, scale = FALSE)

    if(recompute == TRUE) {
      rotation_y <- t(t(yscores) %*% t(t(y)))
      rotation_x <- t(t(xscores) %*% t(t(x)))
    }
  }

  results <- list(yscores = cbind(yscores), yrotation = y_rotation, ycenter = y_center, ytotvar = totvar_y,
                  xscores = cbind(xscores), xrotation = x_rotation, xcenter = x_center, xtotvar = totvar_x,
                  values = values)
  class(results) <- "pls2b"
  return(results)

}

########################################################################################

#' Phylogenetic Two-blocks Partial Least Squares
#'
#' Performs phylogenetic 2B Partial Least Squares allowing for leave-one-out cross-validation. Experimental.
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
#'   equivalent to setting \code{method = "BM"} in [phyl.pca()]).
#'
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
#'   }
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
#' #load data
#' data("tails")
#'
#' #extract mean log sizes and shapes for all species, as well as the phylogenetic tree
#' sp_shapes <- consensus(shapes = tails$shapes, index = tails$data$species)
#' sp_sizes <- tapply(X = log(tails$sizes), INDEX = tails$data$species, FUN = mean)
#' tree <- tails$tree
#'
#' #perform phylogenetic PLS
#' ppls <- phy_pls2b(y = geomorph::two.d.array(sp_shapes), x = sp_sizes, tree = tree)
#'
#' #look at the results
#' names(ppls) #the contents of the resulting object
#' plot(ppls$xscores, ppls$yscores) #ordination
#'
#'
#' #pls_shapes achieves identical results, but it is formatted to be used more
#' #easily when analyzing shape variables and its output is compatible with other
#' #functions from morphospace
#' ppls2 <- pls_shapes(shapes = sp_shapes, X = sp_sizes, tree = tree)
#'
#' names(ppls2) #the contents of the resulting object
#' exp_var(ppls2) #variance explained by each axis
#' plot(ppls2$x2, ppls2$x) #ordination
phy_pls2b <- function(x, y, tree, LOOCV = FALSE, recompute = FALSE) {

  x <- cbind(x)
  y <- cbind(y)

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

  ndims <- min(nrow(x), ncol(x), ncol(y))
  values <- (svd$d^2)[1:ndims]
  rotations <- list(svd$v, svd$u)
  whichy <- which(unlist(lapply(rotations, nrow)) == ncol(y))
  whichx <- which(unlist(lapply(rotations, nrow)) == ncol(x))
  y_rotation <- cbind(rotations[[whichy]][,1:ndims])
  x_rotation <- cbind(rotations[[whichx]][,1:ndims])

  if(LOOCV == FALSE) {

    yscores <- proj_eigen(y, y_rotation, y_anc)[namesy,]
    xscores <- proj_eigen(x, x_rotation, x_anc)[namesx,]

  } else {

    yscores <- c()
    xscores <- c()

    refyaxis <- y_rotation[,1]
    refxaxis <- x_rotation[,1]

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


      if(length(refyaxis) > 1) {
        yscores <- rbind(yscores,
                         proj_eigen(y[i,], suby_rotation, y_anc) *
                           sign(stats::cor(suby_rotation[,1], refyaxis)) )
      } else {
        yscores <- rbind(yscores,
                         proj_eigen(y[i,], suby_rotation, y_anc) *
                           sign(suby_rotation[,1] * refyaxis) )
      }

      if(length(refxaxis) > 1) {
        xscores <- rbind(xscores,
                         proj_eigen(x[i,], subx_rotation, x_anc) *
                           sign(stats::cor(subx_rotation[,1], refxaxis)) )
      } else {
        xscores <- rbind(xscores,
                         proj_eigen(x[i,], subx_rotation, x_anc) *
                           sign(subx_rotation[,1] * refxaxis) )
      }

    }

    yscores <- scale(yscores, center = TRUE, scale = FALSE)
    xscores <- scale(xscores, center = TRUE, scale = FALSE)

    if(recompute == TRUE) {
      y_rotation <- t(t(yscores) %*% t(t(y)))
      x_rotation <- t(t(xscores) %*% t(t(x)))
    }

    rownames(yscores) <- rownames(y)
    rownames(xscores) <- rownames(x)
    yscores <- yscores[namesy,]
    xscores <- xscores[namesx,]

  }

  results <- list(yscores = cbind(yscores), yrotation = y_rotation, ycenter = y_anc, ytotvar = totvar_y,
                  xscores = cbind(xscores), xrotation = y_rotation, xcenter = x_anc, xtotvar = totvar_x,
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
#'   }
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
#' #load data
#' data("tails")
#'
#' #extract shapes and sizes, compute mean shapes and sizes for all species, extract tree
#' shapes <- tails$shapes
#' sizes <- log(tails$sizes)
#' sp_shapes <- consensus(shapes, tails$data$species)
#' sp_sizes <- tapply(X = log(sizes), INDEX = tails$data$species, FUN = mean)
#' tree <- tails$tree
#'
#' #perform PLS between shape and size
#' pls <- pls_shapes(shapes = shapes, X = sizes)
#'
#' #inspect results
#' names(pls) #the contents of the resulting object
#' exp_var(pls) #variance explained by each axis
#' plot(pls$x2, pls$x) #ordination
#'
#' #perform PLS between shape and size with leave-one-out CV
#' pls_cv <- pls_shapes(shapes = shapes, X = sizes, LOOCV = TRUE)
#'
#' #inspect results
#' names(pls_cv) #the contents of the resulting object
#' exp_var(pls_cv) #variance explained by each axis
#' plot(pls_cv$x2, pls_cv$x) #ordination
#'
#' #perform phylogenetic PLS between shape and size (for species means)
#' ppls <- pls_shapes(shapes = sp_shapes, X = sp_sizes, tree = tree)
#'
#' #inspect results
#' names(ppls) #the contents of the resulting object
#' exp_var(ppls) #variance explained by each axis
#' plot(ppls$x2, ppls$x) #ordination
#'
#' #perform phylogenetic PLS between shape and size with leave-one-out CV
#' ppls <- pls_shapes(shapes = sp_shapes, X = sp_sizes, tree = tree, LOOCV = TRUE)
#'
#' #inspect results
#' names(ppls) #the contents of the resulting object
#' exp_var(ppls) #variance explained by each axis
#' plot(ppls$x2, ppls$x) #ordination
pls_shapes <- function(X, shapes, tree = NULL, LOOCV = FALSE, recompute = FALSE) {

  y <- shapes_mat(shapes)$data2d

  if(is.null(tree)) {
    pls <- pls2b(y = y, x = X, LOOCV = LOOCV, recompute = recompute)
  } else {
    pls <- phy_pls2b(y = y, x = X, tree = tree, LOOCV = LOOCV, recompute = recompute)
  }

  results <- list(sdev = apply(pls$yscores, 2, stats::sd),
                  rotation = pls$yrotation,
                  x = pls$yscores,
                  x2 = pls$xscores,
                  center = pls$ycenter,
                  totvar = pls$ytotvar)

  if(class(pls) == "pls2b") {
    class(results) <- "pls_shapes"
  } else {
    class(results) <- "phy_pls_shapes"
  }
  return(results)

}


########################################################################################

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
#' #load data
#' data("shells")
#'
#' #extract shapes, log sizes and species classification
#' species <- shells$data$species
#' shapes <- shells$shapes$coe
#' sizes <- log(shells$sizes)
#'
#' #perform PCA, bgPCA, and PLS
#' pca <- prcomp(shapes)
#' bgpca <- bg_prcomp(shapes, groups = species)
#' pls <- pls_shapes(shapes, X = sizes)
#'
#' #compare percentages of total variation accounted by the first axes when
#' #the analysis is unsupervised (PCA), or supervised by species (bgPCA) or
#' #by log size (PLS)
#' head(exp_var(pca))
#' exp_var(bgpca)
#' exp_var(pls)
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
