################################################################################

#' Phylogenetic Principal Component Analysis
#'
#' @description A wrapper for [phytools::phyl.pca()].
#'
#' @param x A matrix with one or more variables as columns and observations as
#'   rows. Row must be named and match tip labels from the phylogenetic
#'   \code{tree}.
#' @param tree A \code{"phylo"} object containing a phylogenetic tree. Tip
#'   labels should match the row names from \code{x}.
#' @param corr Logical; whether to use correlation instead of covariance matrix
#'   as input.
#' @param ... Further arguments passed to [phytools::phyl.pca()].
#'
#' @details Phylogenetic PCA finds the linear combination of variables (in
#'   the context of \code{morphospace} will generally be a series of shapes
#'   arranged as 2-margin matrix) maximizing the residual variation left after
#'   removing covariation explained by phylogenetic history (i.e., they reflect
#'   the covariance that would correspond to a star phylogeny), assuming a
#'   Brownian model of evolution.
#'
#'   Phylogenetic PCA has some important differences relative to regular PCA.
#'   First, the resulting ordination is centered around the phylogenetic mean
#'   (i.e., the values estimated for the root of the tree) instead of the
#'   overall centroid of the original variables. More importantly, both
#'   phylogenetic PCA's eigenvectors and eigenvalues have been constructed using
#'   covariation that is independent of phylogenetic structure. However, only
#'   the orientation of the scores on those axes, and not their variances,
#'   reflect this adjustment. In other words, orientation of the phylogenetic
#'   PC axes are devoid of phylogenetic structure, but magnitudes measured in
#'   the resulting morphospace (e.g. distances, variances) still retain
#'   phylogenetic information. Because of this, the variances computed using
#'   phylogenetic scores differ from the ones calculated from the phylogenetic
#'   eigenvalues (which represent the total amount of variance among the taxa
#'   after removing covariance accounted by phylogeny, although their magnitude
#'   depends on the units used to measure branch length), and do not necessarily
#'   decrease for subordinate axes. Also, the set of phylogenetic scores are not
#'   uncorrelated, meaning they can contain redundant information. For more
#'   details on the difference between phylogenetic and regular PCA, see Polly
#'   et al. 2013.
#'
#'   Like PCA, phylogenetic PCA does not change the dimensionality of the shape
#'   data set -- the number of resulting pPC axes will be be equal to the
#'   number of original variables (unless this is higher than the number of
#'   observations, in which case the number of resulting pPC will be equal
#'   to the latter).
#'
#'
#' @return A \code{"phy_prcomp"} object formatted following the \code{"prcomp"}
#' class:
#' \itemize{
#'   \item \code{$x:} {a matrix with the scores of observations in the new
#'     ordination axes.}
#'   \item \code{$sdev:} {the standard deviations of the principal components
#'     (i.e., the square roots of the eigenvalues of the covariance/correlation
#'     matrix).}
#'   \item \code{$rotation:} {a matrix of eigenvector coefficients.}
#'   \item \code{$center:} {the phylogenetic mean (i.e., the shape estimated
#'     for the root of the tree).}
#'   \item \code{$totvar:} {the sum of the variances from all the original
#'     variables.}
#'   \item \code{$lambda, $logL:} {fitted value of lambda and log-likelihood
#'     of the model; see \code{\link[phytools]{phyl.pca}}.}
#'   }
#'
#' @seealso \code{\link[phytools]{phyl.pca}}, \code{\link[stats]{prcomp}},
#'   \code{\link{exp_var}}
#'
#' @export
#'
#' @references
#' Revell, L. J. (2009). \emph{Size-correction and principal components for
#'   interspecific comparative studies}. Evolution, 63, 3258-3268.
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
#' #load data and packages
#' library(geomorph)
#' data("tails")
#'
#' #compute mean shapes for all species and extract the phylogenetic tree
#' sp_shapes <- expected_shapes(shapes = tails$shapes, x = tails$data$species)
#' tree <- tails$tree
#'
#' #perform phylogenetic PCA
#' ppca <- phy_prcomp(x = two.d.array(sp_shapes), tree = tree)
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

  phypca <- suppressWarnings(phytools::phyl.pca(Y = x, tree = tree, mode = mode, ...))

  if(is.null(phypca$lambda)) {
    lambda <- 1
  } else {
    lambda <- phypca$lambda
  }

  C <- ape::vcv.phylo(tree)[rownames(x), rownames(x)]
  anc <- c(phytools::phyl.vcv(x, C, lambda)$alpha)

  results <- list(sdev = suppressWarnings(sqrt(diag(phypca$Eval))), rotation = phypca$Evec,
                  x = phypca$S, center = anc, totvar = totvar)
  if(!is.null(phypca$lambda)) results$lambda <- phypca$lambda
  if(!is.null(phypca$logL)) results$logL <- phypca$logL

  class(results) <- "phy_prcomp"
  return(results)

}

################################################################################

#' Phylogenetically Aligned Component Analysis
#'
#' @description Performs Phylogenetically Aligned Component Analysis (PACA).
#'   Experimental.
#'
#' @param x A matrix with one or more variables as columns and observations as
#'   rows. Rows must be named and match tip labels from the phylogenetic
#'   \code{tree}.
#' @param tree A \code{"phylo"} object containing a phylogenetic tree. Tip
#'   labels should match the row names from \code{x}.
#'
#' @details Phylogenetically aligned component analysis is an ordination method
#'   that finds a new set of axes as linear combinations of the original trait
#'   variables that are maximally aligned with the phylogenetic signal (i.e.,
#'   covariation between trait variation and phylogenetic structure) present in
#'   the data.
#'
#'   Phylogenetically aligned component analysis shares some important
#'   properties with phylogenetic PCA that are worth having in mind. Most
#'   importantly, the transformation achieved by Phylogenetically aligned
#'   component analysis to maximize alignment to phylogenetic signal only
#'   concerns the orientation of the new axes (PACs), and not their overall
#'   dispersion (i.e., Euclidean distances among observations are preserved when
#'   measured across the entire set of PACs). Like phylogenetic PCA, scores
#'   projected into the different PACs can be correlated (even though singular
#'   vectors are orthogonal), and their dispersion does not necessarily decrease
#'   for lower axes (even though singular values do progressively decay).
#'   Also, A phylogenetic covariance matrix is computed assuming a Brownian
#'   model of trait evolution, and the ordination is centered around the
#'   phylogenetic mean (i.e., the root of the phylogenetic tree).
#'
#'   This analysis produces the same number of PACs than the number of original
#'   variables. These axes account for increasingly small proportions of
#'   covariation between phylogenetic structure and trait variation -- although
#'   the total amount of phylogenetic signal present in the data is preserved
#'   for the entire set of PACs.
#'
#' @return A \code{"phyalign_comp"} object formatted following the \code{"prcomp"}
#' class:
#' \itemize{
#'   \item \code{$x:} {a matrix with the scores of observations in the new
#'     ordination axes.}
#'   \item \code{$sdev:} {the standard deviations of the principal components
#'     (i.e., the square roots of the eigenvalues of the covariance/correlation
#'     matrix).}
#'   \item \code{$rotation:} {a matrix of eigenvector coefficients.}
#'   \item \code{$center:} {the phylogenetic mean (i.e., the shape estimated
#'     for the root of the tree).}
#'   \item \code{$totvar:} {the sum of the variances from all the original
#'     variables.}
#'   \item \code{$lambda, $logL:} {fitted value of lambda and log-likelihood
#'     of the model; see \code{\link[phytools]{phyl.pca}}.}
#'   }
#'
#' @export
#'
#' @seealso \code{\link[geomorph]{gm.prcomp}}, \code{\link{phy_prcomp}}
#'
#' @references
#' Collyer, M. L., & Adams, D. (2021). \emph{Phylogenetically aligned component
#' analysis}. Methods in Ecology and Evolution, 12(2), 359-372.
#'
#' @examples
#'
phyalign_comp <- function(x, tree, corr = FALSE) {

  C <- ape::vcv.phylo(tree)[rownames(x), rownames(x)]
  anc <- c(phytools::phyl.vcv(x, C, 1)$alpha)
  z <- scale(x = x, center = anc, scale = corr)

  decomp <- svd(C %*% z)
  vectors <- decomp$v
  scores <- z %*% vectors
  values <- decomp$d

  totalRV <- sum(values^2) / (sqrt(sum(diag(t(C) %*% C)) * sum(diag(t(z) %*% z))))
  partialRVs <- values^2 / (sqrt(sum(diag(t(C) %*% C)) * sum(diag(t(z) %*% z))))

  results <- list(sdev = sqrt((decomp$d^2) / (nrow(decomp$u) - 1)), rotation = vectors,
                  x = scores, center = anc, RVs = round(partialRVs / totalRV, 3))

  class(results) <- "phyalign_prcomp"
  return(results)

}


################################################################################

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
#'   ordination axes. This method is preferred here to LDA/CVA as a way to
#'   produce ordinations maximizing separation between groups because it avoids
#'   spherization of shape variation carried out for the former methods.
#'
#'   Recently, it has been pointed out that bgPCA produces spurious separation
#'   between groups when the number of variables exceeds the number of
#'   observations (which is a common situation in geometric morphometrics
#'   analyses). This problem can be alleviated by carrying out a
#'   leave-one-out cross-validation (LOOCV; i.e., each observation is
#'   excluded from the calculation of bgPCA prior to its projection in the
#'   resulting ordination as a way to calculate its score).
#'
#'   The dimensionality of the ordination space generated by bgPCA will be equal
#'   to the number of groups minus one, the number of original variables, or
#'   the number of observations, whichever is lower.
#'
#' @return A \code{"bg_prcomp"} object formatted following the \code{"prcomp"}
#' class:
#' \itemize{
#'   \item \code{$x:} {a matrix with the scores of observations in the new
#'     ordination axes.}
#'   \item \code{$sdev:} {the standard deviations of the principal components
#'     (i.e., the square roots of the eigenvalues of the covariance/correlation
#'     matrix).}
#'   \item \code{$rotation:} {a \code{n x (g - 1)} matrix of eigenvector
#'     coefficients (with \code{g} being the number of groups.}
#'   \item \code{$center:} {the mean values of the original variables for the
#'     entire sample (i.e., the grand mean).}
#'   \item \code{$grcenters:} {the mean values of the original variables for
#'     each group.}
#'   \item \code{$totvar:} {the sum of the variances from all the original
#'     variables.}
#'   }
#'
#' @seealso \code{\link[stats]{prcomp}}, \code{\link{exp_var}}
#'
#' @references
#' Mitteroecker, P., & Bookstein, F. (2011). \emph{Linear discrimination,
#'   ordination, and the visualization of selection gradients in modern
#'   morphometrics}. Evolutionary Biology, 38(1), 100-114.
#'
#' Bookstein, F. L. (2019). \emph{Pathologies of between-groups principal
#'   components analysis in geometric morphometrics}. Evolutionary Biology,
#'   46(4), 271-302.
#'
#' Cardini, A., O'Higgins, P., & Rohlf, F. J. (2019). \emph{Seeing distinct
#'   groups where there are none: spurious patterns from between-group PCA}.
#'   Evolutionary Biology, 46(4), 303-316.
#'
#' Cardini, A., & Polly, P. D. (2020). \emph{Cross-validated between group PCA
#'   scatterplots: A solution to spurious group separation?}. Evolutionary
#'   Biology, 47(1), 85-95.
#'
#' Rohlf, F. J. (2021). \emph{Why clusters and other patterns can seem to be
#'   found in analyses of high-dimensional data}. Evolutionary Biology, 48(1),
#'   1-16.
#'
#' Thioulouse, J., Renaud, S., Dufour, A. B., & Dray, S. (2021).
#'   \emph{Overcoming the spurious groups problem in between-group PCA}.
#'   Evolutionary Biology, 48(4), 458-471.
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

  if(corr == TRUE) {
    vcv_g <- stats::cov.wt(groupmeans, wt = wts, cor = corr)$cor
  } else {
    vcv_g <- stats::cov.wt(groupmeans, wt = wts, cor = corr)$cov
  }
  svd <- svd(x = vcv_g)

  ndims <- min(nrow(x), ncol(x), nlevels(groups) - 1)
  rotation <- cbind(svd$v[,seq_len(ndims)])
  values <- svd$d[seq_len(ndims)]

  if(LOOCV == FALSE) {

    scores <- proj_eigen(x, rotation, colMeans(x))

  } else {

    refaxis <- rotation[, 1]

    scores_l <- lapply(seq_len(nrow(x)), function(i) {
      subx <- x[-i,]
      subgroups <- groups[-i]

      subgroupmeans <- apply(X = subx, MARGIN = 2, FUN = tapply, subgroups, mean)

      subn_g <- tapply(subgroups, subgroups, length)
      if(gweights == TRUE) {
        subwts <- as.numeric(subn_g / sum(subn_g))
      } else {
        subwts <- rep(1, nlevels(subgroups))
      }

      if(corr == TRUE) {
        subvcv_g <- stats::cov.wt(subgroupmeans, wt = subwts,
                                  cor = corr)$cor
      } else {
        subvcv_g <- stats::cov.wt(subgroupmeans, wt = subwts,
                                  cor = corr)$cov
      }
      subsvd <- svd(subvcv_g)

      subrotation <- cbind(subsvd$v[,seq_len(ndims)])

      iscore <- proj_eigen(x[i,], subrotation, center = colMeans(x)) *
        sign(stats::cor(subrotation[,1], refaxis))

    })

    scores <- do.call("rbind", scores_l)
    scores <- scale(scores, center = TRUE, scale = FALSE)

    if(recompute == TRUE) rotation <- matrix(t(t(scores) %*% t(t(x))), ncol = (nlevels(groups) - 1))

  }

  results <- list(sdev = sqrt((values^2) / (nrow(cbind(x)) - 1)),
                  rotation = rotation, x = cbind(scores), center = grandmean,
                  grcenters = groupmeans, totvar = totvar)
  class(results) <- "bg_prcomp"
  return(results)
}


################################################################################

#' Two-blocks Partial Least Squares
#'
#' @description Performs 2B Partial Least Squares allowing for leave-one-out
#'   cross-validation and removal of phylogenetic covariation (experimental).
#'
#' @param x A matrix with one or more variables as columns and observations as
#'   rows, representing the first block.
#' @param y A matrix with one or more variables as columns and observations as
#'   rows, representing the second block.
#' @param tree A \code{"phylo"} object containing a phylogenetic tree. Tip
#'   labels should match the row number and names from \code{x} and \code{y}.
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
#'   carrying out a leave-one-out cross-validation (LOOCV; i.e., each
#'   observation is excluded from the calculation of PLS axes before its
#'   projection in the resulting ordination as a way to calculate its score).
#'
#'   If \code{tree} is supplied, rows of \code{x} and \code{y} are assumed to
#'   be linked by those phylogenetic relationships. In this case, the resulting
#'   axes maximize the residual covariation among blocks left after removing
#'   covariation among blocks accounted by phylogenetic history (assuming a
#'   Brownian model of evolution and 100% phylogenetic signal, which is
#'   equivalent to setting \code{method = "BM"} in [phytools::phyl.pca()]).
#'
#'   The phylogenetic version of [pls2b()] displays the same variational
#'   properties than phylogenetic PCA (i.e., centering on the phylogenetic mean;
#'   orientation of scores reflect non-phylogenetic covariation but their
#'   variance is not scaled and thus contain phylogenetic information; for the
#'   latter reason, variance of scores and eigenvalues differ; scores are
#'   correlated; etc; see Polly et al. 2013).
#'
#'   The number of axes in both blocks will be equal to the number of variables
#'   in the smallest block (i.e., the one with less variables), unless the
#'   number of observations is lower (in which case this will be equal to the
#'   number PLS axes).
#'
#' @return A \code{"pls2b"} or \code{"phy_pls2b"} object, containing:
#' \itemize{
#'   \item \code{$values:} {vector of singular values accounting for the
#'     covariation among blocks explained by each par of axes.}
#'   \item \code{$xscores:} {a matrix with the scores of observations in the
#'     ordination axes of the first block.}
#'   \item \code{$yscores:} {a matrix with the scores of observations in the
#'     ordination axes of the Second block.}
#'   \item \code{$xrotation:} {matrix of vector coefficients for the first
#'     block.}
#'   \item \code{$yrotation:} {matrix of vector coefficients for the second
#'     block.}
#'   \item \code{$xcenter:} {the mean values of the original variables from
#'     the first block (or their phylogenetic mean, if \code{tree} is
#'     provided).}
#'   \item \code{$ycenter:} {the mean values of the original variables from
#'     the second block (or their phylogenetic mean, if \code{tree} is
#'     provided).}
#'   \item \code{$xtotvar:} {the sum of the variances from the original
#'     variables from the second block.}
#'   \item \code{$ytotvar:} {the sum of the variances from the original
#'     variables from the second block.}
#'   }
#'
#' @seealso \code{\link{pls_shapes}}, \code{\link{phy_prcomp}},
#'   \code{\link{exp_var}}
#'
#' @export
#'
#' @references
#' Rohlf, F. J., & Corti, M. (2000). \emph{Use of two-block partial
#'   least-squares to study covariation in shape}. Systematic Biology, 49(4),
#'   740-753.
#'
#' Zelditch, M. L., Swiderski, D. L., & Sheets, H. D. (2012). \emph{Geometric
#'   morphometrics for biologists: A primer}. 2nd ed. Academic Press.
#'
#' Polly, P. D., Lawing, A. M., Fabre, A. C., & Goswami, A. (2013).
#'   \emph{Phylogenetic principal components analysis and geometric
#'   morphometrics}. Hystrix, 24(1), 33.
#'
#' Monteiro, L. R. (2013). \emph{Morphometrics and the comparative method:
#'   studying the evolution of biological shape}. Hystrix, the Italian Journal
#'   of Mammalogy, 24(1), 25-32.
#'
#' Bookstein, F. L. (2019). \emph{Pathologies of between-groups principal
#'   components analysis in geometric morphometrics}. Evolutionary Biology,
#'   46(4), 271-302.
#'
#' @examples
#' #load data and packages
#' library(geomorph)
#' data("tails")
#'
#' #extract mean log sizes and shapes for all species, as well as the
#' #phylogenetic tree
#' shapes <- tails$shapes
#' sizes <- log(tails$sizes)
#' sp_shapes <- expected_shapes(shapes = shapes, x = tails$data$species)
#' sp_sizes <- tapply(X = sizes, INDEX = tails$data$species, FUN = mean)
#' tree <- tails$tree
#'
#' #perform PLS between shape and size
#' pls <- pls2b(y = two.d.array(shapes), x = sizes)
#'
#' #look at the results
#' names(pls) #the contents of the resulting object
#' plot(pls$xscores, pls$yscores) #ordination
#'
#' #pls_shapes achieves identical results, but it is formatted to be used more
#' #easily when analyzing shape variables and its output is compatible with
#' #other functions from morphospace
#' pls2 <- pls_shapes(shapes = shapes, X = sizes)
#'
#' names(pls2) #the contents of the resulting object
#' exp_var(pls2) #variance explained by each axis
#' plot(pls2$x2, pls2$x) #ordination
#'
#'
#'
#' #perform phylogenetic PLS
#' ppls <- pls2b(y = two.d.array(sp_shapes), x = sp_sizes, tree = tree)
#'
#' #look at the results
#' names(ppls) #the contents of the resulting object
#' plot(ppls$xscores, ppls$yscores) #ordination
#'
#'
#' #pls_shapes achieves identical results, but it is formatted to be used more
#' #easily when analyzing shape variables and its output is compatible with
#' #other functions from morphospace
#' ppls2 <- pls_shapes(shapes = sp_shapes, X = sp_sizes, tree = tree)
#'
#' names(ppls2) #the contents of the resulting object
#' exp_var(ppls2) #variance explained by each axis
#' plot(ppls2$x2, ppls2$x) #ordination
pls2b <- function(x, y, tree = NULL, LOOCV = FALSE, recompute = FALSE) {

  x <- cbind(x)
  y <- cbind(y)

  if(is.null(rownames(x))) rownames(x) <- seq_len(nrow(x))
  if(is.null(rownames(y))) rownames(y) <- seq_len(nrow(y))

  namesx <- rownames(x)
  namesy <- rownames(y)

  if(!is.null(tree)){

    if(!all(length(tree$tip.label) == nrow(x), length(tree$tip.label) == nrow(y))) {
      stop("Number of tips in the tree does not match the number of observations in x and/or y data sets")
    }

    if(!all(tree$tip.label %in% namesy, tree$tip.label %in% namesx)) {
      stop("Names in phylogenetic tree does not match names in x and/or y data sets")
    } else {
      x <- cbind(x[tree$tip.label,])
      y <- cbind(y[tree$tip.label,])
    }

    anc <- apply(cbind(x, y), 2, phytools::fastAnc, tree = tree)[1,]
    x_center <- anc[seq_len(ncol(x))]
    y_center <- anc[(ncol(x) + 1):(ncol(x) + ncol(y))]

  } else {
    y_center <- colMeans(y)
    x_center <- colMeans(x)
  }

  totvar_x <- sum(apply(x, 2, stats::var))
  totvar_y <- sum(apply(y, 2, stats::var))

  svd <- svd_block(x, y, tree)

  ndims <- length(svd$d)
  #values <- (svd$d^2)
  values <- svd$d
  y_rotation <- svd$v
  x_rotation <- svd$u


  if(LOOCV == FALSE) {

    yscores <- proj_eigen(y, y_rotation, y_center)[namesy,]
    xscores <- proj_eigen(x, x_rotation, x_center)[namesx,]

  } else {

    yscores <- NULL
    xscores <- NULL


    refyaxis <- y_rotation[,1]
    refxaxis <- x_rotation[,1]

    for(i in seq_len(nrow(y))) {

      subx <- cbind(x[-i,])
      suby <- cbind(y[-i,])

      if(!is.null(tree)) {
        subtree <- ape::drop.tip(phy = tree, tip = i)
      } else {
        subtree <- NULL
      }

      subsvd <- svd_block(subx, suby, subtree)

      suby_rotation <- subsvd$v
      subx_rotation <- subsvd$u


      if(length(refyaxis) > 1) {
        yscores <- rbind(yscores,
                         proj_eigen(y[i,], suby_rotation, y_center) *
                           sign(stats::cor(suby_rotation[,1], refyaxis)) )
      } else {
        yscores <- rbind(yscores,
                         proj_eigen(y[i,], suby_rotation, y_center) *
                           sign(suby_rotation[,1] * refyaxis) )
      }

      if(length(refxaxis) > 1) {
        xscores <- rbind(xscores,
                         proj_eigen(x[i,], subx_rotation, x_center) *
                           sign(stats::cor(subx_rotation[,1], refxaxis)) )
      } else {
        xscores <- rbind(xscores,
                         proj_eigen(x[i,], subx_rotation, x_center) *
                           sign(subx_rotation[,1] * refxaxis) )
      }

    }

    if(is.null(tree)) {
      yscores <- t(t(yscores) - colMeans(yscores))
      xscores <- t(t(xscores) - colMeans(xscores))
    }

    if(recompute == TRUE) {
      y_rotation <- t(t(yscores) %*% t(t(y)))
      x_rotation <- t(t(xscores) %*% t(t(x)))
    }

    rownames(yscores) <- rownames(y)
    rownames(xscores) <- rownames(x)
    yscores <- yscores[namesy,]
    xscores <- xscores[namesx,]

  }


  results <- list(yscores = cbind(yscores), yrotation = cbind(y_rotation),
                  ycenter = y_center, ytotvar = totvar_y,
                  xscores = cbind(xscores), xrotation = cbind(x_rotation),
                  xcenter = x_center, xtotvar = totvar_x,
                  sdev = sqrt((values^2) / (nrow(cbind(xscores)) - 1)))#values = values)
  if(is.null(tree)) {
    class(results) <- "pls2b"
  } else {#
    class(results) <- "phy_pls2b"
  }
  return(results)
}


################################################################################

#' 2B Partial Least Squares for shape data
#'
#' @description A wrapper for [pls2b()] specifically aimed at summarising
#'   covariation between shape data and other external, non-shape variable(s).
#'
#' @param X A matrix with variables as columns and observations as rows,
#'   representing the external variables that supervise ordination
#'   (corresponding to the first block).
#' @param shapes Shape data (corresponding to the second block of
#' \code{pls2b}.
#' @param tree A \code{"phylo"} object containing a phylogenetic tree. Tip
#'   labels should match the row number and names from \code{x} and \code{y}.
#' @param LOOCV Logical; whether to apply leave-one-out cross-validation.
#' @param recompute Logical; whether to re-compute rotation matrix using the
#'   scores resulting from LOOCV.
#'
#' @details This function finds the linear combination maximizing covariation
#'   between a block of variables assumed to be shape variables formatted as
#'   a 2-margins matrix and another block which can be either another set of
#'   shape variables or one or more non-shape variables for supervising the
#'   analysis. If a phylogenetic tree is supplied, observations from \code{x}
#'   and \code{shapes} blocks are treated as data measured for its tips.
#'
#'   It has been reported that PLS (as an algebraic equivalent of bgPCA)
#'   produces spurious covariation between blocks when the number of variables
#'   exceeds the number of observations (which is a common situation in
#'   geometric morphometrics analyses). This problem can be alleviated by
#'   carrying out a leave-one-out cross-validation (LOOCV; i.e., each
#'   observation is excluded from the calculation of PLS axes before its
#'   projection in the resulting ordination as a way to calculate its score).
#'
#'   The number of axes in both blocks will be equal to the number of variables
#'   in the smallest block (i.e., the one with less variables), unless the
#'   number of observations is lower (in which case this will be equal to the
#'   number PLS axes).
#'
#' @return A \code{"pls_shapes"} or \code{"phy_pls_shapes"} object formatted
#' following the \code{"prcomp"} class:
#' \itemize{
#'   \item \code{$x:} {the scores from the shapes block.}
#'   \item \code{$x2:} {the scores from the supervising block (i.e., the
#'     \code{x} scores.}
#'   \item \code{$sdev:} {the standard deviation of the PLS axis (axis of the
#'      shape block).}
#'   \item \code{$rotation:} {a matrix of vector coefficients for the shape
#'     block.}
#'   \item \code{$center:} {the mean values of the original variables from
#'     the shape block.}
#'   \item \code{$totvar:} {the sum of the variances from the original
#'     variables in the shape block.}
#'   }
#'
#' @seealso \code{\link{pls2b}}
#'
#' @export
#'
#' @references
#' Rohlf, F. J., & Corti, M. (2000). \emph{Use of two-block partial
#'   least-squares to study covariation in shape}. Systematic Biology, 49(4),
#'   740-753.
#'
#' Zelditch, M. L., Swiderski, D. L., & Sheets, H. D. (2012). \emph{Geometric
#'   morphometrics for biologists: A primer}. 2nd ed. Academic Press.
#'
#' Bookstein, F. L. (2019). \emph{Pathologies of between-groups principal
#'   components analysis in geometric morphometrics}. Evolutionary Biology,
#'   46(4), 271-302.
#'
#' @examples
#' #load data
#' data("tails")
#'
#' #extract shapes and sizes, compute mean shapes and sizes for all species,
#' #extract tree
#' shapes <- tails$shapes
#' sizes <- log(tails$sizes)
#' sp_shapes <- expected_shapes(shapes, tails$data$species)
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
#' ppls <- pls_shapes(shapes = sp_shapes, X = sp_sizes, tree = tree,
#'                    LOOCV = TRUE)
#'
#' #inspect results
#' names(ppls) #the contents of the resulting object
#' exp_var(ppls) #variance explained by each axis
#' plot(ppls$x2, ppls$x) #ordination
pls_shapes <- function(X, shapes, tree = NULL, LOOCV = FALSE, recompute = FALSE) {

  y <- shapes_mat(shapes)$data2d
  x <- cbind(X)

  pls <- pls2b(y = y, x = x, tree = tree, LOOCV = LOOCV, recompute = recompute)

  results <- list(sdev = sqrt((values^2) / (nrow(cbind(x)) - 1)),
                  #sdev = apply(pls$yscores, 2, stats::sd),
                  rotation = pls$yrotation,
                  x = pls$yscores,
                  x2 = pls$xscores,
                  center = pls$ycenter,
                  totvar = pls$ytotvar)

  if(inherits(pls, "pls2b")) {
    class(results) <- "pls_shapes"
  } else {
    class(results) <- "phy_pls_shapes"
  }
  return(results)

}


#' Burnaby's orthogonal subspace computation
#'
#' @description Performs Burnaby's multivariate approach for standardizing
#'   morphometric data relative to nuisance factors by computing and projecting
#'   data onto an orthogonal subspace.
#'
#' @param x A matrix with one or more variables as columns and observations as
#'   rows, representing the variables the user desires to be rendered
#'   independent (orthogonal) to \code{vars}.
#' @param vars A matrix with one or more variables as columns and observations
#'   as rows, representing the source of variation the user wants to remove from
#'   \code{x}.
#' @param tree A \code{"phylo"} object containing a phylogenetic tree. Tip
#'   labels should match the row names from \code{x} and \code{vars}.
#' @param axmat An optional matrix of axes coefficients to render \code{x}
#'   orthogonal to. If provided, is used instead of \code{vars} to orthogonalize
#'   \code{x}.
#'
#' @details Originally devised by Burnaby (1966) to eliminate allometric
#'   variation from morphometric data sets by computing a subspace that is
#'   orthogonal (i.e., mathematically independent) to an allometric vector,
#'   and then projecting the data onto these new set of axes (However, this
#'   approach can be used to standardize morphometric data for any type of
#'   arbitrary vector(s) -- commonly one or more nuisance factors or covariates
#'   introducing variation the user wishes to exclude from the analysis).
#'
#'   In particular, this function is intended to provide a way for computation
#'   of ordination (sub)spaces orthogonal to other specific sources of
#'   variation, as well as visualization of shape variation associated to their
#'   ordination axes, most commonly as part of the \code{\link{mspace()}} %>%
#'   \code{proj_*} pipeline. If the user wishes to use the Burnaby approach for
#'   standardization of morphometric (shape) data, this is better achieved
#'   through \code{\link{detrend_shapes()}} by setting
#'   \code{method = "orthogonal"}.
#'
#'   The subspace produced in this way has one less dimension for each axis or
#'   variable used for standardization, relative to the original variables.
#'   However, the user should be aware that, when analyzing shape variation, a
#'   rather large portion of variation in \code{x} that is not orthogonal to
#'   \code{vars} will end up 'hidden' in the lowest ordination axes (i.e., those
#'   effectively lost during standardization procedures that remove non-shape
#'   variation, and therefore showing null eigenvalues). Therefore, projection
#'   of the original shape variables into the orthogonal subspace should be
#'   attempted after removing (at least) these null axes.
#'
#' @return A \code{"burnaby"} or \code{"phy_burnaby"} object, containing:
#' \itemize{
#'   \item \code{$x:} {a matrix with the scores of observations in the new
#'     ordination axes.}
#'   \item \code{$sdev:} {the standard deviations of the principal components
#'     (i.e., the square roots of the eigenvalues of the covariance matrix).}
#'   \item \code{$rotation:} {a matrix of eigenvector coefficients.}
#'   \item \code{$center:} {the mean values of the original variables (or their
#'     phylogenetic mean, if \code{tree} is provided).}
#'   \item \code{$totvar:} {the sum of the variances from all the original
#'     variables.}
#'   }
#'
#' @export
#'
#' @seealso \code{\link{detrend_shapes}}, \code{\link{exp_var}}
#'
#' @references
#' Burnaby, T. P. (1966) \emph{Growth-invariant discriminant functions and
#'   generalized distances. Biometrics, 22, 96–110.}
#'
#' Claude, J. (2008). \emph{Morphometrics with R}. Springer Science & Business
#'   Media, 316.
#'
#' Revell, L. J. (2009). \emph{Size-correction and principal components for
#'   interspecific comparative studies}. Evolution, 63, 3258-3268.
#'
#' Klingenberg, C. P. (2016). \emph{Size, shape, and form: concepts of allometry
#'   in geometric morphometrics}. Development Genes and Evolution, 226(3),
#'   113-137.
#'
#' @examples
#'
#'
burnaby <- function(x, vars = NULL, tree = NULL, axmat = NULL) {

  totvar <- sum(apply(x, 2, stats::var))

  if(is.null(vars) & is.null(axmat)) stop("one of vars or axmat must be provided")

  if(!is.null(axmat)) {
    axmat <- cbind(axmat)
    center <- if(is.null(tree)) colMeans(x) else apply(x, 2, phytools::fastAnc, tree = tree)[1,]
    if(!is.null(tree)) warning("\naxmat has been provided directly, it will be assumed its orientation is already phylogenetically corrected; tree will be only used to center x")
  } else {
    vars <- as.data.frame(vars)

    if(is.null(rownames(vars))) rownames(vars) <- seq_len(nrow(vars))
    if(is.null(rownames(x))) rownames(x) <- seq_len(nrow(x))
    namesvars <- rownames(vars)
    namesx <- rownames(x)

    formula <- as.formula(paste("~", paste(colnames(vars), collapse = " + ")))
    designmat <- stats::model.matrix(formula, data = vars)


    if(!is.null(tree)) {

      if(!all(length(tree$tip.label) == nrow(vars), length(tree$tip.label) == nrow(x))) {
        stop("Number of tips in the tree does not match the number of observations in vars and/or x data sets")
      }
      if(!all(tree$tip.label %in% namesx, tree$tip.label %in% namesvars)) {
        stop("Names in phylogenetic tree does not match names in vars and/or x data sets")
      } else {
        vars <- cbind(vars[tree$tip.label,])
        x <- cbind(x[tree$tip.label,])
        designmat <- cbind(designmat[tree$tip.label,])
      }

      center <- apply(x, 2, phytools::fastAnc, tree = tree)[1,]
      C <- ape::vcv.phylo(tree)
      coefs <- solve(t(designmat) %*% solve(C) %*% designmat) %*%
        t(designmat) %*% solve(C) %*% x

    } else {
      center <- colMeans(x)
      coefs <- solve(t(designmat) %*% designmat) %*% t(designmat) %*% x
    }

    axmat <- if(nrow(coefs) < 3) cbind(coefs[-1,]) else cbind(t(coefs[-1,]))

  }

  I <- diag(1, ncol(x))
  orthobasis <- I - axmat %*% MASS::ginv(t(axmat) %*% axmat) %*% t(axmat)

  ndims <- ncol(x) - ncol(axmat)

  z <- t(t(rbind(x)) - center)
  covmat <- cov(z %*% orthobasis)
  values <- eigen(covmat)$values[seq_len(ndims)]
  vectors <- eigen(covmat)$vectors[,seq_len(ndims)]

  scores <- z %*% vectors
  sdev <- suppressWarnings(sqrt(values))
  sdev[is.nan(sdev)] <- 0

  results <- list(x = scores, rotation = vectors, center = center,
                  sdev = sdev, totvar = totvar)
  class(results) <- if(!is.null(tree)) "phy_burnaby" else "burnaby"
  return(results)

}


################################################################################

#' Calculate percentages of variation accounted by ordination axes
#'
#' @description Calculate the percentage of total original variation accounted
#'   by axes generated by different multivariate ordination methods.
#'
#' @param ordination An ordination (i.e., a \code{"prcomp"}, \code{"bg_prcomp"},
#'   \code{"phy_prcomp"}, \code{"pls_shapes"} or \code{"phy_pls_shapes"}
#'   object).
#'
#' @return A table informing the percentages and cumulative percentages of
#'   original variation accounted by each ordination axis.
#'
#' @export
#'
#' @seealso \code{\link[stats]{prcomp}}, \code{\link{bg_prcomp}},
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

  if(inherits(ordination, "prcomp")) {
    totvar <- sum(ax_var)
  } else {
    totvar <- ordination$totvar
  }

  acc_var <- 100 * (ax_var / totvar)
  tab <- as.data.frame(round(x = cbind(variance=acc_var,
                                       cummulative=cumsum(acc_var)),
                             digits = 5))

  if(inherits(ordination, "prcomp")) axname <- "PC"
  if(inherits(ordination, "bg_prcomp")) axname <- "bgPC"
  if(inherits(ordination, "phy_prcomp")) axname <- "phyPC"
  if(inherits(ordination, "pls_shapes")) axname <- "PLS-"
  if(inherits(ordination, "phy_pls_shapes")) axname <- "phyPLS-"

  rownames(tab) <- paste0(axname, seq_len(nrow(tab)))
  return(tab)

}


