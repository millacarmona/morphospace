
###########################################################################

#' Create morphospace
#'
#' @param shapes
#' @param axes
#' @param links
#' @param p
#' @param k
#' @param FUN
#' @param nh
#' @param nv
#' @param mag
#' @param asp
#' @param xlim
#' @param ylim
#' @param xlab
#' @param ylab
#' @param rot.models
#' @param size.models
#' @param asp.models
#' @param col.model
#' @param lwd.model
#' @param points
#' @param cex.ldm
#' @param col.ldm
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
mspace <- function(shapes,
                   axes = c(1,2),
                   links,
                   p = NA,
                   k = NA,
                   FUN = prcomp,
                   nh = 5,
                   nv = 4,
                   mag = 1,
                   asp = NA,
                   xlim = NULL,
                   ylim = NULL,
                   xlab = NULL,
                   ylab = NULL,
                   rot.models = 0,
                   size.models = 1,
                   asp.models = 1,
                   col.model = "#708095",
                   lwd.model = 1,
                   points = FALSE,
                   cex.ldm = 1,
                   col.ldm = "black",
                   ...) {


  if(is.na(p) | is.na(k)) {
    if(length(dim(shapes)) == 3){
      p <- nrow(shapes)
      k <- ncol(shapes)
    } else {stop("Provide values for p and/or k, or deliver shapes as an array")}
  }

  if(length(dim(shapes)) == 3) shapes <- geomorph::two.d.array(shapes)

  FUN <- match.fun(FUN)
  pca <- FUN(shapes, ...)

  scores1 = seq(from = min(pca$x[,axes[1]]),
                to   = max(pca$x[,axes[1]]),
                length.out = nh)
  scores2 = seq(from = min(pca$x[,axes[2]]),
                to   = max(pca$x[,axes[2]]),
                length.out = nv)
  vectors = pca$rotation[,axes]
  center = pca$center


  gridcoords <- data.matrix(expand.grid(scores1, scores2))
  gridcoords_mag <- gridcoords * mag

  #sh_arr <- sample_ax(gridcoords_mag, vectors, center, p, k) * size.models
  sh_mat <- rev_eigen(gridcoords_mag, vectors, center)
  sh_arr <- geomorph::arrayspecs(sh_mat, p = p, k = k) * size.models
  sh_arr[,2,] <- sh_arr[,2,] * asp.models


  if(rot.models!=0) for(i in 1:dim(sh_arr)[3]) {
    sh_arr[,,i] <- spdep::Rotation(sh_arr[,,i], rot.models*0.0174532925199)
  }

  models_coords<-c()
  models_arr<-sh_arr * 0
  for(i in 1:nrow(gridcoords)) {
    descentmat <- matrix(rep(gridcoords[i,], p), p, k, byrow = TRUE)
    models_coords <- rbind(models_coords, (sh_arr[,,i] * 0.07) + descentmat)
    models_arr[,,i] <- (sh_arr[,,i] * 0.07) + descentmat
  }

  plot(models_coords, type="n", asp=asp,
       xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab)

  for(i in 1:nrow(gridcoords)) {
    points(models_arr[,,i], pch = 16, cex = cex.ldm * 0.1, col = col.ldm)
    for(l in 1:length(links)) lines(models_arr[,,i][links[[l]],], col = col.model, lwd = lwd.model)
  }

  if(points == TRUE) points(pca$x)

  results <- list(p = p, k = k, plotax = axes, x = pca$x, rotation = pca$rotation, center = pca$center)
  class(results) <- "mspace"
  return(invisible(results))

}


###############################################################################################

#' Project shapes into morphospace
#'
#' @param shapes
#' @param mspace
#' @param pipe
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
# proj_shapes <- function(shapes, mspace, ...) {
#
#   if(length(dim(shapes)) == 2) {
#     data2d <- shapes
#   } else {
#     data2d <- geomorph::two.d.array(shapes)
#   }
#
#   scores <- proj_eigen(data2d, mspace$rotation)
#
#   if(.Device != "null device") {
#     points(scores[, mspace$plotax], ...)
#     return(invisible(scores))
#   } else {
#     return(scores)
#   }
# }
proj_shapes <- function(shapes, mspace, pipe = TRUE, ...) {

  if(length(dim(shapes)) == 2) {
    data2d <- shapes
  } else {
    data2d <- geomorph::two.d.array(shapes)
  }

  scores <- proj_eigen(data2d, mspace$rotation)

  if(.Device != "null device") points(scores[, mspace$plotax], ...)

  if(pipe == FALSE) return(scores)
  if(pipe == TRUE) return(invisible(mspace))

}

###############################################################################################

#' Project consensus shape(s) into morphospace
#'
#' @param consensus
#' @param mspace
#' @param pipe
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
# proj_consensus <- function(consensus, mspace, add = FALSE, ...) {
#
#   name_mspace <- as.character(match.call()$mspace)
#
#   gr_centroids <- proj_shapes(consensus, mspace)
#
#   if(.Device != "null device") points(gr_centroids[, mspace$plotax], ...)
#   if(add == TRUE) {
#     mspace$gr_centroids <- gr_centroids
#     assign(x = name_mspace, value = mspace, envir = globalenv())
#   } else {
#     return(gr_centroids)
#   }
# }
proj_consensus <- function(consensus, mspace, pipe = TRUE, ...) {

  if(length(dim(consensus)) == 2) {
    data2d <- consensus
  } else {
    data2d <- geomorph::two.d.array(consensus)
  }

  gr_centroids <- proj_eigen(data2d, mspace$rotation)

  # gr_centroids <- proj_shapes(consensus, mspace)
  # if(class(gr_centroids) == "mspace") gr_centroids <- gr_centroids$x
  if(.Device != "null device") points(gr_centroids[, mspace$plotax], ...)

  mspace$gr_centroids <- gr_centroids

  if(pipe == FALSE) return(gr_centroids)
  if(pipe == TRUE) return(invisible(mspace))

}


###############################################################################################

#' Project morphometric axis/axes into morphospace
#'
#' @param neword
#' @param mspace
#' @param ax
#' @param p
#' @param k
#' @param pipe
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
proj_axes <- function(neword, mspace, ax = 1, p = NULL, k = NULL, pipe = TRUE, ...) {

  ext_shapes <- vector(mode = "list", length = length(ax))
  names(ext_shapes) <- paste0("axis_", ax)
  for(i in 1:length(ax)){

    coeffs <- neword$rotation[, ax[i]]
    scores <- apply(matrix(neword$x[, ax[i]]), 2, range)
    center <- neword$center

    ext_shapes2d <- rev_eigen(scores = scores, vector = coeffs, center = center)
    ext_scores <- proj_eigen(ext_shapes2d, mspace$rotation[, mspace$plotax])
    if(pipe == FALSE) ext_shapes[[i]] <- arrayspecs(ext_shapes2d, k, p)

    lines(ext_scores, ...)

  }

  if(pipe == FALSE) return(ext_shapes)
  if(pipe == TRUE) return(invisible(mspace))

}


###############################################################################################

#' Project phylogenetic structure into morphospace
#'
#' @param tree
#' @param mspace
#' @param pipe
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
# proj_phylogeny <- function(tree, mspace, add = TRUE, ...) {
#
#   name_mspace <- as.character(match.call()$mspace)
#
#   if(is.null(mspace$gr_centroids)) stop("Group centroids have not been provided; try with proj_consensus(..., add = TRUE)")
#
#   nodes_scores <- apply(mspace$gr_centroids[tree$tip.label,], 2, phytools::fastAnc, tree = tree)
#   phylo_scores <- rbind(mspace$gr_centroids[tree$tip.label,], nodes_scores)
#
#   if(.Device != "null device") {
#     for(i in 1:nrow(tree$edge)) {
#       lines(rbind(phylo_scores[tree$edge[i, 1], mspace$plotax], phylo_scores[tree$edge[i, 2], mspace$plotax]))
#     }
#   }
#
#   if(add == TRUE) {
#     mspace$phylo_scores <- phylo_scores
#     mspace$phylo <- tree
#     assign(x = name_mspace, value = mspace, envir = globalenv())
#   } else {return(phylo_scores)}
#
# }
proj_phylogeny <- function(tree, mspace, pipe = TRUE, ...) {

  if(is.null(mspace$gr_centroids)) stop("Group centroids have not been provided; add proj_consensus() before")

  nodes_scores <- apply(mspace$gr_centroids[tree$tip.label,], 2, phytools::fastAnc, tree = tree)
  phylo_scores <- rbind(mspace$gr_centroids[tree$tip.label,], nodes_scores)

  if(.Device != "null device") {
    for(i in 1:nrow(tree$edge)) {
      lines(rbind(phylo_scores[tree$edge[i, 1], mspace$plotax], phylo_scores[tree$edge[i, 2], mspace$plotax]), ...)
    }
  }

  mspace$phylo_scores <- phylo_scores
  mspace$phylo <- tree

  if(pipe == FALSE) return(phylo_scores)
  if(pipe == TRUE) return(invisible(mspace))

}


###############################################################################################

#' Delimit groups in morphospace
#'
#' @param mspace
#' @param groups
#' @param pipe
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
proj_groups <- function(mspace, groups, pipe = TRUE, ...) {

  hulls_by_group_2D(mspace$x[, mspace$plotax], fac = groups, ...)

  if(pipe == TRUE) return(invisible(mspace))

  }


###############################################################################################

#' Plot points for mspace objects
#'
#' @param mspace
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
points.mspace <- function(mspace, ...) {points(mspace$x[, mspace$plotax], ...)}
