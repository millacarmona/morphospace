##########################################################################

#' Reverse eigenanalysis-based ordination
#'
#' @param scores
#' @param vectors
#' @param center
#'
#' @return
#' @export
#'
#' @examples
rev_eigen <- function(scores, vectors, center) { t(t(scores %*% t(vectors)) + center) }

###########################################################################

#' Project cases into existing eigenanalysis-based ordination
#'
#' @param x
#' @param vectors
#'
#' @return
#' @export
#'
#' @examples
proj_eigen <- function(x, vectors) { x %*% vectors }


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

