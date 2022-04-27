
###########################################################################

#' Obtain consensus shape(s)
#'
#' @description Compute the mean shape from the entire sample, the mean shape of
#'   a subset of samples, or the mean shape of the levels of a factor (for
#'   landmark configurations).
#'
#' @param shapes A \code{k x p x n} array of superimposed landmarks.
#' @param index Either a numeric vector indicating the configurations to be
#'   averaged or a factor whose levels are used to average groups of
#'   configurations.
#'
#' @return Either a matrix defining a single mean shape \code{p x k} or an array
#'   of shapes, one for each level of \code{index}.
#'
#' @export
#'
#' @examples
consensus <- function(shapes, index = 1:dim(shapes)[3]) {

  p <- nrow(shapes)
  data2d <- geomorph::two.d.array(shapes)

  if(is.numeric(index)) {
    colmeans <- colMeans(data2d[index,])
    consensus <- matrix(colmeans, nrow = p, byrow = TRUE)
  } else {
    if(is.character(index)) index <- factor(index)
    if(is.factor(index)) {
      colmeans <- apply(X = data2d, MARGIN = 2, FUN = tapply, index, mean)
      consensus <- geomorph::arrayspecs(A = colmeans, p = p, k = ncol(shapes))
      dimnames(consensus)[[3]] <- levels(index)
    }
  }

  return(consensus)
}


###########################################################################

#' Remove shape variation associated to external variables
#'
#' @description A function to detrend (i.e. standardize) shape data using the
#'   functional relationship between shape data and some external explanatory
#'   variable(s) (works for both factors and numerics), estimated using a linear
#'   model.
#'
#' @param model A \code{mlm} object created using [lm()].
#' @param xvalue A value (numeric) or level (character) at which shape data is
#'   to be standardized (i.e. centered); If NULL, the mean of the complete
#'   sample is used.
#' @param newx,newy New data to be standardized instead of the used in
#'   \code{model}. Coefficients are taken from the linear model and applied to
#'   the new data to predict shapes at the desired \code{xvalue}.
#'
#' @return
#' @export
#'
#' @examples
detrend_shapes <- function(model, xvalue = NULL, newx = NULL, newy = NULL){

  n <- nrow(model$residuals)
  coefs <- model$coefficients
  resids <- model$resid

  x <- model$model[, ncol(model$model)]

  if(!is.null(newy)) {
    designmat0 <- cbind(1, newx)
    resids <- newy - designmat0 %*% coefs
  } else {
    grandmean <- colMeans(model$resid + model$fitted.values)

    grandmean_vec <- rep(1, n) %*% t(grandmean)
    predicted_mat <- resids + grandmean_vec
  }

  if(!is.null(xvalue)) {
    if(is.numeric(x) == TRUE) {
      designmat <- cbind(1, xvalue)
    }

    if(is.factor(x) == TRUE) {
      designmat <- rep(0, nlevels(x))
      designmat[which(levels(x) == xvalue)] <- 1
      if(isFALSE(designmat[1] == 1)) designmat[1] <- 1
    }

    fitted <- as.numeric(designmat %*% coefs)

    if(!is.null(newy)) n <- nrow(newy)

    fitted_vec <- rep(1, n) %*% t(fitted)
    predicted_mat <- resids + fitted_vec
  }

  return(predicted_mat)

}
