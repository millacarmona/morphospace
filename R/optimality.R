###########################################################################

#' Generate performance space
#'
#' @description Create an empirical performance space using functional
#' performance data. Can be ordinated with a principal component analysis
#' for high dimensional data visualisation.
#'
#' @param performance No. observations-by-No. performance variables matrix of
#' performance data.
#' @param axes Numeric of length 2, indicating the axes to be plotted.
#' @param FUN The function to be used for synthesizing functional
#' variation. the default is NULL, indicating no ordination.
#' @param points Logical indicating whether to plot scatter points of functional
#' variables
#' @param range Logical indicating whether to show the region of explorable
#' performance space
#' @param col.points Scatter plot colour
#' @param col.range Range area colour
#' @param col.border Range border colour
#' @param ... Further arguments passed to \code{FUN}.
#'
#' @details This function is central to the optimality analysis within the
#'   \code{morphospace} workflow. It plots the functional variation within a
#'   morphospace, generated through performance analysis of morphological data,
#'   in order to visualise the range of explorable functional combinations to a
#'   set of morphological data. This can be ordinated with a PCA, or not by
#'   default. The output of \code{pspace} can then be expanded with the
#'   \code{proj_*} family of functions, similar to \code{mspace}. It can also be
#'   used to perform Pareto optimality analyses on the data.
#'
#'
#' @return An object of class \code{"pspace"}, which is a list containing:
#'   \itemize{
#'   \item{$performance:} { raw performance data.}
#'   \item{$x:} { scores of the sample of performances in the synthetic axes.}
#'   \item{$rotation:} { eigenvector's coefficients.}
#'   \item{$center:} { the mean values of the original performance variables.}
#'   \item{$plotinfo:} { a list with the information used to create the plot.}
#'   \item{$ordtype:} { method used for multivariate ordination.}
#'   }
#'


pspace <- function(performance,
                   axes = c(1,2),
                   FUN = NULL,
                   points = TRUE,
                   range = TRUE,
                   col.points = "black",
                   col.range = "#EBEEF1",
                   col.border = "white",
                   ...) {

  # Check performance is the correct size
  if(length(dim(performance)) == 2){
    np <- ncol(performance)
    ns <- nrow(performance)
  } else {stop("Provide performance as a 2D array")}

  # Perform ordination if required
  if (!is.null(FUN)){
    FUN <- match.fun(FUN)
    ordination <- FUN(performance, ...)
    ordtype <- class(ordination)
  } else {
    ordination <- list(sdev = NULL, rotation =  diag(np), x = performance, center = rep(0, times = np))
    ordtype <- NULL
  }

  # Create axes
  plot(ordination$x[,axes], col = "white")

  # Draw the range of explorable performance space if required
  if (range) {
    NE <- pareto_front(performance = ordination$x[,axes], optimality = c(TRUE, TRUE), reverse = TRUE)
    SE <- pareto_front(performance = ordination$x[,axes], optimality = c(TRUE, FALSE), reverse = FALSE)
    SW <- pareto_front(performance = ordination$x[,axes], optimality = c(FALSE, FALSE), reverse = TRUE)
    NW <- pareto_front(performance = ordination$x[,axes], optimality = c(FALSE, TRUE), reverse = FALSE)
    boundary <- c(NE,SE,SW,NW)
    polygon(x = ordination$x[boundary,axes[1]], y = ordination$x[boundary,axes[2]], col = col.range, border = col.border)
  }

  # Draw the data points if required
  if (points)
  {
    graphics::points(ordination$x[,axes], col = col.points, pch = 20)
  }

  # Save plot info
  plotinfo <- list(axes = axes, points = points, range = range,
                   col.points = col.points, col.range = col.range,
                   col.border = col.border)

  # Save results
  results <- list( performance = performance, x = ordination$x,
                   rotation = ordination$rotation,
                   center = ordination$center, ordtype = ordtype,
                   plotinfo = plotinfo)

  class(results) <- "pspace"
  return(results)
}

###########################################################################

#' Project points into performance space
#'
#' @description Project untransformed performance data into a pspace object and
#' plot as individual points
#'
#' @param pspace A \code{"pspace"} object.
#' @param x Performance data to be projected.
#' @param pch.points Point symbol.
#' @param col.points Point colour.
#'

proj_points <- function(pspace,
                        x,
                        pch = 20,
                        col = "black") {
  X <- proj_eigen(x, pspace$rotation, pspace$center)
  graphics::points(X[,pspace$plotinfo$axes], col = col, pch = pch)

}

###########################################################################

#' Project line into performance space
#'
#' @description Project untransformed performance data into a pspace object and
#' plot as a line
#'
#' @param pspace A \code{"pspace"} object.
#' @param x Performance data to be projected.
#' @param lty Line type.
#' @param lwd Line width.
#' @param col Line colour.

proj_line <- function(pspace,
                      x,
                      lty = 1,
                      lwd = 1,
                      col = "black") {
  X <- proj_eigen(x, pspace$rotation, pspace$center)
  graphics::lines(X[,pspace$plotinfo$axes], col = col, lty = lty, lwd = lwd)
}

###########################################################################

#' Find Pareto front
#'
#' @description Find Pareto optimal observations in a set of performance data
#' given the optimality directions of each performance variable
#'
#' @param performance No. observations-by-No. performance variables matrix of
#' performance data.
#' @param optimality Logical vector of length no. performance variables
#' indicating optimality direction of each performance variable.
#' @param reverse Logical indicating whether to reverse the order of output.
#'
#' @details This function finds the Pareto Optimal Subset, or the Pareto front,
#' of a set of performance data. Designed for use with pspace$performance.
#'
#' @return Indices of Pareto optimal observations.
#'
#' @author William J. Deakin
#'
#' @references
#' Deakin, W. J., Anderson, P. S., den Boer, W., Smith, T. J., Hill, J. J.,
#'   Rücklin, M., Donoghue, P. C. J.  & Rayfield, E. J. (2022). \emph{Increasing
#'   morphological disparity and decreasing optimality for jaw speed and
#'   strength during the radiation of jawed vertebrates}. Science Advances,
#'   8(11), eabl3644.
#'
#' @export
pareto_front <- function(performance,
                         optimality,
                         reverse = FALSE) {


  # Find number of variables and observations
  nobs <- nrow(performance)
  nvars <- ncol(performance)

  # Check length of optimality vector matches number of performance variables
  if(length(optimality) != nvars)
    stop("Length of optimality vector doesn't match the number of performance variables")

  #reverse optimality direction of selected variables
  input <- performance
  for(i in seq_len(nvars)) {
    if(!optimality[i]) input[,i] <- input[,i] * -1
  }


  # For every observation...
  opt <- c()
  for (i in seq_len(nobs)) {

    # Assume observation is optimal
    isOpt <- TRUE

    # While observation is optimal, loop through every other observation
    j <- 0
    while(isOpt & j <  nobs) {
      j <- j + 1
      if(i != j) { # Avoid self-comparison
        # If observation i is worse than observation j for all variables...
        subopti <- NULL
        for(k in seq_len(nvars)) subopti[k] <- if(input[i,k] < input[j,k]) TRUE else FALSE

        # ... then observation i is not optimal
        isOpt <- sum(subopti) != nvars
      }
    }


    # If observation i is optimal (i.e., it is better than every other
    # individual observation in at least one performance variable), add it to
    # the optimality vector
    if(isOpt) opt <- c(opt,i)
  }

  # For ease of plotting, sort optimality vector by the size of variable 1
  optNums <- input[opt,1]
  id <- if(reverse){
    sort(optNums, index.return = TRUE)
  } else sort(optNums, decreasing = TRUE, index.return = TRUE)

  return(opt[id$ix])
}


###########################################################################

#' Goldberg Pareto rank
#'
#' @description Find the Goldberg Pareto rank of a set of performance data.
#'
#' @param performance No. observations-by-No. performance variables matrix of
#' performance data.
#' @param optimality Logical vector of length no. performance variables
#' indicating optimality direction of each performance variable.
#'
#' @details This function ranks performance data based on its optimality
#' direction. This is a Goldberg Pareto ranking, which iteratively finds the
#' Pareto front, assigns it an integer rank, then removes it from the dataset
#' and repeats.
#'
#'
#' @return Integer Goldberg Pareto ranks of data.
#'
#' @author William J. Deakin
#'
#' @references
#' Deakin, W. J., Anderson, P. S., den Boer, W., Smith, T. J., Hill, J. J.,
#'   Rücklin, M., Donoghue, P. C. J.  & Rayfield, E. J. (2022). \emph{Increasing
#'   morphological disparity and decreasing optimality for jaw speed and
#'   strength during the radiation of jawed vertebrates}. Science Advances,
#'   8(11), eabl3644.
#'
#' @export
pareto_rank <- function(performance,
                        optimality) {


  # Find number of variables and observations
  nobs <- nrow(performance)
  nvars <- ncol(performance)

  # Assign row names to performance data
  row.names(performance) <- 1:nobs

  # Check length of optimality vector matches number of performance variables
  if(length(optimality) != nvars)
    stop("Length of optimality vector doesn't match the number of performance variables")


  # While the performance array still has data...
  rank <- 0
  output <- vector(mode = "numeric", length = nobs)
  while(nrow(performance) != 0) {

    # Update rank
    rank <- rank + 1

    # Find the current Pareto front
    PF <- pareto_front(performance = performance, optimality = optimality)

    # Assign the current Pareto front its rank
    output[as.numeric(row.names(performance)[PF])] <- rank

    # Remove ranked observations from the performance matrix
    performance <- performance[!rownames(performance) %in% row.names(performance)[PF],]
    if(is.null(nrow(performance))) {
      performance <- rbind(performance)
      rownames(performance) <- "1"
    }
  }

  return(output)
}

###########################################################################

#' Pareto Rank Ratio
#'
#' @description Find the Pareto rank ratio of a set of performance data.
#'
#' @param performance No. observations-by-No. performance variables matrix of
#' performance data.
#' @param optimality Logical vector of length no. performance variables
#' indicating optimality direction of each performance variable.
#'
#' @details This function uses the formula from Deakin et al.
#' (Science Advances, 2022) to calculate the Pareto rank ratio of the
#' performance data.
#'
#' @return Numeric vector containing the Pareto Rank Ratio of the performance
#' data.
#'
#' @author William J. Deakin
#'
#' @references
#' Deakin, W. J., Anderson, P. S., den Boer, W., Smith, T. J., Hill, J. J.,
#'   Rücklin, M., Donoghue, P. C. J.  & Rayfield, E. J. (2022). \emph{Increasing
#'   morphological disparity and decreasing optimality for jaw speed and
#'   strength during the radiation of jawed vertebrates}. Science Advances,
#'   8(11), eabl3644.
#'
#' @export
pareto_rank_ratio <- function(performance,
                              optimality) {

  # Find number of variables and observations
  nobs <- nrow(performance)
  nvars <- ncol(performance)

  # Register original names
  names0 <- rownames(performance)

  # Check length of optimality vector matches number of performance variables
  if(length(optimality) != nvars)
    stop("Length of optimality vector doesn't match the number of performance variables")

  # Get optimal Goldberg ranking
  Ro <- pareto_rank(performance, optimality)

  # Get suboptimal Goldberg ranking
  Rs <- pareto_rank(performance, !optimality)

  # Assign output vector
  output <- vector(mode = "numeric", length = nobs)

  # Calculate the PRR of obs i as: R_i = (Rs_i - 1) / (Ro_i + Rs_i - 2)
  for(i in 1:nobs) {
    if(Ro[i] == 1) {
      output[i] = 1
    } else if(Rs[i] == 1) {
      output[i] = 0
    } else {
      output[i] = (Rs[i] - 1) / (Ro[i] + Rs[i] - 2)
    }
  }

  # Turn PRR into a column and assign original rownames if any.
  output <- cbind(PRR = output)
  rownames(output) <- names0

  return(output)
}

################################################################################
#
# performance <- fdata[,2:3]
# optimality <- c(T,T)
#
# pareto_front(performance, optimality)
# pareto_rank(performance, optimality)
# pareto_rank_ratio(performance, optimality)
# traceback()
#
#
# rev_eigen <- function(scores, vectors, center) { t(t(scores %*% t(vectors)) + center) }
#
# proj_eigen <- function(x, vectors, center) { t(t(rbind(x)) - center) %*% vectors }
