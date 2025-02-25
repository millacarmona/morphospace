###########################################################################

#' Find Pareto front
#'
#' @description Find Pareto optimal observations in a set of performance data
#'   given the optimality directions of each performance variable.
#'
#' @param performance No. observations-by-No. performance variables matrix of
#'   performance data.
#' @param opt.increases Logical vector of length no. performance variables
#'   indicating optimality direction of each performance variable.
#' @param reverse Logical indicating whether to reverse the order of output.
#'
#' @details This function finds the Pareto Optimal Subset, or the Pareto front,
#'   of a set of performance data.
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
                         opt.increases,
                         reverse = FALSE) {


  # Find number of variables and observations
  nobs <- nrow(performance)
  nvars <- ncol(performance)

  # Check length of opt.increases vector matches number of performance variables
  if(length(opt.increases) != nvars)
    stop("Length of opt.increases vector doesn't match the number of performance variables")

  #reverse optimality direction of selected variables
  input <- performance
  for(i in seq_len(nvars)) {
    if(!opt.increases[i]) input[,i] <- input[,i] * -1
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
    # the opt.increases vector
    if(isOpt) opt <- c(opt,i)
  }

  # For ease of plotting, sort opt.increases vector by the size of variable 1
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
#' @param opt.increases Logical vector of length no. performance variables
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
                        opt.increases) {


  # Find number of variables and observations
  nobs <- nrow(performance)
  nvars <- ncol(performance)

  # Assign row names to performance data
  row.names(performance) <- 1:nobs

  # Check length of opt.increases vector matches number of performance variables
  if(length(opt.increases) != nvars)
    stop("Length of opt.increases vector doesn't match the number of performance variables")


  # While the performance array still has data...
  rank <- 0
  output <- vector(mode = "numeric", length = nobs)
  while(nrow(performance) != 0) {

    # Update rank
    rank <- rank + 1

    # Find the current Pareto front
    PF <- pareto_front(performance = performance, opt.increases = opt.increases)

    # Assign the current Pareto front its rank
    output[as.numeric(row.names(performance)[PF])] <- rank

    # Remove ranked observations from the performance matrix
    performance <- performance[!rownames(performance) %in% row.names(performance)[PF],]
    if(is.null(nrow(performance))) {
      performance <- rbind(performance)
      rownames(performance) <- "1"
    }

    # Save first Pareto front
    if(rank == 1) pfront <- PF
  }

  return(list(rank = output, pfront = pfront))
}

###########################################################################

#' Pareto Rank Ratio
#'
#' @description Find the Pareto rank ratio of a set of performance data.
#'
#' @param performance No. observations-by-No. performance variables matrix of
#'   performance data.
#' @param opt.increases Logical vector of length no. performance variables
#'   indicating optimality direction of each performance variable.
#'
#' @details This function uses the formula from Deakin et al. (2022) to
#'   calculate the Pareto rank ratio of the performance data.
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
                              opt.increases) {

  # Find number of variables and observations
  nobs <- nrow(performance)
  nvars <- ncol(performance)

  # Register original names
  names0 <- rownames(performance)

  # Check length of opt.increases vector matches number of performance variables
  if(length(opt.increases) != nvars)
    stop("Length of opt.increases vector doesn't match the number of performance variables")

  # Get optimal Goldberg ranking
  PR <- pareto_rank(performance, opt.increases)
  Ro <- PR$rank
  pfront <- PR$pfront

  # Get suboptimal Goldberg ranking
  Rs <- pareto_rank(performance, !opt.increases)$rank

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

  return(list(PRR = output, pfront = pfront))
}

