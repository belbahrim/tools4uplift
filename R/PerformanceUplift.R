PerformanceUplift <- function(data, treat, outcome, prediction, nb.group = 10, 
                              equal.intervals = TRUE, rank.precision = 2){
  
  # Computes the performance of an uplift estimator.
  #
  # Args:
  #   data: a data frame containing the treatment, the outcome and the predictors.
  #   treat: name of a binary (numeric) vector representing the treatment 
  #          assignment (coded as 0/1).
  #   outcome: name of a binary response (numeric) vector (coded as 0/1).
  #   prediction: a predicted uplift to sorts the observations from highest 
  #               to lowest uplift.
  #   ... and default parameters.
  #
  # Returns:
  #   Descriptive statistics for the performance of an uplift estimator.
  
  # Error handling
  if (nb.group < 2) {
    stop("The number of groups must be greater or equal to 2")
  }
  
  if (rank.precision <= 0) {
    stop("The rank.precision must be 1 or 2")
  }
  
  if (rank.precision > 2) {
    stop("The rank.precision must be 1 or 2")
  }
  
  if (length(unique(data[[prediction]])) < nb.group) {
    equal.intervals <- FALSE
    warning("The number of unique predictions will be used as the number of groups")
  }
  # First, we need to rank and sort the observations
  data$rank <- 0
  data$rank <- round(rank(-data[[prediction]], ties.method = "max") / nrow(data), rank.precision)
  
  if (equal.intervals == TRUE){
    
    dataResults <- data.frame(matrix(rep(0), nb.group, 8))
    colnames(dataResults) <- c("cum_per", "T_Y1", "T_n", "C_Y1", 
                               "C_n", "incremental_Y1", "inc_uplift", "uplift")
    
    data$group <- 0
    for(i in 1:nb.group) {
      data$group[data$rank > (i-1)/nb.group & data$rank <= i/nb.group] <- i
    }
    
    # Incremental observed uplift
    for(i in 1:nb.group){
      subset <- data[data$group <= i, ]
      dataResults[i,1] <- i/nb.group
      dataResults[i,2] <- sum(subset[[treat]] == 1 & subset[[outcome]] == 1)
      dataResults[i,3] <- sum(subset[[treat]] == 1)
      dataResults[i,4] <- sum(subset[[treat]] == 0 & subset[[outcome]] == 1)
      dataResults[i,5] <- sum(subset[[treat]] == 0)
      dataResults[i,6] <- dataResults[i, 2] - dataResults[i, 4]*dataResults[i, 3]/dataResults[i, 5] 
    }
    dataResults[,7] <- dataResults[,6]/dataResults[nb.group,3]*100
    
    
    # Observed uplift in each group
    for (i in 1:nb.group){
      subset <- data[data$group == i, ]
      dataResults[i,8] <- sum(subset[[treat]] == 1 & subset[[outcome]] == 1) / sum(subset[[treat]] == 1) -  
        sum(subset[[treat]] == 0 & subset[[outcome]] == 1) / sum(subset[[treat]] == 0)
    }
    dataResults[,8] <- dataResults[,8]*100
  }
  
  if (equal.intervals == FALSE) {
    
    no_of_groups <- sort(unique(data$rank))
    no_of_groups <- no_of_groups[no_of_groups>0]
    nb.group <- length(no_of_groups)
    
    dataResults <- data.frame(matrix(rep(0), nb.group, 8))
    colnames(dataResults) <- c("cum_per", "T_Y1", "T_n", "C_Y1", 
                               "C_n", "incremental_Y1", "inc_uplift", "uplift")
    
    # Incremental observed uplift
    for (i in 1:nb.group){
      subset <- data[data$rank <= no_of_groups[i], ]
      dataResults[i,1] <- no_of_groups[i]
      dataResults[i,2] <- sum(subset[[treat]] == 1 & subset[[outcome]] == 1)
      dataResults[i,3] <- sum(subset[[treat]] == 1)
      dataResults[i,4] <- sum(subset[[treat]] == 0 & subset[[outcome]] == 1)
      dataResults[i,5] <- sum(subset[[treat]] == 0)
      
      dataResults[i,6] <- dataResults[i, 2] - dataResults[i, 4]*dataResults[i, 3]/dataResults[i, 5] 
    }
    dataResults[,7] <- dataResults[,6]/dataResults[nb.group,3]*100
    
    
    # Observed uplift in each group
    for (i in 1:nb.group){
      if (i == 1) {
        subset <- data[data$rank <= no_of_groups[i], ]
      }
      else {
        subset <- data[no_of_groups[i-1] < data$rank & data$rank <= no_of_groups[i], ]
      }
      dataResults[i,8] <- sum(subset[[treat]] == 1 & subset[[outcome]] == 1) / sum(subset[[treat]] == 1) -  
        sum(subset[[treat]] == 0 & subset[[outcome]] == 1) / sum(subset[[treat]] == 0)
    }
    dataResults[,8] <- dataResults[,8]*100
    
  }
  
  class(dataResults) <- "PerformanceUplift"
  return(dataResults)
}

# END FUN