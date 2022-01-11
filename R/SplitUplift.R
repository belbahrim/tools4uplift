SplitUplift <- function(data, p, group){
  
  # Splits the data with respect to uplift distribution.
  #
  # Args:
  #   data: a data frame containing the treatment, the outcome and the predictors.
  #   p: The desired sample size. p is a value between 0 and 1 expressed as a decimal, 
  #      it is set to be proportional to the number of observations per group.
  #   group: Your grouping variables. Generally, for uplift modelling, this should be 
  #          a vector of treatment and response variables names, e.g. c("treat", "y").
  #
  # Returns:
  #   The training and validation data sets.
  
  data$ID = seq(1:nrow(data))
  
  train <- data %>% group_by(paste(group, collapse = ',')) %>%  sample_frac(p)
  train <- as.data.frame(train[,-ncol(train)])
  valid <- data %>% anti_join(train, by = "ID")
  
  dataSplit <- list(train[,-ncol(train)], valid[,-ncol(valid)])
  class(dataSplit) <- "SplitUplift"
  return(dataSplit)
}

# END FUN