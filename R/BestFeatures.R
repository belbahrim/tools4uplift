BestFeatures <- function (data, treat, outcome, predictors, rank.precision = 2, 
                            equal.intervals = FALSE, nb.group = 10, 
                            validation = TRUE, p = 0.3) {
    # Feature selection for the interaction estimator.
    #
    # Args:
    #   data: a data frame containing the treatment, the outcome and the predictors.
    #   treat: name of a binary (numeric) vector representing the treatment 
    #          assignment (coded as 0/1).
    #   outcome: name of a binary response (numeric) vector (coded as 0/1).
    #   predictors: a vector of names representing the predictors to consider in the model.
    #   ... and default parameters.
    #
    # Returns:
    #   The best features for the interaction estimator.
  
  
    # All variables must be continuous. Before using this function, change 
    # categorical variables to dummies.
    formula <- formulaUplift(treat, outcome, predictors)
      
    # Cross-validation
    # Split data between train (data) and valid using SplitUplift() function
    if (validation == TRUE) {
      split <- SplitUplift(data, 1 - p, c(treat, outcome))
      data <- split[[1]]
      valid <- split[[2]]
    }
    path <- LassoPath(data, formula)
    path <- path[!duplicated(path[,"dimension"]), ]
    # Keep paths of dimension > 0
    path <- path[path[, "dimension"] > 0, ]
    lambda.qini <- c()
    for (k in 1:nrow(path)) {
      features <- path[k, -c(1, 2)]
      # Keep features with non zero estimates only
      features <- features[features != 0]
      
      # Fit the logistic regression model with selected features only
      lambda.model <- InterUplift(data, treat, outcome, names(features), 
                                  input = "best")
      
      # Compute the qini coefficient and add to the list
      if (validation == FALSE) {
        data$lambda.pred <- predict(lambda.model, data, treat)
        lambda.perf <- PerformanceUplift(data, treat, outcome, 
                                         "lambda.pred", 
                                         rank.precision = rank.precision, 
                                         equal.intervals = equal.intervals, 
                                         nb.group = nb.group)
        if (length(lambda.perf[[1]]) == 1) { lambda.qini[k] <- 0}
        else {lambda.qini[k] <- QiniArea(lambda.perf)}
        
      }
      if (validation == TRUE) {
        valid$lambda.pred <- predict(lambda.model, valid, treat)
        lambda.perf <- PerformanceUplift(valid, treat, outcome, 
                                         "lambda.pred", 
                                         rank.precision = rank.precision, 
                                         equal.intervals = equal.intervals, 
                                         nb.group = nb.group)
        if (length(lambda.perf[[1]]) == 1) { lambda.qini[k] <- 0}
        else {lambda.qini[k] <- QiniArea(lambda.perf)}
        
      }
    }
    
    best.model <- cbind(path[, c(1, 2)], lambda.qini)
    
    if (max(best.model[,3]) == 0) { 
      warning("All models result in a Qini coefficient equal to 0. Please check LassoPath().")
    }
    # Take the model that maximizes the qini coefficient
    max.qini <- which.max(best.model[,3])
    best.model <- best.model[max.qini,]
    
    # We also need to know which variables were selected
    best.features <- path[path[, "lambda"] == best.model["lambda"], -c(1, 2, 3)]
    best.features <- names(best.features[best.features != 0])
    class(best.features) <- "BestFeatures"
    return(best.features)
}

# END FUN