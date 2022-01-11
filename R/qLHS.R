qLHS <- function (data, treat, outcome, predictors, lhs_points = 50,
                  lhs_range = 1, adjusted = TRUE, 
                  rank.precision = 2, equal.intervals = FALSE, 
                  nb.group = 10, validation = TRUE, p = 0.3) {
  # Qini-based LHS Uplift Model.
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
  #   A Qini-based LHS optimal uplift interaction model
  
  
  # All variables must be continuous. Before using this function, change 
  # categorical variables to dummies.
  
  #data <- rearrange(data, treat, outcome, predictors)
  data <- standardize(data, treat, outcome)
    
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
  
  # Initialize a new path matrix to collect scores from LHS search
  pathLHS <- vector(mode = "list", length = nrow(path))
  
  for (k in 1:nrow(path)) {
 
    features <- path[k, -c(1, 2)]
    # Keep features with non zero estimates only
    features <- features[features != 0]
    
    # Fit the logistic regression model with selected features only
    lambda.model <- InterUplift(data, treat, outcome, names(features), 
                                input = "best")
    
    ########################TEMP######################################
    ##################################################################
    ##################################################################
    ##################################################################
    # We generate LHS samples from a uniform [-0.5; 0.5]
    betaLHS <- improvedLHS(lhs_points-1, length(features)+1)-0.5
    colnames(betaLHS) <- names(coef(lambda.model))
    # We want to sample points uniformly with mean MLE
    # and variance equal to lhs_range times variance of MLE
    for (l in 1:nrow(betaLHS)){
      betaLHS[l, ] <- sqrt(lhs_range*12)*summary(lambda.model)$coefficient[,2]*betaLHS[l, ] 
      betaLHS[l, ] <- betaLHS[l, ] + lambda.model$coefficients
    }
    
    # We include the MLE as the first estimator
    betaLHS <- cbind(rbind(lambda.model$coefficients, betaLHS), "score" = 0)

    for (m in 1:nrow(betaLHS)){
      lambda.model$coefficients <- betaLHS[m, -ncol(betaLHS)]
      if (validation == FALSE) {
        data$lambda.pred <- predict(lambda.model, data, treat)
        lambda.perf <- PerformanceUplift(data, treat, outcome, 
                                         "lambda.pred", 
                                         rank.precision = rank.precision, 
                                         equal.intervals = equal.intervals, 
                                         nb.group = nb.group)
      }
      if (validation == TRUE) {
        valid$lambda.pred <- predict(lambda.model, valid, treat)
        lambda.perf <- PerformanceUplift(valid, treat, outcome, 
                                         "lambda.pred", 
                                         rank.precision = rank.precision, 
                                         equal.intervals = equal.intervals, 
                                         nb.group = nb.group)
      }
      if (length(lambda.perf[[1]]) == 1) { betaLHS[m, ncol(betaLHS)] <- 0}
      else{betaLHS[m, ncol(betaLHS)] <- QiniArea(lambda.perf, adjusted)}
    }
    
    pathLHS[[k]] <- c(path[k, c(1,2)], betaLHS[which.max(betaLHS[,ncol(betaLHS)]), ])
    
  }

  regularization_matrix <- matrix(nrow=2, ncol=length(pathLHS))
  rownames(regularization_matrix) <- c("lambda", "score")
  bestLHS <- vector(mode = "list", length = length(pathLHS))
    
  for (j in 1:length(pathLHS)) {
    #Save Qini scores with associated regularization constant
    regularization_matrix[1, j] <- pathLHS[[j]][["lambda"]]
    regularization_matrix[2, j] <- pathLHS[[j]][["score"]]
      
    #Save the coefficients found with LHS into an InterUplift model
    #and return it for prediction purpose
    coefLHS <- pathLHS[[j]][-c(1,2,length(pathLHS[[j]]))]  
    predictorsLHS <- names(coefLHS)[-1]
    modelLHS <- InterUplift(data, treat, outcome, predictorsLHS, input = "best")
    #Change the coefficients
    modelLHS$coefficients <- coefLHS
    bestLHS[[j]] <- modelLHS
      
  }
  
  modelLHS <- bestLHS[[which.max(regularization_matrix[2,])]]
  return(modelLHS)
  
}

# END FUN