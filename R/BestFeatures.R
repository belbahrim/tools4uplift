BestFeatures <- function(data, treat, outcome, predictors, nb.lambda = 100, 
                         nb.group = 10, value = FALSE){

  #Here, all variables must be continuous. Before using this function, change 
  #categorical variables to dummies.
  inter.formula <- c()
  for (k in seq(1:length(predictors))) {
    inter.formula <- paste(inter.formula, paste(predictors[k], treat, sep = ":"), sep="+")
  }
  formula <- as.formula(paste(paste(outcome, "~", treat, "+"), paste(predictors, collapse="+"), inter.formula))
  
  path <- LassoPath(data, formula, nb.lambda)
  #Keep paths of dimension > 0
  path <- path[path[, "dimension"] > 0, ]
  
  lambda.qini <- c()
  
  for (k in 1:nrow(path)) {
    features <- path[k, -c(1, 2)]
    #Keep features with positive estimators only
    features <- features[features != 0]
    lambda.formula <- paste(paste(colnames(model.frame(formula, data))[1], "~ "), paste(names(features), collapse=" + "))
    
    #Fit the logistic regression model on selected features only
    lambda.model <- glm(lambda.formula, family=binomial(link="logit"), data)
    
    #Compute the qini coefficient and add to the list
    lambda.qini[k] <- InterPredict(data, treat, colnames(model.frame(formula, data))[1], lambda.model, nb.group)[[2]]
  }
  
  best.model <- cbind(path[, c(1, 2)], lambda.qini)
  #Take the model that maximizes the qini coefficient
  best.model <- best.model[which.max(best.model[, "lambda.qini"]), ]
  
  #We also need to know which variables were selected
  best.features <- path[path[, "lambda"]==best.model["lambda"], -c(1, 2)]
  best.features <- names(best.features[best.features != 0])
  
  if (value == TRUE) {print(best.model)}
  
  return(best.features)   
}