bestFeatures <- function(training, formula, treat, nb_lambda=100, nb_group=10, value=FALSE){

  #Here, all variables must be continuous. Before using this function, change 
  #categorical variables to dummies.
  
  path <- lassoPath(training, formula, nb_lambda)
  #Keep paths of dimension > 0
  path <- path[path[,"dimension"] > 0, ]
  
  lambda_qini <- c()
  
  for (k in 1:nrow(path)) {
    features <- path[k, -c(1,2)]
    #Keep features with positive estimators only
    features <- features[features != 0]
    lambda_formula <- paste(paste(colnames(model.frame(formula,training))[1], "~ "),paste(names(features),collapse=" + "))
    
    #Fit the logistic regression model on selected features only with the training set
    lambda_model <- glm(lambda_formula, family=binomial(link="logit"), training)
    
    #Compute the qini coefficient and add to the list
    lambda_qini[k] <- scoreUplift(training, treat, colnames(model.frame(formula,training))[1], lambda_model, nb_group)
  }
  
  best_model <- cbind(path[,c(1,2)],lambda_qini)
  #Take the model that maximizes the qini coefficient
  best_model <- best_model[which.max(best_model[,"lambda_qini"]),]
  
  #We also need to know which variables were selected
  best_features <- path[path[,"lambda"]==best_model["lambda"], -c(1,2)]
  best_features <- names(best_features[best_features != 0])
  
  if (value == TRUE) {print(best_model)}
  
  return(best_features)   
}