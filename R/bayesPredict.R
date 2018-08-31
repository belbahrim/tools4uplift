bayesPredict <- function(object, newdata){
  
  #Extract information from object bayesUplift
  outcome <- object[[2]]
  treatment <- object[[3]]
  predictors <- object[[4]]
  
  #Check whether or not the outcome, treatment and predictors are included in
  #the newdata
  if (all(predictors %in% colnames(newdata))==FALSE) 
    stop("variables in the training data missing in newdata")
  
  #Check if there are missings for predictors in newdata
  if_missings <- nrow(newdata) - sum(complete.cases(newdata[,which(names(newdata) %in% predictors)]))
  
  if (if_missings > 0){
    warning(paste("The dataset contains",if_missings,"observations with missing values! \n The predicted uplift is performed on complete cases only."))
    newdata <- newdata[complete.cases(newdata[,which(names(newdata) %in% predictors)]),]
  }
  
  
  #Prepare the model formula in order to compute the predictions
  formule <- as.formula(paste(outcome,"~ ."))
  formule0 <- as.formula(paste(treatment,"~ ."))
  
  #Predict the uplift based on the 3 models
  X_data_propensity = model.matrix(formule, newdata[,which(names(newdata) %in% c(predictors,outcome))]) #keep only predictors and outcome
  X_data_conditional = model.matrix(formule0, newdata[,which(names(newdata) %in% c(predictors,treatment))]) #keep only predictors and treatment
  
  newdata$Propensity <- as.vector(predict(object[[5]], newx = X_data_propensity, s=object[[5]]$lambda.1se, type="response"))
  newdata$probability0 <- as.vector(predict(object[[6]], newx = X_data_conditional, s=object[[6]]$lambda.1se, type="response"))
  newdata$probability1 <- as.vector(predict(object[[7]], newx = X_data_conditional, s=object[[7]]$lambda.1se, type="response"))
  
  newdata$uplift <- 2*newdata$probability1*newdata$Propensity + 2*(1-newdata$probability0)*(1-newdata$Propensity) - 1
  
  
  return(newdata)
}

#END FUN