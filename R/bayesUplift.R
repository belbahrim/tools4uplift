bayesUplift <- function(data, outcome, treatment, predictors, folds=10, prior=TRUE){

  #require(glmnet)
  
  #Check if there are missings for predictors
  if_missings <- nrow(data) - sum(complete.cases(data[,which(names(data) %in% predictors)]))
  
  if (if_missings > 0){
    warning(paste("The dataset contains",if_missings,"observations with missing values! \n The fitted uplift is performed on complete cases only."))
    data <- data[complete.cases(data[,which(names(data) %in% predictors)]),]
  }

  
  #create the propensity score model in order to estimate P(Y=1|X)
  #using 10-fold cross validation
  #For the propensity, if prior=TRUE, keep only control group for this model
  #if prior=FALSE, use all observations to estimate P(Y=1|X)
  if (prior==TRUE) {mydata <- data[data[[treatment]] == 0,]}
  if (prior==FALSE) {mydata <- data}

  mydata <- mydata[,which(names(mydata) %in% c(predictors,outcome))] #keep only predictors and outcome
  

  formule <- as.formula(paste(outcome,"~ ."))
  X = model.matrix(formule, mydata)
  cv.fit <- cv.glmnet(x=X, y=mydata[,paste(outcome)], family="binomial", type.measure="auc", nfolds=folds)
  lambda <- cv.fit$lambda.1se
  
  
  #create the propensity score model in order to estimate P(T=1|Y=0)
  #using 10-fold cross validation
  mydata0 <- data[data[[outcome]] == 0,] #keep only Y=0
  mydata0 <- mydata0[,which(names(mydata0) %in% c(predictors,treatment))]  #keep only predictors and treatment
  
  formule0 <- as.formula(paste(treatment,"~ ."))
  X0 = model.matrix(formule0, mydata0)
  cv.fit0 <- cv.glmnet(x=X0, y=mydata0[,paste(treatment)],family="binomial", type.measure="auc", nfolds=folds)
  lambda0 <- cv.fit0$lambda.1se
  
  #create the propensity score model in order to estimate P(T=1|Y=1)
  #using 10-fold cross validation
  mydata1 <- data[data[[outcome]] == 1,] #keep only Y=1
  mydata1 <- mydata1[,which(names(mydata1) %in% c(predictors,treatment))]  #keep only predictors and treatment
  
  X1 = model.matrix(formule0, mydata1)
  cv.fit1 <- cv.glmnet(x=X1, y=mydata1[,paste(treatment)],family="binomial", type.measure="auc", nfolds=folds)
  lambda1 <- cv.fit1$lambda.1se
  
  
  
  #Estimate the uplift based on the 3 models
  X_data_propensity = model.matrix(formule, data[,which(names(data) %in% c(predictors,outcome))]) #keep only predictors and outcome
  X_data_conditional = model.matrix(formule0, data[,which(names(data) %in% c(predictors,treatment))]) #keep only predictors and treatment
  
  data$Propensity <- as.vector(predict(cv.fit, newx = X_data_propensity, s=lambda, type="response"))
  data$probability0 <- as.vector(predict(cv.fit0, newx = X_data_conditional, s=lambda0, type="response"))
  data$probability1 <- as.vector(predict(cv.fit1, newx = X_data_conditional, s=lambda1, type="response"))
  
  data$uplift <- 2*data$probability1*data$Propensity + 2*(1-data$probability0)*(1-data$Propensity) - 1
  
  
  return(list(data, outcome, treatment, predictors, cv.fit, cv.fit0, cv.fit1))
  
}

#END FUN

