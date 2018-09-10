dualUplift <- function(data, treat, outcome, predictors){
  
  #Check if there are missings for predictors
  if_missings <- nrow(data) - sum(complete.cases(data[,which(names(data) %in% predictors)]))
  
  if (if_missings > 0){
    warning(paste("The dataset contains",if_missings,"observations with missing values! \n The fitted uplift is performed on complete cases only."))
    data <- data[complete.cases(data[,which(names(data) %in% predictors)]),]
  }
  
  #Define formula for model fit
  model_formula <- as.formula(paste(outcome,"~ ."))
  
  #create the propensity score model in order to estimate P(Y=1|T=0)
  mydata0 <- data[data[[treat]] == 0,] #keep only T=0
  mydata0 <- mydata0[,which(names(mydata0) %in% c(predictors,outcome))]  #keep only predictors and outcome
  
  model0 <- glm(model_formula, family=binomial(link="logit"), mydata0)
  
  #create the propensity score model in order to estimate P(Y=1|T=1)
  #using 10-fold cross validation
  mydata1 <- data[data[[treat]] == 1,] #keep only T=1
  mydata1 <- mydata1[,which(names(mydata1) %in% c(predictors,outcome))]  #keep only predictors and outcome
  
  model1 <- glm(model_formula, family=binomial(link="logit"), mydata1)  
  
  return(list(model0, model1))
  
}

#END FUN

