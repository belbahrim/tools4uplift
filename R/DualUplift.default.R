DualUplift <- function(data, ...) UseMethod("DualUplift")

DualUplift.default <- function(data, treat, outcome, predictors, ...){
  
  # Two-model estimator.
  #
  # Args:
  #   data: a data frame containing the treatment, the outcome and the predictors.
  #   treat: name of a binary (numeric) vector representing the treatment 
  #          assignment (coded as 0/1).
  #   outcome: name of a binary response (numeric) vector (coded as 0/1).
  #   predictors: a vector of names representing the predictors to consider in the model.
  #
  # Returns:
  #   Fitted model for control group and fitted model for treatment group.
  
  # Error handling
  # Check if there are missings for predictors
  if_missings <- nrow(data) - sum(complete.cases(data[,which(names(data) %in% predictors)]))
  
  if (if_missings > 0){
    warning(paste("The dataset contains",if_missings,"observations with missing values! \n The fitted uplift is performed on complete cases only."))
    data <- data[complete.cases(data[,which(names(data) %in% predictors)]),]
  }
  
  # Define formula for model fit
  model_formula <- as.formula(paste(outcome,"~ ."))
  
  # Create the propensity score model in order to estimate P(Y=1|T=0)
  mydata0 <- data[data[[treat]] == 0,]  # keep only T=0
  mydata0 <- mydata0[,which(names(mydata0) %in% c(predictors,outcome))]  # keep only predictors and outcome
  
  model0 <- glm(model_formula, family=binomial(link="logit"), mydata0, ...)
  
  # Create the propensity score model in order to estimate P(Y=1|T=1)
  mydata1 <- data[data[[treat]] == 1,]  # keep only T=1
  mydata1 <- mydata1[,which(names(mydata1) %in% c(predictors,outcome))]  # keep only predictors and outcome
  
  model1 <- glm(model_formula, family=binomial(link="logit"), mydata1, ...)  
  
  res.dual <- list(model0, model1)
  
  cl <- match.call()
  cl[[1]] <- as.name("DualUplift")
  res.dual[[3]] <- cl
  
  class(res.dual) <- "DualUplift"
  
  return(res.dual)
  
}

# END FUN