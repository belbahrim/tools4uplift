######################################################################
# Functions necessary for other tools4uplift functions
######################################################################

## logistic function
logistic <- function(x){
  return(1/(1+exp(-x)))
}


## formula for interaction model
formulaUplift <- function(treat, outcome, predictors){
  inter.formula <- c()
  for (k in seq(1:length(predictors))) {
    inter.formula <- paste(inter.formula, 
                           paste(predictors[k], treat, sep = ":"), sep="+")
  }
  formula <- as.formula(paste(paste(outcome, "~", treat, "+"), 
                              paste(predictors, collapse="+"), inter.formula))
  return(formula)
}


## rearrange data in order to put outcome and treatment as the first columns
rearrange <- function(data, treat, outcome, predictors){
  
  covariates <- data[,predictors]
  treat <- data[,treat]
  outcome <- data[,outcome]
  df <- data.frame(cbind(outcome, treat, covariates))
  
  return(df)
}


## compute uplift for a given vector of coefficients (Need to rearrange data first)
computeUplift <- function(data, beta){
  
  z = (ncol(data)-2) # number of predictors - the two first columns being Y and T
  
  pt1 <- logistic(beta[[1]] + beta [[2]] 
                  + as.matrix(data[,3:ncol(data)]) %*% t(as.matrix(beta[3:ncol(data)]))
                  + as.matrix(data[,3:ncol(data)]) %*% t(as.matrix(beta[(ncol(data)+1):(ncol(data)+z)])))
  
  pc1 <- logistic(beta[[1]] + as.matrix(data[,3:ncol(data)]) %*% t(as.matrix(beta[3:ncol(data)])))
  
  uplift <-pt1 - pc1
  
  return(uplift)
}


## standardize the data (Need to rearrange data first)
standardize <- function(data, treat, outcome){
  
  df <- cbind(data[,c(outcome,treat)],scale(data[,!(names(data) %in% c(outcome,treat))]))
  
  return(df)
}