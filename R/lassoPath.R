lassoPath <- function(data, formula, nb_lambda=100){
  
  require(glmnet)

  X <- model.matrix(formula, data)
  y <- model.frame(formula,data)[,1]
  
  glmnet.output <- glmnet(X, y, alpha=1, family="binomial", nlambda = nb_lambda, intercept = FALSE)
  
  dimension <- glmnet.output$df
  coeff <- t(as.matrix(glmnet.output$beta))
  lambda <- glmnet.output$lambda
  
  path <- cbind(lambda, dimension, coeff)
  
  
  return(path)  
}