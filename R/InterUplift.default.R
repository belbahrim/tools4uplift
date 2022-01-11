InterUplift <- function(data, ...) UseMethod("InterUplift")

InterUplift.default <- function(data, treat, outcome, predictors, input = "all", ...){
  
  # Interaction estimator.
  #
  # Args:
  #   data: a data frame containing the treatment, the outcome and the predictors.
  #   treat: name of a binary (numeric) vector representing the treatment 
  #          assignment (coded as 0/1).
  #   outcome: name of a binary response (numeric) vector (coded as 0/1).
  #   predictors: a vector of names representing the predictors to consider in the model.
  #   input: controls if we want to use all the predictors ("all")  
  #          or the output of BestFeatures() function as predictors ("best") 
  #
  # Returns:
  #   Fitted interaction estimator model.  
  
  
  # Error handling
  if (input != "all" & input != "best"){
    stop("Choose an appropriate method for input in 'all' or 'best'")
  }
  
  if (input == "all") {
      
    inter_formula <- c()
    for (k in seq(1:length(predictors))) {
      inter_formula <- paste(inter_formula, paste(predictors[k], treat, sep = ":"), sep="+")
    }
    model_formula <- as.formula(paste(paste(outcome, "~", treat, "+"),paste(predictors,collapse="+"),inter_formula))
    model <- glm(model_formula, family=binomial(link="logit"), data, ...)
  }
  
  if (input == "best") {
    model_formula <- as.formula(paste(paste(outcome, "~"),paste(predictors,collapse = "+")))
    model <- glm(model_formula, family=binomial(link="logit"), data, ...)
  }
  
  cl <- match.call()
  cl[[1]] <- as.name("InterUplift")
  model$call <- cl
  class(model) <- c("InterUplift", "glm", "lm")
  return(model)
}

# END FUN