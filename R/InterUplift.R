InterUplift <- function(data, treat, outcome, predictors, input=c("p", "bf")){

    #The parameter input controls if we want to use the predictors ("p") or the output of
    #bestFeatures function as predictors ("bf")
    if (input != "p" & input != "bf") stop("Choose an appropriate method for input in 'p' or 'bf'")
    
    if (input == "p") {
      
      inter_formula <- c()
      for (k in seq(1:length(predictors))) {
        inter_formula <- paste(inter_formula, paste(predictors[k], treat, sep = ":"), sep="+")
      }
      model_formula <- as.formula(paste(paste(outcome, "~", treat, "+"),paste(predictors,collapse="+"),inter_formula))
      
      model <- glm(model_formula, family=binomial(link="logit"), data)
      
    }
  
    if (input == "bf") {
      model_formula <- as.formula(paste(paste(outcome, "~"),paste(predictors,collapse = "+")))
      
      model <- glm(model_formula, family=binomial(link="logit"), data)
    }
  
  return(model)
}