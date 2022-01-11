InterUplift.formula <- function(formula, treat, data, ...){
  
  # Formula interface to InterUplift.
  if (!inherits(formula, "formula"))
    stop("Method is only for formula objects")
  
  mf <- match.call(expand.dots = FALSE)
  args <- match(c("formula", "data"), names(mf), 0)
  mf <- mf[c(1, args)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval.parent(mf)
  
  Terms <- attr(mf, "terms")
  Terms <- names(attr(Terms,"dataClasses")) 
  
  if (length(intersect(treat,colnames(data))) == 0)
    stop("InterUplift: data does not include the control/treatment variable treat).")    
  
  outcome <- Terms[1]
  predictors <- Terms[-1]
  fit <- InterUplift(data=data, treat=treat, outcome=outcome, predictors=predictors, input = "all", ...)

  cl <- match.call()
  cl[[1]] <- as.name("InterUplift")
  fit$call <- cl
  
  return(fit)
}

# END FUN


