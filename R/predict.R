######################################################################
# Predict DualUplift
######################################################################

predict.DualUplift <- function(object, newdata, ...) {
  
  if (!inherits(object, "DualUplift"))
    stop("tools4uplift: object not of class DualUplift")
  if (nrow(newdata) == 0)
    stop("tools4uplift: newdata has 0 rows")
  
  # Extract information from model
  model0 <- object[[1]]
  model1 <- object[[2]]
  
  if (any(!names(model0$coefficients) %in% c("(Intercept)",colnames(newdata))))
    stop("tools4uplift: variables in the training data missing in newdata")
  
  # Predict the 2 probabilities from model0 and model1
  pr.y1_ct1 <- predict.glm(model1, newdata, type = "response", ...)
  pr.y1_ct0 <- predict.glm(model0, newdata, type = "response", ...)
  
  all.res <- pr.y1_ct1 - pr.y1_ct0
  
  return(all.res)
}

# END FUN


######################################################################
# Predict InterUplift
######################################################################

predict.InterUplift <- function(object, newdata, treat, ...) {
  
  if (!inherits(object, "InterUplift"))
    stop("tools4uplift: object not of class InterUplift")
  if (nrow(newdata) == 0)
    stop("tools4uplift: newdata has 0 rows")
  if (is.null(treat) == TRUE)
    stop("tools4uplift: please specify the treatment variable name")
  
  
  #use computeUplift function instead !
  
  data1 <- newdata; data1[treat] <- 1
  pr.y1_ct1 <- predict.glm(object, newdata=data1, type="response", ...)
  
  data0 <- newdata; data0[treat] <- 0
  pr.y1_ct0 <- predict.glm(object, newdata=data0, type="response", ...)
  
  all.res <- pr.y1_ct1 - pr.y1_ct0
  
  return(all.res)
}

# END FUN


######################################################################
# Predict BinUplift
######################################################################

predict.BinUplift <- function(object, newdata, ...){
  
  if (!inherits(object, "BinUplift"))
    stop("tools4uplift: object not of class BinUplift")
  
  if (length(object) == 1)
    stop("tools4uplift: make sure that the variable is quantized")
  
  # Extract information from BinUplift object
  cutoffs <- object$x.cut
  min_var <- min(newdata)
  max_var <- max(newdata)
  
  if (max_var %in% cutoffs) {
    knots <- unique(c(min_var, sort(cutoffs), Inf))
  } else {
    knots <- unique(c(min_var, sort(cutoffs), max_var))
  }
  
  cat_var <- cut(x=newdata, breaks=knots, include.lowest = TRUE, right=FALSE, ...)
  return(cat_var)
}

# END FUN