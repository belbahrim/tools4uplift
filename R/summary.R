summary.DualUplift <- function(object, ...) {
  
  # Summary interface to DualUplift.
  
  res <- list (call         = object[[3]],
               control      = summary.glm(object[[1]]), 
               treatment    = summary.glm(object[[2]]))             
  
  res$control$call <- "Control Group Logistic Model"
  res$treatment$call <- "Treatment Group Logistic Model"
  
  class(res) <- "summary.DualUplift"
  return(res)
}

# END FUN
