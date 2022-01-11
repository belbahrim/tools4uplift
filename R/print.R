######################################################################
# Print DualUplift
######################################################################

print.DualUplift <- function(x, ...) {
  
  # Print interface to DualUplift.
  
  res <- list (call         = x[[3]],
               control      = x[[1]], 
               treatment    = x[[2]])             
  
  res$control$call <- "Control Group Logistic Model"
  res$treatment$call <- "Treatment Group Logistic Model"
  
  class(res) <- "print.DualUplift"
  return(print(res))
}

# END FUN

######################################################################
# Print PerformanceUplift
######################################################################

print.PerformanceUplift <- function(x, ...) {
  
  # Print interface to PerformanceUplift

  res <- cbind("Targeted Population (%)" = x$cum_per, 
              "Incremental Uplift (%)" = x$inc_uplift,
              "Observed Uplift (%)" = x$uplift)  
  
  class(res) <- "print.PerformanceUplift"
  return(print(res))
}

# END FUN



######################################################################
# Print SplitUplift
######################################################################

print.SplitUplift <- function(x, ...) {
  
  # Print interface to SplitUplift
  cat("Call: SplitUplift\n")
  cat("Number of observations in training: ", nrow(x[[1]]), "\n",sep="")
  cat("Number of observations in validation: ", nrow(x[[2]]), "\n",sep="")

}

# END FUN



######################################################################
# Print qLHS
######################################################################

print.qLHS <- function(x, ...) {
  
  res <- x

  # Print interface to qLHS

  for (j in 1:length(res$modelLHS)){
    res$modelLHS[[j]]$call <- paste("LHS Model ", j)
    res$modelLHS[[j]]$aic <- NA
    res$modelLHS[[j]]$null.deviance <- NA
    res$modelLHS[[j]]$deviance <- NA
  }

  class(res) <- "print.qLHS"
  return(print(res))
}

# END FUN