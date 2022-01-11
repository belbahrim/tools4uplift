QiniArea <- function(x, ...)  UseMethod("QiniArea")

QiniArea.default <- function(x, ...)
  stop("No method implemented for this class of object")

QiniArea.PerformanceUplift <- function(x, adjusted=FALSE, ...){
  
  # Computes the area under the Qini curve.
  #
  # Args:
  #   x: a table that must be the output of PerformanceUplift() function.
  #
  # Returns:
  #   The Qini coefficient.
  
  if (!inherits(x, "PerformanceUplift"))
    stop("tools4uplift: object not of class PerformanceUplift")
  
  nb <- length(x$cum_per)
  sum <- x$inc_uplift[1]/2*x$T_n[1]/x$T_n[nb]
  for (i in 2:nb) {
    sum <- sum + (x$inc_uplift[i] + x$inc_uplift[i-1])/2*(x$T_n[i]/x$T_n[nb]-x$T_n[i-1]/x$T_n[nb])
  }
  
  qini_area <- sum-x$inc_uplift[nb]/2
  
  if (adjusted == TRUE){
    #Compute the Kendall's uplift rank correlation
    #Raise a warning if there are missing values in x
    
    x_for_rho <- cbind(x$cum_per, x$uplift)
    complete_x_for_rho <- x_for_rho[complete.cases(x_for_rho),]
    if (nrow(x_for_rho) > nrow(complete_x_for_rho)){
      warning("tools4uplift: there are missing values in your PerformanceUplift object. They will be omitted in the computation of the Kendall's rank correlation.")
    }
    
    kendall_uplift <- cor(seq(nrow(complete_x_for_rho),1), complete_x_for_rho[,2], method="kendall")
    
    
    qini_area <- max(0, qini_area) * kendall_uplift
  }
  
  return(qini_area)
}

# END FUN