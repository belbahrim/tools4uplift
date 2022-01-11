######################################################################
# Plot Qini Curve
######################################################################

plot.PerformanceUplift <- function(x, ...){
  
  # Plots the Qini curve.
  #
  # Args:
  #   x: a table that must be the output of PerformanceUplift() function.
  #
  # Returns:
  #   The Qini curve.
  
  if (!inherits(x, "PerformanceUplift"))
    stop("tools4uplift: object not of class PerformanceUplift")
  
  cGrey   <- rgb(128, 128, 131, 255, maxColorValue = 255)
  
  nb <- length(x$cum_per)
  
  plot(c(0, x$T_n)/x$T_n[nb]*100, c(0, x$inc_uplift),  
       xlab='Proportion of Population Targeted (%)', 
       ylab='Incremental Uplift (%)', 
       ...)
  lines(c(0,100), c(0, x$inc_uplift[nb]), col=cGrey, type="b", lwd=1.5)
  return(QiniArea(x))
}



lines.PerformanceUplift <- function(x, ...){
  
  # Adds an extra Qini curve on an existing Qini curves plot.
  #
  # Args:
  #   x: a table that must be the output of PerformanceUplift() function.
  #
  # Returns:
  #   The Qini curve.
  
  if (!inherits(x, "PerformanceUplift"))
    stop("tools4uplift: object not of class PerformanceUplift")
  

  nb <- length(x$cum_per)
  lines(c(0, x$T_n)/x$T_n[nb]*100, c(0, x$inc_uplift), ...)
  return(QiniArea(x))
}


# END FUN


######################################################################
# Plot Uplift Barplot
######################################################################

barplot.PerformanceUplift <- function(height, ...){
  
  # Plots the Qini barplot.
  #
  # Args:
  #   height: a table that must be the output of PerformanceUplift() function.
  #
  # Returns:
  #   The Qini barplot.
  
  x <- height
  
  if (!inherits(x, "PerformanceUplift"))
    stop("tools4uplift: object not of class PerformanceUplift")
  
  barplot(x$uplift, 
          names.arg = round(x$cum_per*100),
          xlab='Proportion of Population Targeted (%)', 
          ylab='Uplift (%)',
          ...)
  abline(h = 0)
  
  #Compute the Kendall's uplift rank correlation
  #Raise a warning if there are missing values in x
  
  x_for_rho <- cbind(x$cum_per, x$uplift)
  complete_x_for_rho <- x_for_rho[complete.cases(x_for_rho),]
  if (nrow(x_for_rho) > nrow(complete_x_for_rho)){
    warning("tools4uplift: there are missing values in your PerformanceUplift object. They will be omitted in the computation of the Kendall's rank correlation.")
  }
  
  kendall_uplift <- cor(seq(nrow(complete_x_for_rho),1), complete_x_for_rho[,2], method="kendall")
  
  return(kendall_uplift)
  
}

# END FUN
