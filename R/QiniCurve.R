QiniCurve <- function(x, title = "Model Performance: Qini Curve", color = NULL){
  
  
  if (is.null(color)==TRUE) {
    color <- rgb(097, 154, 188, 255, maxColorValue = 255)
  }
  
  cGrey   <- rgb(128, 128, 131, 255, maxColorValue = 255)
  
  plot(c(0, x[, 3])/x[nrow(x), 3]*100, c(0, x[, 7]),
       type = "l",
       lwd = 2.5,
       col = color,
       xlab = "Proportion of Population Targeted (%)", 
       ylab = "Incremental Uplift (%)",
       cex.main = 0.85,
       main = title
  )
  abline(a = 0, b = x[nrow(x), 7]/100, col = cGrey, lwd = 2.0)
}
