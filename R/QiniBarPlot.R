QiniBarPlot <- function(x, title = "Model Performance: Uplift by Group", color = NULL){
  
  if (is.null(color)==TRUE) {
    color <- rgb(097, 154, 188, 255, maxColorValue = 255)
  }
  
  barplot(x[,8]*100, 
          names.arg = x[,1]/nrow(x)*1000, 
          col = color,
          xlab = "Proportion of Population Targeted (%)",
          ylab = "Uplift (%)",
          cex.main = 0.85,
          main = title)
}
