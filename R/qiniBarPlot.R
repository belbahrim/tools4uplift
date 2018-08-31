qiniBarPlot <-
function(dataResults, title="Model Performance: Uplift by Group")
{
  cBlue   <- rgb(097, 154, 188, 255, maxColorValue=255)
  
  barplot(dataResults[,8]*100, 
          names.arg = dataResults[,1]/nrow(dataResults)*1000, 
          col=cBlue,
          xlab = "Proportion of Population Targeted (%)",
          ylab="Uplift (%)",
          cex.main=0.85,
          main=title)
}
