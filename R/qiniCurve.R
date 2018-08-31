qiniCurve <-
function(dataResults, title="Model Performance: Qini Curve")
{
  
  cBlue   <- rgb(097, 154, 188, 255, maxColorValue=255)
  cGrey   <- rgb(128, 128, 131, 255, maxColorValue=255)
  
  plot(c(0, dataResults[,3])/dataResults[nrow(dataResults),3]*100,c(0, dataResults[,7]),
       type="l",
       lwd=2.5,
       col=cBlue,
       xlab="Proportion of Population Targeted (%)", 
       ylab="Incremental Uplift (%)",
       cex.main=0.85,
       main=title
  )
  abline(a=0,b=dataResults[nrow(dataResults),7]/100, col=cGrey,lwd=2.0)
}
