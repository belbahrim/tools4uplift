InterPredict <- function(data, treat, outcome, model, nb.group=10, plotit=FALSE){
  
  data1 <- data; data1[[treat]] <- 1
  pr.y1_ct1 <- predict(model, newdata=data1, type="response")
  
  data0 <- data; data0[[treat]] <- 0
  pr.y1_ct0 <- predict(model, newdata=data0, type="response")
  
  data$uplift_prediction <- pr.y1_ct1 - pr.y1_ct0
  
  perf <- QiniTable(data, treat, outcome, "uplift_prediction", nb.group)
  if (plotit == TRUE) {
    QiniCurve(perf, title = "Qini Curve for interaction model")
    QiniBarPlot(perf, title = "Uplift Bar Plot for interaction model")
  }
  qini <- QiniArea(perf)
  
  return(list(data, qini))
}