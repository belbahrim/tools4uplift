interPredict <- function(data, treat, outcome, model, nb_group=10, plotit=FALSE){
  
  data1 <- data; data1[[treat]] <- 1
  pr.y1_ct1 <- predict(model, newdata=data1, type="response")
  
  data0 <- data; data0[[treat]] <- 0
  pr.y1_ct0 <- predict(model, newdata=data0, type="response")
  
  data$uplift_prediction <- pr.y1_ct1 - pr.y1_ct0
  
  perf <- qiniTable(data, nb_group, treat, outcome, "uplift_prediction")
  if (plotit == TRUE) {
    qiniCurve(perf, title = "Qini Curve for interaction model")
    qiniBarPlot(perf, title = "Uplift Bar Plot for interaction model")
  }
  qini <- qiniArea(perf)
  
  return(list(data, qini))
}