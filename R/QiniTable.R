QiniTable <- function(data, treat, outcome, prediction, nb.group = 10){
  
  data$rank <- 0
  data$rank = rank(-data[[prediction]], ties.method = "min") / nrow(data)
  
  data$group <- 0
  for(i in 1:nb.group) {
      data$group[data$rank > (i-1)/nb.group & data$rank <= i/nb.group] <- i
  }
  
  dataResults <- data.frame(matrix(rep(0), nb.group, 8))
  colnames(dataResults) <- c("cum_per", "T_Y1", "T_n", "C_Y1", 
                             "C_n", "incremental_Y1", "inc_uplift", "uplift")
  
  #incremental observed uplift
  for(i in 1:nb.group){
    subset <- data[data$group <= i, ]
    dataResults[i,1] <- i/10
    dataResults[i,2] <- sum(subset[[treat]] == 1 & subset[[outcome]] == 1)
    dataResults[i,3] <- sum(subset[[treat]] == 1)
    dataResults[i,4] <- sum(subset[[treat]] == 0 & subset[[outcome]] == 1)
    dataResults[i,5] <- sum(subset[[treat]] == 0)
    dataResults[i,6] <- dataResults[i, 2] - dataResults[i, 4]*dataResults[i, 3]/dataResults[i, 5] 
  }
  dataResults[,7] <- dataResults[,6]/dataResults[nb.group,3]*100
  
  
  #observed uplift in each group
  for (i in 1:nb.group){
    subset <- data[data$group == i, ]
    dataResults[i,8] <- sum(subset[[treat]] == 1 & subset[[outcome]] == 1) / sum(subset[[treat]] == 1) -  
                        sum(subset[[treat]] == 0 & subset[[outcome]] == 1) / sum(subset[[treat]] == 0)
  }
  
  return(dataResults)
}
