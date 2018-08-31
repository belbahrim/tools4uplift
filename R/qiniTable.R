qiniTable <-
function(data, nbGroup, treatment, outcome, prediction)
{
  data$rank <- 0
  data$rank = rank(-data[[prediction]], ties.method = "min") / nrow(data)
  
  data$group <- 0
  for(i in 1:nbGroup) {
      data$group[data$rank > (i-1)/nbGroup & data$rank <= i/nbGroup] <- i
  }
  
  dataResults <- data.frame(matrix(rep(0), nbGroup, 8))
  colnames(dataResults) <- c("cum_per", "T_Y1", "T_n", "C_Y1", 
                             "C_n", "incremental_Y1", "inc_uplift", "uplift")
  
  #incremental observed uplift
  for(i in 1:nbGroup){
    subset <- data[data$group <= i, ]
    dataResults[i,1] <- i/10
    dataResults[i,2] <- sum(subset[[treatment]] == 1 & subset[[outcome]] == 1)
    dataResults[i,3] <- sum(subset[[treatment]] == 1)
    dataResults[i,4] <- sum(subset[[treatment]] == 0 & subset[[outcome]] == 1)
    dataResults[i,5] <- sum(subset[[treatment]] == 0)
    dataResults[i,6] <- dataResults[i, 2] - dataResults[i, 4]*dataResults[i, 3]/dataResults[i, 5] 
  }
  dataResults[,7] <- dataResults[,6]/dataResults[nbGroup,3]*100
  
  
  #observed uplift in each group
  for (i in 1:nbGroup){
    subset <- data[data$group == i, ]
    dataResults[i,8] <-  sum(subset[[treatment]] == 1 & subset[[outcome]] == 1) / sum(subset[[treatment]] == 1) -  sum(subset[[treatment]] == 0 & subset[[outcome]] == 1) / sum(subset[[treatment]] == 0)
  }
  
  return(dataResults)
}
