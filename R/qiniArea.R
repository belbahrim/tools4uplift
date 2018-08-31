qiniArea <-
function(dataResults)
{
  nb <- nrow(dataResults)
  sum <- dataResults[1, 7]/2*dataResults[1, 3]/dataResults[nb, 3]
  for (i in 2:nb) {
    sum <- sum + (dataResults[i, 7] + dataResults[i-1, 7])/2*(dataResults[i, 3]/dataResults[nb, 3]-dataResults[i-1, 3]/dataResults[nb, 3])
  }
  
  return(sum-dataResults[nb, 7]/2)
}
