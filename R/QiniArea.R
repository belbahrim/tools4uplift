QiniArea <- function(x){
  
  nb <- nrow(x)
  sum <- x[1, 7]/2*x[1, 3]/x[nb, 3]
  for (i in 2:nb) {
    sum <- sum + (x[i, 7] + x[i-1, 7])/2*(x[i, 3]/x[nb, 3]-x[i-1, 3]/x[nb, 3])
  }
  
  return(sum-x[nb, 7]/2)
}
