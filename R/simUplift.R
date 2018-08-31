simUplift <-
function(n, p, c0, c1, sigma0, sigma1){
  
  #require(optimization)
  #require(MASS)
  
  #Use a maximum of 5 dimensions to generate the probability "distributions"
  d=min(p,5)
  
  #Generate the probability "distributions" from 20 points
  X = data.frame(replicate(d, runif(20,-1,1)))
  
  px <- function(x){
    px = 0
    for (i in 1:nrow(X)){
      px <- px + c0*(2*pi*sigma0^2)^(-d/2)*exp(-0.5*norm(x-X[i,], type="2")/sigma0)
    }
    return(px)
  }
  
  max_px <- optim_nm(px, k=d, maximum = TRUE)$function_value
  
  cx <- function(x){
    px = 0
    vx = 0
    cx = 0
    
    for (i in 1:nrow(X)){
      px <- px + c0*(2*pi*sigma0^2)^(-d/2)*exp(-0.5*norm(x-X[i,], type="2")/sigma0)
      if (i <= nrow(X)/2) vx <- vx - c1*(2*pi*sigma1^2)^(-d/2)*exp(-0.5*norm(x-X[i,], type="2")/sigma1)
      else if (i > nrow(X)/2)  vx <- vx + c1*(2*pi*sigma1^2)^(-d/2)*exp(-0.5*norm(x-X[i,], type="2")/sigma1)
      cx <- px + vx
    }
    return(cx)
  }
  
  max_cx <- optim_nm(cx, k=d, maximum = TRUE)$function_value
  
  #Now, we can use the previous functions to generate a dataframe with uplift value
  Data = data.frame(replicate(p, runif(n,-1,1)))
  y <- c()
  treat <- c()
  p0 <- c()
  p1 <- c()
  for (i in seq(1,n/2,1)){
    p0[i] <- min(max(0,px(Data[i,1:d]))/max(max_px,1),1)
    p1[i] <- min(max(0,cx(Data[i,1:d]))/max(max_cx,1),1)
    y[i] = rbinom(n=1, size=1, prob=p0[i])
    treat[i] = 0
  }
  for (i in seq((n/2)+1,n,1)){
    p0[i] <- min(max(0,px(Data[i,1:d]))/max(max_px,1),1)
    p1[i] <- min(max(0,cx(Data[i,1:d]))/max(max_cx,1),1)
    y[i] = rbinom(n=1, size=1, prob=p1[i])
    treat[i] = 1
  }
  
  ts <- p1-p0
  Sim_Data <- cbind(y,treat,Data,ts)
  
  return(Sim_Data)
  
}
