generateUplift <- function(object, nbobs, dimension){
  
  function_test <- object[[1]]
  function_control <- object[[2]]
  space_dimension <- object[[3]] #needed in order to use the generateScene functions
  
  
  if (dimension < space_dimension) 
    stop(paste("The dimension of your synthetic data should be at least equal to the dimension from your scene:",space_dimension))
  
  #Now, we can use the previous functions to generate a dataframe with uplift value
  Data = data.frame(replicate(dimension, runif(nbobs,-1,1)))
  y <- c()
  treat <- c()
  pt <- c()
  pc <- c()
  
  
  #in order to use the generateScene functions
  #take only 1:space_dimension variables in the dataframe
  
  for (i in seq(1,nbobs/2,1)){
    pt[i] <- function_test(Data[i,1:space_dimension])
    pc[i] <- function_control(Data[i,1:space_dimension])
    y[i] = rbinom(n=1, size=1, prob=pt[i])
    treat[i] = 0
  }
  for (i in seq((nbobs/2)+1,nbobs,1)){
    pt[i] <- function_test(Data[i,1:space_dimension])
    pc[i] <- function_control(Data[i,1:space_dimension])
    y[i] = rbinom(n=1, size=1, prob=pc[i])
    treat[i] = 1
  }
  
  Sim_Data <- cbind(y,treat,Data,pt,pc)
  
  return(Sim_Data)
  
}

#END FUN
