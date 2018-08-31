squareUplift <-
function(train, valid, var1, var2, treat, outcome, step, 
         n.min=1, 
         plotit = TRUE, nbGroup = 10){
  
  require(ggplot2)
  require(gridExtra)
  
  
  #Initalize grid of uplifts for training set
  seq_col = seq(floor(min(train[,paste(var1)])),ceiling(max(train[,paste(var1)])),step)
  seq_row = seq(floor(min(train[,paste(var2)])),ceiling(max(train[,paste(var2)])),step)
  
  n = length(seq_row)-1
  m = length(seq_col)-1
  grid <- matrix(nrow=n, ncol=m)
  train$uplift <- NA
  
  
  for (i in n:1){
    for (j in 1:m){
      
      index <- (train[,paste(var1)] <=  seq_col[j] + step) & (train[,paste(var1)] >  seq_col[j]) & (train[,paste(var2)] <= seq_row[i] + step) & (train[,paste(var2)] > seq_row[i])
      square <- train[index,]
      
      indexT <- square[,paste(treat)]==1
      indexC <- square[,paste(treat)]==0
      sizeT <- sum(indexT)
      sizeC <- sum(indexC)
      
      if (sizeT >= n.min & sizeC >= n.min){
        grid[n-i+1, j] <- sum(square[,paste(treat)] == 1 & square[,paste(outcome)] == 1)/sum(square[,paste(treat)] == 1) - sum(square[,paste(treat)] == 0 & square[,paste(outcome)] == 1)/sum(square[,paste(treat)] == 0)
        train[index, "uplift"] <-  grid[n-i+1, j]
      }
      
    }
  }
  
  #for the observations that were not scored, i.e. uplift=NA, impute the average
  #uplift of all NA
  
  index.naT <- is.na(train$uplift)==TRUE
  train[index.naT,"uplift"] <- sum(train[index.naT,paste(treat)] == 1 & train[index.naT,paste(outcome)] == 1)/sum(train[index.naT,paste(treat)] == 1) - sum(train[index.naT,paste(treat)] == 0 & train[index.naT,paste(outcome)] == 1)/sum(train[index.naT,paste(treat)] == 0)
  
  
  #if length(valid) != 0
  #Initalize grid of uplifts for validation set
  seq_col = seq(floor(min(valid[,paste(var1)])),ceiling(max(valid[,paste(var1)])),step)
  seq_row = seq(floor(min(valid[,paste(var2)])),ceiling(max(valid[,paste(var2)])),step)
  
  n = length(seq_row)-1
  m = length(seq_col)-1
  grid <- matrix(nrow=n, ncol=m)
  valid$uplift <- NA
  
  
  for (i in n:1){
    for (j in 1:m){
      
      index <- (valid[,paste(var1)] <=  seq_col[j] + step) & (valid[,paste(var1)] >  seq_col[j]) & (valid[,paste(var2)] <= seq_row[i] + step) & (valid[,paste(var2)] > seq_row[i])
      square <- valid[index,]
      
      indexT <- square[,paste(treat)]==1
      indexC <- square[,paste(treat)]==0
      sizeT <- sum(indexT)
      sizeC <- sum(indexC)
      
      if (sizeT >= n.min & sizeC >= n.min){
        grid[n-i+1, j] <- sum(square[,paste(treat)] == 1 & square[,paste(outcome)] == 1)/sum(square[,paste(treat)] == 1) - sum(square[,paste(treat)] == 0 & square[,paste(outcome)] == 1)/sum(square[,paste(treat)] == 0)
        valid[index, "uplift"] <-  grid[n-i+1, j]
      }
      
    }
  }
  
  #for the observations that were not scored, i.e. uplift=NA, impute the average
  #uplift of all NA
  
  index.naT <- is.na(valid$uplift)==TRUE
  valid[index.naT,"uplift"] <- sum(valid[index.naT,paste(treat)] == 1 & valid[index.naT,paste(outcome)] == 1)/sum(valid[index.naT,paste(treat)] == 1) - sum(valid[index.naT,paste(treat)] == 0 & valid[index.naT,paste(outcome)] == 1)/sum(valid[index.naT,paste(treat)] == 0)
  
  
  if (plotit) {
    train.rank = rank(-train[["uplift"]], ties.method = "min") / nrow(train)
    valid.rank = rank(-valid[["uplift"]], ties.method = "min") / nrow(valid)
    
    for(i in 1:nbGroup){
      train$cluster[train.rank > (i-1)/nbGroup & train.rank <= i/nbGroup] <- i
      valid$cluster[valid.rank > (i-1)/nbGroup & valid.rank <= i/nbGroup] <- i
      }
  
    train$cluster <- as.factor(train$cluster)
    valid$cluster <- as.factor(valid$cluster)
    
    ptrain <- qplot(x=train[,var1], y=train[,var2], data=train, colour=cluster, 
                    main="Heat Map on Training Set", xlab=paste(var1), ylab=paste(var2)) +scale_fill_gradientn(colours = terrain.colors(10)) 
                    #+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                    #panel.background = element_blank(), axis.line = element_line(colour = "black"))
    pvalid <- qplot(x=valid[,var1], y=valid[,var2], data=valid, colour=cluster, 
                    main="Heat Map on Validation Set", xlab=paste(var1), ylab=paste(var2))+scale_fill_gradientn(colours = terrain.colors(10))
    grid.arrange(ptrain, pvalid, ncol=2)
  }
  
  squareUplift = list(train, valid)
  
  return(squareUplift)
}
