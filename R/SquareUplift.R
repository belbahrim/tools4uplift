SquareUplift <- function(train, valid, var1, var2, treat, outcome, step, 
                         n.min=1, plotit = TRUE, col = heat.colors(100), nb.group = 10){
  
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
  
  #Plot the heat map for the training dataset
  if (plotit) {
    train.rank = rank(-train[["uplift"]], ties.method = "min") / nrow(train)
    
    for(i in 1:nb.group){
      train$cluster[train.rank > (i-1)/nb.group & train.rank <= i/nb.group] <- i
    }
    
    train$cluster <- as.factor(train$cluster)
    
    #reverse and transpose the grid in order to get the actual heat map with image function
    rev_grid <- apply(grid, 2, rev)
    image(seq_row, seq_col, t(rev_grid), xlab = paste(var1), ylab = paste(var2), col = col)
    title(main = "Heat Map for Training Set", font.main = 1)
  }
  
  
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
  
  #Plot the heat map for the validation dataset
  if (plotit) {
    valid.rank = rank(-valid[["uplift"]], ties.method = "min") / nrow(valid)
    
    for(i in 1:nb.group){
      valid$cluster[valid.rank > (i-1)/nb.group & valid.rank <= i/nb.group] <- i
    }
    
    valid$cluster <- as.factor(valid$cluster)
    
    #reverse and transpose the grid in order to get the actual heat map with image function
    rev_grid <- apply(grid, 2, rev)
    image(seq_row, seq_col, t(rev_grid), xlab = paste(var1), ylab = paste(var2), col = col)
    title(main = "Heat Map for Validation Set", font.main = 1)
  }
  
  squareUplift = list(train, valid)
  
  return(squareUplift)
}
