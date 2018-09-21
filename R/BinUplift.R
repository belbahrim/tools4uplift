BinUplift <- function(data, treat, outcome, x, n.split = 10, alpha = 0.05, 
                       n.min = 30, ylim = NULL, ylab = "Uplift",
                       title = "Binning Results", color = NULL){
  
  if (is.null(color)==TRUE) {
    color <- rgb(097, 154, 188, 255, maxColorValue = 255)
  }

  #First, we need to check that the parameters are well filled and we create
  #a data.frame that contains only the variables we need
  
  if (n.split <= 0 ) {
    stop("Number of splits must be positive")
  }
  if (alpha < 0 | alpha > 1 ) {
    stop("alpha must be between 0 and 1")
  }
  
  if (sum(is.na(data[[outcome]])) > 0 ) {
    stop("Dependent variable contains missing values: remove the observations and proceed")
  }
  if (sum(is.na(data[[treat]])) > 0 ) {
    stop("Treatment variable contains missing values: remove the observations and proceed")
  }
  if (sum(is.na(data[[x]])) > 0 ) {
    warning("Independent variable contains missing values: remove the observations and proceed")
  }
  
  #Keep complete cases only for the variable to categorize
  df <- data[is.na(data[[x]])==FALSE,]
  
  data <- data.frame(df[[outcome]], df[[treat]], df[[x]])
  colnames(data) <- c("y", "treat", "x")
  
  #Second, we need to call the subfunctions that create the splits and
  #choose the best split based on z-test
  
  BinUpliftStump <- function(data, y, treat, x, n.split){
    
    y <- deparse(substitute(y))
    treat <- deparse(substitute(treat))
    x <- deparse(substitute(x))
    
    x.cut <- numeric(length = n.split)
    x.step <- (max(data[[x]]) - min(data[[x]]))/n.split
    
    for(i in 1:length(x.cut)){
      x.cut[i] = min(data[[x]]) + i*x.step
    }
    
    splits <- matrix(data = NA, nrow = length(x.cut), ncol = 12)
    colnames(splits) <- c("x.cut", "n.lt", "n.lc", "p.lt", "p.lc", "u.l",
                          "n.rt", "n.rc", "p.rt", "p.rc", "u.r", "diff")

    for(i in 1:length(x.cut)){ 
      
      index.l <- data[[x]] < x.cut[i]
      index.r <- data[[x]] >= x.cut[i]
      splits[i,1] <- x.cut[i]
      
      left <- data[index.l,]
      right <- data[index.r,]
      
      #uplift in left node
      splits[i, 2] <- sum(left[[treat]]==1) #n.lt
      splits[i, 3] <- sum(left[[treat]]==0) #n.lc
      splits[i, 4] <- sum(left[[treat]]==1 & left[[y]]==1)/splits[i, 2] #p.lt
      splits[i, 5] <- sum(left[[treat]]==0 & left[[y]]==1)/splits[i, 3] #p.lc
      
      splits[i, 6] <- splits[i, 4] - splits[i, 5] #u.l
      
      #uplift in right node
      splits[i, 7] <- sum(right[[treat]]==1) #n.rt
      splits[i, 8] <- sum(right[[treat]]==0) #n.rc
      splits[i, 9] <- sum(right[[treat]]==1 & right[[y]]==1)/splits[i, 7] #p.rt
      splits[i, 10] <- sum(right[[treat]]==0 & right[[y]]==1)/splits[i, 8] #p.rc
      
      splits[i, 11] <- splits[i, 9] - splits[i, 10] #u.r
      
      splits[i, 12] <- abs(splits[i, 6] - splits[i, 11]) #diff
      
    }
    
    return(splits)
    
  }
  
  
  BinUpliftTest <- function(splits, alpha, n.min){
    
    z.a <- qnorm(0.5*alpha, mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)
    
    test.splits <- data.frame(splits)
    test.splits$z.num <- (test.splits$p.lt - test.splits$p.lc - test.splits$p.rt + test.splits$p.rc)
    test.splits$z.den <- sqrt((test.splits$p.lt*(1-test.splits$p.lt))/test.splits$n.lt +
                                (test.splits$p.lc*(1-test.splits$p.lc))/test.splits$n.lc +
                                (test.splits$p.rt*(1-test.splits$p.rt))/test.splits$n.rt +
                                (test.splits$p.rc*(1-test.splits$p.rc))/test.splits$n.rc)
    test.splits$z.obs <- abs(test.splits$z.num / test.splits$z.den)
    
    test.splits$sign <- 1*(test.splits$z.obs > z.a)
    
    argmax.size <- test.splits[test.splits$n.lt > n.min & test.splits$n.lc > n.min & test.splits$n.rt > n.min & test.splits$n.rc > n.min,]
    argmax.cut <- argmax.size[which.max(argmax.size$z.obs),]
    
    argmax.cut$x.pos <- which(test.splits$x.cut == argmax.cut$x.cut)	
    
    return(argmax.cut)
    
  }
  
  
  #Now we can call the recursive partitionning subfunction that performs the
  #splits based on the two (2) previous subfunctions
  
  BinUpliftTree <- function(data, y, treat, x, n.split, alpha, n.min){   
    
    y <- deparse(substitute(y))
    treat <- deparse(substitute(treat))
    x <- deparse(substitute(x))
    
    stump <- BinUpliftStump(data, y, treat, x, n.split)
    best <- BinUpliftTest(stump, alpha, n.min)
    
    if (best$sign == 0 || nrow(best) == 0) { 
      
      return("oups..no significant split") 
      
    } else if (best$sign == 1) {
        
      print("The variable has been cut at:")
      print(best$x.cut) #printing the cuts in order
      
      #the x.cuts at the left and at the right of the x.cut optimal 
      l.cuts <- best$x.pos - 1
      r.cuts <- n.split-best$x.pos + 1
      
      l.index <- data[[x]] < best$x.cut
      r.index <- data[[x]] >= best$x.cut
      
      l.create <- data[l.index,]
      r.create <- data[r.index,]
      
      l.stump <- BinUpliftStump(l.create, y, treat, x, l.cuts)
      r.stump <- BinUpliftStump(r.create, y, treat, x, r.cuts)
      
      l.best <- BinUpliftTest(l.stump, alpha, n.min)
      r.best <- BinUpliftTest(r.stump, alpha, n.min)
      
      #possible paths 
      if (l.best$sign == 0 || nrow(l.best) == 0) {
        if (r.best$sign == 0 || nrow(r.best) == 0) {
          return (best)
        } else if (r.best$sign == 1) {
          return (rbind(best, BinUpliftTree(r.create, y, treat, x, r.cuts, alpha, n.min))) 
        }
      } else if (l.best$sign == 1) {
        if (r.best$sign == 0 || nrow(r.best) == 0){
          return (rbind(BinUpliftTree(l.create, y, treat, x, l.cuts, alpha, n.min), best)) 
        } else if (r.best$sign == 1) {
          return (rbind(BinUpliftTree(l.create, y, treat, x, l.cuts, alpha, n.min), best, BinUpliftTree(r.create, y, treat, x, r.cuts, alpha, n.min))) 
        }
      }    
      }
  }
  
  #Now, if the variable we want to bin is not continuous, we need to change 
  #it to a continuous one to be able to split
  #We sort the categories by uplift and create an ordinal variable 1,2,...,K
  #where K is the number of the categories. Then, we can split that variable
  #with n.split=K-1
  
  BinUpliftCatRank <- function(data, y, treat, x){
    
    y <- deparse(substitute(y))
    treat <- deparse(substitute(treat))
    x <- deparse(substitute(x))
    
    splits <- matrix(data = NA, nrow = length(levels(data[[x]])), ncol=6)
    colnames(splits) <- c("cat", "n.t", "n.c", "p.t", "p.c", "u")
    splits <- as.data.frame(splits)
    
    for(i in 1:nrow(splits)){ 
      
      splits[i, 1] <- levels(data[[x]])[[i]]
      
      #uplift in each category
      splits[i, 2] <- sum(data[[treat]]==1 & data[[x]]== splits[i, 1]) #n.t
      splits[i, 3] <- sum(data[[treat]]==0 & data[[x]]== splits[i, 1]) #n.c
      splits[i, 4] <- sum(data[[treat]]==1 & data[[x]]== splits[i, 1] & data[[y]]==1)/splits[i, 2] #p.t
      splits[i, 5] <- sum(data[[treat]]==0 & data[[x]]== splits[i, 1] & data[[y]]==1)/splits[i, 3] #p.c
      
      splits[i, 6] <- splits[i, 4] - splits[i, 5] #u
    }
    
    splits$cat.rank <- rank(splits$u, ties.method = "random")
    
    x.rank <- splits[, c(1,7)]
    x.merge <- merge(x.rank, data, by.x = 'cat', by.y = 'x')
    
    names(x.merge)[names(x.merge) == 'cat.rank'] <- 'x'
    
    x.rank <- x.rank[order(x.rank$cat.rank),]
    x.rank$cat <- paste0("'", x.rank$cat, "'")
    
    res <- list(x.rank, x.merge)   
    return(res)
    
  }
  
  #This part is optional but we could add a plot post-binning
  #Also, a sas code output could help to copy paste and perform the binning directly on sas
  
  BinUpliftPlot <- function(data, out.tree, x, ylim = NULL, ylab = "Uplift"){
    
    y <- deparse(substitute(y))
    treat <- deparse(substitute(treat))
    x <- deparse(substitute(x))
    
    final.cuts <- c(min(data[[x]]), sort(out.tree$x.cut), max(data[[x]]))
    
    data$x.cat <- cut(data[[x]], unique(final.cuts), include.lowest = TRUE, right = FALSE, dig.lab = 5)
    
    ##Create descriptive statistics associated with the binning
    x.cat.stats <- matrix(data = NA, nrow = 6, ncol = length(unique(final.cuts))-1)
    rownames(x.cat.stats) <- c("n", "n.t", "n.c", "p1.t", "p1.c", "uplift")
    colnames(x.cat.stats) <- names(table(data$x.cat))
    
    for(i in 1:length(unique(final.cuts))-1){ 
      
      x.cat.stats[1,] <- as.integer(table(data$x.cat))
      x.cat.stats[2, i] <- as.integer(sum(data[[treat]]==1 & data$x.cat == names(table(data$x.cat))[i] )) #n.t
      x.cat.stats[3, i] <- as.integer(sum(data[[treat]]==0 & data$x.cat == names(table(data$x.cat))[i] )) #n.c
      x.cat.stats[4, i] <- sum(data[[y]]==1 & data[[treat]]==1 & data$x.cat == names(table(data$x.cat))[i] )/ x.cat.stats[2, i] #p1.t
      x.cat.stats[5, i] <- sum(data[[y]]==1 & data[[treat]]==0 & data$x.cat == names(table(data$x.cat))[i] )/ x.cat.stats[3, i] #p1.c
      x.cat.stats[6, i] <- x.cat.stats[4, i] - x.cat.stats[5, i] #uplift
    }
    
    lines <- length(unique(final.cuts)) - 1
    sascode <- as.list(matrix(ncol = 0, nrow = 0))
    
    for (k in 1:lines){
      if (k==1) {sascode <- paste("if .<", x, "<", final.cuts[k+1], " then ", x, "_bin = '", k, "';", sep="")}
      else if (k > 1 &  k < lines) {sascode <- rbind(sascode, paste("else if .<", x, "<", final.cuts[k+1], " then ", x, "_bin = '", k, "';", sep=""))}
      else if (k==lines) {sascode <- rbind(sascode, paste("else if .<", x, " then ", x, "_bin = '", k, "';", sep=""))}
    }
    
    bar.plot <- barplot(x.cat.stats[6,], ylim = ylim, ylab = ylab, main = title, col = color)
    out.plot <- list("sas.code" = sascode, "bar.plot" = bar.plot)
    return(out.plot)
    
  }
  
  BinUpliftPlotC <- function(data, out.tree, out.link, x, ylim = NULL, ylab = "Uplift"){
    
    y <- deparse(substitute(y))
    treat <- deparse(substitute(treat))
    x <- deparse(substitute(x))
    
    final.cuts <- c(min(data[[x]]), sort(out.tree$x.cut), max(data[[x]]))
    
    data$x.cat <- cut(data[[x]], unique(final.cuts), include.lowest = TRUE, right = FALSE, dig.lab = 0)
    
    ##Create descriptive statistics associated with the binning
    
    x.cat.stats <- matrix(data = NA, nrow = 6, ncol = length(unique(final.cuts))-1)
    rownames(x.cat.stats) <- c("n", "n.t", "n.c", "p1.t", "p1.c", "uplift")
    colnames(x.cat.stats) <- names(table(data$x.cat))
    
    for(i in 1:length(unique(final.cuts))-1){ 
      
      x.cat.stats[1,] <- as.integer(table(data$x.cat))
      x.cat.stats[2, i] <- as.integer(sum(data[[treat]]==1 & data$x.cat == names(table(data$x.cat))[i] )) #n.t
      x.cat.stats[3, i] <- as.integer(sum(data[[treat]]==0 & data$x.cat == names(table(data$x.cat))[i] )) #n.c
      x.cat.stats[4, i] <- sum(data[[y]]==1 & data[[treat]]==1 & data$x.cat == names(table(data$x.cat))[i] )/ x.cat.stats[2, i] #p1.t
      x.cat.stats[5, i] <- sum(data[[y]]==1 & data[[treat]]==0 & data$x.cat == names(table(data$x.cat))[i] )/ x.cat.stats[3, i] #p1.c
      x.cat.stats[6, i] <- x.cat.stats[4, i] - x.cat.stats[5, i] #uplift
    }
    
    
    lines <- length(unique(final.cuts)) - 1
    sascode <- as.list(matrix(ncol = 0, nrow = 0))
    cut <- out.tree$x.cut
    
    for (k in 1:lines){
      if (k==1) {sascode <- paste("if ", x, " in (", toString(out.link[out.link$cat.rank < cut[k],1]), ") then ", x, "_bin = '", k, "';", sep="")}
      else if (k > 1 &  k < lines) {sascode <- rbind(sascode, paste("else if ", x, " in (", toString(out.link[out.link$cat.rank < cut[k] & out.link$cat.rank >= cut[k-1],1]), ") then ", x, "_bin = '", k, "';", sep=""))}
      else if (k==lines) {sascode <- rbind(sascode, paste("else if ", x, " in (", toString(out.link[out.link$cat.rank >= cut[k-1],1]), ") then ", x, "_bin = '", k, "';", sep=""))}
    }
    
    bar.plot <- barplot(x.cat.stats[6,], ylim = ylim, ylab = ylab, main = title, col = color)
    out.plot <- list("sas.code" = sascode, "bar.plot" = bar.plot)
    return(out.plot)
    
  }
  
  BinUpliftScore <- function(data, bin.result, y, treat, x){
    
    y <- deparse(substitute(y))
    treat <- deparse(substitute(treat))
    x <- deparse(substitute(x))
    
    Data_To_Score <- data[, c(paste(y), paste(treat), paste(x))]
    Data_To_Score$uplift <- 0
    K <- nrow(bin.result) + 1
    for (k in 1:K) { 
      if (k == 1) {
        index <- Data_To_Score[[x]] < bin.result[["x.cut"]][k]
        category <- Data_To_Score[index,]
        category$uplift <- sum(category[[treat]]==1 & category[[y]]==1)/sum(category[[treat]]==1) - sum(category[[treat]]==0 & category[[y]]==1)/sum(category[[treat]]==0)
        Data_To_Score[index,]$uplift <- category$uplift
        rm(category)
      }
      if (k > 1 & k < K) {
        index <- Data_To_Score[[x]] >= bin.result[["x.cut"]][k-1] & Data_To_Score[[x]] < bin.result[["x.cut"]][k]
        category <- Data_To_Score[index,]
        category$uplift <- sum(category[[treat]]==1 & category[[y]]==1)/sum(category[[treat]]==1) - sum(category[[treat]]==0 & category[[y]]==1)/sum(category[[treat]]==0)
        Data_To_Score[index,]$uplift <- category$uplift
        rm(category)
      }
      if (k == K) {
        index <- Data_To_Score[[x]] >= bin.result[["x.cut"]][k-1]
        category <- Data_To_Score[index,]
        category$uplift <- sum(category[[treat]]==1 & category[[y]]==1)/sum(category[[treat]]==1) - sum(category[[treat]]==0 & category[[y]]==1)/sum(category[[treat]]==0)
        Data_To_Score[index,]$uplift <- category$uplift
        rm(category)
      }
    }
    
    return(Data_To_Score)
  }
  
  
  #Finally, we can grow the tree and return it as out.tree depending on nature of the 
  #variable we want to bin
  
  if (is.factor(data[[x]])==TRUE) {
    out.cat <- BinUpliftCatRank(data, y, treat, x)[[2]]
    out.link <- BinUpliftCatRank(data, y, treat, x)[[1]]
    out.tree <- BinUpliftTree(out.cat, y, treat, x, length(unique(out.cat$x))-1, alpha, n.min)
    sas.code <- NULL
    
    if (length(out.tree) > 1) {
      sas.code <- BinUpliftPlotC(out.cat, out.tree, out.link, x, ylim, ylab, color, title)[[1]]
    }
    
    out.tree.cat <- list("out.tree" = out.tree, "out.link" = out.link, "sas.code" = sas.code)
    return(out.tree.cat)
  }
  
  else if (is.factor(data[[x]])==FALSE) {
    out.tree <- BinUpliftTree(data, y, treat, x, n.split, alpha, n.min)
    sas.code <- NULL
    score.data <- data.frame()
    
    if (length(out.tree) > 1) {
      sas.code <- BinUpliftPlot(data, out.tree, x, ylim, ylab, color, title)[[1]]
      score.data <- BinUpliftScore(data, out.tree, y, treat, x)
    }
    
    out.tree.cat <- list("out.tree" = out.tree, "sas.code" = sas.code, "score" = score.data)
    return(out.tree.cat)
  }
  
}

### END FUN

