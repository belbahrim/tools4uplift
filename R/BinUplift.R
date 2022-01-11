BinUplift <- function(data, treat, outcome, x, n.split = 10, alpha = 0.05, n.min = 30){
  
  # Univariate quantization for uplift.
  #
  # Args:
  #   data: a data frame containing the treatment, the outcome and the predictors.
  #   treat: name of a binary (numeric) vector representing the treatment 
  #          assignment (coded as 0/1).
  #   outcome: name of a binary response (numeric) vector (coded as 0/1).
  #   x: name of the explanatory variable to categorize.
  #   ... and default parameters.
  #
  # Returns:
  #   The cut-offs for x.
  
  # Error handling
  # First, we need to check that the parameters are well filled and we create
  
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
  if (length(unique(data[[x]])) < 3) {
    stop("Independent variable must contain at least 3 different unique values")
  }
  if (sum(is.na(data[[x]])) > 0 ) {
    warning("Independent variable contains missing values: remove the observations and proceed")
  }
  
  # Keep complete cases only for the variable to categorize
  data <- data[is.na(data[[x]])==FALSE,]
  
  # Second, we need to define the subfunctions that create the splits and
  # choose the best split based on z-test
  
  BinUpliftStump <- function(data, outcome, treat, x, n.split){
    
    x.cut <- unique(quantile(data[[x]], seq(0, 1, 1/n.split)))
    
    splits <- matrix(data = NA, nrow = length(x.cut), ncol = 16)
    colnames(splits) <- c("x.cut", "n.lt", "n.lc", "p.lt", "p.lc", "u.l",
                          "n.rt", "n.rc", "p.rt", "p.rc", "u.r", "diff", "p.t", "p.c",
                          "p.t.n.t", "p.c.n.c")
    
    for(i in 1:length(x.cut)){ 
      
      index.l <- data[[x]] < x.cut[i]
      index.r <- data[[x]] >= x.cut[i]
      splits[i,1] <- x.cut[i]
      
      left <- data[index.l,]
      right <- data[index.r,]
      
      # uplift in left node
      splits[i, 2] <- sum(left[[treat]]==1)  # n.lt
      splits[i, 3] <- sum(left[[treat]]==0)  # n.lc
      splits[i, 4] <- sum(left[[treat]]==1 & left[[outcome]]==1)/splits[i, 2]  # p.lt
      splits[i, 5] <- sum(left[[treat]]==0 & left[[outcome]]==1)/splits[i, 3]  # p.lc
      
      splits[i, 6] <- splits[i, 4] - splits[i, 5]  # u.l
      
      # uplift in right node
      splits[i, 7] <- sum(right[[treat]]==1)  # n.rt
      splits[i, 8] <- sum(right[[treat]]==0)  # n.rc
      splits[i, 9] <- sum(right[[treat]]==1 & right[[outcome]]==1)/splits[i, 7]  # p.rt
      splits[i, 10] <- sum(right[[treat]]==0 & right[[outcome]]==1)/splits[i, 8]  # p.rc
      
      splits[i, 11] <- splits[i, 9] - splits[i, 10]  # u.r
      
      splits[i, 12] <- abs(splits[i, 6] - splits[i, 11])  # diff
      
      #probabilities of success in treatment and control groups
      splits[i, 13] <- sum(data[[treat]]==1 & data[[outcome]]==1)/sum(data[[treat]]==1)  # p.t
      splits[i, 14] <- sum(data[[treat]]==0 & data[[outcome]]==1)/sum(data[[treat]]==0)  # p.c
      
      #marginal counts needed for odds ratio
      splits[i, 15] <- sum(data[[treat]]==1 & data[[outcome]]==1)  #p.t.n.t
      splits[i, 16] <- sum(data[[treat]]==0 & data[[outcome]]==1)  #p.c.n.c
      
      
    }
    
    return(splits)
    
  }
  
  
  BinUpliftTest <- function(splits, alpha, n.min){
    
    z.a <- qnorm(0.5*alpha, mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)
    
    test.splits <- data.frame(splits)
    test.splits$sign <- 0
    test.splits <- test.splits[complete.cases(test.splits),]
    
    
    # We need the number of observations per group for the variance
    test.splits$n.t <- test.splits$n.lt + test.splits$n.rt
    test.splits$n.c <- test.splits$n.lc + test.splits$n.rc
    
    # We need the number of responses in the left node per group for the odds ratio
    test.splits$z_t <- test.splits$n.lt * test.splits$p.lt
    test.splits$z_c <- test.splits$n.lc * test.splits$p.lc
    
    test.splits$odds_ratio_t <- (test.splits$z_t / (test.splits$n.lt - test.splits$z_t)) / 
      ((test.splits$p.t * test.splits$n.t - test.splits$z_t)/(test.splits$n.rt - (test.splits$p.t * test.splits$n.t - test.splits$z_t)))
    test.splits$odds_ratio_c <- (test.splits$z_c / (test.splits$n.lc - test.splits$z_c)) / 
      ((test.splits$p.c * test.splits$n.c - test.splits$z_c)/(test.splits$n.rc - (test.splits$p.c * test.splits$n.c - test.splits$z_c)))
    
    # Extra checks in case odds_ratio is Inf of <= 0
    test.splits <- test.splits[is.finite(test.splits$odds_ratio_t) == TRUE,]
    test.splits <- test.splits[is.finite(test.splits$odds_ratio_c) == TRUE,]
    
    test.splits <- test.splits[test.splits$odds_ratio_t > 0,]
    test.splits <- test.splits[test.splits$odds_ratio_c > 0,]
    
    #Error handling in case there are no possible splits
    
    if (nrow(test.splits) == 0) {
      return(test.splits)
    }
    
    # We can now compute the numerator and denominator of the test statistic  
    test.splits$esp_t <- 0
    test.splits$esp_c <- 0
    test.splits$var_t <- 0
    test.splits$var_c <- 0
    
    for (i in 1:nrow(test.splits)){
      test.splits[i,]$esp_t <- BiasedUrn::meanFNCHypergeo(m1=test.splits[i,]$n.lt,
                                                          m2=test.splits[i,]$n.rt,
                                                          n=test.splits[i,]$p.t.n.t,
                                                          odds=test.splits$odds_ratio_t[i])
      test.splits[i,]$esp_c <- BiasedUrn::meanFNCHypergeo(m1=test.splits[i,]$n.lc,
                                                          m2=test.splits[i,]$n.rc,
                                                          n=test.splits[i,]$p.c.n.c,
                                                          odds=test.splits$odds_ratio_c[i])
      test.splits[i,]$var_t <- BiasedUrn::varFNCHypergeo(m1=test.splits[i,]$n.lt,
                                                         m2=test.splits[i,]$n.rt,
                                                         n=test.splits[i,]$p.t.n.t,
                                                         odds=test.splits$odds_ratio_t[i])
      test.splits[i,]$var_c <- BiasedUrn::varFNCHypergeo(m1=test.splits[i,]$n.lc,
                                                         m2=test.splits[i,]$n.rc,
                                                         n=test.splits[i,]$p.c.n.c,
                                                         odds=test.splits$odds_ratio_c[i])
    }
    
    test.splits$z.num <- (test.splits$p.lt - test.splits$p.lc - test.splits$p.rt + test.splits$p.rc) 
    
    test.splits$z.den <- sqrt( (test.splits$n.t^2)*test.splits$var_t / ((test.splits$n.lt^2)*(test.splits$n.rt^2)) +
                                 (test.splits$n.c^2)*test.splits$var_c / ((test.splits$n.lc^2)*(test.splits$n.rc^2)) )
    
    test.splits$z.obs <- abs(test.splits$z.num / test.splits$z.den)
    
    test.splits$sign <- 1*(test.splits$z.obs > z.a)
    
    argmax.size <- test.splits[test.splits$n.lt > n.min & test.splits$n.lc > n.min & test.splits$n.rt > n.min & test.splits$n.rc > n.min,]
    argmax.cut <- argmax.size[which.max(argmax.size$z.obs),]
    
    argmax.cut$x.pos <- which(test.splits$x.cut == argmax.cut$x.cut)	
    
    return(argmax.cut)
    
  }
  
  
  
  # Now we can call the recursive partitionning subfunction that performs the
  # splits based on the two (2) previous subfunctions
  
  BinUpliftTree <- function(data, outcome, treat, x, n.split, alpha, n.min){   
    
    stump <- BinUpliftStump(data, outcome, treat, x, n.split)
    best <- BinUpliftTest(stump, alpha, n.min)
    
    if (best$sign == 0 || nrow(best) == 0) { 
      
      return("oups..no significant split") 
      
    } else if (best$sign == 1) {
      
      # printing the cuts in order
      print(paste("The variable", x, "has been cut at:"))
      print(best$x.cut)
      
      # the x.cuts at the left and at the right of the x.cut optimal 
      l.cuts <- best$x.pos - 1
      r.cuts <- n.split-best$x.pos + 1
      
      l.index <- data[[x]] < best$x.cut
      r.index <- data[[x]] >= best$x.cut
      
      l.create <- data[l.index,]
      r.create <- data[r.index,]
      
      l.stump <- BinUpliftStump(l.create, outcome, treat, x, l.cuts)
      r.stump <- BinUpliftStump(r.create, outcome, treat, x, r.cuts)
      
      l.best <- BinUpliftTest(l.stump, alpha, n.min)  
      r.best <- BinUpliftTest(r.stump, alpha, n.min)
      
      # possible paths 
      if (l.best$sign == 0 || nrow(l.best) == 0) {
        if (r.best$sign == 0 || nrow(r.best) == 0) {
          return (best)
        } else if (r.best$sign == 1) {
          return (rbind(best, BinUpliftTree(r.create, outcome, treat, x, r.cuts, alpha, n.min))) 
        }
      } else if (l.best$sign == 1) {
        if (r.best$sign == 0 || nrow(r.best) == 0){
          return (rbind(BinUpliftTree(l.create, outcome, treat, x, l.cuts, alpha, n.min), best)) 
        } else if (r.best$sign == 1) {
          return (rbind(BinUpliftTree(l.create, outcome, treat, x, l.cuts, alpha, n.min), best, BinUpliftTree(r.create, outcome, treat, x, r.cuts, alpha, n.min))) 
        }
      }    
    }
  }
  
  # Now, if the variable we want to bin is not continuous, we need to change 
  # it to a continuous one to be able to split
  # We sort the categories by uplift and create an ordinal variable 1,2,...,K
  # where K is the number of the categories. Then, we can split that variable
  # with n.split=K-1
  
  BinUpliftCatRank <- function(data, outcome, treat, x){
    
    splits <- matrix(data = NA, nrow = length(levels(data[[x]])), ncol=6)
    colnames(splits) <- c("cat", "n.t", "n.c", "p.t", "p.c", "u")
    splits <- as.data.frame(splits)
    
    for(i in 1:nrow(splits)){ 
      
      splits[i, 1] <- levels(data[[x]])[[i]]
      
      # uplift in each category
      splits[i, 2] <- sum(data[[treat]]==1 & data[[x]]== splits[i, 1])  # n.t
      splits[i, 3] <- sum(data[[treat]]==0 & data[[x]]== splits[i, 1])  # n.c
      splits[i, 4] <- sum(data[[treat]]==1 & data[[x]]== splits[i, 1] & data[[outcome]]==1)/splits[i, 2]  # p.t
      splits[i, 5] <- sum(data[[treat]]==0 & data[[x]]== splits[i, 1] & data[[outcome]]==1)/splits[i, 3]  # p.c
      
      splits[i, 6] <- splits[i, 4] - splits[i, 5]  # u
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
  
  # Finally, we can grow the tree and return it as out.tree depending on nature of the 
  # variable we want to bin
  
  if (is.factor(data[[x]])==TRUE) {
    out.cat <- BinUpliftCatRank(data, outcome, treat, x)[[2]]
    out.link <- BinUpliftCatRank(data, outcome, treat, x)[[1]]
    out.tree <- BinUpliftTree(out.cat, outcome, treat, x, length(unique(out.cat$x))-1, alpha, n.min)
    
    out.tree.cat <- list("out.tree" = out.tree, "out.link" = out.link)
    return(out.tree.cat)
  }  else if (is.factor(data[[x]])==FALSE) {
    out.tree <- BinUpliftTree(data, outcome, treat, x, n.split, alpha, n.min)
    class(out.tree) <- "BinUplift"
    return(out.tree)
  }
  
}

# END FUN