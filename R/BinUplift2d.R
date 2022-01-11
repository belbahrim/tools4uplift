BinUplift2d <- function(data, var1, var2, treat, outcome, valid = NULL, 
                        n.split = 3, n.min = 30, plotit = FALSE, nb.col = 20){

  # Bivariate quantization
  #
  # Args:
  #   data: a data frame containing the treatment, the outcome and the predictors.
  #   var1: name of the first explanatory variable to quantize
  #   var2: name of the second explanatory variable to quantize
  #   treat: name of a binary (numeric) vector representing the treatment 
  #          assignment (coded as 0/1).
  #   outcome: name of a binary response (numeric) vector (coded as 0/1).
  #   valid: a validation dataset containing the treatment and the predictors.
  #   ... and default parameters.
  #
  # Returns:
  #   An augmented dataset with Uplift_var1_var2 variable representing a predicted 
  #   uplift for each observation based on the rectangle it belongs to.
  
  # Error handling
  if (n.split <= 1) stop("n.split must be > 1")
  if (n.split %% 1 != 0) stop("n.split must be an integer")
  
  if (n.min < 0) stop("n.min must be >= 0")
  if (n.min %% 1 != 0) stop("n.min must be an integer")
  
  if (nb.col < n.split) warning("nb.col should be greater than n.split for better visualization")
  if (nb.col %% 1 != 0) stop("nb.col must be an integer")
  
  # Initalize grid of uplifts for training set
  range_var1 <- ceiling(max(data[, paste(var1)])) - floor(min(data[, paste(var1)]))
  step1 <- range_var1 / n.split
  seq_col <- seq(floor(min(data[, paste(var1)])), ceiling(max(data[, paste(var1)])), step1)
  
  range_var2 <- ceiling(max(data[, paste(var2)])) - floor(min(data[, paste(var2)]))
  step2 <- range_var2 / n.split
  seq_row <- seq(floor(min(data[, paste(var2)])), ceiling(max(data[, paste(var2)])), step2)
  
  n <- length(seq_row)-1
  m <- length(seq_col)-1
  grid <- matrix(nrow=n, ncol=m)
  biprediction <- paste0("Uplift_", var1, "_", var2)
  data[, biprediction] <- NA
  
  if (is.null(valid) == FALSE){
    valid[, biprediction] <- NA
  }
  
  for (i in n:1){
    for (j in 1:m){
      
      index <- (data[, paste(var1)] <=  seq_col[j] + step1) & (data[, paste(var1)] >=  seq_col[j]) & (data[, paste(var2)] <= seq_row[n-i+1] + step2) & (data[, paste(var2)] >= seq_row[n-i+1])
      square <- data[index,]
      
      if (is.null(valid) == FALSE){
        index_valid <- (valid[, paste(var1)] <=  seq_col[j] + step1) & (valid[, paste(var1)] >=  seq_col[j]) & (valid[, paste(var2)] <= seq_row[n-i+1] + step2) & (valid[, paste(var2)] >= seq_row[n-i+1])
        square_valid <- valid[index_valid,]
      }
      
      indexT <- square[, paste(treat)]==1
      indexC <- square[, paste(treat)]==0
      sizeT <- sum(indexT)
      sizeC <- sum(indexC)
      
      if (sizeT >= n.min & sizeC >= n.min){
        grid[i, j] <- sum(square[, paste(treat)] == 1 & square[, paste(outcome)] == 1)/
                      sum(square[, paste(treat)] == 1) -
                      sum(square[, paste(treat)] == 0 & square[, paste(outcome)] == 1)/ 
                      sum(square[, paste(treat)] == 0)
        data[index, biprediction] <-  grid[i, j]
        if (is.null(valid) == FALSE){
          valid[index_valid, biprediction] <-  grid[i, j]
        }
        
      }
      
    }
  }
  
  # For the observations that were not scored, i.e. uplift=NA, impute the average
  # uplift of all NA
  
  index.naT <- is.na(data[, biprediction])==TRUE
  
  data[index.naT, biprediction] <- sum(data[index.naT, paste(treat)] == 1 & data[index.naT, paste(outcome)] == 1)/
                                   sum(data[index.naT, paste(treat)] == 1) - 
                                   sum(data[index.naT, paste(treat)] == 0 & data[index.naT, paste(outcome)] == 1)/
                                   sum(data[index.naT, paste(treat)] == 0)
  
  if (is.null(valid) == FALSE){
    index.naT_valid <- is.na(valid[, biprediction])==TRUE
  
    valid[index.naT_valid, biprediction] <- sum(data[index.naT, paste(treat)] == 1 & data[index.naT, paste(outcome)] == 1)/
      sum(data[index.naT, paste(treat)] == 1) - 
      sum(data[index.naT, paste(treat)] == 0 & data[index.naT, paste(outcome)] == 1)/
      sum(data[index.naT, paste(treat)] == 0)
    
  }
  
  main_train = ""
  
  if (is.null(valid) == FALSE){
    main_train = "Training Heatmap"
    main_valid = "Validation Heatmap"
  }
  
  
  if (plotit == TRUE){
    color.ramp.length <- nb.col
    negative.length <- round(abs(min(data[, biprediction]) - mean(data[, biprediction])) / 
                               diff(range(data[, biprediction])) * color.ramp.length)
    
    positive.length <- color.ramp.length - negative.length
    cols <- c(colorRampPalette(c("firebrick", "white"))(negative.length),
              colorRampPalette(c("white", "seagreen"))(positive.length))
    
    
    print(lattice::levelplot(as.formula(noquote(paste(biprediction, "~", var1,"*", var2))),
              xlab = var1,
              ylab = var2,
              data = data, 
              panel = panel.levelplot.points,
              col.regions = cols,
              main = main_train,
              colorkey = list(col = cols, 
                              at = lattice::do.breaks(range(data[, biprediction]), 
                                            color.ramp.length)),
              cex = 1.5) + layer_(panel.2dsmoother(method="lm", n = 100)))

    
    if (is.null(valid) == FALSE){
      print(lattice::levelplot(as.formula(noquote(paste(biprediction, "~", var1,"*", var2))),
                xlab = var1,
                ylab = var2,
                data = valid, 
                panel = panel.levelplot.points,
                col.regions = cols,
                main = main_valid,
                colorkey = list(col = cols, 
                                at = lattice::do.breaks(range(valid[, biprediction]),
                                              color.ramp.length)),
                cex = 1.5) + layer_(panel.2dsmoother(method="lm", n = 100)))
    }
    
  }
 
  if (is.null(valid) == FALSE){
    return(list(data, valid))
  } else {
    return(data)
  }
}

# END FUN