UpliftPerCat <- function(data, treat, outcome, x, ...){
  
  # Compute uplift per category and plot it.
  #
  # Args:
  #   data: a data frame containing the treatment, the outcome and the variable
  #         of interest
  #   treat: name of a binary (numeric) vector representing the treatment 
  #          assignment (coded as 0/1).
  #   outcome: name of a binary response (numeric) vector (coded as 0/1).
  #   x: name of the explanatory variable to study.
  #   ... extra parameters for barplot.
  #
  # Returns:
  #   A barplot representing the uplift per category.
  
  
  # Create descriptive statistics associated with the binning
  x.cat.stats <- matrix(data = NA, nrow = 6, ncol = length(levels(data[[x]])))
  
  rownames(x.cat.stats) <- c("n", "n.t", "n.c", "p1.t", "p1.c", "uplift")
  colnames(x.cat.stats) <- levels(data[[x]])
  
  for(i in 1:ncol(x.cat.stats)){ 
    x.cat.stats[1,] <- as.integer(table(data[[x]]))
    x.cat.stats[2, i] <- as.integer(sum(data[[treat]]==1 & data[[x]] == levels(data[[x]])[i] ))  # n.t
    x.cat.stats[3, i] <- as.integer(sum(data[[treat]]==0 & data[[x]] == levels(data[[x]])[i] ))  # n.c
    x.cat.stats[4, i] <- sum(data[[outcome]]==1 & data[[treat]]==1 & data[[x]] == levels(data[[x]])[i] )/ x.cat.stats[2, i]  # p1.t
    x.cat.stats[5, i] <- sum(data[[outcome]]==1 & data[[treat]]==0 & data[[x]] == levels(data[[x]])[i] )/ x.cat.stats[3, i]  # p1.c
    x.cat.stats[6, i] <- x.cat.stats[4, i] - x.cat.stats[5, i]  # uplift
  }
  
  bar.plot <- barplot(x.cat.stats[6,], ...)
  abline(h = 0)
  
}

# END FUN