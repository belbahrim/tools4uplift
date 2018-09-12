splitUplift <- function(data, p, group){
  
  data$ID = seq(1:nrow(data))
  
  train <- data %>% group_by(paste(group, collapse = ',')) %>%  sample_frac(p)
  train <- as.data.frame(train[,-ncol(train)])
  valid <- data %>% anti_join(train, by = "ID")
  
  dataSplit <- list(train[,-ncol(train)], valid[,-ncol(valid)])
  return(dataSplit)
}