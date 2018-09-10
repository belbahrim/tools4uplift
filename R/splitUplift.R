splitUplift <- function(data, p, group){
  
  data$ID = seq(1:nrow(data))
  
  train <- stratified(df = data, group=paste(group) , size = p)
  valid <- data %>% anti_join(train, by = "ID")
  
  dataSplit <- list(train[,-ncol(data)], valid[,-ncol(data)])
  return(dataSplit)
}
