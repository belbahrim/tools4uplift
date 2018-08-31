splitUplift <-
function(data, p, group){
  
  require(fifer, quietly = TRUE)
  require(dplyr, quietly = TRUE)
  
  data$ID = seq(1:nrow(data))
  
  train <- stratified(df = data, group=paste(group) , size = p)
  valid <- data %>% anti_join(train, by = "ID")
  
  dataSplit <- list(train[,-ncol(data)], valid[,-ncol(data)])
  return(dataSplit)
}
