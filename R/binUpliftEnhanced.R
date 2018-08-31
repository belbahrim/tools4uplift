binUpliftEnhanced <- function(data, 
                               var_list,
                               response_var_name, 
                               group_var_name, 
                               n.split=10, 
                               alpha=0.05, 
                               n.min=30, 
                               ylim=NULL,
                               ylab="Uplift",
                               title="Binning Results",...){
  
  for (var_name in var_list){
    # For now no gestion of categorical variables
      print(var_name)
      current_binning  <- binUplift(data,response_var_name,group_var_name,var_name,n.split,alpha,n.min,ylim,ylab)
      if(current_binning[[1]]== "oups..no significant split"){ 
        print(paste("Building: The variable",var_name,"will not be enhanced",sep =" "))
      }else{
        current_enhanced_var <- paste0("BinU_",var_name)
        data[,current_enhanced_var] = data[,var_name]
        data[,current_enhanced_var] = 0
        for (threshold in current_binning[[1]]$x.cut){
            data[,current_enhanced_var] = as.double(data[,var_name]>threshold)+ data[,current_enhanced_var]
        }
        
      }
  }
  return(data)
}
