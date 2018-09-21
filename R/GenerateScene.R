GenerateScene <- function(space_dimension=2,
                        Principal_effect=c(0.1,1),
                        Principal_std=c(1,4),
                        Control_effect_positive=c(0.2,0.3),
                        Control_std_positive=c(1.5,2),
                        Control_effect_negative=c(0.3,0.4),
                        Control_std_negative=c(1.5,2),
                        Nb_point=c(20,10,20),
                        plotit=TRUE){
  
  #Main effet definition
  Principal_effect_max = Principal_effect[2]
  Principal_effect_min = Principal_effect[1]
  Principal_std_max = Principal_std[[2]]
  Principal_std_min = Principal_std[[1]]
  
  # Control the effect positive
  Control_effect_positive_max = Control_effect_positive[[2]]
  Control_effect_positive_min = Control_effect_positive[[1]]
  Control_std_positive_max = Control_std_positive[[2]]
  Control_std_positive_min = Control_std_positive[[1]]
  
  # Control the effect
  Control_effect_negative_max = Control_effect_negative[[2]]
  Control_effect_negative_min = Control_effect_negative[[1]]
  Control_std_negative_max = Control_std_negative[[2]]
  Control_std_negative_min = Control_std_negative[[1]]
  
  #Define the number of point for each generating set
  Nb_point_principale = Nb_point[[1]]
  Nb_point_control_positive = Nb_point[[2]]
  Nb_point_control_negative = Nb_point[[3]]
  
  #Genrating point uniformly in the compact square [-1,1]^dimension
  X_Principal = data.frame(replicate(space_dimension, runif(Nb_point_principale,-1,1)))
  X_Control_positive   = data.frame(replicate(space_dimension, runif(Nb_point_control_positive,-1,1)))
  X_Control_negative   = data.frame(replicate(space_dimension, runif(Nb_point_control_negative,-1,1)))
  
  Vect_effect_principale = runif(Nb_point_principale,Principal_effect_min,Principal_effect_max)
  Vect_princ_std         = runif(Nb_point_principale, Principal_std_min,Principal_std_max)
  
  Vect_effect_control_positive = runif(Nb_point_control_positive,Control_effect_positive_min,Control_effect_positive_max)
  Vect_control_positive_std    = runif(Nb_point_control_positive, Control_std_positive_min,Control_std_positive_max)
  
  Vect_effect_control_negative = runif(Nb_point_control_negative,Control_effect_negative_min,Control_effect_negative_max)
  Vect_control_negative_std    = runif(Nb_point_control_negative, Control_std_negative_min,Control_std_negative_max)
  
  
  #parameter meant to control for each point the weight of influence of each variable
  distortion_principal = data.frame(replicate(space_dimension,runif(Nb_point_principale,0.001,1)))
  distortion_positive = data.frame(replicate(space_dimension,runif(Nb_point_principale,0.001,1)))
  distortion_negative = data.frame(replicate(space_dimension,runif(Nb_point_principale,0.001,1)))
  
  
  Test_propensity <- function(x){
    Test_propensity = 0
    for (i in 1:nrow(X_Principal)){
      Test_propensity <- max(Vect_effect_principale[i]*exp(-norm((x-X_Principal[i,])/distortion_principal[i,], type="2")/Vect_princ_std[i]),
                             Test_propensity)
    }
    return(Test_propensity)
  }
  
  
  perturbation_positive <- function(x){
    perturbation_positive = 0
    for (i in 1:nrow(X_Control_positive)){
      perturbation_positive = max(
        Vect_effect_control_positive[i]*
          exp(-norm((x-X_Control_positive[i,])/distortion_positive[i,], type="2")/Vect_control_positive_std[i] ),
        perturbation_positive)
    }
    return(perturbation_positive)
  }
  
  perturbation_negative <- function(x){
    perturbation_negative = 0
    for (i in 1:nrow(X_Control_negative)){
      perturbation_negative = max(
        Vect_effect_control_negative[i]*
          exp(-norm((x-X_Control_negative[i,])/distortion_negative[i,], type="2")/Vect_control_negative_std[i] ),perturbation_negative)
    }
    return(perturbation_negative)
  }
  

  Control_propensity <- function(x) {
    return(min(  Test_propensity(x) - perturbation_negative(x) +perturbation_positive(x),1))
  }
  
  
  if (plotit == TRUE){
    a=seq(-2,2,0.05)
    resultp <- c()
    resultc <- c()
    null_vector = replicate(space_dimension-1,0)
    
    
    for (i in 1:length(a)){
      x=cbind(a[i],null_vector)
      resultp[i] <-  Test_propensity(x)
      resultc[i] <-  Control_propensity(x)
    }
    
    resultv <- resultc - resultp
    
    plot(a,resultc, col="red", type="l",xlab="dimension 1",ylab="distribution",ylim=c(-1,1))
    lines(a,resultp, col="blue")
    lines(a,resultv, col="green")
    legend(-1.5,-0.2,legend=c("Treatment Propensity","Control Propensity", "Uplift"),pch = "-",col=c("blue","red","green"))
    abline(0,0)
    
    for (i in 1:length(a)){
      x=cbind(null_vector,a[i])
      resultp[i] <-  Test_propensity(x)
      resultc[i] <-  Control_propensity(x)
    }
    
    resultv <- resultc - resultp
    
    plot(a,resultc, col="red", type="l",xlab="dimension 2",ylab="distribution",ylim=c(-1,1))
    lines(a,resultp, col="blue")
    lines(a,resultv, col="green")
    legend(-1.5,-0.2,legend=c("Treatment Propensity","Control Propensity", "Uplift"),pch = "-",col=c("blue","red","green"))
    abline(0,0)
    
  }
  
  return(list(Test_propensity, Control_propensity, space_dimension))
  
}

#END FUN
