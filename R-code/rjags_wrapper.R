# Wrapper for running rJAGS and extracting credible intervals for parameters

library(rjags)
parameterEstimates = function(jagsdata=jagsdata,jagsmodel,variable.names=variable.names){
  library(jagsUI)
  jagsout = autojags(data=jagsdata, inits=NULL, 
                     parameters.to.save=variable.names, model.file=textConnection(jagsmodel),
                     n.chains= 4, n.adapt=NULL, iter.increment=1000, n.burnin=0, n.thin=1)
  maxIndex = as.integer(sapply(jagsout$mean, function(x) length(x)))
  Var = Index = X50. = X2.5. = X97.5. = mean = sd = c()
  for(i in 1:length(maxIndex)){
    for(j in 1:maxIndex[i]){
      Vari = names(jagsout$q50)[i]
      Var = c(Var,Vari)
      Index = c(Index,j)
      X50. = c(X50.,unlist(jagsout$q50[Vari])[j])
      X2.5. = c(X2.5.,unlist(jagsout$q2.5[Vari])[j])
      X97.5. = c(X97.5.,unlist(jagsout$q97.5[Vari])[j])
      mean = c(mean,unlist(jagsout$mean[Vari])[j])
      sd = c(sd,unlist(jagsout$sd[Vari])[j])
    }
  }
  df = data.frame(Var,Index,X50.,X2.5.,X97.5.,mean,sd)  
  return(df)
}
