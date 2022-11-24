######################### RJAGS WRAPPER AND FUNCTIONS FOR POST-PROCESSING MODEL ESTIMATES
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

# Parse rownames in quantile dataframe
extract_varnames=function(x) {
  Var = unlist(strsplit(x=x, split='[[]'))[1]
  Var = gsub("[[:digit:]]", "", Var)
  return(Var)
}
extract_index=function(x) {
  tempvar = unlist(strsplit(x=x, split='[[]'))[2]
  Index = unlist(strsplit(x=tempvar, split='[]]'))[1]
  return(Index)
}

# Post process results to separate estimates of emergence rates from model parameters
post_process_quantiles = function(file=all_quantiles){
  Var=apply(X=as.data.frame(rownames(file)),1,FUN=extract_varnames)
  Index=apply(X=as.data.frame(rownames(file)),1,FUN=extract_index)
  estimates=cbind(data.frame(Var),data.frame(Index),file)
  emergence_rates = data.frame(estimates[estimates$Var=='emergence',])
  parameters = data.frame(estimates[!(estimates$Var=='emergence'),])
  results=list(emergence_rates=emergence_rates, parameters=parameters)
  return(results)
}



extract_quantiles = function(jags.samples.out=jags.samples.out){
  quantiles = data.frame(summary(jags.samples.out)$quantiles)
  statistics = data.frame(summary(jags.samples.out)$statistics)
  quantiles$Mean = statistics$Mean
  quantiles$SD = statistics$SD
  return(quantiles)
}

