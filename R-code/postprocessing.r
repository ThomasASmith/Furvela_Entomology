######################### FUNCTIONS FOR POST-PROCESSING AND PLOTTING MODEL ESTIMATES
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

