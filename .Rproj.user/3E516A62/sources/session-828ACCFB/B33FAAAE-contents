################# ANALYSIS OF TIME SERIES OF EXIT TRAP DATA 

analyseExitTrap = function(firstquarter=15,lastquarter=16){   
  data=selectData(firstquarter=firstquarter,lastquarter=lastquarter,A0=proportions$A0,P=0.528) 
  quantiles = analyseExitTrap_allmodels(all_quantiles=NULL,jagsdata=data$jagsdata,firstquarter=firstquarter)
  # create a vector of emergence rates noting that these are mostly missing
  emergence = as.vector(quantiles[quantiles$Var=='emergence' 
                                  & quantiles$jagsModel=='P = est.','X50.'])
  emergence_full = rep(NA,nrow(data$plotdata))
  emergence_full[data$jagsdata$ptr]=emergence
  scaled = data.frame(emergence = emergence_full/mean(emergence_full,na.rm = TRUE))
  scaled$unfed1 = with(data$plotdata,mean.Af.unfed1/mean(mean.Af.unfed1,na.rm=TRUE))
  scaled$gravid = with(data$plotdata,mean.Af.gravid/mean(mean.Af.gravid,na.rm=TRUE))
  scaled$male = with(data$plotdata,mean.Af.maleEx/mean(mean.Af.maleEx,na.rm=TRUE))
  value = with(data$plotdata,c(emergence_full,mean.Af.unfed1,mean.Af.gravid,mean.Af.maleEx)) 
  variable = rep(c('Emergent','Unfed','Gravid','Males'),each=nrow(data$plotdata))
  levels(variable)=c('Emergent','Unfed','Gravid','Males')
  rescaled = with(scaled,c(emergence,unfed1,gravid,male))
  plotdate = rep(data$plotdata$date,times=4)
  toPlot = data.frame(plotdate,variable,value,rescaled)
  library(ggplot2)
  get_mp_recomCol = function() {return(c("#D55E00", "#009E73", "#0072A7","#C879C8"))}
  plt1 =ggplot(data = toPlot,aes(x = plotdate,group=variable)) + theme_bw() +
    theme(text = element_text(size=12)) +
    scale_colour_manual(values=get_mp_recomCol(),name='Mosquito category') +
    geom_point(aes(y=value,colour=variable),size=2, show.legend=TRUE) +
    theme(legend.position = c(.75,.8)) +
    theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid')) +
    #scale_x_continuous(name = 'Day',limits=c(0,80)) +
    scale_y_continuous(name = 'Mosquitoes per exit trap')
  plt2 =ggplot(data = toPlot,aes(x = plotdate,group=variable)) + theme_bw() +
    theme(text = element_text(size=12)) +
    scale_colour_manual(values=get_mp_recomCol(),name='Mosquito category') +
    geom_point(aes(y=rescaled,colour=variable),size=2, show.legend=FALSE) +
    scale_x_continuous(name = 'Day',limits=c(0,80)) +
    scale_y_log10(name = 'Relative numbers',limits=c(0.1,10))
  library(cowplot)
  plt=plot_grid(plt1, plt2, labels = c('A', 'B'), label_size = 12,ncol=2)
  result = list(quantiles=quantiles,plt=plt,jagsdata=data$jagsdata,jagsmodel=exittrap1)
  return(result)}

analyseExitTrap_allmodels = function(all_quantiles,jagsdata,firstquarter,repeats=1) {
  jagsdata$P = rep(jagsdata$P,repeats)
  quantiles = parameterEstimates(jagsdata=jagsdata,jagsmodel=exittrap1,variable.names=variable.names)
  quantiles$firstquarter = firstquarter
  quantiles$jagsModel='Reference'
  all_quantiles = rbind(all_quantiles,quantiles)
  
  jagsdata$P = rep(0.75,repeats)
  quantiles = parameterEstimates(jagsdata=jagsdata,jagsmodel=exittrap1,variable.names=variable.names)
  quantiles$firstquarter = firstquarter
  quantiles$jagsModel='P = 0.75'
  all_quantiles = rbind(all_quantiles,quantiles)
  
  jagsdata$P = NULL
  quantiles = parameterEstimates(jagsdata=jagsdata,jagsmodel=exittrap1,variable.names=variable.names)
  quantiles$firstquarter = firstquarter
  quantiles$jagsModel='P = est.'
  all_quantiles = rbind(all_quantiles,quantiles)
  
  jagsdata$Teu = 1
  quantiles = parameterEstimates(jagsdata=jagsdata,jagsmodel=exittrap1,variable.names=variable.names)
  quantiles$firstquarter = firstquarter
  quantiles$jagsModel='Teu = 1'
  all_quantiles = rbind(all_quantiles,quantiles)
  
  jagsdata$Tem = 1
  quantiles = parameterEstimates(jagsdata=jagsdata,jagsmodel=exittrap1,variable.names=variable.names)
  quantiles$firstquarter = firstquarter
  quantiles$jagsModel='Tem = 1'
  all_quantiles = rbind(all_quantiles,quantiles)
return(all_quantiles)}


# Estimate survival, cycle length and sac rate from exit-trap data by year
analyseExitTrapByYear = function(){
  all_quantiles = NULL
  fq = seq(1:6)*4 + 4
  for(firstquarter in fq){
    lastquarter=firstquarter + 3
    jagsdata = selectData(firstquarter=firstquarter,lastquarter=lastquarter,
                          A0=proportions$A0,P=0.528)$jagsdata
    all_quantiles=analyseExitTrap_allmodels(all_quantiles=all_quantiles,jagsdata=jagsdata,firstquarter=firstquarter)
  }
  # include only results from models that are relevant: females
  results = all_quantiles[which(all_quantiles$jagsModel %in% c('Reference','P = est.','Teu = 1')),]
  p4a=plotEstimatesByTime(param='Teu', paramLabel= 'trapping efficiency (gravids)',results=results,requirelegend=TRUE)
  p4b=plotEstimatesByTime(param='cycle', paramLabel= 'cycle duration (days)',results=results)
  p4c=plotEstimatesByTime(param='P', paramLabel= 'survival per cycle',results=results)
  p4d=plotEstimatesByTime(param='p', paramLabel= 'daily survival of females',results=results)
  
  library(cowplot)
  plt1=plot_grid(p4a, p4b, p4c, p4d, labels = c('A', 'B', 'C', 'D'), label_size = 12,ncol=2)
  
  # include only results from models that are relevant: males
  results = all_quantiles[which(all_quantiles$jagsModel %in% c('Reference','Tem = 1')),]
  p4e=plotEstimatesByTime(param='Tem', paramLabel= 'trapping efficiency (males)',results=results,requirelegend=TRUE)
  p4f=plotEstimatesByTime(param='Pm', paramLabel= 'daily survival of males',results=results)
  plt2=plot_grid(p4e, p4f, labels = c('A', 'B'), label_size = 12,ncol=2)
  
  all_results = list(quantiles=all_quantiles,plt1=plt1,plt2=plt2,p4e=p4e,jagsmodel=exittrap1)
  return(all_results)}

analyseExitTrapByTemp = function(){
  data=selectData(firstquarter=1,lastquarter=32,model = 'temperature',A0=proportions$A0,P=0.528) 
  data$jagsdata$tcat = cut(data$jagsdata$mean_temp, breaks=c(18, 23, 25, 27, 29, 32), labels = FALSE)
  data$jagsdata$mean_temp = NULL
  all_quantiles = NULL
  all_quantiles=analyseExitTrap_allmodels(all_quantiles=all_quantiles,
                                          jagsdata=data$jagsdata,
                                          firstquarter=1,repeats=5)
  
  results = all_quantiles[which(all_quantiles$jagsModel %in% c('Reference','P = est.','Teu = 1')),]
  p5a=plotEstimatesByTemp(param='cycle', paramLabel= 'cycle duration (days)',results=results)
  p5b=plotEstimatesByTemp(param='P', paramLabel= 'survival per cycle',results=results)
  p5c=plotEstimatesByTemp(param='p', paramLabel= 'daily survival of females',results=results,requirelegend=TRUE)
  
  # include only results from models that are relevant: males
  results = all_quantiles[which(all_quantiles$jagsModel %in% c('Reference','Tem = 1')),]
  p5d=plotEstimatesByTemp(param='Pm', paramLabel= 'daily survival of males',results=results,requirelegend=TRUE)
  
  library(cowplot)
  plt1=plot_grid(p5a, p5b, p5c, p5d, labels = c('A', 'B', 'C', 'D'), label_size = 12,ncol=2)
  results=list(quantiles=all_quantiles,plt1=plt1,jagsmodel=exittrap1)
  return(results)}

analyseExitTrapByTempTrend = function(){
  data=selectData(firstquarter=1,lastquarter=32,model = 'temperature',A0=proportions$A0,P=0.528) 
  data$jagsdata$tcat = cut(data$jagsdata$mean_temp, breaks=c(18, 23, 25, 27, 29, 32), labels = FALSE)
  all_quantiles = NULL
  all_quantiles=analyseExitTrap_allmodels(all_quantiles=all_quantiles,
                                          jagsdata=data$jagsdata,
                                          firstquarter=1,repeats=5)
  
  results=list(quantiles=all_quantiles,jagsmodel=exittrap1)
  return(results)}

# JAGS code for estimation of emergence rates and both male and female survival from exit trap data 
# Unfed mosquitoes and emergent mosquitoes are equivalent in this implementation
exittrap_base ="
model{
  for(d in 1:sampleddates){
    # numbers of unfed (+ part fed) per trap (equivalent to emergence rate)
    T0[d] <- unfed1[d]/traps[d]    
    # gravid per trap
    eg[d] <- gravid[d]/traps[d]
  }
  for(t in 1:days_with_complete_data){
    # expected number gravid per trap (Teu is relative trapping efficiency of gravid vs unfed)

    # expected gravids in first cycle
    eg0[t] <- (T0[lag1[t]]*pr[1] + T0[lag2[t]]*pr[2] + T0[lag3[t]]*pr[3] + T0[lag4[t]]*pr[4])* Teu 
    # expected gravids in subsequent cycles
    eg1[t] <- eg[lag1[t]]*pr[1] + eg[lag2[t]]*pr[2] + eg[lag3[t]]*pr[3] + eg[lag4[t]]*pr[4]
    eg2[t] <- eg[lag2[t]]*pr[1] + eg[lag3[t]]*pr[2] + eg[lag4[t]]*pr[3] + eg[lag5[t]]*pr[4]
    G[t] <- (eg0[t] + A0*eg1[t] + (1-A0)*eg2[t])*P
    pg[t] <- rg/(G[t]*traps[ptr[t]] + rg)  # Negative binomial p 
    gravid[ptr[t]] ~ dnegbin(pg[t],rg)
    
    
    #emergence is equivalent to unfed
    emergence[t] <- T0[ptr[t]]
    
    # expected number of males:
    #     Tem is relative trapping efficiency of males vs unfed assuming 50:50 sex ratio on emergence
    #     Pm is daily survival of males
    
    # Expected males 
    Emales[t] <- (Tem*emergence[t] + Pm * males[lag1[t]]/traps[lag1[t]])*traps[ptr[t]]
    
    pm[t] <- rm/(Emales[t] + rm)  # Negative binomial p 
    males[ptr[t]] ~ dnegbin(pm[t],rm)
    
  }

  
  ################# PRIORS ########################

  # discretised normal kernel for the cycle duration
 
  pr[1] <- pnorm(1.5,resting,tau.b)
  pr[2] <- pnorm(2.5,resting,tau.b) - pr[1]
  pr[3] <- pnorm(3.5,resting,tau.b) - pr[1] - pr[2]
  pr[4] <- 1 - pr[1] - pr[2] - pr[3]
  
  tau.b ~ dgamma(1,1)    # precision of distribution of cycle duration 
  rg ~ dgamma(1,1)       # Negative binomial overdispersion for gravid females
  rm ~ dgamma(1,1)       # Negative binomial overdispersion for males

  resting ~ dunif(0.5,4.5) # duration of resting period
  cycle <- resting + (1 - A0)
  Teu ~ dgamma(1,1)       # relative trapping efficiency of gravids relative to unfeds 
  Tem ~ dgamma(1,1)       # relative trapping efficiency of males vs unfed assuming 50:50 sex ratio on emergence
  Pm ~ dunif(0,1)         # daily survival of males
  P ~ dunif(0,1)          # overall survival of females through oviposition cycle
  A0 ~ dunif(0,1)         # proportion returning to feed within same night
  p <- pow(P,1/cycle)   # daily survival of females
}  
"

# Estimation of temperature dependent emergence rates and male survival from exit trap data 
# Unfed mosquitoes and emergent mosquitoes are equivalent
exittrap_temp ="
model{
  for(d in 1:sampleddates){
    # numbers of unfed (+ part fed) per trap (equivalent to emergence rate)
    T0[d] <- unfed1[d]/traps[d]    
    # gravid per trap
    eg[d] <- gravid[d]/traps[d]
  }
  for(t in 1:days_with_complete_data){
    # expected number gravid per trap (Teu is relative trapping efficiency of gravid vs unfed)

    # expected gravids in first cycle
    eg0[t] <- (T0[lag1[t]]*pr[1,tcat[t]] + T0[lag2[t]]*pr[2,tcat[t]] + T0[lag3[t]]*pr[3,tcat[t]] + T0[lag4[t]]*pr[4,tcat[t]])* Teu 
    # expected gravids in subsequent cycles
    eg1[t] <- eg[lag1[t]]*pr[1,tcat[t]] + eg[lag2[t]]*pr[2,tcat[t]] + eg[lag3[t]]*pr[3,tcat[t]] + eg[lag4[t]]*pr[4,tcat[t]]
    eg2[t] <- eg[lag2[t]]*pr[1,tcat[t]] + eg[lag3[t]]*pr[2,tcat[t]] + eg[lag4[t]]*pr[3,tcat[t]] + eg[lag5[t]]*pr[4,tcat[t]]
    G[t] <- (eg0[t] + A0*eg1[t] + (1-A0)*eg2[t])*P[tcat[t]]
    pg[t] <- rg/(G[t]*traps[ptr[t]] + rg)  # Negative binomial p 
    gravid[ptr[t]] ~ dnegbin(pg[t],rg)
    
    
    #emergence is equivalent to unfed
    emergence[t] <- T0[ptr[t]]
    
    # expected number of males:
    #     Tem is relative trapping efficiency of males vs unfed assuming 50:50 sex ratio on emergence
    #     Pm is daily survival of males
    
    # Expected males 
    Emales[t] <- (Tem*emergence[t] + Pm[tcat[t]] * males[lag1[t]]/traps[lag1[t]])*traps[ptr[t]]
    
    pm[t] <- rm/(Emales[t] + rm)  # Negative binomial p 
    males[ptr[t]] ~ dnegbin(pm[t],rm)
    
  }

  
  ################# PRIORS ########################
  
  # discretised normal kernel for the cycle duration
  for(h in 1:5){
    pr[1,h] <- pnorm(1.5,resting[h],tau.b)
    pr[2,h] <- pnorm(2.5,resting[h],tau.b) - pr[1,h]
    pr[3,h] <- pnorm(3.5,resting[h],tau.b) - pr[1,h] - pr[2,h]
    pr[4,h] <- 1 - pr[1,h] - pr[2,h] - pr[3,h]
    resting[h] ~ dunif(0.5,4.5) # duration of resting period
    cycle[h] <- resting[h] + (1 - A0)
    Pm[h] ~ dunif(0,1)         # daily survival of males
    P[h] ~ dunif(0,1)          # overall survival through oviposition cycle
    p[h] <- pow(P[h],1/cycle[h])   # daily survival of females
  }
  
  tau.b ~ dgamma(1,1)    # precision of distribution of cycle duration 
  rg ~ dgamma(1,1)       # Negative binomial overdispersion for gravid females
  rm ~ dgamma(1,1)       # Negative binomial overdispersion for males
  Teu ~ dgamma(1,1)       # relative trapping efficiency of gravids relative to unfeds 
  Tem ~ dgamma(1,1)       # relative trapping efficiency of males vs unfed assuming 50:50 sex ratio on emergence
  A0 ~ dunif(0,1)         # proportion returning to feed within same night
}  
"

# Estimation of temperature dependent emergence rates and male survival from exit trap data 
# Unfed mosquitoes and emergent mosquitoes are equivalent
exittrap_temptrend ="
model{
  for(d in 1:sampleddates){
    # numbers of unfed (+ part fed) per trap (equivalent to emergence rate)
    T0[d] <- unfed1[d]/traps[d]    
    # gravid per trap
    eg[d] <- gravid[d]/traps[d]
  }
  for(t in 1:days_with_complete_data){
    Pm[t] <- 1/(exp(Pm0 - (mean_temp[t]-25)*Pm1))        # daily survival of males
    Pt[t] <- 1/(exp(P0 - (mean_temp[t]-25)*P1))          # overall survival through oviposition cycle

    # expected number gravid per trap (Teu is relative trapping efficiency of gravid vs unfed)

    # expected gravids in first cycle
    eg0[t] <- (T0[lag1[t]]*pr[1,tcat[t]] + T0[lag2[t]]*pr[2,tcat[t]] + T0[lag3[t]]*pr[3,tcat[t]] + T0[lag4[t]]*pr[4,tcat[t]])* Teu 
    # expected gravids in subsequent cycles
    eg1[t] <- eg[lag1[t]]*pr[1,tcat[t]] + eg[lag2[t]]*pr[2,tcat[t]] + eg[lag3[t]]*pr[3,tcat[t]] + eg[lag4[t]]*pr[4,tcat[t]]
    eg2[t] <- eg[lag2[t]]*pr[1,tcat[t]] + eg[lag3[t]]*pr[2,tcat[t]] + eg[lag4[t]]*pr[3,tcat[t]] + eg[lag5[t]]*pr[4,tcat[t]]
    G[t] <- (eg0[t] + A0*eg1[t] + (1-A0)*eg2[t])*Pt[t]
    pg[t] <- rg/(G[t]*traps[ptr[t]] + rg)  # Negative binomial p 
    gravid[ptr[t]] ~ dnegbin(pg[t],rg)
    
    
    #emergence is equivalent to unfed
    emergence[t] <- T0[ptr[t]]
    
    # expected number of males:
    #     Tem is relative trapping efficiency of males vs unfed assuming 50:50 sex ratio on emergence
    #     Pm is daily survival of males
    
    # Expected males 
    Emales[t] <- (Tem*emergence[t] + Pm[t] * males[lag1[t]]/traps[lag1[t]])*traps[ptr[t]]
    
    pm[t] <- rm/(Emales[t] + rm)  # Negative binomial p 
    males[ptr[t]] ~ dnegbin(pm[t],rm)
    
  }

  
  ################# PRIORS ########################
  
  # discretised normal kernel for the cycle duration
  for(h in 1:5){
    pr[1,h] <- pnorm(1.5,resting[h],tau.b)
    pr[2,h] <- pnorm(2.5,resting[h],tau.b) - pr[1,h]
    pr[3,h] <- pnorm(3.5,resting[h],tau.b) - pr[1,h] - pr[2,h]
    pr[4,h] <- 1 - pr[1,h] - pr[2,h] - pr[3,h]
    resting[h] ~ dunif(0.5,4.5) # duration of resting period
    cycle[h] <- resting[h] + (1 - A0)
    P[h] ~ dunif(0,1) # Redundant assignment included to avert warning
  }
  
  P0 ~ dnorm(0,1.0E-3)
  P1 ~ dnorm(0,1.0E-3)
  Pm0 ~ dnorm(0,1.0E-3)
  Pm1 ~ dnorm(0,1.0E-3)
  
  tau.b ~ dgamma(1,1)    # precision of distribution of cycle duration 
  rg ~ dgamma(1,1)       # Negative binomial overdispersion for gravid females
  rm ~ dgamma(1,1)       # Negative binomial overdispersion for males
  Teu ~ dgamma(1,1)       # relative trapping efficiency of gravids relative to unfeds 
  Tem ~ dgamma(1,1)       # relative trapping efficiency of males vs unfed assuming 50:50 sex ratio on emergence
  A0 ~ dunif(0,1)         # proportion returning to feed within same night
}  
"

SimulationsExitTrapsFixedSampleSize = function(nsimulations){
  # sampling quarters 8-19 gives 207 days with complete data (precedes introduction of LLINs)
  jagsdata = selectData(firstquarter=8,lastquarter=19,A0=proportions$A0,P=0.528)$jagsdata
  # Analyse simulated datasets
  all_simulated_quantiles = NULL
  inputs_to_simulations = NULL
  results = simulationsExitTraps(nsimulations = nsimulations,jagsdata=jagsdata)
  return(results)  }

consistencyAnalysis = function(nsimulations){
  # Simulate and analyse datasets of varying sizes
  for(i in 1:nsimulations){
    firstquarter = sample(8:27, size=1)
    lastquarter = sample(firstquarter:28, size=1)
    jagsdata = selectData(firstquarter=firstquarter,lastquarter=lastquarter,A0=proportions$A0,P=0.528)$jagsdata
    results = simulationsExitTraps(nsimulations = 1,jagsdata=jagsdata)
    results$inputs_to_simulations$days_with_complete_data=jagsdata$days_with_complete_data
    if(i == 1){
      all_simulated_quantiles = results$simresults$parameters
      inputs_to_simulations = results$inputs_to_simulations
    } else {
      all_simulated_quantiles = rbind(all_simulated_quantiles,results$simresults$parameters)
      inputs_to_simulations = rbind(inputs_to_simulations,results$inputs_to_simulations)
    }
  }
  all_results = list(all_simulated_quantiles, inputs_to_simulations)
  simresults = post_process_quantiles(file=all_simulated_quantiles)
  
  all_results = list(simresults=simresults,inputs_to_simulations=inputs_to_simulations,jagsmodel=exittrap1)
  return(all_results)}

simulationsExitTraps = function(nsimulations = nsimulations,jagsdata=jagsdata){  
  all_simulated_quantiles = NULL
  inputs_to_simulations = NULL
  variable.names = c("tau.b","rg","Teu","resting","pr","Tem","Pm","rm","P")
  for(i in 1:nsimulations){
    # parameters that might be recoverable
    tau.b = runif(1, min = 0.5, max = 1.5)
    rg = runif(1, min = 1, max = 4)
    te = runif(1, min = 0.5, max = 2.0)
    resting = runif(1, min = 1.0, max = 4.0)
    Teu= runif(1, min = 1, max = 5)
    Tem= runif(1, min = 1, max = 5)
    rm = runif(1, min = 1, max = 4)
    Pm = runif(1, min = 0, max = 1)
    P = runif(1, min = 0, max = 1)
    data0 = simulateExitTrapData(tau.b = tau.b, rg=rg, resting=resting,Teu=Teu,
                                 Tem = Tem, Pm = Pm,rm = rm,A0=proportions$A0,P=P,jagsdata=jagsdata)
    inputs=data.frame(data0$inputs)
    inputs$days_with_complete_data=jagsdata$days_with_complete_data
    # P known, sac rate as input
    simulated_quantiles = parameterEstimates(jagsdata=data0$jagsdata,jagsmodel=exittrap1,variable.names=variable.names)
    simulated_quantiles$simulation=i
    simulated_quantiles$jagsModel='exittrap1'
    all_simulated_quantiles = rbind(all_simulated_quantiles,simulated_quantiles)
    # P fixed at 0.75
    data0$jagsdata$P = 0.75
    simulated_quantiles = parameterEstimates(jagsdata=data0$jagsdata,jagsmodel=exittrap1,variable.names=variable.names)
    simulated_quantiles$simulation=i
    simulated_quantiles$jagsModel='P_075'
    all_simulated_quantiles = rbind(all_simulated_quantiles,simulated_quantiles)
    # P known, A0 fixed at 0.3
    data0$jagsdata$P = inputs$P
    data0$jagsdata$A0 = 0.3 
    simulated_quantiles = parameterEstimates(jagsdata=data0$jagsdata,jagsmodel=exittrap1,variable.names=variable.names)
    simulated_quantiles$simulation=i
    simulated_quantiles$jagsModel='A_03'
    all_simulated_quantiles = rbind(all_simulated_quantiles,simulated_quantiles)
    # P estimated
    data0$jagsdata$P = NULL
    data0$jagsdata$A0 = inputs$A0 
    simulated_quantiles = parameterEstimates(jagsdata=data0$jagsdata,jagsmodel=exittrap1,variable.names=variable.names)
    simulated_quantiles$simulation=i
    simulated_quantiles$jagsModel='P_est'
    all_simulated_quantiles = rbind(all_simulated_quantiles,simulated_quantiles)
    # Teu known, P Estimated
    data0$jagsdata$Teu = inputs$Teu
    simulated_quantiles = parameterEstimates(jagsdata=data0$jagsdata,jagsmodel=exittrap1,variable.names=variable.names)
    simulated_quantiles$simulation=i
    simulated_quantiles$jagsModel='Teu_known'
    all_simulated_quantiles = rbind(all_simulated_quantiles,simulated_quantiles)
    # Teu, Tem known, P Estimated
    data0$jagsdata$Tem = inputs$Tem
    simulated_quantiles = parameterEstimates(jagsdata=data0$jagsdata,jagsmodel=exittrap1,variable.names=variable.names)
    simulated_quantiles$simulation=i
    simulated_quantiles$jagsModel='Tem_known'
    all_simulated_quantiles = rbind(all_simulated_quantiles,simulated_quantiles)
    inputs_to_simulations = rbind(inputs_to_simulations,inputs)
  }  
  simresults = post_process_quantiles(file=all_simulated_quantiles)
  results=list(simresults=simresults,inputs_to_simulations=inputs_to_simulations,jagsmodel=exittrap1)
  return(results)}

# Create simulated exit trap dataset for analysis of survival and cycle duration
simulateExitTrapData = function(tau.b = 2.0, 
                                rg = 3,
                                resting = 2,
                                Teu = 3,
                                model='default',
                                Tem = 4,
                                Pm = 0.5,
                                rm = 3,
                                A0=A0,
                                P=M,
                                jagsdata=jagsdata,
                                requirePlot=FALSE){
  # discretise the normal kernel for the cycle duration
  pr = rep(NA,4)
  pr[1] = pnorm(1.5,resting,tau.b)
  pr[2] = pnorm(2.5,resting,tau.b) - pr[1]
  pr[3] = pnorm(3.5,resting,tau.b) - pr[1] - pr[2]
  pr[4] = 1 - pr[1] - pr[2] - pr[3]
  
  # generate simulated numbers of unfed per trap with autocorrelation
  U=jagsdata$unfed1/jagsdata$traps
  acf.obs = as.vector(acf(U,lag.max=1)[[1]][]) 
  # draw an iid sample extended by xtr extra dates so that a 'burn-in' period for the number of gravids can be included
  xtr=150
  T0_iid = sample(U,size=length(U)+xtr,replace=TRUE)
  T0a <- stats::filter(T0_iid, filter=acf.obs, circular = T)
  T0b = T0a*mean(U)/mean(T0a)
  sdratio = sd(U)/sd(T0b)
  T0x = mean(U) + (T0b - mean(U))*sdratio
  T0x = ifelse(T0x < 1, 1,T0x)
  emergence = em = eg = rep(NA,length(T0x))
  emergenceRate = 'unfed'
  if(emergenceRate == 'unfed'){
    # if unfeds measure emergence then numbers of gravid per trap depend on previous gravids as well as on unfeds
    # eg is the simulated expected number of gravids per trap
    # the values at the start of the burn-in are irrelevant, as are mosquitoes emerging prior to the burn-in (these have negligible effect on the main sequence)
    # so the sequence can be initialised at an arbitrary value
    for (t in 1:5){
      #expected gravid
      eg[t] <- 1
      emergence[t] <- T0x[t]
    }  
    for (t in 6:length(T0x)){
      emergence[t] <- T0x[t]
      # expected gravids in first cycle
      eg0 <- (T0x[t-1]*pr[1] + T0x[t-2]*pr[2] + T0x[t-3]*pr[3] + T0x[t-4]*pr[4])* Teu 
      # expected gravids in subsequent cycles
      eg1 <- eg[t-1]*pr[1] + eg[t-2]*pr[2] + eg[t-3]*pr[3] + eg[t-4]*pr[4]
      eg2 <- eg[t-2]*pr[1] + eg[t-3]*pr[2] + eg[t-4]*pr[3] + eg[t-5]*pr[4]
      eg[t] <- (eg0 + A0*eg1 + (1-A0)*eg2)*P
      eg[t] <- max(eg[t],1) # constrain the expectation to be positive
    }
  } else {
    # if unfeds include survivors from previous cycle then ignore previous gravids and Teu includes the survival effect
    eg[t] <- (T0x[t-1]*pr[1] + T0x[t-2]*pr[2] + T0x[t-3]*pr[3] + T0x[t-4]*pr[4])* Teu * P 
    eg[t] <- max(eg[t],1) # constrain the expectation to be positive
    
    #emergence calculated as the difference between total unfed and survivors from the previous cycle
    # - for the case where there is no delay from oviposition to host seeking
    e1 <- T0x[ptr[t]] - (T0x[t-1]*pr[1] + T0x[t-2]*pr[2] + T0x[t-3]*pr[3] + T0x[t-4]*pr[4])*P 
    # - for the case where there is a one night delay
    e2 <- T0x[ptr[t]] - (T0x[t-2]*pr[1] + T0x[t-3]*pr[2] + T0x[t-4]*pr[3] + T0x[t-5]*pr[4])*P  
    # overall average
    emergence[t] <- max(0.001,A0*e1 + (1-A0)*e2)
  }
  #expected males
  em[1] <- Tem*emergence[1] 
  for (t in 2:length(T0x)){
    em[t] <- Tem*emergence[t] + Pm * em[t-1]
  }
  # remove the 'burn-in'
  T0 = T0x[(xtr+1):length(T0x)]
  em = em[(xtr+1):length(T0x)]
  emergence = emergence[(xtr+1):length(T0x)]
  eg = eg[(xtr+1):length(T0x)]
  
  # Multiply by number of traps and add negative binomial variation to the simulated values
  traps = jagsdata$traps
  males = gravid = rep(NA,length(traps))
  for(t in 1:length(traps)){
    males[t] = rnbinom(n=1, size=rm, mu=traps[t]*em[t])
    gravid[t] = rnbinom(n=1, size=rg, mu=traps[t] * eg[t])
  }
  unfed1 = round(as.numeric(traps*T0),0)
  if(requirePlot) {
    plt= plotSimulatedData() 
  } else { plt = NULL }
  inputs = list(tau.b = tau.b, rg = rg, resting = resting, Teu = Teu, Tem = Tem, rm = rm, Pm = Pm, A0 = A0, P=P)
  jagsdata$unfed1= unfed1
  jagsdata$gravid= gravid
  # duplicate the vector of numbers of gravids to allow fitting of two different sets of parameters
  jagsdata$log.unfed= log(unfed1)
  jagsdata$males = males
  jagsdata$gravid1 = gravid
  jagsdata$males1 = males
  jagsdata$A0 = A0
  jagsdata$P = P
  simulation=list(inputs=inputs,jagsdata=jagsdata, plt=plt)
  return(simulation)}


# Function calls start here
# The analysis of exit trap data requires a valid workfile. 

if( !exists('workfile') )
{
  source('summarize_data.r')
}

# Creation of JAGS input dataset for time-series analysis of Exit trap data
source('rjags_wrapper.R')
source('postprocessing.r')
source('plotting.r')

variable.names = c("cycle","tau.b","Teu","resting","pr","Tem","Pm","P","p","emergence")

switch (option, 
    {
      # Estimate survival, cycle length and sac rate from exit-trap data overall 
      exittrap1 = exittrap_base
      ExitTrapAnalysis = analyseExitTrap(firstquarter=8,lastquarter=28) 
      resting_sample = rnorm(1000, mean=ExitTrapAnalysis$quantiles$mean[ExitTrapAnalysis$quantiles$Var == 'resting'], 
                             sd=ExitTrapAnalysis$quantiles$sd[ExitTrapAnalysis$quantiles$Var == 'resting'])
      theta_exit = resting_sample + (1-proportions$A0_sample)/proportions$A0_sample
      ExitTrapAnalysis$theta_quantiles =quantile(x = theta_exit, probs = c(0.025,0.5,0.975))
    },
    {
      exittrap1 = exittrap_base  
      SimulationsExitTraps = SimulationsExitTrapsFixedSampleSize(nsimulations=5)   #change to 100 
    },
    {
      exittrap1 = exittrap_base
      ConsistencyAnalysis = consistencyAnalysis(nsimulations=5)  #change to 50 
    },
    {
      # Estimate survival, cycle length and sac rate from exit-trap data by year
      exittrap1 = exittrap_base
      ExitTrapAnalysisByYear = analyseExitTrapByYear()
    },
    {
      # Estimate survival, cycle length and sac rate from exit-trap data by temperature 
      variable.names = c("resting","cycle","Teu","Tem","Pm","P","p")
      exittrap1 = exittrap_temp
      ExitTrapAnalysisByTemp = analyseExitTrapByTemp()
      variable.names = c("resting","cycle","Teu","Tem","Pm0","Pm1","P0","P1")
      exittrap1 = exittrap_temptrend
      ExitTrapAnalysisByTempTrend = analyseExitTrapByTempTrend()
    }
)

# remove functions and variables that are required only in this script
suppressWarnings(
remove(list=c("analyseExitTrap","resting_sample","analyseExitTrap_allmodels","exittrap_base",
              "exittrap_temp","exittrap_temptrend","exittrap1", "option", "parameterEstimates", 
              "theta_exit", "variable.names","analyseExitTrapByTemp",
              "analyseExitTrapByTempTrend","analyseExitTrapByYear","consistencyAnalysis",
              "extract_index","extract_quantiles","extract_varnames","post_process_quantiles",
              "simulateExitTrapData","simulationsExitTraps",
              "SimulationsExitTrapsFixedSampleSize","summarise_by_method","summarise_by_quarter",
              "get_mp_recomCol","plotConsistency","plotCorrelations","plotEstimatesByTemp","errbarPlot",
              "plotEstimatesByTime","plotSimulatedData","plotSimulations_vs_Inputs")))

