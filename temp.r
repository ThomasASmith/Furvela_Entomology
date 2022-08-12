


# Create simulated light trap dataset for analysis of survival and cycle duration
simulateBirleyData = function(tau.b = 0.75,
                              r = 3,
                              cycle = 4,
                              P = 0.6){
  # discretise the normal kernel for the cycle duration
  pr = rep(NA,3)
  pr[1] = pnorm(2.5,cycle,tau.b)
  pr[2] = pnorm(3.5,cycle,tau.b) - pr[1]
  pr[3] = 1 - pr[1] - pr[2]
  # 
  
  EX = mean(BirleyEmergenceRates)
  VarX = var(BirleyEmergenceRates)
  # extract autocorrelation function 
  acf.obs = as.vector(acf(BirleyEmergenceRates)[[1]][]) 
  muLogN = log(EX*EX/sqrt(VarX + EX*EX))
  sigmaLogN = sqrt(log(VarX/(EX*EX) + 1))
  emergence_iid = rlnorm(513, meanlog = muLogN, sdlog = sigmaLogN)
  # introduce autocorrelation
  emergence <- stats::filter(emergence_iid, filter=acf.obs, circular = T)
  # set negative values to zero
  emergence[emergence < 0] = 0
  
  meanFemales = rep(NA,513)
  for (t in 1:513){
    meanFemales[MBData$ptrs$ptr[t]]= emergence[MBData$ptrs$ptr[t]]
    + (emergence[MBData$ptrs$lag2[t]] *pr[1]
       + emergence[MBData$ptrs$lag3[t]] *pr[2]
       + emergence[MBData$ptrs$lag4[t]] *pr[3])*P
  }
  # trap rare near-zero females and replace with 1
  meanFemales[meanFemales < 1] = 1
  
  T0=T_1=T_2=T_3=T_4=M=eParous=parous=rep(NA,87)
  jagsdata = MBData$jagsdata
  for(t in 1:87){
    T0[t] = meanFemales[jagsdata$ptr[t]]
    T_1[t] = meanFemales[jagsdata$lag1[t]]
    T_2[t] = meanFemales[jagsdata$lag2[t]]
    T_3[t] = meanFemales[jagsdata$lag3[t]]
    T_4[t] = meanFemales[jagsdata$lag4[t]]
    # expected number parous per trap
    M[t] = (T_2[t]*pr[1] + T_3[t]*pr[2] + T_4[t]*pr[3])*P
    # expected total parous among those tested
    eParous[t] = M[t]*jagsdata$dissected[t]/T0[t]
    parous[t] = rnbinom(n=1, size=r, mu=eParous[t])
  }
  simulatedData = jagsdata
  simulatedData$meanFemalesRecorded=meanFemales[1:296]
  simulatedData$parous = parous
  return(simulatedData)}



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



################ ANALYSIS OF SIMULATED DATASETS ########################
SimulationsBirleyModel = function(nsimulations,BirleyEmergenceRates){
  all_simulated_quantiles = NULL
  inputs_to_simulations = NULL
  variable.names = c("tau.b","r","cycle","P")
  for(i in 1:nsimulations){
    # parameters that might be recoverable
    tau.b = runif(1, min = 0.5, max = 1)
    r = runif(1, min = 1, max = 5)
    cycle = runif(1, min = 1.9, max = 4.1)
    P = runif(1, min = 0, max = 1)
    inputs = data.frame(tau.b = tau.b,r = r,cycle = cycle,P = P)
    inputs_to_simulations = rbind(inputs_to_simulations,inputs)
    simulatedData=simulateBirleyData(tau.b = tau.b,r = r,cycle = cycle,P = P)
    simulated_quantiles = parameterEstimates(jagsdata= simulatedData,jagsmodel=modifiedBirleyModel,variable.names=variable.names)
    simulated_quantiles$simulation=i
    simulated_quantiles$jagsModel='modifiedBirleyModel'
    all_simulated_quantiles = rbind(all_simulated_quantiles,simulated_quantiles)
  }
  simresults = post_process_quantiles(file=all_simulated_quantiles)
  toPlot6a = data.frame(simresults$parameters[simresults$parameters$Var=='P',],
                        input=inputs_to_simulations$P)
  plt1=plotSimulations_vs_Inputs(data= toPlot6a, textLabel='survival per cycle', nmodels=1, xlim= c(0,1))
  toPlot6b = data.frame(simresults$parameters[simresults$parameters$Var=='cycle',],
                        input=inputs_to_simulations$cycle)
  plt2=plotSimulations_vs_Inputs(data= toPlot6b, textLabel='cycle length (days)', nmodels=1, xlim= c(1.5,4.5))
  library(cowplot)
  plt = plot_grid(plt1, plt2, labels = c('A', 'B'), label_size = 12,ncol=2)
  results=list(simresults=simresults,inputs_to_simulations=inputs_to_simulations,plt=plt)
  return(results)}

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

createTables_Plots = function(input,plottype='bias'){  
  # Plots of simulation results for females
  simresults=input$simresults
  inputs_to_simulations=input$inputs_to_simulations
  
  plt1 = NULL
  plt2 = NULL
  CCCTable = NULL
  library(tidyr)
  library(qwraps2)
  
  
  calculateBias = function(df,Var){ 
    df1 = df %>% 
      group_by(jagsModel) %>%  
      summarise(N=n(),
                mean.ci = list(mean_ci((X50.-input)/input))) %>% 
      unnest_wider(mean.ci)
    df1$Var=Var
    return(df1)}
  
  resultsFemales = simresults$parameters[which(simresults$parameters$jagsModel %in% 
                                                 c('exittrap1','P_075','A_03','P_est','Teu_known')),]
  nmodels=5
  resultsFemales$jagsModel[resultsFemales$jagsModel=='exittrap1']= 'a'
  resultsFemales$jagsModel[resultsFemales$jagsModel=='P_est']= 'b'
  resultsFemales$jagsModel[resultsFemales$jagsModel=='Teu_known']= 'c'
  resultsFemales$jagsModel[resultsFemales$jagsModel=='P_075']= 'd'
  resultsFemales$jagsModel[resultsFemales$jagsModel=='A_03']= 'e'
  toPlot1a = data.frame(resultsFemales[resultsFemales$Var=='Teu',],
                        input=rep(inputs_to_simulations$Teu,each=nmodels),
                        days_with_complete_data=rep(inputs_to_simulations$days_with_complete_data,each=nmodels))
  
  toPlot1b = data.frame(resultsFemales[resultsFemales$Var=='P',],
                        input=rep(inputs_to_simulations$P,each=nmodels),
                        days_with_complete_data=rep(inputs_to_simulations$days_with_complete_data,each=nmodels))
  toPlot1c = data.frame(resultsFemales[resultsFemales$Var=='resting',],
                        input=rep(inputs_to_simulations$resting,each=nmodels),
                        days_with_complete_data=rep(inputs_to_simulations$days_with_complete_data,each=nmodels))
  
  # Plots of simulation results for males
  resultsMales = simresults$parameters[which(simresults$parameters$jagsModel %in% 
                                               c('Teu_known','Tem_known')),]
  nmodels=2
  resultsMales$jagsModel[resultsMales$jagsModel=='Teu_known']= 'a'
  resultsMales$jagsModel[resultsMales$jagsModel=='Tem_known']= 'b'
  toPlot2a = data.frame(resultsMales[resultsMales$Var=='Pm',], input=rep(inputs_to_simulations$Pm,each=nmodels),
                        days_with_complete_data=rep(inputs_to_simulations$days_with_complete_data,each=nmodels))
  toPlot2b = data.frame(resultsMales[resultsMales$Var=='Tem',],input=rep(inputs_to_simulations$Tem,each=nmodels),
                        days_with_complete_data=rep(inputs_to_simulations$days_with_complete_data,each=nmodels))
  bias = calculateBias(df=toPlot1a,Var='Teu')
  bias = rbind(bias, calculateBias(df=toPlot1b,Var='P'))
  bias = rbind(bias, calculateBias(df=toPlot1c,Var='resting'))
  bias = rbind(bias, calculateBias(df=toPlot2a,Var='Pm'))
  bias = rbind(bias, calculateBias(df=toPlot2b,Var='Tem'))
  library(ggplot2)  
  library(cowplot)
  if (plottype == 'bias'){
    p1a=plotSimulations_vs_Inputs(data= toPlot1a, textLabel='trapping efficiency: gravids', xlim= c(0,8),requirelegend=TRUE)
    p1b=plotSimulations_vs_Inputs(data= toPlot1b, textLabel='survival per cycle: females', xlim= c(0,1))
    p1c=plotSimulations_vs_Inputs(data= toPlot1c, textLabel='resting period (days)', xlim= c(1,5))
    toPlot1c$Var = 'cycle'
    toPlot1c$X50. = toPlot1c$X50. + 1 - rep(inputs_to_simulations$A0,each=5)
    toPlot1c$input = toPlot1c$input + 1 - rep(inputs_to_simulations$A0,each=5)
    p1d = plotSimulations_vs_Inputs(data= toPlot1c, textLabel='cycle duration (days)', xlim= c(1,5))
    plt1 = plot_grid(p1a, p1b, p1c, p1d, labels = c('A', 'B', 'C', 'D'), label_size = 12,ncol=2)
    p2a=plotSimulations_vs_Inputs(data= toPlot2a, textLabel='daily survival: males',nmodels=2, xlim= c(0,1))
    p2b=plotSimulations_vs_Inputs(data= toPlot2b, textLabel='trapping efficiency: males',nmodels=2, xlim= c(1,10),requirelegend=TRUE)
    plt2 = plot_grid(p2a, p2b, labels = c('A', 'B'), label_size = 12,ncol=2)
  }
  if (plottype == 'consistency'){
    p1a=plotConsistency(data=toPlot1a, textLabel='trapping efficiency (gravids)')
    p1b=plotConsistency(data=toPlot1b, textLabel='survival per cycle: females')
    p1c=plotConsistency(data=toPlot1c, textLabel='resting period (days)',requirelegend=TRUE)
    plt1 = plot_grid(p1a, p1b, p1c, labels = c('A', 'B', 'C'), label_size = 12,ncol=2)
    p2a=plotConsistency(data=toPlot2a, textLabel='daily survival: males',nmodels=2,requirelegend=TRUE)
    p2b=plotConsistency(data=toPlot2b, textLabel='trapping efficiency: males',nmodels=2)
    plt2 = plot_grid(p2a, p2b, labels = c('A', 'B'), label_size = 12,ncol=2)
  }
  
  # calculation of concordance correlation coefficients
  CCCTable = function(inputdf){
    library(epiR)    
    df = data.frame(Var=c(),model=c(),CCC=c(),lower=c(),upper=c())
    # use only variables that are found in the input data frame
    for(Var in levels(as.factor(as.character(inputdf$Var)))){
      for(model in levels(as.factor((inputdf$jagsModel)))){
        dfsub = inputdf[inputdf$Var==Var & inputdf$jagsModel==model,]
        rho.c = with(dfsub,epi.ccc(X50., input, ci = "z-transform", conf.level = 0.95, rep.measure = FALSE))$rho.c
        df1 = data.frame(Var=Var,model=model,CCC=rho.c[1],lower=rho.c[2],upper=rho.c[3])
        df=rbind(df,df1)
      }  
    }
    return(df)
  }
  CCCTable = CCCTable(rbind(toPlot1a,toPlot1b,toPlot1c,toPlot2a,toPlot2b))
  results=list(plt1=plt1,plt2=plt2,CCC=CCCTable,bias=bias)
  return(results)}

# saving plots
library(grid)
savePlot <- function(plot,Plotname,vertical_panels=2){
  print(Plotname)
  grid.newpage()
  png(Plotname,width=18.5,height=18.5*vertical_panels/2,units="cm",res=600)
  grid.draw(plot)
  dev.off()
}
