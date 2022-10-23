


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

