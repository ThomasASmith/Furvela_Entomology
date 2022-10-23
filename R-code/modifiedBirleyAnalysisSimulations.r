####### Use simulations to evaluate performance of modified Birley model ##################

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

####### CALLING SCRIPT

BirleyEmergenceRates = BirleyResults$X50.[BirleyResults$Var == 'emergence']
BirleyModelSimulations=SimulationsBirleyModel(nsimulations=20,
                                              BirleyEmergenceRates = BirleyEmergenceRates)

