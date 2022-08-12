##### TIME SERIES ANALYSIS OF FURVELA ENTOMOLOGY ################

requireTablesPlots = FALSE
requireResimulation = FALSE


################ MAIN PROGRAM STARTS HERE BY READING IN DATA #######################
# if necessary set working directory
setwd("../Furvela_Entomology")

# Creation of data summaries and workfile (essential for subsequent analysis)
source('summarize_data.r')

# Analysis of resting collections (equations # in Charlwood et al,....)
source('resting_collection_analysis.r')

# Analysis via rootsolving (equations # in Charlwood et al,....)
source('rootsolving.r')

# Analysis of time series of light trap data (equations # in Charlwood et al,....)
source('modifiedBirleyAnalysis.r')

source('plotting.r')


source('temp.r')
source('JAGS_code.R')



# Use simulations to evaluate performance of modified Birley model 
BirleyEmergenceRates = BirleyResults$X50.[BirleyResults$Var == 'emergence']
BirleyModelSimulations=SimulationsBirleyModel(nsimulations=20,
                                              BirleyEmergenceRates = BirleyEmergenceRates)


############### Analysis of Exit trap data
# Use simulations to evaluate performance of models for exit traps 

  exittrap1 = createListOfModels()$H
  variable.names = c("cycle","tau.b","Teu","resting","pr","Tem","Pm","P","p","emergence")
  SimulationsExitTraps = SimulationsExitTrapsFixedSampleSize(nsimulations=100)
  ConsistencyAnalysis = consistencyAnalysis(nsimulations=50)
  
  # Estimate survival, cycle length and sac rate from exit-trap data overall and by year

  ExitTrapAnalysis = analyseExitTrap(firstquarter=8,lastquarter=28) 
  resting_sample = rnorm(1000, mean=ExitTrapAnalysis$quantiles$mean[ExitTrapAnalysis$quantiles$Var == 'resting'], 
                         sd=ExitTrapAnalysis$quantiles$sd[ExitTrapAnalysis$quantiles$Var == 'resting'])
  theta_exit = resting_sample + (1-proportions$A0_sample)/proportions$A0_sample
  ExitTrapAnalysis$theta_quantiles =quantile(x = theta_exit, probs = c(0.025,0.5,0.975))

  ExitTrapAnalysisByYear = analyseExitTrapByYear()

  # Analysis by temperature  
  variable.names = c("resting","cycle","Teu","Tem","Pm","P","p")
  exittrap1 = createListOfModels()$P
  ExitTrapAnalysisByTemp = analyseExitTrapByTemp()
  savePlot(ExitTrapAnalysisByTemp$plt1,'ExitTrapAnalysisByTemp.png',vertical_panels=2)
  variable.names = c("resting","cycle","Teu","Tem","Pm0","Pm1","P0","P1")
  exittrap1 = createListOfModels()$Q
  ExitTrapAnalysisByTempTrend = analyseExitTrapByTempTrend()

} else {   
  load("C:/git_repos/charlwood/CBMC/fullenvironment.RData")
}

# To remove easily reproducible objects
# rm(list= ls()[! (ls() %in% c('SimulationsExitTraps','ConsistencyAnalysis',
#                              'BirleyResults','BirleyModelSimulations',
#                              'ExitTrapAnalysis','ExitTrapAnalysisByYear',
#                              'ExitTrapAnalysisByTemp','ExitTrapAnalysisByTempTrend'))])

if(requireTablesPlots){
  BiasPlots = createTables_Plots(input=SimulationsExitTraps,plottype='bias')
  ConsistencyPlots = createTables_Plots(input=ConsistencyAnalysis,plottype='consistency')
  savePlot(ExitTrapAnalysisByYear$plt1,'ExitTrapAnalysisByYearFemales.png',vertical_panels=2)
  savePlot(ExitTrapAnalysisByYear$plt2,'ExitTrapAnalysisByYearMales.png',vertical_panels=1)
  savePlot(BiasPlots$plt1,'ExitTrapSimulationsFemales.png',vertical_panels=2)
  savePlot(BiasPlots$plt2,'ExitTrapSimulationsMales.png',vertical_panels=1)
  savePlot(ConsistencyPlots$plt1,'ExitTrapConsistencyFemales.png',vertical_panels=2)
  savePlot(ConsistencyPlots$plt2,'ExitTrapConsistencyMales.png',vertical_panels=1)
}

save.image("C:/git_repos/charlwood/CBMC/fullenvironment.RData") 


