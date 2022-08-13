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

# Estimate survival, cycle length and sac rate from exit-trap data overall 
option = 1
source('ExitTrapTimeSeriesAnalysis.r')

# Evaluate bias and precision of estimates from time-series of exit traps using simulations
option = 2
source('ExitTrapTimeSeriesAnalysis.r')

# Evaluate consistency of estimates from time-series of exit traps 
option = 3
source('ExitTrapTimeSeriesAnalysis.r')

# Estimate survival, cycle length and sac rate from exit-trap data by year
option = 4
source('ExitTrapTimeSeriesAnalysis.r')

# Estimate survival, cycle length and sac rate from exit-trap data by temperature
option = 5
source('ExitTrapTimeSeriesAnalysis.r')



CorrelationPlot=plotCorrelations()
remove('plotCorrelations')

source('temp.r')
source('JAGS_code.R')



# Use simulations to evaluate performance of modified Birley model 
BirleyEmergenceRates = BirleyResults$X50.[BirleyResults$Var == 'emergence']
BirleyModelSimulations=SimulationsBirleyModel(nsimulations=20,
                                              BirleyEmergenceRates = BirleyEmergenceRates)


  
  load("C:/git_repos/charlwood/CBMC/fullenvironment.RData")

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
  savePlot(ExitTrapAnalysisByTemp$plt1,'ExitTrapAnalysisByTemp.png',vertical_panels=2)
  savePlot(BiasPlots$plt1,'ExitTrapSimulationsFemales.png',vertical_panels=2)
  savePlot(BiasPlots$plt2,'ExitTrapSimulationsMales.png',vertical_panels=1)
  savePlot(ConsistencyPlots$plt1,'ExitTrapConsistencyFemales.png',vertical_panels=2)
  savePlot(ConsistencyPlots$plt2,'ExitTrapConsistencyMales.png',vertical_panels=1)
}

save.image("C:/git_repos/charlwood/CBMC/fullenvironment.RData") 


