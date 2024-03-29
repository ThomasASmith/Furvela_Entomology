##### TIME SERIES ANALYSIS OF FURVELA ENTOMOLOGY ################

requireResimulation <- FALSE

################ MAIN PROGRAM STARTS HERE BY READING IN DATA #####################
# if necessary set working directory
setwd("../Furvela_Entomology")

############## DATA DESCRIPTION #############
# Creation of data summaries and workfile (essential for subsequent analysis)
source("R-code/summarize_data.r")

############## DATA ANALYSIS #############
# To load pre-existing analysis results
load("timeseriesAnalyses.RData")

# To carry out new analyses First replace the pre-existing data with the data to be analysed If you need to run
# new analyses make sure you have JAGS correctly installed and that rjags can find the JAGS installation on your
# computer

# Analysis of resting collections (equations # in Charlwood et al,....)
source("R-code/resting_collection_analysis.r")

# Analysis via rootsolving (equations # in Charlwood et al,....)
source("R-code/rootsolving.r")

# Analysis of time series of light trap data (equations # in Charlwood et al,....)
source("R-code/modifiedBirleyAnalysis.r")

# Estimate survival, cycle length and sac rate from exit-trap data overall
option <- 1
source("R-code/ExitTrapTimeSeriesAnalysis.r")

# Estimate survival, cycle length and sac rate from exit-trap data by year
option <- 4
source("R-code/ExitTrapTimeSeriesAnalysis.r")

# Estimate survival, cycle length and sac rate from exit-trap data by temperature
option <- 5
source("R-code/ExitTrapTimeSeriesAnalysis.r")

############## EVALUATION OF ESTIMATORS USING SIMULATIONS #####
# Evaluate bias and precision of estimates from time-series using simulations

# Use simulations to evaluate performance of modified Birley model
source("R-code/modifiedBirleyAnalysisSimulations.r")

# Use simulations to evaluate performance of time series analysis of exit traps
option <- 2
source("R-code/ExitTrapTimeSeriesAnalysis.r")

# Evaluate consistency of estimates from time-series of exit traps using simulations
option <- 3
source("R-code/ExitTrapTimeSeriesAnalysis.r")

############ PLOTTING OF RESULTS #############
requirePlotoutput <- TRUE
source("R-code/plotting.r")

# To remove easily reproducible objects
rm(list = ls()[!(ls() %in% c("SimulationsExitTraps", "ConsistencyAnalysis", "BirleyAnalysis", "BirleyModelSimulations",
  "ExitTrapAnalysis", "ExitTrapAnalysisByYear", "ExitTrapAnalysisByTemp", "ExitTrapAnalysisByTempTrend"))])


