scale_colour_manual(values=get_mp_recomCol(),name='Model') +
geom_errorbar(aes(ymin=X2.5.,
ymax=X97.5., width = 0.5, colour = jagsModel), show.legend=FALSE) +
scale_fill_manual(values=get_mp_recomCol(), guide = 'none') +
geom_point(aes(y=X50., colour = jagsModel), shape=21, size=3, stroke = 1.5, show.legend=requirelegend) +
theme(legend.position = c(.8,.75))+
scale_y_continuous(name = paramLabel) +
scale_x_continuous(name = 'Year',limits=c(8,32), breaks=c(8.5,12.5,16.5,20.5,24.5,28.5),labels=ticklabels) +
geom_vline(xintercept=22.5)
return(plt)}
plotSimulations_vs_Inputs = function(data, textLabel,xlim, requirelegend= FALSE){
plt= ggplot(data = data,aes(x = input, y=X50., group = jagsModel)) + theme_bw(base_size = 9) +
scale_colour_manual(values=get_mp_recomCol(),name='Model') +
#geom_errorbar(aes(ymin=X2.5.,
#                  ymax=X97.5., width = 0.02, colour = jagsModel), show.legend=FALSE) +
theme(legend.position = c(.2,.75))+
scale_fill_manual(values=get_mp_recomCol(), guide = 'none') +
geom_point(aes(y=X50.,colour = jagsModel) , shape = 21, size=3, stroke = 1.5, show.legend=requirelegend) +
geom_line(aes(x=input,y=input)) +
scale_x_continuous(name = paste('Input',textLabel), limits=xlim) +
scale_y_continuous(name = paste('Modelled',textLabel), limits = xlim)
return(plt)}
modifiedBirleyModel = "
model{
for(t in 1:87){
parous[t] ~ dnegbin(ppar[t],r)
T0[t] <- meanFemales[ptr[t]]
T_1[t] <- meanFemales[lag1[t]]
T_2[t] <- meanFemales[lag2[t]]
T_3[t] <- meanFemales[lag3[t]]
T_4[t] <- meanFemales[lag4[t]]
# expected number parous per trap
M[t] <- (T_2[t]*pr[1] + T_3[t]*pr[2] + T_4[t]*pr[3])*P
# expected total parous among those tested
eParous[t] <- M[t]*dissected[t]/T0[t]
ppar[t] <- r/(eParous[t] + r)  # Negative binomial p
# emergence rate (scaled to light trap collections) n.b. this may need to be constrained to be strictly positive
emergence[t] <- T0[t] - M[t]
}
for(i in 1:296){
meanFemales[i] <- meanFemalesRecorded[i]
logMeanFemales[i] <- log(meanFemales[i])
}
for(i in 297:513){
variation[i] ~ dnorm(0,tau.t)
meanFemales[i] <- exp(logtbar+variation[i])
}
logtbar <- mean(logMeanFemales)
# discretise the normal kernel for the cycle duration
pr[1] <- pnorm(2.5,cycle,tau.b)
pr[2] <- pnorm(3.5,cycle,tau.b) - pr[1]
pr[3] <- 1 - pr[1] - pr[2]
# priors
tau.t ~ dgamma(1,1)
sigma2 <- 1/tau.t
tau.b ~ dgamma(1,1)
r ~ dgamma(1,1)
cycle ~ dunif(1.5,4.5)
P ~dunif(0,1)
p <- exp(log(P)/cycle)
}
"
# Estimation of resting duration from exit traps: no imputation of missing data
exittrapE ="
model{
for(t in 1:days_with_complete_data){
# expected number gravid per trap (Teu is relative trapping efficiency of gravid vs unfed)
G[t] <- (T0[lag1[t]]*pr[1] + T0[lag2[t]]*pr[2] + T0[lag3[t]]*pr[3] + T0[lag4[t]]*pr[4]) * Teu
pg[t] <- rg/(G[t]*traps[ptr[t]] + rg)  # Negative binomial p
gravid[t] ~ dnegbin(pg[t],rg)
}
for(d in 1:sampleddates){
# numbers of unfed (+ part fed)
T0[d] <- unfed1[d]/traps[d]
}
################# PRIORS ########################
# discretised normal kernel for the cycle duration
pr[1] <- pnorm(1.5,resting,tau.b)
pr[2] <- pnorm(2.5,resting,tau.b) - pr[1]
pr[3] <- pnorm(3.5,resting,tau.b) - pr[1] - pr[2]
pr[4] <- 1 - pr[1] - pr[2] - pr[3]
tau.b ~ dgamma(1,1)    # precision of distribution of cycle duration
rg ~ dgamma(1,1)       # Negative binomial overdispersion for gravid females
resting ~ dunif(0.5,4.5) # Oviposition interval/ duration of resting period
Teu ~ dgamma(1,1)       # relative trapping efficiency of gravids relative to unfeds
}
"
# Estimation of emergence rates and male survival from exit trap data with parous and sac rates known: no imputation of missing data
exittrapF ="
model{
for(t in 1:days_with_complete_data){
# expected number gravid per trap (Teu is relative trapping efficiency of gravid vs unfed)
G[t] <- (T0[lag1[t]]*pr[1] + T0[lag2[t]]*pr[2] + T0[lag3[t]]*pr[3] + T0[lag4[t]]*pr[4]) * Teu
pg[t] <- rg/(G[t]*traps[ptr[t]] + rg)  # Negative binomial p
gravid[t] ~ dnegbin(pg[t],rg)
#emergence calculated as the difference between total unfed and survivors from the previous cycle
# - for the case where there is no delay from oviposition to host seeking
E0[t] <- T0[ptr[t]] - (T0[lag1[t]]*pr[1] + T0[lag2[t]]*pr[2] + T0[lag3[t]]*pr[3] + T0[lag4[t]]*pr[4])*P
# - for the case where there is a one night delay
E1[t] <- T0[ptr[t]] - (T0[lag2[t]]*pr[1] + T0[lag3[t]]*pr[2] + T0[lag4[t]]*pr[3] + T0[lag5[t]]*pr[4])*P
# overall average
emergence[t] <- max(0.001,A0*E0[t] + (1-A0)*E1[t])
# expected number of males:
#     Tem is relative trapping efficiency of males vs unfed assuming 50:50 sex ratio on emergence
#     Pm is daily survival of males
# Expected males per trap
Emales[t] <- (Tem*emergence[t] + Pm * males[lag1[t]]/traps[lag1[t]])*traps[ptr[t]]
pm[t] <- rm/(Emales[t] + rm)  # Negative binomial p
males[ptr[t]] ~ dnegbin(pm[t],rm)
}
for(d in 1:sampleddates){
# numbers of unfed (+ part fed)
T0[d] <- unfed1[d]/traps[d]
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
Teu ~ dgamma(1,1)       # relative trapping efficiency of gravids relative to unfeds
Tem ~ dgamma(1,1)       # relative trapping efficiency of males vs unfed assuming 50:50 sex ratio on emergence
Pm ~ dunif(0,1)         # daily survival of males
}
"
# Estimation of emergence rates and male survival from exit trap data with sac rates known: imputation of P
exittrapG ="
model{
for(t in 1:days_with_complete_data){
# expected number gravid per trap (Teu is relative trapping efficiency of gravid vs unfed)
G[t] <- (T0[lag1[t]]*pr[1] + T0[lag2[t]]*pr[2] + T0[lag3[t]]*pr[3] + T0[lag4[t]]*pr[4]) * Teu
pg[t] <- rg/(G[t]*traps[ptr[t]] + rg)  # Negative binomial p
gravid[t] ~ dnegbin(pg[t],rg)
#emergence calculated as the difference between total unfed and survivors from the previous cycle
# - for the case where there is no delay from oviposition to host seeking
E0[t] <- T0[ptr[t]] - (T0[lag1[t]]*pr[1] + T0[lag2[t]]*pr[2] + T0[lag3[t]]*pr[3] + T0[lag4[t]]*pr[4])*P
# - for the case where there is a one night delay
E1[t] <- T0[ptr[t]] - (T0[lag2[t]]*pr[1] + T0[lag3[t]]*pr[2] + T0[lag4[t]]*pr[3] + T0[lag5[t]]*pr[4])*P
# overall average
emergence[t] <- max(0.001,A0*E0[t] + (1-A0)*E1[t])
# expected number of males:
#     Tem is relative trapping efficiency of males vs unfed assuming 50:50 sex ratio on emergence
#     Pm is daily survival of males
# Expected males per trap
Emales[t] <- (Tem*emergence[t] + Pm * males[lag1[t]]/traps[lag1[t]])
pm[t] <- rm/(Emales[t]*traps[ptr[t]] + rm)  # Negative binomial p
males[ptr[t]] ~ dnegbin(pm[t],rm)
}
for(d in 1:sampleddates){
# numbers of unfed (+ part fed)
T0[d] <- unfed1[d]/traps[d]
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
Teu ~ dgamma(1,1)       # relative trapping efficiency of gravids relative to unfeds
Tem ~ dgamma(1,1)       # relative trapping efficiency of males vs unfed assuming 50:50 sex ratio on emergence
Pm ~ dunif(0,1)         # daily survival of males
P ~ dunif(0,1)         # overall survival through oviposition cycle
}
"
library(rjags)
extract_quantiles = function(jags.samples.out=jags.samples.out){
quantiles = data.frame(summary(jags.samples.out)$quantiles)
statistics = data.frame(summary(jags.samples.out)$statistics)
quantiles$Mean = statistics$Mean
quantiles$SD = statistics$SD
return(quantiles)
}
parameterEstimates = function(jagsdata=jagsdata,jagsmodel){
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
################ DATA DESCRIPTION ################################
data_description = function(){
drops = c("Casa","Data","data","Sheet","ISO_Week","Tempo_aqui","Nome_de_fora","Cumprimento","Largura","Preco","period_lado","Name",
"Observacoes","Age","preco..mil.","PIN","Ano_Nasc","Data_Nasc","number","plate")
keeps = c("Tempo_aqui","Cumprimento","Largura","Preco","Age","Ano_Nasc")
animais_summary = mapply(summary,animais[,!(names(animais) %in% drops)])
casa_table      = mapply(table, casa[,!(names(casa) %in% drops)])
casa_summary    = mapply(summary, casa[,names(casa) %in% keeps])
preco_summary   = summary(as.numeric(casa$Preco[casa$Preco != 'variavel']))
morte_table     = mapply(table, morte[,!(names(morte) %in% drops)])
mosquiteiro_summary = mapply(summary, mosquiteiro[,!(names(mosquiteiro) %in% drops)])
collection_table = mapply(table, mosquito_collection[,!(names(mosquito_collection) %in% drops)])
collection_summary = mapply(summary, mosquito_collection[,!(names(mosquito_collection) %in% drops)])
pessao_table = mapply(table, pessao[,!(names(pessao) %in% drops)])
pessao_summary = mapply(summary, pessao[,names(pessao) %in% keeps])
sporozoites_table =mapply(table, sporozoites[,!(names(sporozoites) %in% drops)])
description = list(
animais_summary	= animais_summary,
casa_table	= casa_table,
casa_summary	= casa_summary,
preco_summary	= preco_summary,
morte_table	= morte_table,
mosquiteiro_summary	= mosquiteiro_summary,
collection_table	= collection_table,
collection_summary	= collection_summary,
pessao_table	= pessao_table,
pessao_summary	= pessao_summary,
sporozoites_table	= sporozoites_table)
return(description)
}
################# Interval estimates of proportions #####################
analysis_of_proportions = function(){
parous = sum(dissections$Sac + dissections$No.sac)
dissected = sum(dissections$total_dissected)
M = parous/dissected
M_interval = prop.test(x=parous, n=dissected, conf.level=.95, correct=FALSE)
sac = sum(dissections$Sac)
A0 = sac/parous
A0_sample = rbinom(1000,parous,A0)/parous
A0_interval = prop.test(x=sac, n=parous, conf.level=.95, correct=FALSE)
proportions = list(M=M, M_interval=M_interval, A0=A0, A0_sample=A0_sample, A0_interval=A0_interval)
return(proportions)}
################## Creation of dataset for extension of Birley model
createModifiedBirleyData = function(){
library(dplyr)
dissDate = as.numeric(as.Date(dissections[,1], format = "%d-%b-%Y"))+716827
df4 = data.frame(cbind(date=dissDate, dissected=dissections$total_dissected,parous = dissections$Sac + dissections$No.sac))
df5 = data.frame(df4 %>%
group_by(date) %>%
dplyr::summarize(parous = sum(parous, na.rm=TRUE),
total = sum(dissected, na.rm=TRUE),
n = n()))
# fields lag1-lag4 were added manually to df5 to give the jags input
jagsDataB = read.csv(file='jags_input_for_Birley_method.csv')
df9 = data.frame(read.csv(file='Moonlight collections dissections Furvela.csv'))
df9$date = as.numeric(as.Date(as.character(df9[,1]), format = "%d-%b-%Y"))-13658
# multiply numbers captured by 3 for timed traps (which are for 1/3 of night)
df9$Af_tot_fem[df9$Time > 0]= df9$Af_tot_fem[df9$Time > 0]*3
df10 = data.frame(df9 %>%
group_by(date) %>%
dplyr::summarize(meanAf.female = mean(Af_tot_fem, na.rm=TRUE)))
df8 = data.frame(with(df0, df0[Collection=='Light' & date > 2155 & date <= (2159+509),] %>%
group_by(date) %>%
dplyr::summarize(meanAf.female = mean(Af.unfed+Af.part+Af.fed+Af.semi+Af_gravid, na.rm=TRUE))))
df8$date = df8$date-2159
df11= rbind(df8,df10)
df12= merge(df11,data.frame(date=seq(1:509)),by='date',all=TRUE)
# create a pointer for finding lagged numbers caught
i1=0
i2=length(which(!is.na(df12$meanAf.female)))
df12$ptr = NA
for(i in 1:nrow(df12)){
if(!is.na(df12$meanAf.female[i])) {
i1=i1+1
df12$ptr[i] = i1
} else {
i2=i2+1
df12$ptr[i] = i2
}
}
df12$lag1=c(NA,df12$ptr[1:512])
df12$lag2=c(NA,df12$lag1[1:512])
df12$lag3=c(NA,df12$lag2[1:512])
df12$lag4=c(NA,df12$lag3[1:512])
df13 = merge(df12,df5,by='date',all.y=TRUE)
jagsdata = list(meanFemalesRecorded=df11$meanAf.female,
dissected =df13$total,
parous=df13$parous,
ptr = df13$ptr,
lag1 = df13$lag1,
lag2 = df13$lag2,
lag3 = df13$lag3,
lag4 = df13$lag4)
MBData= list(jagsdata=jagsdata,ptrs=df12)
return(MBData)
}
################# Analysis of resting collections ################
resting_collection_analysis = function(){
fdenom = 6379+791+5121
f = prop.test(x=6379, n=fdenom, conf.level=.95, correct=FALSE)$estimate
f_sample = rbinom(1000,fdenom,f)/fdenom
f_quantiles = quantile(x = f_sample, probs = c(0.025,0.5,0.975))
# duration of the resting period as an explicit function of f
# using P from the Birley model and a normal approximation (sd of P is 0.03770995)
P = 0.528
P_sample = rnorm(1000, mean=P, sd=0.03770995)
theta_prime = log(P_sample)/log(1-f_sample*(1-P_sample))
theta_prime_quantiles = quantile(x = theta_prime, probs = c(0.025,0.5,0.975))
theta = log(P_sample)/log(1-f_sample*(1-P_sample)) + (1-proportions$A0_sample)/proportions$A0_sample
theta_quantiles =quantile(x = theta, probs = c(0.025,0.5,0.975))
results = list(f_quantiles=f_quantiles,
theta_prime_quantiles = theta_prime_quantiles,
theta_quantiles = theta_quantiles)
return(results)}
################# Analysis of Exit Trap collections ################
analyseExitTrap = function(firstquarter=15,lastquarter=16){
data=selectData(firstquarter=firstquarter,lastquarter=lastquarter,model='maleSurvival',A0=proportions$A0,P=0.528)
quantiles = parameterEstimates(jagsdata=data$jagsdata,jagsmodel=exittrapF)
# create a vector of emergence rates noting that these are mostly missing
emergence = as.vector(quantiles[quantiles$Var=='emergence','X50.'])
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
plt1 =ggplot(data = toPlot,aes(x = plotdate,group=variable)) + theme_bw() +
theme(text = element_text(size=12)) +
scale_colour_manual(values=get_mp_recomCol(),name='Mosquito category') +
geom_point(aes(y=value,colour=variable),size=2, show.legend=TRUE) +
theme(legend.position = c(.75,.8)) +
theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid')) +
scale_x_continuous(name = 'Day',limits=c(0,80)) +
scale_y_continuous(name = 'Mosquitoes per exit trap')
plt2 =ggplot(data = toPlot,aes(x = plotdate,group=variable)) + theme_bw() +
theme(text = element_text(size=12)) +
scale_colour_manual(values=get_mp_recomCol(),name='Mosquito category') +
geom_point(aes(y=rescaled,colour=variable),size=2, show.legend=FALSE) +
scale_x_continuous(name = 'Day',limits=c(0,80)) +
scale_y_log10(name = 'Relative numbers',limits=c(0.1,10))
library(cowplot)
plt=plot_grid(plt1, plt2, labels = c('A', 'B'), label_size = 12,ncol=2)
result = list(quantiles=quantiles,plt=plt)
return(result)}
# Estimate survival, cycle length and sac rate from exit-trap data by year
analyseExitTrapByYear = function(){
all_quantiles = NULL
fq = seq(1:6)*4 + 4
for(firstquarter in fq){
lastquarter=firstquarter + 3
jagsdata = selectData(firstquarter=firstquarter,lastquarter=lastquarter,
model='maleSurvival',A0=proportions$A0,P=0.528)$jagsdata
quantiles = parameterEstimates(jagsdata=jagsdata,jagsmodel=exittrapF)
quantiles$firstquarter = firstquarter
quantiles$jagsModel='Reference'
all_quantiles = rbind(all_quantiles,quantiles)
jagsdata$P = 0.75
quantiles = parameterEstimates(jagsdata=jagsdata,jagsmodel=exittrapF)
quantiles$firstquarter = firstquarter
quantiles$jagsModel='P = 0.75'
all_quantiles = rbind(all_quantiles,quantiles)
jagsdata$P = NULL
quantiles = parameterEstimates(jagsdata=jagsdata,jagsmodel=exittrapG)
quantiles$firstquarter = firstquarter
quantiles$jagsModel='P = est.'
all_quantiles = rbind(all_quantiles,quantiles)
}
results = all_quantiles
# plot results of fitting
p4a=plotEstimatesByTime(param='resting', paramLabel= 'Resting period (days)')
p4b=plotEstimatesByTime(param='Pm', paramLabel= 'daily survival of males')
p4c=plotEstimatesByTime(param='Teu', paramLabel= 'trapping efficiency (gravids)')
p4d=plotEstimatesByTime(param='Tem', paramLabel= 'trapping efficiency (males)',requirelegend=TRUE)
p4e=plotEstimatesByTime(param='P', paramLabel= 'survival per cycle')
library(cowplot)
plt=plot_grid(p4a, p4b, p4c, p4d, labels = c('A', 'B', 'C', 'D'), label_size = 12,ncol=2)
all_results = list(quantiles=results,plt=plt,p4e=p4e)
return(all_results)}
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
simulated_quantiles = parameterEstimates(jagsdata= simulatedData,jagsmodel=modifiedBirleyModel)
simulated_quantiles$simulation=i
simulated_quantiles$jagsModel='modifiedBirleyModel'
all_simulated_quantiles = rbind(all_simulated_quantiles,simulated_quantiles)
}
simresults = post_process_quantiles(file=all_simulated_quantiles)
toPlot6a = data.frame(simresults$parameters[simresults$parameters$Var=='P',],
input=inputs_to_simulations$P)
plt1=plotSimulations_vs_Inputs(data= toPlot6a, textLabel='survival per cycle', xlim= c(0,1))
toPlot6b = data.frame(simresults$parameters[simresults$parameters$Var=='cycle',],
input=inputs_to_simulations$cycle)
plt2=plotSimulations_vs_Inputs(data= toPlot6b, textLabel='cycle length (days)', xlim= c(1.5,4.5))
library(cowplot)
plt = plot_grid(plt1, plt2, labels = c('A', 'B'), label_size = 12,ncol=2)
results=list(simresults=simresults,inputs_to_simulations=inputs_to_simulations,plt=plt)
return(results)}
SimulationsExitTraps= function(nsimulations){
jagsdata = selectData(firstquarter=15,lastquarter=16,model='maleSurvival',A0=proportions$A0,P=0.528)$jagsdata
# Analyse simulated datasets
all_simulated_quantiles = NULL
inputs_to_simulations = NULL
variable.names = c("tau.b","rg","Teu","resting","pr","Tem","Pm","rm")
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
data0 = simulateExitTrapData(tau.b = tau.b, rg=rg, resting=resting,Teu=Teu,model='maleSurvival',Tem = Tem, Pm = Pm,rm = rm,A0=proportions$A0,P=0.528)
inputs=data.frame(data0$inputs)
# Model F, parous and sac rates as input
simulated_quantiles = parameterEstimates(jagsdata=data0$jagsdata,jagsmodel=exittrapF)
simulated_quantiles$simulation=i
simulated_quantiles$jagsModel='exittrapF'
all_simulated_quantiles = rbind(all_simulated_quantiles,simulated_quantiles)
# P fixed at 0.75
data0$jagsdata$P = 0.75
simulated_quantiles = parameterEstimates(jagsdata=data0$jagsdata,jagsmodel=exittrapF)
simulated_quantiles$simulation=i
simulated_quantiles$jagsModel='P_075'
all_simulated_quantiles = rbind(all_simulated_quantiles,simulated_quantiles)
# A0 fixed at 0.3
data0$jagsdata$P = inputs$P
data0$jagsdata$A0 = 0.3
simulated_quantiles = parameterEstimates(jagsdata=data0$jagsdata,jagsmodel=exittrapF)
simulated_quantiles$simulation=i
simulated_quantiles$jagsModel='A_03'
inputs_to_simulations = rbind(inputs_to_simulations,inputs)
all_simulated_quantiles = rbind(all_simulated_quantiles,simulated_quantiles)
}
simresults = post_process_quantiles(file=all_simulated_quantiles)
simresults$parameters$jagsModel[simresults$parameters$jagsModel=='exittrapF']= 'Reference'
simresults$parameters$jagsModel[simresults$parameters$jagsModel=='P_075']= 'P = 0.75'
simresults$parameters$jagsModel[simresults$parameters$jagsModel=='A_03']= 'A0 = 0.3'
toPlot5a = data.frame(simresults$parameters[simresults$parameters$Var=='resting',],
input=rep(inputs_to_simulations$resting,each=3))
p5a=plotSimulations_vs_Inputs(data= toPlot5a, textLabel='resting period (days)', xlim= c(1,4), requirelegend=TRUE)
toPlot5b = data.frame(simresults$parameters[simresults$parameters$Var=='Teu',],
input=rep(inputs_to_simulations$Teu,each=3))
p5b=plotSimulations_vs_Inputs(data= toPlot5b, textLabel='trapping efficiency (gravids)', xlim= c(1,5))
toPlot5c = data.frame(simresults$parameters[simresults$parameters$Var=='Pm',],
input=rep(inputs_to_simulations$Pm,each=3))
p5c=plotSimulations_vs_Inputs(data= toPlot5c, textLabel='daily survival of males', xlim= c(0,1))
toPlot5d = data.frame(simresults$parameters[simresults$parameters$Var=='Tem',],
input=rep(inputs_to_simulations$Tem,each=3))
p5d=plotSimulations_vs_Inputs(data= toPlot5d, textLabel='trapping efficiency (males)', xlim= c(1,5))
library(cowplot)
plt = plot_grid(p5a, p5b, p5c, p5d, labels = c('A', 'B', 'C', 'D'), label_size = 12,ncol=2)
results=list(simresults=simresults,inputs_to_simulations=inputs_to_simulations,plt=plt)
return(results)
}
################ MAIN PROGRAM STARTS HERE BY READING IN DATA #######################
# if necessary set working directory
setwd("C:/git_repos/charlwood/CBMC")
animais = read.csv(file='Animais.csv')
casa = read.csv(file='Casa.csv')
morte = read.csv(file='Morte.csv')
mosquiteiro = read.csv(file='Mosquiteiro.csv')
mosquito_collection = read.csv2(file='Mosquito_Collection.csv',sep = ",")
pessao = read.csv(file='Pessao.csv')
sporozoites=read.csv(file='Sporozoites.csv')
dissections=read.csv(file='dissections.csv')
################ Data descriptions
description = data_description()
proportions = analysis_of_proportions()
################ Summarise time series #####################
#mosquito_collection$date = as.numeric(as.Date(mosquito_collection$Data, format = "%d/%m/%Y"))-as.numeric(as.Date("26/06/2001", format = "%d/%m/%Y"))
# ensure that there is at least one record for each date
alldates = data.frame(date=as.numeric(as.Date(mosquito_collection$date,'%d/%m/%Y')))-as.numeric(as.Date("26/06/2001", format = "%d/%m/%Y"))
df0 = mosquito_collection
df0$date = alldates$date
CorrelationPlot=plotCorrelations()
################# Analysis of resting collections
resting_analysis = resting_collection_analysis()
requireModifiedBirley
################# Analysis of Light trap data using extension of Birley model
MBData = createModifiedBirleyData()
library(dplyr)
dissDate = as.numeric(as.Date(dissections[,1], format = "%d-%b-%Y"))+716827
df4 = data.frame(cbind(date=dissDate, dissected=dissections$total_dissected,parous = dissections$Sac + dissections$No.sac))
df5 = data.frame(df4 %>%
group_by(date) %>%
dplyr::summarize(parous = sum(parous, na.rm=TRUE),
total = sum(dissected, na.rm=TRUE),
n = n()))
# fields lag1-lag4 were added manually to df5 to give the jags input
jagsDataB = read.csv(file='jags_input_for_Birley_method.csv')
savehistory("C:/git_repos/charlwood/CBMC/history.r")
