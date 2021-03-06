##### TIME SERIES ANALYSIS OF FURVELA ENTOMOLOGY ################

requireTablesPlots = FALSE
requireResimulation = FALSE

############ functions for summarising data
selectData = function(firstquarter=1,lastquarter=31,model='default', A0=0.6037, P= 0.4812){
  threemonths=365.25/4
  mindate = round((firstquarter-1)*threemonths)
  maxdate = round(lastquarter*threemonths)
  date=df0$date
  # for models of temperature dependence link in temperature data
  if(model == 'temperature'){
    library("readxl")
    meteo = read_excel('../some_like_it_hot/Mozambique_temperature_RAW_data.xlsx')
    keeps = c('Year','Week','Mean Temp')
    meteo_week = meteo[grep('Average', meteo$Day), keeps]
    meteo_week = meteo_week[!is.na(meteo_week[['Mean Temp']]) & meteo_week[['Year']]== round(meteo_week[['Year']]),]
    # compute week number for use as lookup
    meteo_week$weekno=meteo_week[['Year']]*52 + meteo_week[['Week']] - 52*2001
    
    lookup_avg_temp = function(x){
      weekno = trunc((176 + x)/7)
      return(meteo_week[['Mean Temp']][which(meteo_week$weekno==weekno)[1]])
    }
    avg_temp = sapply(date,lookup_avg_temp)
  } else { avg_temp = NA }
  
  # For exit trap analysis restrict to period with exit traps
  start = min(date[!is.na(df0$Af.unfed)])
  date_offset = date[!is.na(df0$Af.unfed)] - start + 1
  maxdate = min(maxdate,max(date_offset))
  Af.unfed1 = with(df0,ifelse(Collection=='exit',Af.unfed+Af.part,NA))
  Af.old = with(df0,ifelse(Collection=='exit',Af.fed+Af_gravid,NA))
  Af.gravid = with(df0,ifelse(Collection=='exit',Af_gravid,NA))
  Af.maleEx = with(df0,ifelse(Collection=='exit',Af.male,NA))
  Af.female = Af.unfed1 + Af.old + df0$Af.semi
  Af.femaleLT = with(df0,ifelse(Collection=='Light',Af.female,NA))
  Af.unfed1[date > maxdate | date < mindate] = NA
  Af.maleEx[date > maxdate | date < mindate] = NA
  Af.old[date > maxdate | date < mindate] = NA 
  df = data.frame(cbind(Af.unfed1,Af.maleEx,Af.gravid,Af.old,date,avg_temp))
  df = as.data.frame(sapply(df[,1:6],function(x) as.numeric(as.character(x))))
  df = df[!is.na(df$date) & df$date <= maxdate & df$date >= mindate & !is.na(df$Af.gravid),]
  df1 = data.frame(df %>%
    group_by(date) %>%
    dplyr::summarize(traps = n(),
                      sum.Af.unfed1 = sum(Af.unfed1, na.rm=TRUE),
                      mean.Af.unfed1 = mean(Af.unfed1, na.rm=TRUE),
                      mean.Af.maleEx = mean(Af.maleEx, na.rm=TRUE),
                      sum.Af.maleEx = sum(Af.maleEx, na.rm=TRUE),
                      mean.Af.gravid = mean(Af.gravid, na.rm=TRUE),
                      sum.Af.gravid = sum(Af.gravid, na.rm=TRUE),
                      mean.Af.old = mean(Af.old, na.rm=TRUE),
                      mean_temp = mean(avg_temp)))
  df1$date=df1$date-min(df1$date) + 1
  ndates=max(df1$date)
  
  # create a pointer for finding lagged numbers of unfed 
  df1b= merge(df1,data.frame(date=seq(1:max(df1$date))),by='date',all=TRUE)
  i1=0
  i20=length(which(!is.na(df1$sum.Af.unfed1)))
  i2=i20+5
  df1b$ptr = NA
  for(i in 1:nrow(df1b)){
    if(!is.na(df1b$sum.Af.unfed1[i])) {
      i1=i1+1
      df1b$ptr[i] = i1
    } else {
      i2=i2+1
      df1b$ptr[i] = i2
    }
  }
  df1b$lag1=c(i20+5,df1b$ptr[1:(ndates-1)])
  df1b$lag2=c(i20+4,df1b$lag1[1:(ndates-1)])
  df1b$lag3=c(i20+3,df1b$lag2[1:(ndates-1)])
  df1b$lag4=c(i20+2,df1b$lag3[1:(ndates-1)])
  df1b$lag5=c(i20+1,df1b$lag4[1:(ndates-1)])
  
  # df1b - all data
  # df1c - all sampled dates
  # df1d - all sampled dates with complete data for lags
  df1c = with(df1b, df1b[!is.na(sum.Af.unfed1),]) 
  df1d = with(df1c, df1c[lag1 <= i20  & lag2 <= i20  & lag3 <= i20  & lag4 <= i20,])
  
  # lag5 data are required for the model of male survival
  df1d = df1d[df1d$lag5 <= i20,]
  jagsdata = list(sampleddates = i20,
                  unfed1=df1c$sum.Af.unfed1[1:i20],
                  traps = df1c$traps[1:i20],
                  days_with_complete_data = length(df1d$ptr),
                  gravid=df1c$sum.Af.gravid[1:i20],
                  ptr = df1d$ptr,
                  lag1 = df1d$lag1,
                  lag2 = df1d$lag2,
                  lag3 = df1d$lag3,
                  lag4 = df1d$lag4)
  
      jagsdata$lag5 = df1d$lag5
      jagsdata$P = P
      jagsdata$A0 = A0
      jagsdata$males = df1c$sum.Af.maleEx
  if(model == 'temperature'){ jagsdata$mean_temp = df1c$mean_temp}    
return(list(df1=df1,jagsdata=jagsdata,plotdata=df1b))}

####### Summarise by trapping method for tabulation 

summarise_by_method = function(){
  # houses observed
  df3=data.frame(table(df0$Casa,df0$Collection))
  df3$Freq[df3$Freq > 0] = 1
  df4 = as.data.frame(table(df3$Var2,df3$Freq))
  # nights observed
  df3=data.frame(table(df0$Data,df0$Collection))
  df3$Freq[df3$Freq > 0] = 1
  df5 = data.frame(table(df3$Var2,df3$Freq))
  df0$l_Af.unfed = log(df0$Af.unfed+1)
  df0$l_Af.part = log(df0$Af.part+1)
  df0$l_Af.fed = log(df0$Af.fed+1)
  df0$l_Af.semi = log(df0$Af.semi+1)
  df0$l_Af.gravid= log(df0$Af_gravid+1)
  df0$l_Af.male = log(df0$Af.male+1)
  df2 = data.frame(df0 %>%
    group_by(Collection) %>%
    dplyr::summarize(mu.Af.unfed = mean(l_Af.unfed, na.rm=TRUE),
                     mu.Af.part = mean(l_Af.part, na.rm=TRUE),
                     mu.Af.fed = mean(l_Af.fed, na.rm=TRUE),
                     mu.Af.semi = mean(l_Af.semi, na.rm=TRUE),
                     mu.Af.gravid = mean(l_Af.gravid, na.rm=TRUE),
                     mu.Af.male = mean(l_Af.male, na.rm=TRUE),
                     Af.unfed = sum(Af.unfed, na.rm=TRUE),
                     Af.part = sum(Af.part, na.rm=TRUE),
                     Af.fed = sum(Af.fed, na.rm=TRUE),
                     Af.semi = sum(Af.semi, na.rm=TRUE),
                     Af.gravid = sum(Af_gravid, na.rm=TRUE),
                     Af.male = sum(Af.male, na.rm=TRUE)))
 df2$w_Af.unfed=exp(as.numeric(df2$mu.Af.unfed))-1
 df2$w_Af.part=exp(df2$mu.Af.part)-1
 df2$w_Af.fed=exp(df2$mu.Af.fed)-1
 df2$w_Af.semi=exp(df2$mu.Af.semi)-1
 df2$w_Af.gravid=exp(df2$mu.Af.gravid)-1
 df2$w_Af.male=exp(df2$mu.Af.male)-1
 df2=df2[!is.na(df2$Collection),c(1,(seq(1:12)+7))]
 df2$houses = df4$Freq[df4$Var2 == 1]
 df2$dates = df5$Freq[df5$Var2 == 1]
 df2$trapnights = data.frame(table(df0$Collection))$Freq
return(df2)
}

####### Summarise time series by quarter for plotting
summarise_by_quarter = function(){
  threemonths=365.25/4
  quarter = ceiling(df0$date/threemonths)
  Af.unfed1 = with(df0,ifelse(Collection=='exit',log(Af.unfed+Af.part+1),NA))
  Af.gravid = with(df0,ifelse(Collection=='exit',log(Af_gravid+1),NA))
  Af.maleEx = with(df0,ifelse(Collection=='Light',log(Af.male+1),NA))
  Af.femaleLT = with(df0,ifelse(Collection=='Light',log(Af.unfed+Af.part+Af.fed+Af.semi+Af_gravid+1),NA))
  df = data.frame(cbind(Af.unfed1,Af.maleEx,Af.femaleLT,Af.gravid,quarter))
  df2 = df %>%
    group_by(quarter) %>%
    dplyr::summarize(mu.Af.unfed1 = mean(Af.unfed1, na.rm=TRUE),
                     mu.Af.maleEx = mean(Af.maleEx, na.rm=TRUE),
                     mu.Af.femaleLT = mean(Af.femaleLT, na.rm=TRUE),
                     mu.Af.gravid = mean(Af.gravid, na.rm=TRUE))
  df2 = df2[!is.na(df2$quarter),]
  library(data.table)
  df3 <- melt(df2, id.vars = "quarter")
  df3$value = exp(df3$value)-1
  ticklabels = c(rep('',3),'2002',rep('',3),'2003',rep('',3),'2004',rep('',3),'2005',rep('',3),'2006',rep('',3),'2007',rep('',3),'2008',rep('',3))
  ggplot(data = df3,aes(x = as.factor(quarter),group=variable)) + theme_bw() +
    theme(text = element_text(size=12)) +
    scale_colour_manual(values=get_mp_recomCol(),name='Mosquito category',
                        labels=c('An. funestus unfed or part fed',
                                 'An. funestus males (light trap)',
                                 'An. funestus females (light trap)',
                                 'An. funestus fed or gravid')) +
    geom_line(aes(y=value,colour=variable),size=2) +
    geom_vline(xintercept=22.5) +
    scale_y_log10(name  = 'Williams mean (per trap night)')+
    scale_x_discrete(name = 'Year',labels=ticklabels)
return(df2)}

################# ANALYSIS OF CORRELATIONS
# Function to Calculate Spearman cross-correlations for different lags
calculate_correlations <- function(df1=df1) {
  corr = spearman_CI(x=df1$mean.Af.unfed1,y=df1$mean.Af.maleEx,yname="male")
  corr1 = spearman_CI(x=df1$mean.Af.unfed1,y=df1$mean.Af.gravid,yname="gravid")
  corr$lag=-0.1
  corr1$lag=0.1
  corr = rbind(corr,corr1)
  for (lag in 1:10){
    na.vector= rep(NA,lag)
    mean.Af.unfed1_lag = c(na.vector,df1$mean.Af.unfed1)
    corr_lag = spearman_CI(x=mean.Af.unfed1_lag,y=df1$mean.Af.maleEx,yname="male")
    corr1_lag = spearman_CI(x=mean.Af.unfed1_lag,y=df1$mean.Af.gravid,yname="gravid")
    corr_lag$lag=lag-0.1
    corr1_lag$lag=lag+0.1
    corr = rbind(corr,corr_lag,corr1_lag)
  }
return(corr)}

# Function to calculate Spearman autocorrelations for different lags
calculate_autocorrelations <- function(df1=df1) {
  corr = spearman_CI(x=df1$mean.Af.maleEx,y=df1$mean.Af.maleEx,yname="male")
  corr1 = spearman_CI(x=df1$mean.Af.unfed1,y=df1$mean.Af.unfed1,yname="unfed1")
  corr2 = spearman_CI(x=df1$mean.Af.gravid,y=df1$mean.Af.gravid,yname="gravid")
  corr$lag=-0.2
  corr1$lag=0.0
  corr2$lag=0.2
  corr = rbind(corr,corr1,corr2)
  for (lag in 1:10){
    na.vector= rep(NA,lag)
    mean.Af.maleEx_lag = c(na.vector,df1$mean.Af.maleEx)
    mean.Af.unfed1_lag = c(na.vector,df1$mean.Af.unfed1)
    mean.Af.gravid_lag = c(na.vector,df1$mean.Af.gravid)
    corr_lag = spearman_CI(x=mean.Af.maleEx_lag,y=df1$mean.Af.maleEx,yname="male")
    corr1_lag = spearman_CI(x=mean.Af.unfed1_lag,y=df1$mean.Af.unfed1,yname="unfed1")
    corr2_lag = spearman_CI(x=mean.Af.gravid_lag,y=df1$mean.Af.gravid,yname="gravid")
    corr_lag$lag=lag-0.2
    corr1_lag$lag = lag
    corr2_lag$lag=lag+0.2
    corr = rbind(corr,corr_lag,corr1_lag,corr2_lag)
  }
  return(corr)}

##### Calculate correlations for plotting ##############
spearman_CI <- function(x, y, yname, alpha = 0.05){
  # approximate 95% CI following Bonett and Wright (2000) https://doi.org/10.1007/BF02294183 
  y = y[1:length(x)]
  rs = cor(x, y, method = "spearman", use = "complete.obs")
  n = sum(complete.cases(x, y))
  CI = sort(tanh(atanh(rs) + c(-1,1)*sqrt((1+rs^2/2)/(n-3))*qnorm(p = alpha/2)))
  result = data.frame(yname=yname,nobs=n, rs=rs, lower=CI[1], upper=CI[2])
  return(result)}

plotCorrelations = function(){
  data= selectData(firstquarter=1,lastquarter=32)
  df1=data$df1
  corr= calculate_correlations(df1=df1)
  levels(corr$yname) = list('Males'='male' , 'Gravid' = 'gravid')
  plt1= with(corr,errbarPlot(data=corr,low=lower,med=rs,high=upper,xname='lag (days)',yname=yname,ylimits=c(-0.3,0.7),logscale=FALSE))
  autocorr= calculate_autocorrelations(df1=df1)
  levels(autocorr$yname) = list('Male'='male' , 'Gravid' = 'gravid', 'Unfed or part fed' = 'unfed1' )
  plt2= with(autocorr,errbarPlot(data=autocorr,low=lower,med=rs,high=upper,xname='lag (days)',yname=yname,ylimits=c(0.2,0.7),logscale=FALSE,requireLegend=TRUE))
  library(cowplot)
  plt=plot_grid(plt1, plt2, labels = c('A', 'B'), label_size = 12,ncol=2)
return(plt)  }

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

# Plot Simulated Exit Trap Data
plotSimulatedData = function() {
  toPlot = data.frame(value=c(emergence,em,eg))
  toPlot$var = as.factor(rep(c('Emergence','Males','Gravid'),each=length(em)))
  toPlot$time = rep(seq(1:length(T0x)),3)
  ggplot(data=toPlot,aes(x=time,y=value, colour=var)) +
  geom_point() + 
  theme_bw() +
  theme(text = element_text(size=12)) 
return()
}

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

#################### PLOTTING FUNCTIONS ######################
library(ggplot2)
get_mp_recomCol = function() {return(c("#D55E00", "#009E73", "#0072A7","#C879C8"))}
errbarPlot = function(data,low,med,high,xname,yname,ylimits,logscale=TRUE,requireLegend=FALSE){
  data$low=low
  data$med=med
  data$high=high
  plt = ggplot(data = data,aes(x = lag, group = yname)) + theme_bw() +
    theme(text = element_text(size=12)) +
    
    scale_colour_manual(values=get_mp_recomCol(),name='Mosquito category') +
    geom_errorbar(aes(ymin=low,
                      ymax=high, width = 0.5, colour = yname), show.legend=requireLegend) +
    theme(legend.position = c(.75,.82))+ 
    scale_fill_manual(values=get_mp_recomCol(), guide = 'none') +
    geom_point(aes(y=rs, colour = yname), size=5,show.legend=FALSE) +
    scale_x_continuous(name = xname,limits=c(-1,11), breaks=c(0,2,4,6,8,10))
  if(logscale) plt=plt + scale_y_log10(name = yname, limits=ylimits)
  if(!logscale) plt=plt + scale_y_continuous(name = 'Spearman r', limits=ylimits)
  return(plt)
}

plotEstimatesByTime = function(param,paramLabel,results=results,requirelegend=FALSE){
  toPlot = results[results$Var==param,]
  toPlot$jagsModel[toPlot$jagsModel=='Reference']='a'
  toPlot$jagsModel[toPlot$jagsModel=='P = est.']='b'
  toPlot$jagsModel[toPlot$jagsModel=='Teu = 1']='c'
  toPlot$jagsModel[toPlot$jagsModel=='P = 0.75']='d'
  toPlot$jagsModel[toPlot$jagsModel=='Tem = 1']='b' # appears only in plots for males
  toPlot$firstquarter[toPlot$jagsModel=='b'] = toPlot$firstquarter[toPlot$jagsModel=='b'] - 0.3
  toPlot$firstquarter[toPlot$jagsModel=='c'] = toPlot$firstquarter[toPlot$jagsModel=='c'] + 0.3
  ticklabels = c('2003','2004','2005','2006','2007','2008')
  plt = ggplot(data = toPlot,aes(x = (firstquarter+1), group = jagsModel)) + theme_bw() +
    theme(text = element_text(size=10)) +
    scale_colour_manual(values=get_mp_recomCol(),name='Model') +
    geom_errorbar(aes(ymin=X2.5.,
                      ymax=X97.5., width = 0.5, colour = jagsModel), show.legend=FALSE) +
    scale_fill_manual(values=get_mp_recomCol(), guide = 'none') +
    geom_point(aes(y=X50., colour = jagsModel), shape=21, size=3, stroke = 1.5, show.legend=requirelegend) +
    theme(legend.position = c(.6,.75))+  
    theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid')) +
    scale_y_continuous(name = paramLabel) +
    scale_x_continuous(name = 'Year',limits=c(8,32), breaks=c(8.5,12.5,16.5,20.5,24.5,28.5),labels=ticklabels) +
    geom_vline(xintercept=22.5) 
  return(plt)}

plotEstimatesByTemp = function(param,paramLabel,results=results,requirelegend=FALSE){
  toPlot = results[results$Var==param,]
  toPlot$jagsModel[toPlot$jagsModel=='Reference']='a'
  toPlot$jagsModel[toPlot$jagsModel=='P = est.']='b'
  toPlot$jagsModel[toPlot$jagsModel=='Teu = 1']='c'
  toPlot$jagsModel[toPlot$jagsModel=='P = 0.75']='x'
  toPlot$jagsModel[toPlot$jagsModel=='Tem = 1']='d' # appears only in plots for males
  ticklabels = c('<23?','23-25?','25-27?','27-29?','>29?')
  toPlot$Index = ifelse(toPlot$jagsModel =='a',toPlot$Index - 0.2, toPlot$Index)
  toPlot$Index = ifelse(toPlot$jagsModel =='c' | toPlot$jagsModel =='d',toPlot$Index + 0.2, toPlot$Index)
  plt = ggplot(data = toPlot,aes(x = Index, group = jagsModel)) + theme_bw() +
    theme(text = element_text(size=10)) +
    scale_colour_manual(values=get_mp_recomCol(),name=NULL) +
    geom_errorbar(aes(ymin=X2.5.,
                      ymax=X97.5., width = 0.2, colour = jagsModel), show.legend=FALSE) +
    scale_fill_manual(values=get_mp_recomCol()) +
    geom_point(aes(y=X50., colour = jagsModel), shape=21, size=3, stroke = 1.5, show.legend=requirelegend) +
    theme(legend.position = c(.35,.09))+ 
    guides(colour = guide_legend(nrow = 1,label.position = "left"))+
    theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'),
          legend.title=element_blank(),
          legend.key.height = unit(0.1, 'cm')) +
    scale_y_continuous(name = paramLabel) +
    scale_x_continuous(name = 'Weekly Mean Temperature',limits=c(0.5,5.5), breaks=c(1,2,3,4,5),labels=ticklabels)
  return(plt)}



plotSimulations_vs_Inputs = function(data, textLabel,xlim, nmodels=5, requirelegend= FALSE){
  plt= ggplot(data = data,aes(x = input, y=X50., group = jagsModel)) + theme_bw(base_size = 9) +
  theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid')) +
            theme(legend.position = c(.8,.5)) +    
  geom_point(aes(y=X50.,colour = jagsModel, shape = jagsModel), size=2, stroke = 1.5, show.legend=requirelegend) +
  scale_colour_discrete(name='Constraints') +
  scale_shape_manual(name='Constraints',values=(seq(1:nmodels)-1)) +
  geom_line(aes(x=input,y=input)) +
  scale_x_continuous(name = paste('Input',textLabel), limits=xlim) +
  scale_y_continuous(name = paste('Modelled',textLabel), limits = xlim)
return(plt)}

plotConsistency = function(data, textLabel,nmodels=5, requirelegend= FALSE){
  data$ratio = data$X50./data$input
  data$sqe = (data$X50.- data$input)^2
  plt= ggplot(data = data,aes(x = days_with_complete_data, y=ratio, group = jagsModel)) + theme_bw(base_size = 9) +
    theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid')) +
    theme(legend.position = c(.1,.32)) +    
    #geom_point(aes(y=sqe,colour = jagsModel, shape = jagsModel), size=2, stroke = 1.5) +
    geom_smooth(aes(y=sqe,colour = jagsModel), size=1, show.legend=requirelegend) +
    scale_colour_discrete(name=element_blank()) +
    scale_shape_manual(name=element_blank(),values=(seq(1:nmodels)-1)) +
    scale_x_continuous(name = 'Days with complete data') +
    scale_y_log10(name = paste('SqErr:' ,textLabel))
  return(plt)}

source('JAGS_code.R')

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
    quantiles = parameterEstimates(jagsdata=jagsdata,jagsmodel=exittrap1)
    quantiles$firstquarter = firstquarter
    quantiles$jagsModel='Reference'
    all_quantiles = rbind(all_quantiles,quantiles)
    
    jagsdata$P = rep(0.75,repeats)
    quantiles = parameterEstimates(jagsdata=jagsdata,jagsmodel=exittrap1)
    quantiles$firstquarter = firstquarter
    quantiles$jagsModel='P = 0.75'
    all_quantiles = rbind(all_quantiles,quantiles)
    
    jagsdata$P = NULL
    quantiles = parameterEstimates(jagsdata=jagsdata,jagsmodel=exittrap1)
    quantiles$firstquarter = firstquarter
    quantiles$jagsModel='P = est.'
    all_quantiles = rbind(all_quantiles,quantiles)
    
    jagsdata$Teu = 1
    quantiles = parameterEstimates(jagsdata=jagsdata,jagsmodel=exittrap1)
    quantiles$firstquarter = firstquarter
    quantiles$jagsModel='Teu = 1'
    all_quantiles = rbind(all_quantiles,quantiles)
    
    jagsdata$Tem = 1
    quantiles = parameterEstimates(jagsdata=jagsdata,jagsmodel=exittrap1)
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
    simulated_quantiles = parameterEstimates(jagsdata=data0$jagsdata,jagsmodel=exittrap1)
    simulated_quantiles$simulation=i
    simulated_quantiles$jagsModel='exittrap1'
    all_simulated_quantiles = rbind(all_simulated_quantiles,simulated_quantiles)
    # P fixed at 0.75
    data0$jagsdata$P = 0.75
    simulated_quantiles = parameterEstimates(jagsdata=data0$jagsdata,jagsmodel=exittrap1)
    simulated_quantiles$simulation=i
    simulated_quantiles$jagsModel='P_075'
    all_simulated_quantiles = rbind(all_simulated_quantiles,simulated_quantiles)
    # P known, A0 fixed at 0.3
    data0$jagsdata$P = inputs$P
    data0$jagsdata$A0 = 0.3 
    simulated_quantiles = parameterEstimates(jagsdata=data0$jagsdata,jagsmodel=exittrap1)
    simulated_quantiles$simulation=i
    simulated_quantiles$jagsModel='A_03'
    all_simulated_quantiles = rbind(all_simulated_quantiles,simulated_quantiles)
    # P estimated
    data0$jagsdata$P = NULL
    data0$jagsdata$A0 = inputs$A0 
    simulated_quantiles = parameterEstimates(jagsdata=data0$jagsdata,jagsmodel=exittrap1)
    simulated_quantiles$simulation=i
    simulated_quantiles$jagsModel='P_est'
    all_simulated_quantiles = rbind(all_simulated_quantiles,simulated_quantiles)
    # Teu known, P Estimated
    data0$jagsdata$Teu = inputs$Teu
    simulated_quantiles = parameterEstimates(jagsdata=data0$jagsdata,jagsmodel=exittrap1)
    simulated_quantiles$simulation=i
    simulated_quantiles$jagsModel='Teu_known'
    all_simulated_quantiles = rbind(all_simulated_quantiles,simulated_quantiles)
    # Teu, Tem known, P Estimated
    data0$jagsdata$Tem = inputs$Tem
    simulated_quantiles = parameterEstimates(jagsdata=data0$jagsdata,jagsmodel=exittrap1)
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

################ MAIN PROGRAM STARTS HERE BY READING IN DATA #######################
# if necessary set working directory
setwd("C:/git_repos/charlwood/CBMC")
library(dplyr)

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

################# Analysis of Light trap data using extension of Birley model

if (requireResimulation){
  MBData = createModifiedBirleyData()
  variable.names = c("tau.t","tau.b","r","cycle","P","p","sigma2","emergence")
  BirleyResults=parameterEstimates(jagsdata=MBData$jagsdata,jagsmodel=modifiedBirleyModel)
  cycle_sample = rnorm(1000, mean=BirleyResults$mean[BirleyResults$Var == 'cycle'], 
                               sd=BirleyResults$sd[BirleyResults$Var == 'cycle'])
  theta_prime_parous = cycle_sample - (1-proportions$A0_sample)/proportions$A0_sample
  theta_prime_quantiles =quantile(x = theta_prime_parous, probs = c(0.025,0.5,0.975))
  theta_prime_record=data.frame(Var='theta_prime',Index=1,
                                X50.=theta_prime_quantiles[2],
                                X2.5.=theta_prime_quantiles[1],
                                X97.5.=theta_prime_quantiles[3],
                                mean=mean(theta_prime_parous),
                                sd=sd(theta_prime_parous))   
  BirleyResults = rbind(BirleyResults,theta_prime_record)
  
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


