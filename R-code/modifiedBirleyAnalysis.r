################# ANALYSIS OF TIME SERIES OF LIGHT TRAP DATA

# The analysis of light trap data requires a valid workfile. Some other information is hardcoded (see comments below)

if( !exists('workfile') )
{
  source('R-code/summarize_data.r')
}

# Creation of  JAGS input dataset for analysis via extended Birley model
modifiedBirleyAnalysis = function(){

  # The creation of the JAGS input dataset requires datset specific processing of dates.
  createModifiedBirleyData = function(){
    dissDate = as.numeric(as.Date(Furvela_data$dissections[,1], format = "%d-%b-%Y"))+716827
    df4 = data.frame(cbind(date=dissDate,
                           dissected=Furvela_data$dissections$total_dissected,
                           parous = Furvela_data$dissections$Sac + Furvela_data$dissections$No.sac))
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

    df8 = data.frame(with(workfile, workfile[Collection=='Light' & date > 2155 & date <= (2159+509),] %>%
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

  MBData = createModifiedBirleyData()
  variable.names = c("tau.t","tau.b","r","cycle","P","p","sigma2","emergence")
  BirleyResults=parameterEstimates(jagsdata=MBData$jagsdata,
                                   jagsmodel=modifiedBirleyModel,
                                   variable.names=variable.names)
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
  BirleyAnalysis = list(BirleyResults=BirleyResults,MBData=MBData)
  return(BirleyAnalysis)
}


# JAGS code
modifiedBirleyModel = "
  model{
    for(t in 1:87){
      parous[t] ~ dnegbin(ppar[t],r)
      T_n[t] <- meanFemales[ptr[t]]
      T_n1[t] <- meanFemales[lag1[t]]
      T_n2[t] <- meanFemales[lag2[t]]
      T_n3[t] <- meanFemales[lag3[t]]
      T_n4[t] <- meanFemales[lag4[t]]

      # expected number parous per trap
      m[t] <- (T_n2[t]*pr[1] + T_n3[t]*pr[2] + T_n4[t]*pr[3])*P

      # expected total parous among those tested
      eParous[t] <- m[t]*dissected[t]/T_n[t]
      ppar[t] <- r/(eParous[t] + r)  # Negative binomial p

      # emergence rate (scaled to light trap collections) n.b. this may need to be constrained to be strictly positive
      emergence[t] <- T_n[t] - m[t]

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
source('R-code/rjags_wrapper.R')

BirleyAnalysis = modifiedBirleyAnalysis()



