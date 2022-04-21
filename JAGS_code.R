################ JAGS CODE STORED AS STRINGS ############################

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

createListOfModels = function(){
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
  # Unfed mosquitoes are assumed to include survivors from previous cycles (contrary to Derek's understanding)
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
  # Unfed mosquitoes are assumed to include survivors from previous cycles (contrary to Derek's understanding)
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
  
  # Estimation of emergence rates and male survival from exit trap data 
  # Unfed mosquitoes and emergent mosquitoes are equivalent
  exittrapH ="
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
  P ~ dunif(0,1)          # overall survival through oviposition cycle
  A0 ~ dunif(0,1)         # proportion returning to feed within same night
  p <- pow(P,1/cycle)   # daily survival of females
}  
"
  exittrapI ="
model{
  for(d in 1:sampleddates){
    # numbers of unfed (+ part fed) per trap (equivalent to emergence rate)
    em[d] ~  dgamma(re, lambda_e)
    T0[d] <- em[d]*traps[d]   
    unfed1[d] ~ dpois(T0[d])
    # gravid per trap
    eg[d] <- gravid[d]/traps[d]
  }
  re ~ dunif(0, 3)
  lambda_e ~dunif(0, 3)
  
  for(t in 1:days_with_complete_data){
    # expected number gravid per trap (Teu is relative trapping efficiency of gravid vs unfed)

    # expected gravids in first cycle
    eg0[t] <- (em[lag1[t]]*pr[1] + em[lag2[t]]*pr[2] + em[lag3[t]]*pr[3] + em[lag4[t]]*pr[4])* Teu 
    # expected gravids in subsequent cycles
    eg1[t] <- eg[lag1[t]]*pr[1] + eg[lag2[t]]*pr[2] + eg[lag3[t]]*pr[3] + eg[lag4[t]]*pr[4]
    eg2[t] <- eg[lag2[t]]*pr[1] + eg[lag3[t]]*pr[2] + eg[lag4[t]]*pr[3] + eg[lag5[t]]*pr[4]
    G[t] <- (eg0[t] + A0*eg1[t] + (1-A0)*eg2[t])*P
    pg[t] <- rg/(G[t]*traps[ptr[t]] + rg)  # Negative binomial p 
    gravid[ptr[t]] ~ dnegbin(pg[t],rg)
    
    
    #emergence is equivalent to unfed
    emergence[t] <- em[ptr[t]]
    
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
  P ~ dunif(0,1)          # overall survival through oviposition cycle
  A0 ~ dunif(0,1)         # proportion returning to feed within same night
  p <- pow(P,1/cycle)   # daily survival of females
}  
"
  
  exittrapJ ="
model{
  for(d in 1:sampleddates){
    # numbers of unfed (+ part fed) per trap (equivalent to emergence rate)
    em[d] ~  dgamma(r_e, lambda_e)
    T0[d] <- em[d]*traps[d]   
    unfed1[d] ~ dpois(T0[d])
    # gravid per trap
    eg[d] ~ dgamma(r_g, lambda_g)
    gm[d] <- eg[d]*traps[d]
    gravid1[d] ~ dpois(gm[d])
    
        # gravid per trap
    eg[d] <- gravid[d]/traps[d]
  }
  r_e ~ dgamma(1, 1)
  lambda_e ~ dgamma(1, 1)
  r_g ~ dgamma(1, 1)
  lambda_g ~ dgamma(1, 1)

  for(t in 1:days_with_complete_data){
    # expected number gravid per trap (Teu is relative trapping efficiency of gravid vs unfed)

    # expected gravids in first cycle
    eg0[t] <- (em[lag1[t]]*pr[1] + em[lag2[t]]*pr[2] + em[lag3[t]]*pr[3] + em[lag4[t]]*pr[4])* Teu 
    # expected gravids in subsequent cycles
    eg1[t] <- eg[lag1[t]]*pr[1] + eg[lag2[t]]*pr[2] + eg[lag3[t]]*pr[3] + eg[lag4[t]]*pr[4]
    eg2[t] <- eg[lag2[t]]*pr[1] + eg[lag3[t]]*pr[2] + eg[lag4[t]]*pr[3] + eg[lag5[t]]*pr[4]
    G[t] <- (eg0[t] + A0*eg1[t] + (1-A0)*eg2[t])*P
    pg[t] <- rg/(G[t]*traps[ptr[t]] + rg)  # Negative binomial p 
    gravid[ptr[t]] ~ dnegbin(pg[t],rg)
    
    
    #emergence is equivalent to unfed
    emergence[t] <- em[ptr[t]]
    
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
  P ~ dunif(0,1)          # overall survival through oviposition cycle
  A0 ~ dunif(0,1)         # proportion returning to feed within same night
  p <- pow(P,1/cycle)   # daily survival of females
}  
"
  exittrapK ="
model{

  # distribution of emergence rates
  # expected value should be mean(unfed1[])
  
  mean.u <- mean(unfed1)
  mean.y <- log(mean.u) - var.y/2
  
  for(d in 1:sampleddates){
    # numbers of unfed (+ part fed) per trap (equivalent to emergence rate)
    log.em[d] ~ dnorm(mean.y,tau.y)
    # em[d] <- exp(log.em[d])
    log.unfed[d] ~ dnorm(log.em[d],tau.u)
    em[d] <- exp(log.unfed[d])
    # gravid per trap
    eg[d] <- gravid[d]/traps[d]
  }

  tau.y ~ dgamma(1, 1)
  var.y <- 1/tau.y
  tau.u ~ dgamma(1, 1)


  for(t in 1:days_with_complete_data){
    # expected number gravid per trap (Teu is relative trapping efficiency of gravid vs unfed)

    # expected gravids in first cycle
    eg0[t] <- (em[lag1[t]]*pr[1] + em[lag2[t]]*pr[2] + em[lag3[t]]*pr[3] + em[lag4[t]]*pr[4])* Teu 
    # expected gravids in subsequent cycles
    eg1[t] <- eg[lag1[t]]*pr[1] + eg[lag2[t]]*pr[2] + eg[lag3[t]]*pr[3] + eg[lag4[t]]*pr[4]
    eg2[t] <- eg[lag2[t]]*pr[1] + eg[lag3[t]]*pr[2] + eg[lag4[t]]*pr[3] + eg[lag5[t]]*pr[4]
    G[t] <- (eg0[t] + A0*eg1[t] + (1-A0)*eg2[t])*P
    pg[t] <- rg/(G[t]*traps[ptr[t]] + rg)  # Negative binomial p 
    gravid[ptr[t]] ~ dnegbin(pg[t],rg)
    
    
    #emergence is equivalent to unfed
    emergence[t] <- em[ptr[t]]
    
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
  P ~ dunif(0,1)          # overall survival through oviposition cycle
  A0 ~ dunif(0,1)         # proportion returning to feed within same night
  p <- pow(P,1/cycle)   # daily survival of females
}  
"
  exittrapL ="
model{
  for(d in 1:sampleddates){
    # numbers of unfed (+ part fed) per trap (equivalent to emergence rate)
    em[d] ~  dgamma(r_e, lambda_e)
    T0[d] <- em[d]*traps[d]   
    unfed1[d] ~ dpois(T0[d])
    # gravid per trap
    eg[d] <- gravid[d]/traps[d]
  }
  r_e ~ dgamma(1, 1)
  lambda_e ~ dgamma(1, 1)

  for(t in 1:days_with_complete_data){
    # expected number gravid per trap (Teu is relative trapping efficiency of gravid vs unfed)

    # expected gravids in first cycle
    eg0[t] <- (em[lag1[t]]*pr[1] + em[lag2[t]]*pr[2] + em[lag3[t]]*pr[3] + em[lag4[t]]*pr[4])* Teu 
    # expected gravids in subsequent cycles
    eg1[t] <- eg[lag1[t]]*pr[1] + eg[lag2[t]]*pr[2] + eg[lag3[t]]*pr[3] + eg[lag4[t]]*pr[4]
    eg2[t] <- eg[lag2[t]]*pr[1] + eg[lag3[t]]*pr[2] + eg[lag4[t]]*pr[3] + eg[lag5[t]]*pr[4]
    G[t] <- (eg0[t] + A0*eg1[t] + (1-A0)*eg2[t])*P
    pg[t] <- rg/(G[t]*traps[ptr[t]] + rg)  # Negative binomial p 
    gravid[ptr[t]] ~ dnegbin(pg[t],rg)
    
    
    #emergence is equivalent to unfed
    emergence[t] <- em[ptr[t]]
    
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
  P ~ dunif(0,1)          # overall survival through oviposition cycle
  A0 ~ dunif(0,1)         # proportion returning to feed within same night
  p <- pow(P,1/cycle)   # daily survival of females
}  
"
  
  # Binomial fits to proportions
  exittrapM ="
model{
  for(d in 1:sampleddates){
    # numbers of unfed (+ part fed) per trap (equivalent to emergence rate)
    em[d] ~  dgamma(r_e, lambda_e)
    T0[d] <- em[d]*traps[d]   
    unfed1[d] ~ dpois(T0[d])
    # gravid per trap
    eg[d] <- gravid[d]/traps[d]
    # denominators for binomial models
    tg[d] <- unfed1[d] + gravid[d]
    tm[d] <- unfed1[d] + males[d]
  }
  r_e ~ dgamma(1, 1)
  lambda_e ~ dgamma(1, 1)

  for(t in 1:days_with_complete_data){
  
    #emergence is equivalent to unfed per trap
    emergence[t] <- em[ptr[t]]
    
    # expected number gravid per trap (Teu is relative trapping efficiency of gravid vs unfed)
    # expected gravids in first cycle (with 100% survival)
    eg0[t] <- (em[lag1[t]]*pr[1] + em[lag2[t]]*pr[2] + em[lag3[t]]*pr[3] + em[lag4[t]]*pr[4])* Teu 
    # expected gravids in subsequent cycles (with 100% survival)
    eg1[t] <- eg[lag1[t]]*pr[1] + eg[lag2[t]]*pr[2] + eg[lag3[t]]*pr[3] + eg[lag4[t]]*pr[4]
    eg2[t] <- eg[lag2[t]]*pr[1] + eg[lag3[t]]*pr[2] + eg[lag4[t]]*pr[3] + eg[lag5[t]]*pr[4]

    # overall expected gravids per trap
    G[t] <- (eg0[t] + A0*eg1[t] + (1-A0)*eg2[t])*P
    
    pg[t] <- G[t]/(emergence[t] + G[t])
    gravid1[ptr[t]] ~ dbin(pg[t],tg[ptr[t]])
    
    
    # expected number of males:
    #     Tem is relative trapping efficiency of males vs unfed assuming 50:50 sex ratio on emergence
    #     Pm is daily survival of males
    
    # Expected males per trap 
    Emales[t] <- (Tem*emergence[t] + Pm * males[lag1[t]]/traps[lag1[t]])
    
    pm[t] <- Emales[t]/(emergence[t] + Emales[t])
    males1[ptr[t]] ~ dbin(pm[t],tm[ptr[t]])
    
  }

  
  ################# PRIORS ########################

  # discretised normal kernel for the cycle duration
 
  pr[1] <- pnorm(1.5,resting,tau.b)
  pr[2] <- pnorm(2.5,resting,tau.b) - pr[1]
  pr[3] <- pnorm(3.5,resting,tau.b) - pr[1] - pr[2]
  pr[4] <- 1 - pr[1] - pr[2] - pr[3]
  
  tau.b ~ dgamma(1,1)    # precision of distribution of cycle duration 

  resting ~ dunif(0.5,4.5) # duration of resting period
  cycle <- resting + (1 - A0)
  Teu ~ dgamma(1,1)       # relative trapping efficiency of gravids relative to unfeds 
  Tem ~ dgamma(1,1)       # relative trapping efficiency of males vs unfed assuming 50:50 sex ratio on emergence
  Pm ~ dunif(0,1)         # daily survival of males
  P ~ dunif(0,1)          # overall survival through oviposition cycle
  A0 ~ dunif(0,1)         # proportion returning to feed within same night
  p <- pow(P,1/cycle)   # daily survival of females
}  
"
  
  exittrapN ="
model{
  for(d in 1:sampleddates){
    # numbers of unfed (+ part fed) per trap (equivalent to emergence rate)
    em[d] ~  dgamma(re, lambda_e)
    T0[d] <- em[d]*traps[d]   
    unfed1[d] ~ dpois(T0[d])
    # gravid per trap
    eg[d] <- gravid[d]/traps[d]
  }
  

    # gravid per trap
    eg[d] <- gravid[d]/traps[d]
    # denominators for binomial models
    tg[d] <- unfed1[d] + gravid[d]
    tm[d] <- unfed1[d] + males[d]
  }
  tau.y ~ dgamma(1, 1)
  var.y <- 1/tau.y
  tau.u ~ dgamma(1, 1)


  for(t in 1:days_with_complete_data){
  
    #emergence is equivalent to unfed per trap
    emergence[t] <- em[ptr[t]]
    
    # expected number gravid per trap (Teu is relative trapping efficiency of gravid vs unfed)
    # expected gravids in first cycle (with 100% survival)
    eg0[t] <- (em[lag1[t]]*pr[1] + em[lag2[t]]*pr[2] + em[lag3[t]]*pr[3] + em[lag4[t]]*pr[4])* Teu 
    # expected gravids in subsequent cycles (with 100% survival)
    eg1[t] <- eg[lag1[t]]*pr[1] + eg[lag2[t]]*pr[2] + eg[lag3[t]]*pr[3] + eg[lag4[t]]*pr[4]
    eg2[t] <- eg[lag2[t]]*pr[1] + eg[lag3[t]]*pr[2] + eg[lag4[t]]*pr[3] + eg[lag5[t]]*pr[4]

    # overall expected gravids per trap
    G[t] <- (eg0[t] + A0*eg1[t] + (1-A0)*eg2[t])*P
    
    pg[t] <- G[t]/(emergence[t] + G[t])
    gravid1[ptr[t]] ~ dbin(pg[t],tg[ptr[t]])
    
    
    # expected number of males:
    #     Tem is relative trapping efficiency of males vs unfed assuming 50:50 sex ratio on emergence
    #     Pm is daily survival of males
    
    # Expected males per trap 
    Emales[t] <- (Tem*emergence[t] + Pm * males[lag1[t]]/traps[lag1[t]])
    
    pm[t] <- Emales[t]/(emergence[t] + Emales[t])
    males1[ptr[t]] ~ dbin(pm[t],tm[ptr[t]])
    
  }

  
  ################# PRIORS ########################

  # discretised normal kernel for the cycle duration
 
  pr[1] <- pnorm(1.5,resting,tau.b)
  pr[2] <- pnorm(2.5,resting,tau.b) - pr[1]
  pr[3] <- pnorm(3.5,resting,tau.b) - pr[1] - pr[2]
  pr[4] <- 1 - pr[1] - pr[2] - pr[3]
  
  tau.b ~ dgamma(1,1)    # precision of distribution of cycle duration 

  resting ~ dunif(0.5,4.5) # duration of resting period
  cycle <- resting + (1 - A0)
  Teu ~ dgamma(1,1)       # relative trapping efficiency of gravids relative to unfeds 
  Tem ~ dgamma(1,1)       # relative trapping efficiency of males vs unfed assuming 50:50 sex ratio on emergence
  Pm ~ dunif(0,1)         # daily survival of males
  P ~ dunif(0,1)          # overall survival through oviposition cycle
  A0 ~ dunif(0,1)         # proportion returning to feed within same night
  p <- pow(P,1/cycle)   # daily survival of females
}  
"
  
  exittrapO ="
model{
  for(d in 1:sampleddates){
    # numbers of unfed (+ part fed) per trap (equivalent to emergence rate)
    em[d] <- unfed1[d]/traps[d]    
    # gravid per trap
    eg[d] <- gravid[d]/traps[d]
    # denominators for binomial models
    tg[d] <- unfed1[d] + gravid[d]
    tm[d] <- unfed1[d] + males[d]
  }

  for(t in 1:days_with_complete_data){
  
    #emergence is equivalent to unfed per trap
    emergence[t] <- em[ptr[t]]
    
    # expected number gravid per trap (Teu is relative trapping efficiency of gravid vs unfed)
    # expected gravids in first cycle (with 100% survival)
    eg0[t] <- (em[lag1[t]]*pr[1] + em[lag2[t]]*pr[2] + em[lag3[t]]*pr[3] + em[lag4[t]]*pr[4])* Teu 
    # expected gravids in subsequent cycles (with 100% survival)
    eg1[t] <- eg[lag1[t]]*pr[1] + eg[lag2[t]]*pr[2] + eg[lag3[t]]*pr[3] + eg[lag4[t]]*pr[4]
    eg2[t] <- eg[lag2[t]]*pr[1] + eg[lag3[t]]*pr[2] + eg[lag4[t]]*pr[3] + eg[lag5[t]]*pr[4]

    # overall expected gravids per trap
    G[t] <- (eg0[t] + A0*eg1[t] + (1-A0)*eg2[t])*P
    
    pg[t] <- G[t]/(emergence[t] + G[t])
    gravid1[ptr[t]] ~ dbin(pg[t],tg[ptr[t]])
    
    
    # expected number of males:
    #     Tem is relative trapping efficiency of males vs unfed assuming 50:50 sex ratio on emergence
    #     Pm is daily survival of males
    
    # Expected males per trap 
    Emales[t] <- (Tem*emergence[t] + Pm * males[lag1[t]]/traps[lag1[t]])
    
    pm[t] <- Emales[t]/(emergence[t] + Emales[t])
    males1[ptr[t]] ~ dbin(pm[t],tm[ptr[t]])
    
  }

  
  ################# PRIORS ########################

  # discretised normal kernel for the cycle duration
 
  pr[1] <- pnorm(1.5,resting,tau.b)
  pr[2] <- pnorm(2.5,resting,tau.b) - pr[1]
  pr[3] <- pnorm(3.5,resting,tau.b) - pr[1] - pr[2]
  pr[4] <- 1 - pr[1] - pr[2] - pr[3]
  
  tau.b ~ dgamma(1,1)    # precision of distribution of cycle duration 

  resting ~ dunif(0.5,4.5) # duration of resting period
  cycle <- resting + (1 - A0)
  Teu ~ dgamma(1,1)       # relative trapping efficiency of gravids relative to unfeds 
  Tem ~ dgamma(1,1)       # relative trapping efficiency of males vs unfed assuming 50:50 sex ratio on emergence
  Pm ~ dunif(0,1)         # daily survival of males
  P ~ dunif(0,1)          # overall survival through oviposition cycle
  A0 ~ dunif(0,1)         # proportion returning to feed within same night
  p <- pow(P,1/cycle)   # daily survival of females
}  
"
  
  # Estimation of temperature dependent emergence rates and male survival from exit trap data 
  # Unfed mosquitoes and emergent mosquitoes are equivalent
  exittrapT ="
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
  exittrapU ="
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
  
  listOfModels = list(E = exittrapE, F= exittrapF, G= exittrapG, H= exittrapH,
                      I = exittrapE, J= exittrapF, K= exittrapG, L= exittrapH,
                      M = exittrapE, N= exittrapF, O= exittrapG, P= exittrapT,
                      Q = exittrapU)
  return(listOfModels)
}