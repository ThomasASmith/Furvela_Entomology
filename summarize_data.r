library(dplyr)

get_Furvela_data = function(){

  data_description = function(datalist){
    ################ generates simple descriptive statistics for the Furvela entomology database ################################
    drops = c("Casa","Data","data","Sheet","ISO_Week","Tempo_aqui","Nome_de_fora","Cumprimento","Largura","Preco","period_lado","Name",
              "Observacoes","Age","preco..mil.","PIN","Ano_Nasc","Data_Nasc","number","plate")
    keeps = c("Tempo_aqui","Cumprimento","Largura","Preco","Age","Ano_Nasc")
    animais_summary = mapply(summary,datalist$animais[,!(names(datalist$animais) %in% drops)])
    casa_table      = mapply(table, datalist$casa[,!(names(datalist$casa) %in% drops)])
    casa_summary    = mapply(summary, datalist$casa[,names(datalist$casa) %in% keeps])
    preco_summary   = summary(as.numeric(datalist$casa$Preco[datalist$casa$Preco != 'variavel']))
    morte_table     = mapply(table, datalist$morte[,!(names(datalist$morte) %in% drops)])
    mosquiteiro_summary = mapply(summary, datalist$mosquiteiro[,!(names(datalist$mosquiteiro) %in% drops)])
    collection_table = mapply(table, datalist$mosquito_collection[,!(names(datalist$mosquito_collection) %in% drops)])
    collection_summary = mapply(summary, datalist$mosquito_collection[,!(names(datalist$mosquito_collection) %in% drops)])
    pessao_table = mapply(table, datalist$pessao[,!(names(datalist$pessao) %in% drops)])
    pessao_summary = mapply(summary, datalist$pessao[,names(datalist$pessao) %in% keeps])
    sporozoites_table =mapply(table, datalist$sporozoites[,!(names(datalist$sporozoites) %in% drops)])
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
  
  analysis_of_proportions = function(dissections){
  ################# Interval estimates of proportions #####################
    parous = sum(dissections$Sac + dissections$No.sac)
    dissected = sum(dissections$total_dissected)
    M = parous/dissected
    M_interval = prop.test(x=parous, n=dissected, conf.level=.95, correct=FALSE)
    M_sample = rbinom(1000,dissected,M)/dissected
    sac = sum(dissections$Sac)
    A0 = sac/parous
    A0_sample = rbinom(1000,parous,A0)/parous
    A0_interval = prop.test(x=sac, n=parous, conf.level=.95, correct=FALSE)
    fdenom = 6379+791+5121
    fvalue = 6379/fdenom
    f_sample = rbinom(1000,fdenom,fvalue)/fdenom
    proportions = list(M=M, M_interval=M_interval, M_sample=M_sample, 
                       A0=A0, A0_sample=A0_sample, A0_interval=A0_interval,
                       fvalue = fvalue, f_sample=f_sample)
    return(proportions)}
  
  get_timeseries = function(datalist) {
    # Create a time series with integer dates ensuring that there is at least one record for each date
    alldates = data.frame(date=as.numeric(as.Date(
      datalist$mosquito_collection$date,'%d/%m/%Y')))-as.numeric(as.Date("26/06/2001", format = "%d/%m/%Y")) 
    workfile = datalist$mosquito_collection
    workfile$date = alldates$date
    return(workfile)
  }


  animais = read.csv(file='Animais.csv')
  casa = read.csv(file='Casa.csv')
  morte = read.csv(file='Morte.csv')
  mosquiteiro = read.csv(file='Mosquiteiro.csv')
  mosquito_collection = read.csv2(file='Mosquito_Collection.csv',sep = ",")
  pessao = read.csv(file='Pessao.csv')
  sporozoites=read.csv(file='Sporozoites.csv')
  dissections=read.csv(file='dissections.csv')
  
  ################ Data descriptions
  Furvela_data = list(animais=animais,
                      casa=casa,
                      morte=morte,
                      mosquiteiro=mosquiteiro,
                      mosquito_collection=mosquito_collection,
                      pessao=pessao,
                      sporozoites = sporozoites,
                      dissections = dissections)
  
  Furvela_data$description = data_description(datalist=Furvela_data)
  Furvela_data$proportions = analysis_of_proportions(dissections=Furvela_data$dissections)
  Furvela_data$workfile = get_timeseries(datalist=Furvela_data)
  
  return(Furvela_data)
}

############ functions for summarising data
selectData = function(firstquarter=1,lastquarter=31,model='default', A0=0.6037, P= 0.4812){
  threemonths=365.25/4
  mindate = round((firstquarter-1)*threemonths)
  maxdate = round(lastquarter*threemonths)
  date=workfile$date
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
  start = min(date[!is.na(workfile$Af.unfed)])
  date_offset = date[!is.na(workfile$Af.unfed)] - start + 1
  maxdate = min(maxdate,max(date_offset))
  Af.unfed1 = with(workfile,ifelse(Collection=='exit',Af.unfed+Af.part,NA))
  Af.old = with(workfile,ifelse(Collection=='exit',Af.fed+Af_gravid,NA))
  Af.gravid = with(workfile,ifelse(Collection=='exit',Af_gravid,NA))
  Af.maleEx = with(workfile,ifelse(Collection=='exit',Af.male,NA))
  Af.female = Af.unfed1 + Af.old + workfile$Af.semi
  Af.femaleLT = with(workfile,ifelse(Collection=='Light',Af.female,NA))
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
  df3=data.frame(table(workfile$Casa,workfile$Collection))
  df3$Freq[df3$Freq > 0] = 1
  df4 = as.data.frame(table(df3$Var2,df3$Freq))
  # nights observed
  df3=data.frame(table(workfile$Data,workfile$Collection))
  df3$Freq[df3$Freq > 0] = 1
  df5 = data.frame(table(df3$Var2,df3$Freq))
  workfile$l_Af.unfed = log(workfile$Af.unfed+1)
  workfile$l_Af.part = log(workfile$Af.part+1)
  workfile$l_Af.fed = log(workfile$Af.fed+1)
  workfile$l_Af.semi = log(workfile$Af.semi+1)
  workfile$l_Af.gravid= log(workfile$Af_gravid+1)
  workfile$l_Af.male = log(workfile$Af.male+1)
  df2 = data.frame(workfile %>%
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
  df2$trapnights = data.frame(table(workfile$Collection))$Freq
  return(df2)
}

####### Summarise time series by quarter for plotting
summarise_by_quarter = function(){
  threemonths=365.25/4
  quarter = ceiling(workfile$date/threemonths)
  Af.unfed1 = with(workfile,ifelse(Collection=='exit',log(Af.unfed+Af.part+1),NA))
  Af.gravid = with(workfile,ifelse(Collection=='exit',log(Af_gravid+1),NA))
  Af.maleEx = with(workfile,ifelse(Collection=='Light',log(Af.male+1),NA))
  Af.femaleLT = with(workfile,ifelse(Collection=='Light',log(Af.unfed+Af.part+Af.fed+Af.semi+Af_gravid+1),NA))
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

Furvela_data =  get_Furvela_data()
description = Furvela_data$description
proportions = Furvela_data$proportions
workfile = Furvela_data$workfile

remove('get_Furvela_data')
