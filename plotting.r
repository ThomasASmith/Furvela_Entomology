
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

#Functions for analysis and plotting of correlations
################# ANALYSIS OF CORRELATIONS
plotCorrelations = function(){
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


CorrelationPlot=plotCorrelations()
remove('plotCorrelations')