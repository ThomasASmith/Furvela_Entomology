################# Analysis of resting collections ################


# The analysis of resting collections requires the list of proportions
if( !exists('proportions') )
{
  source('summarize_data.r')
}

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

resting_analysis = resting_collection_analysis()
remove(resting_collection_analysis)
