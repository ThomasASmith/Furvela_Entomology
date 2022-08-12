



# The solutions require the list of proportions
if( !exists('proportions') )
{
  source('summarize_data.r')
}

solveroots = function(proportions=proportions){
  # fvalue=0.519
  # M=0.528
  # A0=0.603
  
  library(stats)
  
  theta_solution = function(p, Po, fvalue=proportions$fvalue, A0=proportions$A0){
    theta_o = log(Po)/log(p)
    theta_r = theta_o - (1-A0)/A0
    result = (1-Po^(1/theta_o))/(1-Po^(theta_r/theta_o)) - fvalue
    return(result)}
  
  p_2sample = p_1sample = theta_r_2= as.numeric(rep(NA,1000))
  for(i in 1:1000){
    fvalue = proportions$f_sample[i]
    M = proportions$M_sample[i]
    A0 = proportions$A0_sample[i]
    theta_r_1 = 1/proportions$fvalue
    Po_1 = M
    interval = c(M,0.999)
    p_2 = uniroot(f=theta_solution, interval=interval, Po=Po_1, fvalue=fvalue, A0=A0)$root[1]
    theta_r_2[i]=(log(fvalue-1+p_2)-log(fvalue))/log(p_2)
    p_2sample[i] = p_2
  }
  
    # Parous rate (from Davidson, 1954)
  Po_1sample = proportions$M_sample
  CI_Po_1=quantile(Po_1sample, probs = c(0.5, 0.025, 0.975))
  
  # from equation 1
  theta_r_1 = 1/ proportions$f_sample
  CI_theta_r_1=quantile(theta_r_1, probs = c(0.5, 0.025, 0.975))

  
  # sac rate (from Charlwood, 1985)
  CI_A0=quantile( proportions$A0_sample, probs = c(0.5, 0.025, 0.975))
  
  # from equation 2
  
  theta_o_1 = theta_r_1 + (1 -  proportions$A0_sample)/ proportions$A0_sample
  CI_theta_o_1=quantile(theta_o_1, probs = c(0.5, 0.025, 0.975))
  
  # from equation 3
  p_1sample = exp(log(proportions$M_sample)/theta_o_1)
  p_1 = mean(p_1sample)
  CI_p_1=quantile(p_1sample, probs = c(0.5, 0.025, 0.975))
  
  # from solution of equation 5
  
  CI_theta_r_2=quantile(theta_r_2, probs = c(0.5, 0.025, 0.975))
  CI_p_2=quantile(p_2sample, probs = c(0.5, 0.025, 0.975))
  
  # from solution of equation 5 and equation 2
  
  theta_o_2 = theta_r_2 + (1 - proportions$A0_sample)/proportions$A0_sample
  CI_theta_o_2=quantile(theta_o_2, probs = c(0.5, 0.025, 0.975))
  
  roots = list(p_1 = p_1, 
               p_2=p_2, 
               theta_r_1=mean(theta_r_1),
               theta_o_1=mean(theta_o_1), 
               theta_r_2=mean(theta_r_2),
               theta_o_2=mean(theta_o_2), 
               Po_1 = Po_1, 
               A0 = A0,
               CI_p_1 = CI_p_1, 
               CI_p_2=CI_p_2, 
               CI_theta_r_1=CI_theta_r_1,
               CI_theta_o_1=CI_theta_o_1,
               CI_theta_r_2=CI_theta_r_2,
               CI_theta_o_2=CI_theta_o_2,
               CI_Po_1 = CI_Po_1, 
               CI_A0 = CI_A0)
  return(roots)
}

roots = solveroots(proportions=proportions)
remove('solveroots')

