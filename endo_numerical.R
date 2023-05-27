library(tidyverse)

coho<-readxl::read_xlsx(here::here("data","OC Coho Abundance.xlsx"),skip=1,col_names=TRUE) %>% 
  janitor::clean_names()

coho_abun<-coho %>% 
  column_to_rownames(var="year")



library(nloptr)

means=c(2,3,.5)

sig1=0.5
sig2=1.2
sig3=0.1

var=c(sig1,sig2,sig3)
cov=matrix(data=c(sig1^2,sig1*sig2*0.5,sig1*sig3*0.2,sig1*sig2*0.5,sig2^2,sig2*sig3*0.7,sig1*sig3*0.2,sig2*sig3*0.7,sig3^2),ncol=3,nrow=3)

objective<-function(w,mu=means,sigma=cov,var=var,gamma=1,budget=1){
  port_mu=sum(w*mu)
  
  var1=sum(w^2*var^2)
  
  j_sum=rep(0,length(w))
  
  for(j in 1:length(w)){
    k_sum=rep(0,length(w))
    for(k in 1:length(w)){
      
      
      
      k_sum[k]=w[j]*w[k]*sigma[j,k]
      
      if(k==j){
        k_sum[k]=0
      }
      
    }
    j_sum[j]=sum(k_sum)
  }
  var2=sum(j_sum)
  
  port_var=var1+var2
  utility=port_mu-gamma*port_var
  
  return(-utility)
}


constraints<-function(w,mu,sigma,var,gamma,budget=1){
  const=sum(w)-budget
  
  return(const)
}



local_opts<-list("algorithm"="NLOPT_LN_COBYLA",xtol_rel=1e-15)
options=list("algorithm"="NLOPT_LN_AUGLAG",xtol_rel=1e-15,maxeval=16000,"local_opts"=local_opts)

out=nloptr(x0=c(0.3,0.3,0.4),
       eval_f=objective,
       eval_g_eq = constraints,
       lb=c(0,0,0),
       ub=c(1,1,1),
       opts=options,
       mu=means,
       sigma=cov,
       var=var,
       gamma=1,
       budget=1)


sum(out$solution)

##### Endogenize just in means ####

means=c(2,2,2)

alpha=c(.001,.005,.001)

sig1=0.5
sig2=1.2
sig3=0.1

var=c(sig1,sig2,sig3)
cov=matrix(data=c(sig1^2,sig1*sig2*0.5,sig1*sig3*0.2,sig1*sig2*0.5,sig2^2,sig2*sig3*0.7,sig1*sig3*0.2,sig2*sig3*0.7,sig3^2),ncol=3,nrow=3)

# gamma is risk aversion, budget is how much money we can spend, alpha is the vector of weight improvements
objective_endo<-function(w,mu=means,sigma=cov,var=var,gamma=1,budget=1,alpha){
  port_mu=sum(w*mu*alpha*log(w+1))
  
  var1=sum(w^2*var^2)
  
  j_sum=rep(0,length(w))
  
  for(j in 1:length(w)){
    k_sum=rep(0,length(w))
    for(k in 1:length(w)){
      
      
      
      k_sum[k]=w[j]*w[k]*sigma[j,k]
      
      if(k==j){
        k_sum[k]=0
      }
      
    }
    j_sum[j]=sum(k_sum)
  }
  var2=sum(j_sum)
  
  port_var=var1+var2
  utility=port_mu-gamma*port_var
  
  return(-utility)
}


constraints_endo<-function(w,mu,sigma,var,gamma,budget=1,alpha){
  const=sum(w)-budget
  
  return(const)
}



local_opts<-list("algorithm"="NLOPT_LN_COBYLA",xtol_rel=1e-15)
options=list("algorithm"="NLOPT_LN_AUGLAG",xtol_rel=1e-15,maxeval=16000,"local_opts"=local_opts)

out=nloptr(x0=temp,
           eval_f=objective_endo,
           eval_g_eq = constraints_endo,
           lb=c(0,0,0),
           ub=c(1,1,1),
           opts=options,
           mu=means,
           sigma=cov,
           var=var,
           gamma=0,
           budget=1,
           alpha=alpha)

out$solution
out$objective

## Convergence is dependent on starting values

#### grid search ####
# Make all possible starting combinations
a<-expand_grid(x=seq(0.05,1,by=.05),y=seq(0.05,1,by=.05),z=seq(0.05,1,by=.05)) %>% 
  mutate(sum=rowSums(across(everything()))) %>% 
  filter(sum==1) %>% 
  select(-sum) %>% 
  filter(if_all(everything())>=1)


grid_list<-split(a,seq(nrow(a)))

max_fcn<-function(x){
  temp=slice(x) %>% unlist()
  
  out=nloptr(x0=temp,
             eval_f=objective_endo,
             eval_g_eq = constraints_endo,
             lb=c(0,0,0),
             ub=c(1,1,1),
             opts=options,
             mu=means,
             sigma=cov,
             var=var,
             gamma=0.05,
             budget=1,
             alpha=alpha)
  
  tempsol=out$solution
  
  
  return(round(data.frame(x=out$solution[1],y=out$solution[2],z=out$solution[3],obj=out$objective),5))
}
try=map_df(.x=grid_list,~max_fcn(.x))


c=filter(try,obj==min(try$obj))
