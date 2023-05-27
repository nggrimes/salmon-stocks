## PT optimization

library(tidyverse)
library(quadprog)


## Make up data

means=c(2,3,.5)

sig1=0.5
sig2=1.2
sig3=0.1

cov=matrix(data=c(sig1^2,sig1*sig2*0.5,sig1*sig3*0.2,sig1*sig2*0.5,sig2^2,sig2*sig3*0.7,sig1*sig3*0.2,sig2*sig3*0.7,sig3^2),ncol=3,nrow=3)
dvec=rep(0,3)



## Construct constraint matrix
Amat <- means %>% rbind(c(1,1,1)) %>% rbind(diag(1,nrow=3)) %>% t()

eff_front=function(x){
  bvec<-c(x,1,0,0,0)

  out=solve.QP(cov,dvec,Amat,bvec,meq=2)
  
  return(round(data.frame(x=out$solution[1],y=out$solution[2],z=out$solution[3],var=out$value),5))
}


mu_list<-seq(0.5,3,by=0.01)

front=map_df(mu_list,~eff_front(.x))


### plot

var_fcn<-function(weights,cov){
  
  temp=as.matrix(cov) %*% as.matrix(weights)
  out=t(as.matrix(weights))%*%temp
  return(out)
}

mu_fcn<-function(weights,means){
  out=weights %*% means
  return(out)
}

var_fcn(weight1,cov)


#### Make a bunch of random plots
n=1000
rand<-matrix(0,nrow=n,ncol=5)
for(i in 1:n){
  pull=runif(3)
  rand[i,1:3]<-pull/sum(pull)
  
 
  rand[i,4]<-var_fcn(rand[i,1:3],cov)
  
  rand[i,5]<-mu_fcn(rand[i,1:3],means)
  
}

colnames(rand)<-c("weight_x","weight_y","weight_z","var","mu")

ggplot(as.data.frame(rand),aes(x=var,y=mu))+
  geom_point()+
  geom_point(aes(x=cov[1,1],y=means[1]),color="blue",size=3)+
  geom_point(aes(x=cov[2,2],y=means[2]),color="green",size=3)+
  geom_point(aes(x=cov[3,3],y=means[3]),color="red",size=3)+
  theme_classic()

