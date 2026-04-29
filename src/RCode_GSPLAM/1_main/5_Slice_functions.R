
#####################################################################################
#####################################################################################
################## Slice functions ##################################################
#####################################################################################
#####################################################################################

### This script contains slice functions which are required for efficient MCMC updates
### of the variance parameters tau2.

####################################################### general functions

### Slice sampling with stepping out (see Neal 2003)

stepping_out<-function(x,theta,logy,w=1,m=30,log_target){
  
  U=runif(1,0,1)
  L=x-w*U
  R=L+w
  
  V=runif(1,0,1)
  J=floor(m*V)
  K=(m-1)-J
  
  while((J>0)&(logy<log_target(L,theta))){
    L=L-w
    J=J-1
  }
  while((K>0)&(logy<log_target(R,theta))){
    R=R+w
    K=K-1
  }
  
  return(c(L,R))
}

## here we preform a single slice step

perform_slice_step<-function(x0,theta,log_target){
  
  logy=log_target(x0,theta)-rexp(n=1,rate=1) ## here we need to evaluate the log-kernel
  
  #d=sqrt(-2*log(u))
  
  fit=stepping_out(x0,theta,logy,log_target=log_target)
  L=fit[1]
  R=fit[2]
  
  paste(L,R)
  
  xtemp=runif(1,L,R)
  while(log_target(xtemp,theta)<logy){ ## sample uniformly from the horizontal slice
    xtemp=runif(1,L,R)  
  }
  
  return(xtemp)
}

slice_sampler<-function(x0,theta,T=10^3){
  
  x=rep(0,T)
  x[1]=x0
  
  for(t in 2:T){
    
    #u=runif(n=1,0,exp(-x[t-1]^2/2)) ## draw from vertical slice
    
    logy=log_target(x[t-1],theta)-rexp(n=1,rate=1) ## here we need to evaluate the log-kernel
    
    #d=sqrt(-2*log(u))
    
    fit=stepping_out(x[t-1],theta,logy)
    L=fit[1]
    R=fit[2]
    
    xtemp=runif(1,L,R)
    while(log_target(xtemp,theta)<logy){ ## sample uniformly from the horizontal slice
      xtemp=runif(1,L,R)  
    }
    
    x[t]=xtemp
    
    #d=-2*log(sqrt(2*pi))-2*logu
    #d=sqrt(-2*logy)
    #L#=-d
    #R=d
    #x[t]=runif(n=1,min=L,max=R) ## draw from horizontal slice
    
    #print(x[t])
    
  }
  
  return(x)
  
}

#################################################### here we need to specify the log-target

## this is the density of the scaled beta prime distribution
dscaled_beta_prime<-function(x,a,b,c){
  return((x/c)^(a-1)*(1+x/c)^(-(a+b))/beta(a,b)*1/c)
}

## this is the log density of the FCP of log tau2
## when tau2 follows a Ga(shape=a,scale=c) prior

log_target_ga<-function(x,theta){
  d=theta[1];q=theta[2];a=theta[3];c=theta[4];
  return(-d/2*x-1/2*q/exp(x)+dgamma(exp(x),shape=a,scale=c,log=TRUE)+x)
}

log_target_sbp<-function(x,theta){
  d=theta[1];q=theta[2];a=theta[3];b=theta[4];c=theta[5]
  return(-d/2*x-1/2*q/exp(x)+dbetapr(exp(x),a,b,c,log=TRUE)+x)
}

######################################################################################
######################################################################################
######################################################################################
