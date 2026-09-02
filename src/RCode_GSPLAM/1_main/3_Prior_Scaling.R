
################################################################################
################################################################################
#################### Prior scaling #############################################
################################################################################
################################################################################

## In this script we implement prior scaling routines for 
## both, Gamma and scaled beta prime (SBP) hyperpriors.

## this is required to perform prior scaling for NBPSS
prior_scaling_sbp <-function(X,shape1=1/2,shape2=5,th=c(0.1,1),alpha=c(0.1,0.1),T=10^3)
{
  
  X=as.matrix(X)
  d=ncol(X)
  
  tau2=rbetapr(n=T,shape1=shape1,shape2=shape2,scale=1)

  beta=matrix(rnorm(d*T),T,d) 
  
  beta=Diagonal(x=sqrt(tau2))%*%beta ## the rows of beta contain prior coeff draws
  
  F=beta%*%t(X) ## the rows of F contain prior function draws
  
  fsup=apply(abs(F),MARGIN=1,FUN=max)
  
  q=quantile(fsup,1-alpha)
  
  c=(th/q)^2
  c=as.numeric(c)
  
  return(c)
}

## this is required to perform prior scaling for the SSGL prior
## note: this functions returns scale parameters (not rate parameters)
prior_scaling_ga <-function(X,th=c(0.1,1),alpha=c(0.1,0.1),T=10^3)
{
  
  X=as.matrix(X)
  d=ncol(X)
  
  tau2=rgamma(n=T,shape=(d+1)/2,scale=1)
  
  beta=matrix(rnorm(d*T),T,d) 
  
  beta=Diagonal(x=sqrt(tau2))%*%beta ## the rows of beta contain prior coeff draws
  
  F=beta%*%t(X) ## the rows of F contain prior function draws
  
  fsup=apply(abs(F),MARGIN=1,FUN=max)
  
  q=quantile(fsup,1-alpha)
  
  c=(th/q)^2
  c=as.numeric(c)
  
  return(c)
}

