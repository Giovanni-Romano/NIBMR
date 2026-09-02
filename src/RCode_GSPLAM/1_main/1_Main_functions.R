
################################################################################
################################################################################
################# Main functions ###############################################
################################################################################
################################################################################

## This script contains the main functions.

## The key function is "nbpsscox", which implements a MCMC sampler for NBPSS in the geoadditive Cox model.

## The other functions in this script are plotting routines required 
## to visualize the estimated log baseline hazard rate,
## the effects of the continuous covariates or
## the spatial effect (if present).

## => Run the whole script.

################################################################################
################################################################################
################################################################################

## load required packages

library(mgcv)
library(geigen)
library(posterior)
library(Matrix)
library(mvtnorm)
library(acid)
library(bayesplot)
library(spBayesSurv)
library(extraDistr)
library(survival)
library(geometry)
library(xtable)
library(mboost)
library(GIGrvg)
library(latex2exp)

## load other scripts

source("1_main/2_Design_functions.R")
source("1_main/3_Prior_Scaling.R")
source("1_main/4_IWLS_functions.R")
source("1_main/5_Slice_functions.R")

################################################################################
################################################################################
################################################################################

## This function implements Bayesian effect selection for the geoadditive
## Cox regression model

## y are the observed time points (note: they need to be sorted increasingly!)
## delta are the event indicators (delta=1 means that the event has occurred)
## Xdummy contains dummy covariates
## Xcont contains continuous covariates
## Xspat contains continuous spatial information
## nonlin indicates the basis expansion used for the nonlinear effects of continuous covariates
## nonlin can take three values: "DR", "MM", "Eilers". For other values (e.g. "none") there is no nonlinear effect.
## prior can take two values: "SBP" corresponds to NBPSS and "GA" corresponds to the SSGL

nbpsscox<-function(y,delta,Xdummy=NULL,Xcont=NULL,Xspat=NULL,nonlin="DR",blockwise=FALSE,globalw=FALSE,
                 prior="SBP",a=1/2,b=5,T=10^3,burn=100,optim=100,d0=10,d1=9,ngrid=10^3,th=c(0.1,1),w1=1,w2=1,gamma.ini=NULL)
{
  
  ## 1. the preparator #########################################################
  
  yorig=y
  rescale=sum(delta)/sum(y)
  y=y*rescale ## rescale time points to avoid numerical issues
  
  ## 2. the designer ###########################################################
  
  ## we use P-splines for g0 (log baseline hazard) and DR splines for f(x)
  
  #d0=10
  
  ytemp=c(0,y)
  fit=smoothCon(s(ytemp,bs="ps",k=d0),data=data.frame(ytemp))
  #fit=smoothCon(s(ytemp,bs=bs,k=d0),data=data.frame(ytemp))
  X0=fit[[1]]$X[-1,]
  K0=fit[[1]]$S[[1]]
  
  ## the following stuff is needed to obtain a good approximation for the integrals
  ## that appear in the likelihood and its derivatives & for nice predictions of the baseline hazard
  
  fine=fineGrid(y,length=ngrid,d0=d0,bs="ps")
  X0fine=fine[[1]];ynew=fine[[4]];
  
  # beta0=rnorm(d0)
  # plot(y,X0%*%beta0)
  # lines(ynew,X0fine%*%beta0,col="red")
  
  #### set up covariate design matrix X
  
  if(is.null(Xdummy)){
    X=scale(Xcont)
  } else {
    X=cbind(Xdummy,scale(Xcont))
  }
  
  indices=as.list(1:ncol(X)) ## indices of the groups
  J=length(indices) ## number of groups
  
  if(nonlin=="DR"){ ### nonlinear effects of continuous covariates
    Trafo=list() 
    for(j in 1:ncol(Xcont)){
      DR=DemmlerReinsch(Xcont[,j],d=d1)
      Trafo[[j]]=DR$Trafo
      X=cbind(X,DR$X)
      indices=c(indices,list(J+(j-1)*d1+1:d1))
    }
  }
  
  if(nonlin=="MM"){ ### nonlinear effects of continuous covariates
    Trafo=list() 
    for(j in 1:ncol(Xcont)){
      MM=MixedModel(Xcont[,j],d=d1)
      Trafo[[j]]=MM$Trafo
      X=cbind(X,MM$X)
      indices=c(indices,list(J+(j-1)*d1+1:d1))
    }
  }
  
  if(nonlin=="Eilers"){ ### nonlinear effects of continuous covariates
    Trafo=list() 
    for(j in 1:ncol(Xcont)){
      MM=Eilers(Xcont[,j],d=d1)
      Trafo[[j]]=MM$Trafo
      X=cbind(X,MM$X)
      indices=c(indices,list(J+(j-1)*d1+1:d1))
    }
  }
  
  Trafos=list()
  if(!is.null(Xspat)){
    DRS=DemmlerReinschSpatial(Xspat[,1],Xspat[,2],k1=6,k2=6)
    Trafos=DRS$Trafo
    X=cbind(X,DRS$X)
    ds=ncol(DRS$X)
    indices=c(indices,list(max(unlist(indices))+1:ds))
  }
  
  J=length(indices) ## number of groups
  d=ncol(X) ## dimension of overall design matrix
  dim=unlist(lapply(indices,length))
  
  #### 3. the initializer ########################################################
  
  eps=10^-3
  beta0.ini=rep(0,d0); tau20.ini=100; ## start with large values (otherwise we may be stuck)
  beta.ini=rep(0,d); tau2.ini=rep(100,J);
  
 Beta.ini=matrix(0,optim,d)
  
  for(t in 1:optim){
    fit=evaluate_log_FCP_beta0(beta0.ini,beta.ini,X0,X,y,delta,K0,tau20.ini,fine) ## use some Newton-Raphson steps first
    beta0.ini=beta0.ini-solve(fit[[3]])%*%fit[[2]]
    
    tau20.ini=as.numeric((eps+t(beta0.ini)%*%K0%*%beta0.ini/2)/(eps+(d0-2)/2))
    
    fit=evaluate_log_FCP_beta(beta0.ini,beta.ini,X0,X,y,delta,Omega=diag(rep(1/tau2.ini,dim)),fine) 
    beta.ini=beta.ini-solve(fit[[3]])%*%fit[[2]]
    
    Beta.ini[t,]=beta.ini
    
    for(j in 1:J){
    tau2.ini[j]=(eps+sum(beta.ini[indices[[j]]]^2)/2)/(eps+dim[j]/2)
    }
  }
  
  #### 4. the prior scaler #######################################################
  
  #a=1/2;b=5;#5; ## shape parameters of tau2j
  w1=w1; w2=w2; ## shape parameters of theta/thetaj
  
  #c=c(0.1,1); ## scale parameters of tau2j
  
  #c=matrix(0,J,2); ## scale parameters of tau2j
  #c[,1]=c0
  #c[,2]=c1
  
  c=matrix(0,J,2)
  for(j in 1:J){
    if(prior=="SBP"){c[j,]=prior_scaling_sbp(X[,indices[[j]]],th=th)}
    if(prior=="GA"){c[j,]=prior_scaling_ga(X[,indices[[j]]],th=th)}
  }
  
  #### 5. the MCMC sampler #######################################################
  
  #T=2*10^2
  beta0=matrix(0,T,d0); beta0[1,]=beta0.ini
  beta=matrix(0,T,d); beta[1,]=beta.ini
  
  tau20=rep(0,T);tau20[1]=tau20.ini; 
  tau2=matrix(0,T,J);tau2[1,]=tau2.ini;
  
  #gamma=matrix(1,T,J)
  gamma=matrix(0,T,J)
  
  if(!is.null(gamma.ini)){gamma[1,]=gamma.ini}
    
  probi=matrix(0,T,J) ## store transition probabilities for Rao-Blackwellized
  
  if(globalw==FALSE){theta=matrix(0.5,T,J)} ## local wjs
  if(globalw==TRUE){theta=rep(0.5,T)} ## global w
  
  t0=Sys.time()
  
  for(t in 1:(T-1)){
    
    print(t)
    
    ## update beta0
    
    beta0[t+1,]=perform_iwls_step_beta0(beta0[t,],beta[t,],X0,X,y,delta,K0,tau20[t],fine,step=1)
    #print(beta0[t+1,])
    
    ## update tau20
    
    tau20[t+1]=exp(perform_slice_step(log(tau20[t]),theta=c(d0-2,t(beta0[t+1,])%*%K0%*%beta0[t+1,],1/2,1/2,1),log_target_sbp)) ## smoothing variance of log baseline hazard
    
    ## update beta
    
    Omega=diag(rep(1/tau2[t,],dim)) ## current prior precision matrix
    beta[t+1,]=perform_iwls_step_beta(beta0[t+1,],beta[t,],X0,X,y,delta,Omega,fine,step=1)
    #print(beta[t+1,])
    
    ## update tau2j and the selection indicators gammaj
    
    for(j in 1:J){
      
      if(blockwise==TRUE){ ## update tau2j and gammaj in a block
        spike=dmargNSBP(beta[t+1,indices[[j]]],a,b,c[j,1],dim[j])
        slab=dmargNSBP(beta[t+1,indices[[j]]],a,b,c[j,2],dim[j])
        ratio=spike/slab
        print(ratio)
        probi[t+1,j]=1/(1+ratio*(1-theta[t,j])/theta[t,j])
        gamma[t+1,j]=rbinom(n=1,size=1,prob=probi[t+1,j]) 
        tau2[t+1,j]=exp(perform_slice_step(log(tau2[t,j]),c(dim[j],sum(beta[t+1,indices[[j]]]^2),a,b,c[j,gamma[t+1,j]+1]),log_target_sbp)) ## smoothing variances of covariate effects
        
      } else {
        
      if(prior=="SBP"){
        tau2[t+1,j]=exp(perform_slice_step(log(tau2[t,j]),c(dim[j],sum(beta[t+1,indices[[j]]]^2),a,b,c[j,gamma[t,j]+1]),log_target_sbp)) ## smoothing variances of covariate effects
        spike=dbetapr(tau2[t+1,j],a,b,c[j,1],log=TRUE)
        slab=dbetapr(tau2[t+1,j],a,b,c[j,2],log=TRUE)
        ratio=exp(spike-slab)
      }
      if(prior=="GA"){
        #param=c(dim[j],sum(beta[t+1,indices[[j]]]^2),(dim[j]+1)/2,c[j,gamma[t,j]+1])
        #print(param)
        #tau2[t+1,j]=exp(perform_slice_step(log(tau2[t,j]),c(dim[j],sum(beta[t+1,indices[[j]]]^2),(dim[j]+1)/2,c[j,gamma[t,j]+1]),log_target_ga)) ## smoothing variances of covariate effects 
        tau2[t+1,j]=rgig(n=1,lambda=1/2,chi=sum(beta[t+1,indices[[j]]]^2),psi=2/c[j,gamma[t,j]+1]) ## smoothing variances of covariate effects 
        spike=dgamma(tau2[t+1,j],shape=(dim[j]+1)/2,scale=c[j,1],log=TRUE)
        slab=dgamma(tau2[t+1,j],shape=(dim[j]+1)/2,scale=c[j,2],log=TRUE)
        ratio=exp(spike-slab)
      }
      
      if(globalw==FALSE){probi[t+1,j]=1/(1+ratio*(1-theta[t,j])/theta[t,j])} ## local wjs  
      if(globalw==TRUE){probi[t+1,j]=1/(1+ratio*(1-theta[t])/theta[t])} ## global w
      
      gamma[t+1,j]=rbinom(n=1,size=1,prob=probi[t+1,j]) 
      }
      
      if(globalw==FALSE){theta[t+1,j]=rbeta(n=1,shape1=w1+gamma[t+1,j],shape2=w2+1-gamma[t+1,j])} ## local wjs
      #theta[t+1,j]=1/2
    }
    
    ## update the global complexity parameter theta
    
    if(globalw==TRUE){theta[t+1]=rbeta(n=1,shape1=w1+sum(gamma[t+1,]),shape2=w2+J-sum(gamma[t+1,]))}
    #theta[t+1]=0.5
    #theta[t+1]=1 ## no selection
  }
  
  t1=Sys.time()
  t1-t0
  
  ## 6. the postprocessor ########################################################
  
  if(burn>0){
    beta0=beta0[-c(1:burn),] ### throw away burn-in samples
   # beta0=beta0+log(rescale) ### 
    tau20=tau20[-c(1:burn)]
    beta=beta[-c(1:burn),]
    tau2=tau2[-c(1:burn),]
    gamma=gamma[-c(1:burn),]
    probi=probi[-c(1:burn),]
    
    if(globalw==FALSE){theta=theta[-c(1:burn),]} ## local wjs
    if(globalw==TRUE){theta=theta[-c(1:burn)]} ## global
  }
  
  ## 7. the returner ###########################################################
  
  ret=list()
  
  samples=list(beta0,tau20,beta,tau2,gamma,theta,probi,Beta.ini)
  names(samples)=c("beta0","tau20","beta","tau2","gamma","theta","probi","Beta.ini")
  ret[[1]]=samples
  
  design0=list(X0,X0fine,yorig,y,ynew)
  names(design0)=c("X0","X0fine","yorig","y","ynew")
  ret[[2]]=design0
  
  design=list(X,J,indices,dim,Xdummy,Xcont,Trafo,Trafos,c)
  names(design)=c("X","J","indices","dim","Xdummy","Xcont","Trafo","Trafos","c")
  ret[[3]]=design
  
  names(ret)=c("samples","design0","design")
  
  return(ret)
}

###################################################################################
###################################################################################
###################################################################################

################ Plotting routines ################################################


################ this function plots the estimated log baseline hazard rate

plot_hazard_rate<-function(fit,ylim=c(-5,5))
{
  
  X0=fit$design0$X0;
  X0fine=fit$design0$X0fine;
  ynew=fit$design0$ynew;
  beta0=fit$samples$beta0;
  
  F=beta0%*%t(X0fine)
  plot(ynew*mean(y)/mean(delta),apply(F,2,FUN=mean),type="l",ylim=ylim,main="log baseline hazard rate",ylab="",xlab="")
  lines(ynew*mean(y)/mean(delta),apply(F,2,FUN=quantile,prob=0.025),lty="dashed")
  lines(ynew*mean(y)/mean(delta),apply(F,2,FUN=quantile,prob=1-0.025),lty="dashed")
  
  rug(fit$design0$yorig)
}

################ this function plots the estimated effect for a single covariate

plot_effect<-function(fit,ylim=c(-1.2,1.2))
{
  
  # x=fit$design$Xcont
  # X=fit$design$X;
  # 
  # F=beta%*%t(X)
  # 
  # f.hat=apply(F,2,FUN=mean)
  # low=apply(F,2,FUN=quantile,prob=0.025)
  # up=apply(F,2,FUN=quantile,prob=1-0.025)
  # 
  # plot(x,f.hat,type="p",main="covariate effect")
  # lines(x,low,lty="dashed",type="p")
  # lines(x,up,lty="dashed",type="p")

  ###################### now on fine grid
  
  beta=fit$samples$beta;
  x=fit$design$Xcont 
  t=seq(min(x),max(x),length=200)
  BT=smoothCon(s(t,bs="ps",k=9+2),data=data.frame(t))[[1]]$X
  
  m=mean(x)
  s=sd(x)
  
  ###################################### overall effect
  
  XT=BT%*%fit$design$Trafo[[1]]
  XT=cbind((t-m)/s,XT)
  
  F=beta%*%t(XT) ## the rows contain posterior function draws 
  
  f.hat=apply(F,2,FUN=mean)
  low=apply(F,2,FUN=quantile,prob=0.025)
  up=apply(F,2,FUN=quantile,prob=1-0.025)
  
  plot(t,f.hat,type="l",main="overall effect",ylim=ylim)
  lines(t,low,lty="dashed")
  lines(t,up,lty="dashed")
  
  
  ##################################### nonlinear effect
  
  XT=BT%*%fit$design$Trafo[[1]]
  
  F=beta[,-1]%*%t(XT) ## the rows contain posterior function draws 
  
  f.hat=apply(F,2,FUN=mean)
  low=apply(F,2,FUN=quantile,prob=0.025)
  up=apply(F,2,FUN=quantile,prob=1-0.025)
  
  plot(t,f.hat,type="l",main="nonlinear effect",ylim=ylim)
  lines(t,low,lty="dashed")
  lines(t,up,lty="dashed")
  
  ##################################### linear effect
  
  XT=(t-m)/s
  
  F=beta[,1]%*%t(XT) ## the rows contain posterior function draws 
  
  f.hat=apply(F,2,FUN=mean)
  low=apply(F,2,FUN=quantile,prob=0.025)
  up=apply(F,2,FUN=quantile,prob=1-0.025)
  
  plot(t,f.hat,type="l",main="linear effect",ylim=ylim)
  lines(t,low,lty="dashed")
  lines(t,up,lty="dashed")
  
}

plot_effect_2<-function(x,beta,Trafo,X=NULL,ylim=c(-1.2,1.2),name="",xlab="",ylab="")
{
  
  # x=fit$design$Xcont
  # X=fit$design$X;
  # 
  # F=beta%*%t(X)
  # 
  # f.hat=apply(F,2,FUN=mean)
  # low=apply(F,2,FUN=quantile,prob=0.025)
  # up=apply(F,2,FUN=quantile,prob=1-0.025)
  # 
  # plot(x,f.hat,type="p",main="covariate effect")
  # lines(x,low,lty="dashed",type="p")
  # lines(x,up,lty="dashed",type="p")
  
  ###################### now on fine grid
  
  #beta=fit$samples$beta;
  #x=fit$design$Xcont 
  t=seq(min(x),max(x),length=200)
  BT=smoothCon(s(t,bs="ps",k=9+2),data=data.frame(t))[[1]]$X
  
  m=mean(x)
  s=sd(x)
  
  
  ##################################### linear effect
  
  XT=(t-m)/s
  
  F=beta[,1]%*%t(XT) ## the rows contain posterior function draws 
  
  f.hat=apply(F,2,FUN=mean)
  low=apply(F,2,FUN=quantile,prob=0.025)
  up=apply(F,2,FUN=quantile,prob=1-0.025)
  
  plot(t,f.hat,type="l",main=paste("linear effect",name,sep=""),ylim=ylim,xlab="",ylab="")
  lines(t,low,lty="dashed")
  lines(t,up,lty="dashed")
  rug(x)
  
  title(ylab="linear effect", line=2.3)
  title(xlab=xlab, line=2.6)
 
  ##################################### nonlinear effect
  
  XT=BT%*%Trafo
  
  F=beta[,-1]%*%t(XT) ## the rows contain posterior function draws 
  
  f.hat=apply(F,2,FUN=mean)
  low=apply(F,2,FUN=quantile,prob=0.025)
  up=apply(F,2,FUN=quantile,prob=1-0.025)
  
  plot(t,f.hat,type="l",main=paste("nonlinear effect",name,sep=""),ylim=ylim,xlab="",ylab="")
  lines(t,low,lty="dashed")
  lines(t,up,lty="dashed")
  rug(x)
  
  title(ylab="nonlinear effect", line=2.3)
  title(xlab=xlab, line=2.6)
  
  ###################################### overall effect
  
  XT=BT%*%Trafo
  XT=cbind((t-m)/s,XT)
  
  F=beta%*%t(XT) ## the rows contain posterior function draws 
  
  f.hat=apply(F,2,FUN=mean)
  low=apply(F,2,FUN=quantile,prob=0.025)
  up=apply(F,2,FUN=quantile,prob=1-0.025)
  
  plot(t,f.hat,type="l",main=paste("overall effect",name,sep=""),ylim=ylim,xlab="",ylab="")
  lines(t,low,lty="dashed")
  lines(t,up,lty="dashed")
  rug(x)
  
  title(ylab="overall effect", line=2.3)
  title(xlab=xlab, line=2.6)
  
  if(!is.null(X)){
    abline(lm(X%*%colMeans(beta)~x),col="red",lwd=2,lty="dashed") 
  }
}

plot_effect_spat<-function(Xspat,beta,Trafo,X=NULL,xlim=c(-1.2,1.2),ylim=c(-1.2,1.2),name="",ncol=30,ngrid=100,cex=cex)
{
  
  ###################### now on fine grid
  
  #beta=fit$samples$beta;
  #x=fit$design$Xcont 
  
  x1=Xspat[,1]
  x2=Xspat[,2]
  
  t1=seq(min(x1),max(x1),length=ngrid)
  t2=seq(min(x2),max(x2),length=ngrid)
  T=expand.grid(t1,t2)
  t1n=T[,1]
  t2n=T[,2]
  
  BT1=smoothCon(s(t1n,bs="ps",k=6),data=data.frame(t1n))[[1]]$X
  BT2=smoothCon(s(t2n,bs="ps",k=6),data=data.frame(t2n))[[1]]$X
  
  BT=tensor.prod.model.matrix(list(BT1,BT2))

  ##################################### overall effect
  
  XT=BT%*%Trafo

  #F=beta%*%t(XT) ## the rows contain posterior function draws 
  #f.hat=apply(F,2,FUN=mean)
  
  f.hat=XT%*%colMeans(beta)
  #low=apply(F,2,FUN=quantile,prob=0.025)
  #up=apply(F,2,FUN=quantile,prob=1-0.025)
  
  #scatter3D(t1n,t2n,f.hat,type="l",main="overall effect")
  fhat=matrix(f.hat,ngrid,ngrid)
  
  ch=convhulln(Xspat)
  inside=inhulln(ch,p=as.matrix(T))
  fhat[!inside]=NA
  image2D(t1,t2,z=fhat,contour=FALSE,main="Spatial effect",col=hcl.colors(ncol, "Temps"),xlab="Easting",ylab="Northing",xlim=xlim,ylim=ylim)
  points(x1,x2,pch=16,cex=cex)
 
  
  #persp3D(t1,t2,z=fhat)
}

################ this function plots the estimated overall effects for all continuous covariates

plot_continuous_effects<-function(fit,ylim=c(-1,1))
{
  
  X=fit$design$X
  indices=fit$design$indices
  dim=fit$design$dim
  beta=fit$samples$beta
  Xcont=fit$design$Xcont
  
  if(is.null(ncol(fit$design$Xdummy))){
    nd=0
  } else {
    nd=ncol(fit$design$Xdummy) ## number of dummy covariates  
  }
  
  nc=ncol(fit$design$Xcont) ## number of continuous covariates
  
  for(j in 1:nc){
    
    # xtemp=Xcont[,j]
    # betatemp=cbind(beta[,j+nd],beta[,indices[[j+(nd+nc)]]])
    # Xtemp=cbind(X[,j+nd],X[,indices[[j+(nd+nc)]]])
    # 
    # F=betatemp%*%t(Xtemp)
    # 
    # f.hat=apply(F,2,FUN=mean)
    # low=apply(F,2,FUN=quantile,prob=0.025)
    # up=apply(F,2,FUN=quantile,prob=1-0.025)
    # 
    # plot(xtemp,f.hat,ylim=ylim,type="p")
    # points(xtemp,low,lty="dashed")
    # points(xtemp,up,lty="dashed")
    
    x=Xcont[,j]
    t=seq(min(x),max(x),length=200)
    
    BT=smoothCon(s(t,bs="ps",k=dim[j+nd+nc]+2),data=data.frame(t))[[1]]$X
    
    m=mean(x)
    s=sd(x)
    
    XT=BT%*%fit$design$Trafo[[j]]
    XT=cbind((t-m)/s,XT)
    
    betatemp=cbind(beta[,j+nd],beta[,indices[[j+(nd+nc)]]])
    
    F=betatemp%*%t(XT)
    
    f.hat=apply(F,2,FUN=mean)
    low=apply(F,2,FUN=quantile,prob=0.025)
    up=apply(F,2,FUN=quantile,prob=1-0.025)
    
    plot(t,f.hat,ylim=ylim,type="l")
    lines(t,low,lty="dashed")
    lines(t,up,lty="dashed")
    
  }
}

###################################################################################
###################################################################################
###################################################################################
