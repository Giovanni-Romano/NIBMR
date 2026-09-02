
################################################################################
################################################################################
################################################################################

## Basis comparison (DRB vs. MMR vs. Eilers) for Gaussian scatterplot smoothing.

## In this script we create Figure 3 in the manuscript.

## this function returns the design matrix for the centered truncated
## power series basis functions used by Hu et al. (2015).
truncated_power_series<-function(x,d=20,center=TRUE){

  Z=cbind(x^2,x^3)  
  
  knots=seq(min(x),max(x),length=d)[-c(1,d)]
  
  for(j in 1:length(knots)){
    Z=cbind(Z,(pmax(x-knots[j],0))^3)
  }
  
  if(center==TRUE){
    for(j in 1:d){
      Z[,j]=Z[,j]-mean(Z[,j])
    }  
  }
  return(Z)
}

## generate data

n=2000

set.seed(123)

#x=sort(rbeta(n,2,2))
x=sort(runif(n));

f=sin(2*pi*x); ## true function
f=x^2

flin=-6/pi*(x-1/2) ## true linear effect (compute L2-projection of f onto span(1,x))
flin=x-1/6

fnonlin=sin(2*pi*x)+6/pi*(x-1/2) ## true nonlinear effect (compute residual from projection)
fnonlin=f-flin

y=f+rnorm(n,sd=0.1)

# par(mfrow=c(1,1))
# plot(x,y)

## set up design matrices for estimation

## y=X*beta+Z*u+eps

X=cbind(rep(1,n),scale(x)) ## linear effect

d=20
Z=DemmlerReinsch(x,d=d)$X ## nonlinear effect
#Z=MixedModel(x,d=d)$X ## nonlinear
#Z=truncated_power_series(x,d=d) ## nonlinear
#Z=smoothCon(s(x,bs="ps",k=d+2),data=data.frame(x))[[1]]$X[,-c(1,d+2)]
#Z=Eilers(x,d=d)$X ## nonlinear

C=cbind(X,Z) ## complete design matrix 

## estimate f and sigma2 using a Gibbs sampler #################################

a=0.001 ## IG hyperparameters
b=0.001

T=10^3

G=t(C)%*%C
K=diag(c(0,0,rep(1,d)))

rk=sum(K)

sigma2=rep(1,T)
tau2=rep(1,T)
rho2=rep(1,T)
beta=matrix(0,T,d+2)

for(t in 1:(T-1)){
  print(t)
  
  Sigma=solve(G/sigma2[t]+K/tau2[t])
  
  mu=Sigma%*%t(C)%*%y/sigma2[t]
  
  beta[t+1,]=rmvnorm(n=1,mu,Sigma)
  
  tau2[t+1]=rinvgamma(n=1,a+rk/2,b+t(beta[t+1,])%*%K%*%beta[t+1,]/2)
  
  RSS=sum((y-C%*%beta[t+1,])^2)
  
  sigma2[t+1]=rinvgamma(n=1,n/2,RSS/2)
  
}

beta=beta[-c(1:100),] ## remove burn-in

##### overall estimate

par(mfrow=c(1,1))
#par(mar=c(2,2,2,2))

##### linear effect

plot(x,y,main="linear effect")
Flin=beta[,c(1,2)]%*%t(X)
fhatlin=colMeans(Flin)
lines(x,fhatlin)
lines(x,apply(Flin,2,quantile,probs=0.975))
lines(x,apply(Flin,2,quantile,probs=0.025))
lines(x,flin,col="red",lwd=2)

##### nonlinear effect

plot(x,y,main="Nonlinear effect")
Fnonlin=beta[,-c(1,2)]%*%t(Z)
fhatnonlin=colMeans(Fnonlin)
lines(x,fhatnonlin)
lines(x,apply(Fnonlin,2,quantile,probs=0.975))
lines(x,apply(Fnonlin,2,quantile,probs=0.025))
lines(x,fnonlin,col="red",lwd=2)

### overall effect

plot(x,y,main="Overall effect")
F=beta%*%t(C)
fhat=colMeans(F)
lines(x,fhat)
lines(x,apply(F,2,quantile,probs=0.975))
lines(x,apply(F,2,quantile,probs=0.025))
lines(x,f,col="red",lwd=2)

## compute RMSEs

sqrt(mean((fhatlin-flin)^2))
sqrt(mean((fhatnonlin-fnonlin)^2))
sqrt(mean((fhat-f)^2))

################################################################################
################################################################################
#### repeat many times and store results #######################################
################################################################################
################################################################################

set.seed(123)

R=50 ## number of replicate data sets
T=5*10^2 ## chain length
n=500 ## sample size
#n=5000

storef=matrix(0,n,R) ## truth 
storeflin=matrix(0,n,R)
storefnonlin=matrix(0,n,R)

storefhat=matrix(0,n,R) ## estimates
storefhatlin=matrix(0,n,R)
storefhatnonlin=matrix(0,n,R)

storefhatMMR=matrix(0,n,R) ## estimates
storefhatlinMMR=matrix(0,n,R)
storefhatnonlinMMR=matrix(0,n,R)

storefhatEil=matrix(0,n,R) ## estimates
storefhatlinEil=matrix(0,n,R)
storefhatnonlinEil=matrix(0,n,R)

storefhatTrunc=matrix(0,n,R) ## estimates
storefhatlinTrunc=matrix(0,n,R)
storefhatnonlinTrunc=matrix(0,n,R)

storefhatB=matrix(0,n,R) ## estimates
storefhatlinB=matrix(0,n,R)
storefhatnonlinB=matrix(0,n,R)

storeX=matrix(0,n,R)

for (r in 1:R){
  
  print(r)
  
  ## generate data
  
  x=sort(runif(n))
  #x=sort(rbeta(n,2,2))
  
  storeX[,r]=x
  
  f=sin(2*pi*x); ## true function
  #f=sin(2*pi*x)+6/pi*(x-1/2)
  #f=x
  #f=x^2
  
  flin=(-6/pi)*(x-1/2) ## true linear effect (compute L2-projection of f onto span(1,x))
  #flin=0*x
  #flin=x
  #flin=(-90/pi^3)*(x-1/2)
  #flin=x-1/6
  
  fnonlin=f-flin  ## true nonlinear effect (compute residual from projection)
  
  y=f+rnorm(n,sd=1)
  
  ### store relevant quantities
  
  storef[,r]=f
  storeflin[,r]=flin
  storefnonlin[,r]=fnonlin
  
  ## set up design matrices for estimation (use DRB first)
  
  ## y=X*beta+Z*u+eps
  
  X=cbind(rep(1,n),scale(x)) ## linear effect
  
  d=20
  Z=DemmlerReinsch(x,d=d)$X ## nonlinear effect
  #Z=MixedModel(x,d=d)$X ## nonlinear
  #Z=Eilers(x,d=d)$X ## nonlinear
  
  C=cbind(X,Z) ## complete design matrix 
  
  ## estimate f and sigma2 using a Gibbs sampler #################################
  
  a=0.001 ## IG hyperparameters
  b=0.001
  
  G=t(C)%*%C
  K=diag(c(0,0,rep(1,d)))
  
  rk=sum(K)
  
  sigma2=rep(1,T)
  tau2=rep(1,T)
  rho2=rep(1,T)
  beta=matrix(0,T,d+2)
  
  for(t in 1:(T-1)){
    #print(t)
    
    Sigma=solve(G/sigma2[t]+K/tau2[t])
    
    mu=Sigma%*%t(C)%*%y/sigma2[t]
    
    beta[t+1,]=rmvnorm(n=1,mu,Sigma)
    
    tau2[t+1]=rinvgamma(n=1,a+rk/2,b+t(beta[t+1,])%*%K%*%beta[t+1,]/2)
    
    RSS=sum((y-C%*%beta[t+1,])^2)
    
    sigma2[t+1]=rinvgamma(n=1,n/2,RSS/2)
    
  }
  
  beta=beta[-c(1:100),] ## remove burn-in
  betahat=colMeans(beta)
  
  fhatlin=X%*%betahat[c(1,2)] ##### linear 
  fhatnonlin=Z%*%betahat[-c(1,2)] ## nonlinear 
  fhat=C%*%betahat ## overall
  
  storefhat[,r]=fhat
  storefhatlin[,r]=fhatlin
  storefhatnonlin[,r]=fhatnonlin
  
  ########################### repeat using MMR for nonlinear effect
  
  X=cbind(rep(1,n),scale(x)) ## linear effect
  
  d=20
  #Z=DemmlerReinsch(x,d=d)$X ## nonlinear effect
  Z=MixedModel(x,d=d)$X ## nonlinear
  #Z=Eilers(x,d=d)$X ## nonlinear
  
  C=cbind(X,Z) ## complete design matrix 
  
  ## estimate f and sigma2 using a Gibbs sampler #################################
  
  a=0.001 ## IG hyperparameters
  b=0.001
  
  G=t(C)%*%C
  K=diag(c(0,0,rep(1,d)))
  
  rk=sum(K)
  
  sigma2=rep(1,T)
  tau2=rep(1,T)
  rho2=rep(1,T)
  beta=matrix(0,T,d+2)
  
  for(t in 1:(T-1)){
    #print(t)
    
    Sigma=solve(G/sigma2[t]+K/tau2[t])
    
    mu=Sigma%*%t(C)%*%y/sigma2[t]
    
    beta[t+1,]=rmvnorm(n=1,mu,Sigma)
    
    tau2[t+1]=rinvgamma(n=1,a+rk/2,b+t(beta[t+1,])%*%K%*%beta[t+1,]/2)
    
    RSS=sum((y-C%*%beta[t+1,])^2)
    
    sigma2[t+1]=rinvgamma(n=1,n/2,RSS/2)
    
  }
  
  beta=beta[-c(1:100),] ## remove burn-in
  betahat=colMeans(beta)
  
  fhatlin=X%*%betahat[c(1,2)] ##### linear 
  fhatnonlin=Z%*%betahat[-c(1,2)] ## nonlinear 
  fhat=C%*%betahat ## overall
  
  storefhatMMR[,r]=fhat
  storefhatlinMMR[,r]=fhatlin
  storefhatnonlinMMR[,r]=fhatnonlin
  
  ########################### repeat using Eilers for nonlinear effect
  
  X=cbind(rep(1,n),scale(x)) ## linear effect
  
  d=20
  #Z=DemmlerReinsch(x,d=d)$X ## nonlinear effect
  #Z=MixedModel(x,d=d)$X ## nonlinear
  Z=Eilers(x,d=d)$X ## nonlinear
  #Z=truncated_power_series(x,d=d)
  
  C=cbind(X,Z) ## complete design matrix 
  
  ## estimate f and sigma2 using a Gibbs sampler #################################
  
  a=0.001 ## IG hyperparameters
  b=0.001
  
  G=t(C)%*%C
  K=diag(c(0,0,rep(1,d)))
  
  rk=sum(K)
  
  sigma2=rep(1,T)
  tau2=rep(1,T)
  rho2=rep(1,T)
  beta=matrix(0,T,d+2)
  
  for(t in 1:(T-1)){
    #print(t)
    
    Sigma=solve(G/sigma2[t]+K/tau2[t])
    
    mu=Sigma%*%t(C)%*%y/sigma2[t]
    
    beta[t+1,]=rmvnorm(n=1,mu,Sigma)
    
    tau2[t+1]=rinvgamma(n=1,a+rk/2,b+t(beta[t+1,])%*%K%*%beta[t+1,]/2)
    
    RSS=sum((y-C%*%beta[t+1,])^2)
    
    sigma2[t+1]=rinvgamma(n=1,n/2,RSS/2)
    
  }
  
  beta=beta[-c(1:100),] ## remove burn-in
  betahat=colMeans(beta)
  
  fhatlin=X%*%betahat[c(1,2)] ##### linear 
  fhatnonlin=Z%*%betahat[-c(1,2)] ## nonlinear 
  fhat=C%*%betahat ## overall
  
  storefhatEil[,r]=fhat
  storefhatlinEil[,r]=fhatlin
  storefhatnonlinEil[,r]=fhatnonlin
  
  ########################### repeat using truncated power series for nonlinear effect
  
  X=cbind(rep(1,n),scale(x)) ## linear effect
  
  d=20
  #Z=DemmlerReinsch(x,d=d)$X ## nonlinear effect
  #Z=MixedModel(x,d=d)$X ## nonlinear
  #Z=Eilers(x,d=d)$X ## nonlinear
  Z=truncated_power_series(x,d=d,center=TRUE)
  
  C=cbind(X,Z) ## complete design matrix 
  
  ## estimate f and sigma2 using a Gibbs sampler #################################
  
  a=0.001 ## IG hyperparameters
  b=0.001
  
  G=t(C)%*%C
  #K=diag(c(0,0,0,0,rep(1,d-2)))
  K=diag(c(0,0,rep(1,d)))
  
  rk=sum(K)
  
  sigma2=rep(1,T)
  tau2=rep(1,T)
  rho2=rep(1,T)
  beta=matrix(0,T,d+2)
  
  for(t in 1:(T-1)){
    #print(t)
    
    Sigma=solve(G/sigma2[t]+K/tau2[t])
    
    mu=Sigma%*%t(C)%*%y/sigma2[t]
    
    beta[t+1,]=rmvnorm(n=1,mu,Sigma)
    
    tau2[t+1]=rinvgamma(n=1,a+rk/2,b+t(beta[t+1,])%*%K%*%beta[t+1,]/2)
    
    RSS=sum((y-C%*%beta[t+1,])^2)
    
    sigma2[t+1]=rinvgamma(n=1,n/2,RSS/2)
    
  }
  
  beta=beta[-c(1:100),] ## remove burn-in
  betahat=colMeans(beta)
  
  fhatlin=X%*%betahat[c(1,2)] ##### linear 
  fhatnonlin=Z%*%betahat[-c(1,2)] ## nonlinear 
  fhat=C%*%betahat ## overall
  
  storefhatTrunc[,r]=fhat
  storefhatlinTrunc[,r]=fhatlin
  storefhatnonlinTrunc[,r]=fhatnonlin
  
}


###################################################### evaluation ##################

t=seq(0,1,0.01)
ft=sin(2*pi*t)
#ft=t^2
#ft=sin(2*pi*t)+6/pi*(t-1/2)
#ft=t
#flint=t
flint=(-6/pi)*(t-1/2)
#flint=t-1/6
#flint=0*t
#flint=-90/pi^3*(t-1/2)
fnonlint=ft-flint

par(mfrow=c(5,4))
par(mar=c(3.5, 4, 3, 1) + 0.2)
par(mar=c(3, 3.5, 0, 0.5) + 0.2)
col="darkgray"
col2="black"
par(cex.lab=1.5)
par(cex.axis=1.5)
par(cex.main=2.1)
cex=1.4
lwd=3

ind=round(seq(1,n,length=100))
#ind=seq(1,5000,10) ## this reduces the size of the pdf-file for n=5,000

z=1.4

layout(matrix(1:24,6,4,byrow=TRUE),widths=c(1.5,3,3,3),heights=c(1,2,2,2,2,2))
layout(matrix(1:20,5,4,byrow=TRUE),widths=c(2,3,3,3),heights=c(1,2,2,2,2))


plot(x=0,y=0,col="white",ylim=c(-z,z),xlim=c(0,1),xaxt="n",xlab="",ylab="",yaxt="n",bty="n")
#rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "lightgray")

plot(x=0,y=0,col="white",ylim=c(-z,z),xlim=c(0,1),xaxt="n",xlab="",ylab="",yaxt="n",bty="n")
#rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "lightgray")
title(main="Linear effect",line=-5)
plot(x=0,y=0,col="white",ylim=c(-z,z),xlim=c(0,1),xaxt="n",xlab="",ylab="",yaxt="n",bty="n")
#rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "lightgray")
title(main="Nonlinear effect",line=-5)
plot(x=0,y=0,col="white",ylim=c(-z,z),xlim=c(0,1),xaxt="n",xlab="",ylab="",yaxt="n",bty="n")
#rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "lightgray")
title(main="Overall effect",line=-5)


############## DR basis ###########################################

plot(x=0,y=0,col="white",ylim=c(-z,z),xlim=c(0,1),xaxt="n",xlab="",ylab="",yaxt="n",bty="n")
#rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "lightgray")
title(main="DR basis",line=-5)

plot(storeX[ind,1],storefhatlin[ind,1],type="l",ylim=c(-z,z),main="",col=col,xlab="",ylab="")
for(r in 2:R){
  lines(storeX[ind,r],storefhatlin[ind,r],col=col)
}
lines(t,flint,col=col2,lwd=lwd)

title(ylab="y", line=2.3)
title(xlab="x", line=1.5)
legend(x="topright",legend=c("true effect","estimates"),pch=c(16,16),col=c(col2,col),bty="n",cex=cex)

plot(storeX[ind,1],storefhatnonlin[ind,1],type="l",ylim=c(-z,z),main="",col=col,xlab="",ylab="")
for(r in 2:R){
  lines(storeX[ind,r],storefhatnonlin[ind,r],col=col)
}
lines(t,fnonlint,col=col2,lwd=lwd)

title(ylab="y", line=2.3)
title(xlab="x", line=1.5)
#legend(x=0.49,y=1.05,legend=c("true effect","estimates"),pch=c(16,16),col=c(col2,col),bty="n",cex=cex)


plot(storeX[ind,1],storefhat[ind,1],type="l",ylim=c(-z,z),main="",col=col,xlab="",ylab="")
for(r in 2:R){
  lines(storeX[ind,r],storefhat[ind,r],col=col)
}
lines(t,ft,col=col2,lwd=lwd)

title(ylab="y", line=2.3)
title(xlab="x", line=1.5)
#legend(x=0.49,y=1.05,legend=c("true effect","estimates"),pch=c(16,16),col=c(col2,col),bty="n",cex=cex)

########## MMR #############################################################

plot(x=0,y=0,col="white",ylim=c(-z,z),xlim=c(0,1),xaxt="n",xlab="",ylab="",yaxt="n",bty="n")
#rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "lightgray")
title(main="MMR",line=-5)

plot(storeX[ind,1],storefhatlinMMR[ind,1],type="l",ylim=c(-z,z),main="",col=col,xlab="",ylab="")
for(r in 2:R){
  lines(storeX[ind,r],storefhatlinMMR[ind,r],col=col)
}
lines(t,flint,col=col2,lwd=lwd)

title(ylab="y", line=2.3)
title(xlab="x", line=1.5)
#legend(x=0.49,y=1.05,legend=c("true effect","estimates"),pch=c(16,16),col=c(col2,col),bty="n",cex=cex)

plot(storeX[ind,1],storefhatnonlinMMR[ind,1],type="l",ylim=c(-z,z),main="",col=col,xlab="",ylab="")
for(r in 2:R){
  lines(storeX[ind,r],storefhatnonlinMMR[ind,r],col=col)
}
lines(t,fnonlint,col=col2,lwd=lwd)

title(ylab="y", line=2.3)
title(xlab="x", line=1.5)
#legend(x=0.49,y=1.05,legend=c("true effect","estimates"),pch=c(16,16),col=c(col2,col),bty="n",cex=cex)


plot(storeX[ind,1],storefhatMMR[ind,1],type="l",ylim=c(-z,z),main="",col=col,xlab="",ylab="")
for(r in 2:R){
  lines(storeX[ind,r],storefhatMMR[ind,r],col=col)
}
lines(t,ft,col=col2,lwd=lwd)

title(ylab="y", line=2.3)
title(xlab="x", line=1.5)
#legend(x=0.49,y=1.05,legend=c("true effect","estimates"),pch=c(16,16),col=c(col2,col),bty="n",cex=cex)

############# Eilers #############################################

plot(x=0,y=0,col="white",ylim=c(-z,z),xlim=c(0,1),xaxt="n",xlab="",ylab="",yaxt="n",bty="n")
title(main="Eilers",line=-5)

plot(storeX[ind,1],storefhatlinEil[ind,1],type="l",ylim=c(-z,z),main="",col=col,xlab="",ylab="")
for(r in 2:R){
  lines(storeX[ind,r],storefhatlinEil[ind,r],col=col)
}
lines(t,flint,col=col2,lwd=lwd)

title(ylab="y", line=2.3)
title(xlab="x", line=1.5)
#legend(x=0.49,y=1.05,legend=c("true effect","estimates"),pch=c(16,16),col=c(col2,col),bty="n",cex=cex)

plot(storeX[ind,1],storefhatnonlinEil[ind,1],type="l",ylim=c(-z,z),main="",col=col,xlab="",ylab="")
for(r in 2:R){
  lines(storeX[ind,r],storefhatnonlinEil[ind,r],col=col)
}
lines(t,fnonlint,col=col2,lwd=lwd)

title(ylab="y", line=2.3)
title(xlab="x", line=1.5)
#legend(x=0.49,y=1.05,legend=c("true effect","estimates"),pch=c(16,16),col=c(col2,col),bty="n",cex=cex)


plot(storeX[ind,1],storefhatEil[ind,1],type="l",ylim=c(-z,z),main="",col=col,xlab="",ylab="")
for(r in 2:R){
  lines(storeX[ind,r],storefhatEil[ind,r],col=col)
}
lines(t,ft,col=col2,lwd=lwd)

title(ylab="y", line=2.3)
title(xlab="x", line=1.5)
#legend(x=0.49,y=1.05,legend=c("true effect","estimates"),pch=c(16,16),col=c(col2,col),bty="n",cex=cex)
# 
# ############# Trunc #############################################
# 
# 
# plot(x=0,y=0,col="white",ylim=c(-z,z),xlim=c(0,1),xaxt="n",xlab="",ylab="",yaxt="n",bty="n")
# title(main="cTPS",line=-5)
# 
# a=5
# plot(storeX[ind,1],storefhatlinTrunc[ind,1],type="l",ylim=c(-a*z,a*z),main="",col=col,xlab="",ylab="")
# for(r in 2:R){
#   lines(storeX[ind,r],storefhatlinTrunc[ind,r],col=col)
# }
# lines(t,flint,col=col2,lwd=lwd)
# 
# title(ylab="y", line=2.3)
# title(xlab="x", line=1.5)
# #legend(x=0.49,y=-0.5,legend=c("true effect","estimates"),pch=c(16,16),col=c(col2,col),bty="n",cex=cex)
# 
# plot(storeX[ind,1],storefhatnonlinTrunc[ind,1],type="l",ylim=c(-a*z,a*z),main="",col=col,xlab="",ylab="")
# for(r in 2:R){
#   lines(storeX[ind,r],storefhatnonlinTrunc[ind,r],col=col)
# }
# lines(t,fnonlint,col=col2,lwd=lwd)
# 
# title(ylab="y", line=2.3)
# title(xlab="x", line=1.5)
# #legend(x=0.49,y=5.55,legend=c("true effect","estimates"),pch=c(16,16),col=c(col2,col),bty="n",cex=cex)
# 
# plot(storeX[ind,1],storefhatTrunc[ind,1],type="l",ylim=c(-z,z),main="",col=col,xlab="",ylab="")
# for(r in 2:R){
#   lines(storeX[ind,r],storefhatTrunc[ind,r],col=col)
# }
# lines(t,ft,col=col2,lwd=lwd)
# 
# title(ylab="y", line=2.3)
# title(xlab="x", line=1.5)
#legend(x=0.49,y=1.05,legend=c("true effect","estimates"),pch=c(16,16),col=c(col2,col),bty="n",cex=cex)

############# B-splines #############################################

# plot(x=0,y=0,col="white",ylim=c(-z,z),xlim=c(0,1),xaxt="n",xlab="",ylab="",yaxt="n",bty="n")
# title(main="Bsplines",line=-5)
# 
# a=1
# plot(storeX[ind,1],storefhatlinB[ind,1],type="l",ylim=c(-a*z,a*z),main="",col=col,xlab="",ylab="")
# for(r in 2:R){
#   lines(storeX[ind,r],storefhatlinB[ind,r],col=col)
# }
# lines(t,flint,col=col2,lwd=lwd)
# 
# title(ylab="y", line=2.3)
# title(xlab="x", line=2)
# #legend(x=0.49,y=1.05,legend=c("true effect","estimates"),pch=c(16,16),col=c(col2,col),bty="n",cex=cex)
# 
# plot(storeX[ind,1],storefhatnonlinB[ind,1],type="l",ylim=c(-a*z,a*z),main="",col=col,xlab="",ylab="")
# for(r in 2:R){
#   lines(storeX[ind,r],storefhatnonlinB[ind,r],col=col)
# }
# lines(t,fnonlint,col=col2,lwd=lwd)
# 
# title(ylab="y", line=2.3)
# title(xlab="x", line=2)
# #legend(x=0.49,y=1.05,legend=c("true effect","estimates"),pch=c(16,16),col=c(col2,col),bty="n",cex=cex)
# 
# 
# plot(storeX[ind,1],storefhatB[ind,1],type="l",ylim=c(-z,z),main="",col=col,xlab="",ylab="")
# for(r in 2:R){
#   lines(storeX[ind,r],storefhatB[ind,r],col=col)
# }
# lines(t,ft,col=col2,lwd=lwd)
# 
# title(ylab="y", line=2.3)
# title(xlab="x", line=2)
#legend(x=0.49,y=1.05,legend=c("true effect","estimates"),pch=c(16,16),col=c(col2,col),bty="n",cex=cex)


### RMSE #############################################################

plot(x=0,y=0,col="white",ylim=c(-z,z),xlim=c(0,1),xaxt="n",xlab="",ylab="",yaxt="n",bty="n")
title(main="RMSE",line=-5)

boxplot(main="",(sqrt(colMeans((storeflin-storefhatlin)^2))),(sqrt(colMeans((storeflin-storefhatlinMMR)^2))),(sqrt(colMeans((storeflin-storefhatlinEil)^2))),names=c("DR basis","MMR","Eilers"),ylim=c(0,0.4))
title(ylab="", line=2.3)
#mtext("logRMSE", line=2.3,las=1,side=2)
#title(xlab="basis", line=2)

boxplot(main="",(sqrt(colMeans((storefnonlin-storefhatnonlin)^2))),(sqrt(colMeans((storefnonlin-storefhatnonlinMMR)^2))),(sqrt(colMeans((storefnonlin-storefhatnonlinEil)^2))),names=c("DR basis","MMR","Eilers"),ylim=c(0,0.4))
title(ylab="", line=2.3)
#title(xlab="basis", line=2)

boxplot(main="",(sqrt(colMeans((storef-storefhat)^2))),(sqrt(colMeans((storef-storefhatMMR)^2))),(sqrt(colMeans((storef-storefhatEil)^2))),names=c("DR basis","MMR","Eilers"),ylim=c(0,0.4))
title(ylab="", line=2.3)
#title(xlab="basis", line=2)

 
