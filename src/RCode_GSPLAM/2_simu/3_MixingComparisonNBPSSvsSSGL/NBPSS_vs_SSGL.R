
################################################################################
################################################################################
#################### NBPSS vs. SSGL ############################################
################################################################################
################################################################################

## In this script we compare the performance of our method when
## using the NBPSS prior for the group coefficients and when using
## the SSGL prior for the group coefficients.
## We are particularly interested in the MCMC mixing of the binary inclusion 
## indicators of the nonlinear effects.
## The results of the mixing comparison are contained in the supplement.


myfunc<-function(x,nr=1){
  
  if(nr==1){f=x-1/2}
  if(nr==2){f=-2*x+1}
  if(nr==3){f=4*(x-1/2)^2-1/3}
  if(nr==4){f=sin(2*pi*x)+6/pi*(x-1/2)}
  if(nr==5){f=2*exp(-3*x)-2*(1/3-1/(3*exp(3)))}
  if(nr==6){f=sin(2*pi*x)}
  if(nr>6){f=rep(0,length(x))}
  
  return(f)
  
}

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

### 1. local weights wj in (0,1) with uniform prior

#### Generate data 

set.seed(123) ## here we set the seed

n=10^3 ## sample size
p=8 ## number of continuous covariates
J=19 ## overall number of effects

Xcont=data.frame(matrix(runif(n*p),n,p)) ## continuous covariates

Fc=matrix(0,n,p)
for(j in 1:8){
  Fc[,j]=myfunc(Xcont[,j],nr=j)
}

Xdummy=data.frame(cbind(rbinom(n,1,0.5),rbinom(n,1,0.5))) ## dummies
colnames(Xdummy)=c("d1","d2")

Fd=cbind(0*Xdummy[,1],0.5*Xdummy[,2])

lon=runif(n) ## spatial effect
lat=runif(n)
Xspat=data.frame(lon,lat)

Fs=rep(0,n)
#Fs=sin(2*pi*lon)*cos(2*pi*lat)

eta=rowSums(cbind(Fd,Fc,Fs))  ################# true predictor

k=3/2
y=rweibull(n,shape=k,scale=exp(-eta/k)) ########## time points


delta=rep(1,n) ############################ censoring 
cens=rep(2,n)

delta[y>=cens]=0
y=pmin(y,cens)

data=data.frame(y,delta,eta,Xdummy,Xcont,Xspat) ########### sort time points increasingly
data=data[order(y),] 
y=data$y;delta=data$delta;eta=data$eta;
Xdummy=cbind(data$d1,data$d2);
Xcont=data[,5+1:p];
Xspat=cbind(data$lon,data$lat)


R=20 ### number of chains
J=19 ## number of effects

phatNBPSS=matrix(0,R,J)
phatSSGL=matrix(0,R,J)

rhatNBPSS=matrix(0,R,J)
rhatSSGL=matrix(0,R,J)

essNBPSS=matrix(0,R,J)
essSSGL=matrix(0,R,J)

gammaStore=matrix(0,R,J)

StoreNBPSS=list()
StoreSSGL=list()

T=10^4
t0=Sys.time()
for(r in 1:R){

print(paste("Replicate",r,"######################################"))
  
gamma.ini=rbinom(n=J,1,0.5) ### this is the initial model 
gammaStore[r,]=gamma.ini

## apply NBPSS
fit=nbpsscox(y=y,delta=delta,Xdummy=as.matrix(Xdummy),Xcont=as.matrix(Xcont),Xspat=as.matrix(Xspat),nonlin="DR",d1=9,ngrid=1000,prior="SBP",T=T,w1=1,w2=1,gamma.ini=gamma.ini,burn=10^3)
phatNBPSS[r,]=colMeans(fit$samples$probi)
rhatNBPSS[r,]=summarize_draws(as_draws(fit$samples$probi))$rhat
essNBPSS[r,]=summarize_draws(as_draws(fit$samples$probi))$ess_bulk

StoreNBPSS[[r]]=fit$samples$probi

## apply SSGL
fit=nbpsscox(y=y,delta=delta,Xdummy=as.matrix(Xdummy),Xcont=as.matrix(Xcont),Xspat=as.matrix(Xspat),nonlin="DR",d1=9,ngrid=1000,prior="GA",T=T,w1=1,w2=1,gamma.ini=gamma.ini,burn=10^3)
phatSSGL[r,]=colMeans(fit$samples$probi)
rhatSSGL[r,]=summarize_draws(as_draws(fit$samples$probi))$rhat
essSSGL[r,]=summarize_draws(as_draws(fit$samples$probi))$ess_bulk

StoreSSGL[[r]]=fit$samples$probi

}

t1=Sys.time()
t1-t0


## MCMC approximate posterior inclusion probabilities across the R chains

par(mfrow=c(1,2))
boxplot(phatNBPSS,ylim=c(0,1),main="Posterior inclusion probabilities (NBPSS)")
boxplot(phatSSGL,ylim=c(0,1),main="Posterior inclusion probabilities (SSGL)")

## Rhat across the R chains

par(mfrow=c(1,2))
boxplot(rhatNBPSS,ylim=c(0.9,2),main="Rhat (NBPSS)")
boxplot(rhatSSGL,ylim=c(0.9,2),main="Rhat (SSGL)")

## trace plots for one of the chains (here we select chain number 3)

par(mfrow=c(4,5))
plot.new()
text(0.5,0.5,"NBPSS",cex=2,font=2)

for(i in 1:19){
  plot(StoreNBPSS[[3]][,i],type="l",ylim=c(0,1),main=paste("Effect",i))  
}

par(mfrow=c(4,5))
plot.new()
text(0.5,0.5,"SSGL",cex=2,font=2)

for(i in 1:19){
plot(StoreSSGL[[3]][,i],type="l",ylim=c(0,1),main=paste("Effect",i))  
}

### Store all results

# Store=list()
# Store[[1]]=StoreNBPSS
# Store[[2]]=StoreSSGL
# names(Store)=c("NBPSS","SSGL")
#  
# save(Store,file="C:\\Users\\Nutzer\\Documents\\R Programme\\NBPSSCOX\\2_simu\\3_MixingComparisonNBPSSvsSSGL\\Results\\Localwj\\Store.RData")
# ## load results
# load("C:\\Users\\Nutzer\\Documents\\R Programme\\NBPSSCOX\\2_simu\\3_MixingComparisonNBPSSvsSSGL\\Results\\Localwj\\Store.RData")
# 
# StoreNBPSS=Store$NBPSS
# StoreSSGL=Store$SSGL
# for(r in 1:R){
#   temp=StoreNBPSS[[r]]
#   phatNBPSS[r,]=colMeans(temp)
#   rhatNBPSS[r,]=summarize_draws(as_draws(temp))$rhat
#   essNBPSS[r,]=summarize_draws(as_draws(temp))$ess_bulk
# 
#   temp=StoreSSGL[[r]]
#   phatSSGL[r,]=colMeans(temp)
#   rhatSSGL[r,]=summarize_draws(as_draws(temp))$rhat
#   essSSGL[r,]=summarize_draws(as_draws(temp))$ess_bulk
# }
# StoreNBPSSloc=StoreNBPSS
# StoreSSGLloc=StoreSSGL
# phatNBPSSloc=phatNBPSS
# phatSSGLloc=phatSSGL
# rhatNBPSSloc=rhatNBPSS
# rhatSSGLloc=rhatSSGL

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

### 2. global w in (0,1) with Beta(1,J) prior

#### Generate data 

set.seed(123) ## here we set the seed

n=10^3 ## sample size
p=8 ## number of continuous covariates
J=19 ## overall number of effects

Xcont=data.frame(matrix(runif(n*p),n,p)) ## continuous covariates

Fc=matrix(0,n,p)
for(j in 1:8){
  Fc[,j]=myfunc(Xcont[,j],nr=j)
}

Xdummy=data.frame(cbind(rbinom(n,1,0.5),rbinom(n,1,0.5))) ## dummies
colnames(Xdummy)=c("d1","d2")

Fd=cbind(0*Xdummy[,1],0.5*Xdummy[,2])

lon=runif(n) ## spatial effect
lat=runif(n)
Xspat=data.frame(lon,lat)

Fs=rep(0,n)
#Fs=sin(2*pi*lon)*cos(2*pi*lat)

eta=rowSums(cbind(Fd,Fc,Fs))  ################# true predictor

k=3/2
y=rweibull(n,shape=k,scale=exp(-eta/k)) ########## time points


delta=rep(1,n) ############################ censoring 
cens=rep(2,n)

delta[y>=cens]=0
y=pmin(y,cens)

data=data.frame(y,delta,eta,Xdummy,Xcont,Xspat) ########### sort time points increasingly
data=data[order(y),] 
y=data$y;delta=data$delta;eta=data$eta;
Xdummy=cbind(data$d1,data$d2);
Xcont=data[,5+1:p];
Xspat=cbind(data$lon,data$lat)

R=20 ### number of chains

phatNBPSS=matrix(0,R,J)
phatSSGL=matrix(0,R,J)

rhatNBPSS=matrix(0,R,J)
rhatSSGL=matrix(0,R,J)

essNBPSS=matrix(0,R,J)
essSSGL=matrix(0,R,J)

gammaStore=matrix(0,R,J)

StoreNBPSS=list()
StoreSSGL=list()

T=10^4
burn=10^3
w1=1; w2=J

t0=Sys.time()
for(r in 1:R){
  
  print(paste("Replicate",r,"######################################"))
  
  gamma.ini=rbinom(n=J,1,0.5) ### this is the initial model 
  gammaStore[r,]=gamma.ini
  
  ## apply NBPSS
  fit=nbpsscox(y=y,delta=delta,Xdummy=as.matrix(Xdummy),Xcont=as.matrix(Xcont),Xspat=as.matrix(Xspat),nonlin="DR",d1=9,ngrid=1000,prior="SBP",T=T,w1=w1,w2=w2,gamma.ini=gamma.ini,burn=burn,globalw=TRUE)
  phatNBPSS[r,]=colMeans(fit$samples$probi)
  rhatNBPSS[r,]=summarize_draws(as_draws(fit$samples$probi))$rhat
  essNBPSS[r,]=summarize_draws(as_draws(fit$samples$probi))$ess_bulk
  
  StoreNBPSS[[r]]=fit$samples$probi
  
  ## apply SSGL
  fit=nbpsscox(y=y,delta=delta,Xdummy=as.matrix(Xdummy),Xcont=as.matrix(Xcont),Xspat=as.matrix(Xspat),nonlin="DR",d1=9,ngrid=1000,prior="GA",T=T,w1=w1,w2=w2,gamma.ini=gamma.ini,burn=burn,globalw=TRUE)
  phatSSGL[r,]=colMeans(fit$samples$probi)
  rhatSSGL[r,]=summarize_draws(as_draws(fit$samples$probi))$rhat
  essSSGL[r,]=summarize_draws(as_draws(fit$samples$probi))$ess_bulk
  
  StoreSSGL[[r]]=fit$samples$probi
  
}

t1=Sys.time()
t1-t0


## MCMC approximate posterior inclusion probabilities across the R chains

par(mfrow=c(1,2))
boxplot(phatNBPSS,ylim=c(0,1),main="Posterior inclusion probabilities (NBPSS)")
boxplot(phatSSGL,ylim=c(0,1),main="Posterior inclusion probabilities (SSGL)")

## Rhat across the R chains

par(mfrow=c(1,2))
boxplot(rhatNBPSS,ylim=c(0.9,2),main="Rhat (NBPSS)",names=names,las=1)
boxplot(rhatSSGL,ylim=c(0.9,2),main="Rhat (SSGL)",names=names,las=1)

## trace plots for one of the chains (here we select chain number 13)

par(mfrow=c(4,5))
plot.new()
text(0.5,0.5,"NBPSS",cex=2,font=2)

j=13
for(i in 1:19){
  plot(StoreNBPSS[[j]][,i],type="l",ylim=c(0,1),main=paste("Effect",i))  
}

par(mfrow=c(4,5))
plot.new()
text(0.5,0.5,"SSGL",cex=2,font=2)

for(i in 1:19){
  plot(StoreSSGL[[j]][,i],type="l",ylim=c(0,1),main=paste("Effect",i))  
}

### Store all results

# Store=list()
# Store[[1]]=StoreNBPSS
# Store[[2]]=StoreSSGL
# names(Store)=c("NBPSS","SSGL")
#  
# save(Store,file="C:\\Users\\Nutzer\\Documents\\R Programme\\NBPSSCOX\\2_simu\\3_MixingComparisonNBPSSvsSSGL\\Results\\Globalw\\Store.RData")

# ## load results
# load("C:\\Users\\Nutzer\\Documents\\R Programme\\NBPSSCOX\\2_simu\\3_MixingComparisonNBPSSvsSSGL\\Results\\Globalw\\Store.RData")
# 
# StoreNBPSS=Store$NBPSS
# StoreSSGL=Store$SSGL
# for(r in 1:R){
#   temp=StoreNBPSS[[r]]
#   phatNBPSS[r,]=colMeans(temp)
#   rhatNBPSS[r,]=summarize_draws(as_draws(temp))$rhat
#   essNBPSS[r,]=summarize_draws(as_draws(temp))$ess_bulk
# 
#   temp=StoreSSGL[[r]]
#   phatSSGL[r,]=colMeans(temp)
#   rhatSSGL[r,]=summarize_draws(as_draws(temp))$rhat
#   essSSGL[r,]=summarize_draws(as_draws(temp))$ess_bulk
# }
# 
# StoreNBPSSglob=StoreNBPSS
# StoreSSGLglob=StoreSSGL
# phatNBPSSglob=phatNBPSS
# phatSSGLglob=phatSSGL
# rhatNBPSSglob=rhatNBPSS
# rhatSSGLglob=rhatSSGL

#################################################### joint evaluation ##########

## here we create the Figures shown in the supplement

#### MCMC approximate posterior inclusion probabilities

names=c(expression(D[1]),
        expression(D[2]),
        expression(L[1]),
        expression(L[2]),
        expression(L[3]),
        expression(L[4]),
        expression(L[5]),
        expression(L[6]),
        expression(L[7]),
        expression(L[8]),
        expression(N[1]),
        expression(N[2]),
        expression(N[3]),
        expression(N[4]),
        expression(N[5]),
        expression(N[6]),
        expression(N[7]),
        expression(N[8]),
        "S")

layout(matrix(c(1,1,2,3,4,5),ncol=2,byrow=TRUE),heights=c(0.5,3,3))
par(mar=c(0,0,0,0))
plot.new()
text(0.5,0.3,"MCMC approximate posterior effect inclusion probabilities",cex=2.5,font=2)
par(mar=c(3,4,4,1))
par(cex.main=2)
par(cex.axis=1.5)
par(cex.lab=1.5)
ylab="Posterior effect inclusion probabilities"
boxplot(phatNBPSSloc,ylim=c(0,1),main="NBPSS (local wj)",names=names)
title(ylab=ylab, line=2.3)
boxplot(phatSSGLloc,ylim=c(0,1),main="SSGL (local wj)",names=names)
title(ylab=ylab, line=2.3)
boxplot(phatNBPSSglob,ylim=c(0,1),main="NBPSS (global w)",names=names)
title(ylab=ylab, line=2.3)
boxplot(phatSSGLglob,ylim=c(0,1),main="SSGL (global w)",names=names)
title(ylab=ylab, line=2.3)
#dev.off()

#### rhat

layout(matrix(c(1,1,2,3,4,5),ncol=2,byrow=TRUE),heights=c(0.5,3,3))
par(mar=c(0,0,0,0))
plot.new()
text(0.5,0.3,"Rhat for the conditional posterior effect inclusion probabilities",cex=2.5,font=2)
par(mar=c(3,4,4,1))
par(cex.main=2)
par(cex.axis=1.5)
par(cex.lab=1.5)
ylab="Rhat"
boxplot(rhatNBPSSloc,ylim=c(0.9,2),main="NBPSS (local wj)",names=names)
title(ylab=ylab, line=2.3)
boxplot(rhatSSGLloc,ylim=c(0.9,2),main="SSGL (local wj)",names=names)
title(ylab=ylab, line=2.3)
boxplot(rhatNBPSSglob,ylim=c(0.9,2),main="NBPSS (global w)",names=names)
title(ylab=ylab, line=2.3)
boxplot(rhatSSGLglob,ylim=c(0.9,2),main="SSGL (global w)",names=names)
title(ylab=ylab, line=2.3)
#dev.off()

## trace plots for one of the chains (here we select chain number 13)

par(mfrow=c(4,5))
plot.new()
text(0.5,0.5,"NBPSS",cex=2,font=2)
ylab="cond. post. inclu."

par(mar=c(2,2,2,1))
j=13
for(i in 1:19){
  plot(StoreNBPSSloc[[j]][,i],type="l",ylim=c(0,1.1),main=names[i],ylab="") 
}

#par(mfrow=c(4,5))
plot.new()
text(0.5,0.5,"SSGL",cex=2,font=2)

for(i in 1:19){
  plot(StoreSSGLloc[[j]][,i],type="l",ylim=c(0,1.1),main=names[i],ylab="")  
}

## trace plots for one of the chains (here we select chain number 13)

par(mfrow=c(4,5))
plot.new()
text(0.5,0.5,"NBPSS",cex=2,font=2)

j=13
for(i in 1:19){
  plot(StoreNBPSSglob[[j]][,i],type="l",ylim=c(0,1.1),main=names[i],ylab="")  
}

par(mfrow=c(4,5))
plot.new()
text(0.5,0.5,"SSGL",cex=2,font=2)

for(i in 1:19){
  plot(StoreSSGLglob[[j]][,i],type="l",ylim=c(0,1.1),main=names[i],ylab="")  
}

