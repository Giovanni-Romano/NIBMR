
################################################################################
################################################################################
########## Simulation study ####################################################
################################################################################
################################################################################

########### In this script we conduct the main simulation study
########### from the manuscript where we compare the performance of NBPSS and mboost
########### in the geoadditive Cox regression model.

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
################################################################################
################################################################################

### This is the code for the simulations in the manuscript

### Repeat a lot of times and store results

set.seed(123)

n=1000 ## sample size
p=8 ## number of continuous covariates
J=19 ## overall number of effects

dim1=9 ## number of basis functions for the nonlinear effects

df=4 ## degrees of freedom for boosting
mstop=10^3 ## number of stopping it for boosting (we use bootstrapping to optimize afterwards)

R=30 ## number of replicate data sets

rho=0.5
Sigma=matrix(0,8,8) ## for correlated covariates
for(j in 1:8){
  for(k in 1:8){
    Sigma[j,k]=rho^{abs(j-k)}
  }
}

etaStore=matrix(0,R,n) ## store true predictor
deltaStore=rep(0,R) ## store event rate (to compute average censoring rate)
  
etahatNBPSS1=matrix(0,R,n) ########### store estimated predictors
etahatNBPSS2=matrix(0,R,n) 
etahatNBPSS3=matrix(0,R,n) 

etahatmboost1=matrix(0,R,n)  
etahatmboost2=matrix(0,R,n) 
etahatmboost3=matrix(0,R,n) 

phatNBPSS1=matrix(0,R,J) ### store estimated inclusion probabilities / model
phatNBPSS2=matrix(0,R,J)
phatNBPSS3=matrix(0,R,J)

phatmboost1=matrix(0,R,J)
phatmboost2=matrix(0,R,J)
phatmboost3=matrix(0,R,J)

mstop1=rep(0,R)
mstop2=rep(0,R)
mstop3=rep(0,R)

t0=Sys.time()
for(r in 1:R){
  
  print(paste("######## Replicate",r))
  
  ############################ generate data
  
  Xcont=data.frame(matrix(runif(n*p),n,p)) ## continuous covariates
  
  #Xcont=pnorm(rmvnorm(n=n,sigma=Sigma)) with correlation
  
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
  
  #Fs=rep(0,n)
  Fs=sin(2*pi*lon)*cos(2*pi*lat) ## true spatial effect
  
  eta=rowSums(cbind(Fd,Fc,Fs))  ################# true predictor
  
  k=3/2
  y=rweibull(n,shape=k,scale=exp(-eta/k)) ########## time points
  
  
  delta=rep(1,n) ############################ censoring 
  #cens=rep(quantile(y,probs=0.85),n)
  #cens=rexp(n,rate=1/5) 
  cens=rep(2,n)
  
  delta[y>=cens]=0
  y=pmin(y,cens)
  
  data=data.frame(y,delta,eta,Xdummy,Xcont,Xspat) ########### sort time points increasingly
  data=data[order(y),] 
  y=data$y;delta=data$delta;eta=data$eta;
  Xdummy=cbind(data$d1,data$d2);
  Xcont=data[,5+1:p];
  Xspat=cbind(data$lon,data$lat)
  
  ##############################################################################
  
  deltaStore[r]=mean(delta) ## store event rate
  etaStore[r,]=eta ### store true predictor
  
  ##############################################################################
  ##############################################################################
  ##############################################################################
  ######################################################################### apply NBPSS with Eilers
  
  fit=nbpsscox(y=y,delta=delta,Xdummy=as.matrix(Xdummy),Xcont=as.matrix(Xcont),Xspat=as.matrix(Xspat),nonlin="Eilers",d1=dim1,ngrid=1000,prior="SBP",T=10^3,w1=1,w2=1)
  
  phat=colMeans(fit$samples$probi)
  etahat=fit$design$X%*%colMeans(fit$samples$beta)
  
  etahatNBPSS1[r,]=etahat
  phatNBPSS1[r,]=phat
  
  ######################################################################### NBPSS with MMR
  
  fit=nbpsscox(y=y,delta=delta,Xdummy=as.matrix(Xdummy),Xcont=as.matrix(Xcont),Xspat=as.matrix(Xspat),nonlin="MM",d1=dim1,ngrid=1000,prior="SBP",T=10^3,w1=1,w2=1)
  
  phat=colMeans(fit$samples$probi)
  etahat=fit$design$X%*%colMeans(fit$samples$beta)
  
  etahatNBPSS2[r,]=etahat
  phatNBPSS2[r,]=phat
  
  
  ######################################################################### NBPSS with DR
  
  fit=nbpsscox(y=y,delta=delta,Xdummy=as.matrix(Xdummy),Xcont=as.matrix(Xcont),Xspat=as.matrix(Xspat),nonlin="DR",d1=dim1,ngrid=1000,prior="SBP",T=10^3,w1=1,w2=1)
  
  phat=colMeans(fit$samples$probi)
  etahat=fit$design$X%*%colMeans(fit$samples$beta)
  
  etahatNBPSS3[r,]=etahat
  phatNBPSS3[r,]=phat
  
  #####################################################################################################
  #####################################################################################################
  #####################################################################################################
  
  ################################################################################## apply mboost
  
  ### prepare dataset #####################################################
  
  data$X1S=scale(data$X1)
  data$X2S=scale(data$X2)
  data$X3S=scale(data$X3)
  data$X4S=scale(data$X4)
  data$X5S=scale(data$X5)
  data$X6S=scale(data$X6)
  data$X7S=scale(data$X7)
  data$X8S=scale(data$X8)
  
  data$d1=as.factor(data$d1)
  data$d2=as.factor(data$d2)
  
  ### run mboost for Cox model ###################################################### mboost with Eilers (default)
  
  fit=mboost(Surv(y,delta)~
               bols(d1,intercept=FALSE)+ ## dummies
               bols(d2,intercept=FALSE)+
               bols(X1S,intercept=FALSE)+ ## linear effects
               bols(X2S,intercept=FALSE)+
               bols(X3S,intercept=FALSE)+
               bols(X4S,intercept=FALSE)+
               bols(X5S,intercept=FALSE)+
               bols(X6S,intercept=FALSE)+
               bols(X7S,intercept=FALSE)+
               bols(X8S,intercept=FALSE)+
               bbs(X1S,center=TRUE,df=df,knots=dim1-2)+ ## nonlinear effects
               bbs(X2S,center=TRUE,df=df,knots=dim1-2)+
               bbs(X3S,center=TRUE,df=df,knots=dim1-2)+
               bbs(X4S,center=TRUE,df=df,knots=dim1-2)+
               bbs(X5S,center=TRUE,df=df,knots=dim1-2)+
               bbs(X6S,center=TRUE,df=df,knots=dim1-2)+
               bbs(X7S,center=TRUE,df=df,knots=dim1-2)+
               bbs(X8S,center=TRUE,df=df,knots=dim1-2)+
               bspatial(x1,x2,knots=2,differences=1,center=TRUE),           
             family=CoxPH(),data=data,control=boost_control(mstop=mstop,trace=TRUE))
  
  #brad(lon,lat,knots=35,df=df), ## spatial effect
  
  cv=cvrisk(fit)
  mstop(cv)
  
  mstop1[r]=mstop(cv)
  
  fit=fit[mstop(cv)]
  
  #summary(fit)
  
  ## extract selected model ########################################################
  
  pos=sort(unique(selected(fit)))
  
  gamma.hat=rep(0,19)
  gamma.hat[pos]=1
  gamma.hat
  
  phatmboost1[r,]=gamma.hat
  
  ## extract predictor eta ##########################################################
  
  etahat=fit$fitted()
  
  etahatmboost1[r,]=etahat
  
  
  #################################################################################### mboost with MMR
  
  fit=mboost(Surv(y,delta)~
               bols(d1,intercept=FALSE)+
               bols(d2,intercept=FALSE)+
               bols(X1S,intercept=FALSE)+
               bols(X2S,intercept=FALSE)+
               bols(X3S,intercept=FALSE)+
               bols(X4S,intercept=FALSE)+
               bols(X5S,intercept=FALSE)+
               bols(X6S,intercept=FALSE)+
               bols(X7S,intercept=FALSE)+
               bols(X8S,intercept=FALSE)+
               buser(X=MixedModel(data$X1S,d=dim1)$X,K=diag(dim1),df=df)+
               buser(X=MixedModel(data$X2S,d=dim1)$X,K=diag(dim1),df=df)+
               buser(X=MixedModel(data$X3S,d=dim1)$X,K=diag(dim1),df=df)+
               buser(X=MixedModel(data$X4S,d=dim1)$X,K=diag(dim1),df=df)+
               buser(X=MixedModel(data$X5S,d=dim1)$X,K=diag(dim1),df=df)+
               buser(X=MixedModel(data$X6S,d=dim1)$X,K=diag(dim1),df=df)+
               buser(X=MixedModel(data$X7S,d=dim1)$X,K=diag(dim1),df=df)+
               buser(X=MixedModel(data$X8S,d=dim1)$X,K=diag(dim1),df=df)+
               bspatial(x1,x2,knots=2,differences=1,center=TRUE),
             family=CoxPH(),data=data,control=boost_control(mstop=mstop,trace=TRUE))
  
  # brad(lon,lat,knots=35,df=df),
  
  cv=cvrisk(fit)
  mstop(cv)
  
  mstop2[r]=mstop(cv)
  
  fit=fit[mstop(cv)]
  
  #summary(fit)
  
  ## extract selected model ########################################################
  
  pos=sort(unique(selected(fit)))
  
  gamma.hat=rep(0,19)
  gamma.hat[pos]=1
  gamma.hat
  
  phatmboost2[r,]=gamma.hat
  
  ## extract predictor eta ##########################################################
  
  etahat=fit$fitted()
  
  etahatmboost2[r,]=etahat
  
  #################################################################################### mboost with DRB
  
  fit=mboost(Surv(y,delta)~
               bols(d1,intercept=FALSE)+
               bols(d2,intercept=FALSE)+
               bols(X1S,intercept=FALSE)+
               bols(X2S,intercept=FALSE)+
               bols(X3S,intercept=FALSE)+
               bols(X4S,intercept=FALSE)+
               bols(X5S,intercept=FALSE)+
               bols(X6S,intercept=FALSE)+
               bols(X7S,intercept=FALSE)+
               bols(X8S,intercept=FALSE)+
               buser(X=DemmlerReinsch(data$X1S,d=dim1)$X,K=diag(dim1),df=df)+
               buser(X=DemmlerReinsch(data$X2S,d=dim1)$X,K=diag(dim1),df=df)+
               buser(X=DemmlerReinsch(data$X3S,d=dim1)$X,K=diag(dim1),df=df)+
               buser(X=DemmlerReinsch(data$X4S,d=dim1)$X,K=diag(dim1),df=df)+
               buser(X=DemmlerReinsch(data$X5S,d=dim1)$X,K=diag(dim1),df=df)+
               buser(X=DemmlerReinsch(data$X6S,d=dim1)$X,K=diag(dim1),df=df)+
               buser(X=DemmlerReinsch(data$X7S,d=dim1)$X,K=diag(dim1),df=df)+
               buser(X=DemmlerReinsch(data$X8S,d=dim1)$X,K=diag(dim1),df=df)+
               bspatial(x1,x2,knots=2,differences=1,center=TRUE),
             family=CoxPH(),data=data,control=boost_control(mstop=mstop,trace=TRUE))
  
  #brad(lon,lat,knots=35,df=df),
  
  cv=cvrisk(fit)
  mstop(cv)
  
  mstop3[r]=mstop(cv)
  
  fit=fit[mstop(cv)]
  
  #summary(fit)
  
  ## extract selected model ########################################################
  
  pos=sort(unique(selected(fit)))
  
  gamma.hat=rep(0,19)
  gamma.hat[pos]=1
  gamma.hat
  
  phatmboost3[r,]=gamma.hat
  
  ## extract predictor eta ##########################################################
  
  etahat=fit$fitted()
  
  etahatmboost3[r,]=etahat
  
  ####################################################################################
  ####################################################################################
  ####################################################################################
  
  ## Store all results
  
  Store=list()
  Store[[1]]=list(etahatNBPSS1,etahatNBPSS2,etahatNBPSS3,etahatmboost1,etahatmboost2,etahatmboost3)
  Store[[2]]=list(phatNBPSS1,phatNBPSS2,phatNBPSS3,phatmboost1,phatmboost2,phatmboost3)
  Store[[3]]=list(etaStore,deltaStore)
  save(Store,file="Store.RData")
  
  storeMstop=list(mstop1,mstop2,mstop3)
  save(storeMstop,file="mstop.RData")
  
}
t1=Sys.time()
t1-t0

######################################################################################
######################################################################################
######################################################################################

#load("C:\\Users\\Nutzer\\Documents\\R Programme\\NBPSSCOX\\2_simu\\2_NBPSSvsMboost\\Results\\20231216\\Store.RData")
# etahatNBPSS1=Store[[1]][[1]]
# etahatNBPSS2=Store[[1]][[2]]
# etahatNBPSS3=Store[[1]][[3]]
# etahatmboost1=Store[[1]][[4]]
# etahatmboost2=Store[[1]][[5]]
# etahatmboost3=Store[[1]][[6]]
# phatNBPSS1=Store[[2]][[1]]
# phatNBPSS2=Store[[2]][[2]]
# phatNBPSS3=Store[[2]][[3]]
# phatmboost1=Store[[2]][[4]]
# phatmboost2=Store[[2]][[5]]
# phatmboost3=Store[[2]][[6]]
# etaStore=Store[[3]][[1]]
# R=40

### Compute MSEs and selection accuracies

gamma.ast=c(0,1,1,1,0,0,1,1,0,0,0,0,1,1,1,1,0,0,0) ## true model

#gamma.ast=c(0,1,1,1,0,0,1,1,0,0,0,0,1,1,1,1,0,0,1) ## true model with active spatial effect

RMSE_NBPSS1=rep(0,R)
RMSE_NBPSS2=rep(0,R)
RMSE_NBPSS3=rep(0,R)
RMSE_mboost1=rep(0,R)
RMSE_mboost2=rep(0,R)
RMSE_mboost3=rep(0,R)

missclass_NBPSS1=matrix(0,R,J)
missclass_NBPSS2=matrix(0,R,J)
missclass_NBPSS3=matrix(0,R,J)
missclass_mboost1=matrix(0,R,J)
missclass_mboost2=matrix(0,R,J)
missclass_mboost3=matrix(0,R,J)

for(r in 1:R){
  
  RMSE_NBPSS1[r]=sqrt(mean((etaStore[r,]-etahatNBPSS1[r,])^2))
  RMSE_NBPSS2[r]=sqrt(mean((etaStore[r,]-etahatNBPSS2[r,])^2))
  RMSE_NBPSS3[r]=sqrt(mean((etaStore[r,]-etahatNBPSS3[r,])^2))
  
  RMSE_mboost1[r]=sqrt(mean((etaStore[r,]-etahatmboost1[r,])^2))
  RMSE_mboost2[r]=sqrt(mean((etaStore[r,]-etahatmboost2[r,])^2))
  RMSE_mboost3[r]=sqrt(mean((etaStore[r,]-etahatmboost3[r,])^2))
  
  missclass_NBPSS1[r,]=abs(1.0*(phatNBPSS1[r,]>0.5)-gamma.ast)
  missclass_NBPSS2[r,]=abs(1.0*(phatNBPSS2[r,]>0.5)-gamma.ast)
  missclass_NBPSS3[r,]=abs(1.0*(phatNBPSS3[r,]>0.5)-gamma.ast)
  
  missclass_mboost1[r,]=abs(phatmboost1[r,]-gamma.ast)
  missclass_mboost2[r,]=abs(phatmboost2[r,]-gamma.ast)
  missclass_mboost3[r,]=abs(phatmboost3[r,]-gamma.ast)
}

RMSE_NBPSS1
RMSE_NBPSS2
RMSE_NBPSS3
RMSE_mboost1
RMSE_mboost2
RMSE_mboost3

colMeans(missclass_NBPSS1)
colMeans(missclass_mboost1)

RMSE=c(RMSE_mboost1,RMSE_mboost2,RMSE_mboost3,RMSE_NBPSS1,RMSE_NBPSS2,RMSE_NBPSS3)
g=c(rep("mboost (E)",R),rep("mboost (M)",R),rep("mboost (D)",R),rep("NBPSS (E)",R),rep("NBPSS (M)",R),rep("NBPSS (D)",R))

par(mar=c(4,4,2,2))
#boxplot(RMSE~g,main="Overall RMSE",xlab="")

names=c("mboost (Eilers)","mboost (MMR)","mboost (DRB)","NBPSS (Eilers)","NBPSS (MMR)","NBPSS (DRB)")
boxplot(RMSE_mboost1,RMSE_mboost2,RMSE_mboost3,RMSE_NBPSS1,RMSE_NBPSS2,RMSE_NBPSS3,names=rep("",6),main="Overall RMSE",xlab="")

# Draw the x-axis labels.
text(x = 1:6+0.3,
     ## Move labels to just below bottom of chart.
     y = par("usr")[3] -0.03,
     ## Use names from the data list.
     labels = names,
     ## Change the clipping region.
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 25,
     ## Adjust the labels to almost 100% right-justified.
     adj = 1,
     ## Increase label size.
     cex = 1)

# plot(rep(1,100),RMSE_mboost1,xlim=c(0,7),ylim=c(0,1))
# points(rep(2,100),RMSE_mboost2)
# points(rep(3,100),RMSE_mboost3)
# points(rep(4,100),RMSE_NBPSS1)
# points(rep(5,100),RMSE_NBPSS2)
# points(rep(6,100),RMSE_NBPSS3)
# 
# A=cbind(RMSE_mboost1,RMSE_mboost2,RMSE_mboost3,RMSE_NBPSS1,RMSE_NBPSS2,RMSE_NBPSS3)
# 
# for(i in 1:100){
#   lines(1:6,A[i,],col="lightgray")
# }

#Store=list(RMSE_mboost,RMSE_NBPSS,phatmboost,phatNBPSS)

M=rbind(colMeans(missclass_mboost1),colMeans(missclass_mboost2),colMeans(missclass_mboost3),colMeans(missclass_NBPSS1),colMeans(missclass_NBPSS2),colMeans(missclass_NBPSS3))*100
M=cbind(M,rowMeans(M))
M=xtable(M,digits=0)
M


#storeMstop=list(mstop1,mstop2,mstop3)
#save(storeMstop,file="C:/Users/Nutzer/Documents/R Programme/NBPSSCOX/6_SimulationResults/20231211/Mstop.Rdata")
# load("C:\\Users\\Nutzer\\Documents\\R Programme\\NBPSSCOX\\2_simu\\2_NBPSSvsMboost\\Results\\20231216\\Mstop.RData")
# mstop1=storeMstop[[1]][1:R]
# mstop2=storeMstop[[2]][1:R]
# mstop3=storeMstop[[3]][1:R]

boxplot(mstop1,mstop2,mstop3,ylim=c(0,10^3),main=paste("Optimal stopping iteration for mboost with df=",df,sep=""),names=c("Eilers","MMR","DRB"))

