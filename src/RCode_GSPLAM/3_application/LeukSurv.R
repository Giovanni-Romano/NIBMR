
################################################################################
################################################################################
################# LeukSurv #####################################################
################################################################################
################################################################################

## This script contains the code for the real data application (leukemia survival).

data("LeukSurv")

########################################### Example 3 LeukSurv (multivariate analysis without spatial effect)

set.seed(123)

y=LeukSurv$time
n=length(y)
delta=LeukSurv$cens

Xdummy=as.matrix(LeukSurv$sex)
Xcont=cbind(LeukSurv$age,LeukSurv$wbc,LeukSurv$tpi)
Xspat=cbind(LeukSurv$xcoord,LeukSurv$ycoord)

## posterior sampling

fit=nbpsscox(y=y,delta=delta,Xdummy=Xdummy,Xcont=Xcont,T=10^4,th=c(0.1,1),w2=1)

## evaluation

par(mfrow=c(2,2))

plot_hazard_rate(fit,ylim=c(-20,10))
plot_continuous_effects(fit,ylim=c(-1.5,2))


plot(colMeans(fit$samples$probi),ylim=c(0,1),pch=16,main="Posterior effect inclusion probabilities")
abline(h=0.5,lty="dashed")
plot(colMeans(fit$samples$gamma),ylim=c(0,1))

abline(h=0.5,col="red")
plot(fit$samples$theta[,1],type="l")

mcmc_trace(as_draws(fit$samples$probi))
summarize_draws(fit$samples$probi)

##################################################################################

### plot effect estimates

nd=1
nc=3

par(mfrow=c(3,3))

j=1
l=2
plot_effect_2(fit$design$Xcont[,j],fit$samples$beta[,c(j+nd,fit$design$indices[[j+(nd+nc)]])],fit$design$Trafo[[j]],ylim=c(-l,l),name=" of age")
j=2
plot_effect_2(fit$design$Xcont[,j],fit$samples$beta[,c(j+nd,fit$design$indices[[j+(nd+nc)]])],fit$design$Trafo[[j]],ylim=c(-l,l),name=" of wbc")
j=3
plot_effect_2(fit$design$Xcont[,j],fit$samples$beta[,c(j+nd,fit$design$indices[[j+(nd+nc)]])],fit$design$Trafo[[j]],name=" of tpi",ylim=c(-l,l))

### plot baseline hazard

par(mfrow=c(1,1))

plot_hazard_rate(fit,ylim=c(-15,0))

##################################################################################
##################################################################################

############ compare with bamlss

f=list(Surv(time,event=delta)~s(time),gamma~s(age)+s(wbc)+s(tpi))
fit2=bamlss(f,family="cox",data=LeukSurv)
plot(fit2)

plot(x,fit$fitted.values$gamma-mean(fit$fitted.values$gamma))
points(x,f,col="red")

########################################### Example 4 LeukSurv (multivariate analysis augmented with noise)

y=LeukSurv$time
delta=LeukSurv$cens

n=length(y);p=5
Xnoise=matrix(runif(n*p),n,p)

Xdummy=as.matrix(LeukSurv$sex)
Xcont=cbind(LeukSurv$age,LeukSurv$wbc,LeukSurv$tpi,Xnoise)
Xspat=cbind(LeukSurv$xcoord,LeukSurv$ycoord)

## posterior sampling

fit=besgacox(y=y,delta=delta,Xdummy=Xdummy,Xcont=Xcont,T=10^3)

## evaluation

par(mfrow=c(3,2))

plot_hazard_rate(fit,ylim=c(-20,10))
plot_continuous_effects(fit,ylim=c(-1.5,2))

plot(colMeans(fit$samples$gamma),ylim=c(0,1))
abline(h=0.5,col="red")

plot(fit$samples$theta,ylim=c(0,1))

mcmc_trace(as_draws(fit$samples$probi))


#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################

########################################### Example 5 LeukSurv (multivariate analysis with continuous spatial effect)

### this is the analysis shown in the manuscript

set.seed(123)

y=LeukSurv$time
n=length(y)
delta=LeukSurv$cens

Xdummy=as.matrix(LeukSurv$sex)
Xcont=cbind(LeukSurv$age,LeukSurv$wbc,LeukSurv$tpi)
Xspat=cbind(LeukSurv$xcoord,LeukSurv$ycoord)

## posterior sampling

time0=Sys.time()
fit=nbpsscox(y=y,delta=delta,Xdummy=Xdummy,Xcont=Xcont,Xspat=Xspat,T=10^3,th=c(0.1,1),w2=1,burn=10^2)
#fit=nbpsscox(y=y,delta=delta,Xdummy=Xdummy,Xcont=Xcont,Xspat=Xspat,T=5*10^4,th=c(0.1,1),w1=1,w2=1,burn=10^4,blockwise=FALSE,prior="SBP")
time1=Sys.time()
dtime=time1-time0

#save(fit,file="C:\\Users\\Nutzer\\Documents\\R Programme\\NBPSSCOX\\3_application\\fit.RData") ## store results
#load(file="C:\\Users\\Nutzer\\Documents\\R Programme\\NBPSSCOX\\3_application\\fit.RData") ## load fit with T=5*10^4

## evaluation

par(mfrow=c(2,2))

plot_hazard_rate(fit,ylim=c(-20,10))
plot_continuous_effects(fit,ylim=c(-1.5,2))

########################### posterior effect inclusion probabilities

par(mfrow=c(1,1))
par(mar=c(3,4,2,2))
plot(colMeans(fit$samples$probi),ylim=c(0,1),pch=16,main="Posterior effect inclusion probabilities",ylab="",xlab="",xaxt="n")
abline(h=0.5,lty="dashed")
axis(1,labels=FALSE)
names=c("sex(Dummy)","age(Linear)","wbc(Linear)","tpi(Linear)","age(Nonlin)","wbc(Nonlin)","tpi(Nonlin)","Space")

# Draw the x-axis labels.
text(x = 1:8+0.5,
     ## Move labels to just below bottom of chart.
     y = par("usr")[3] -0.05,
     ## Use names from the data list.
     labels = names,
     ## Change the clipping region.
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 35,
     ## Adjust the labels to almost 100% right-justified.
     adj = 1,
     ## Increase label size.
     cex = 1)
#points(colMeans(fit$samples$gamma),ylim=c(0,1))

ylab=TeX(r"($P(\gamma_j=1 |data)$)")
title(ylab=ylab, line=2.3, cex.lab=1)

#abline(h=0.5,col="red")
plot(fit$samples$theta[,1],type="l")

mcmc_trace(as_draws(fit$samples$probi))

xtable(x=summarize_draws(fit$samples$probi))

############################## plot covariate effect estimates

nd=1
nc=3

par(mfrow=c(3,3))
par(mar=c(5,4,2,1))
j=1
l=2
par(cex.main=1.5)
par(cex.lab=1.5)
par(cex.axis=1.5)

plot_effect_2(fit$design$Xcont[,j],fit$samples$beta[,c(j+nd,fit$design$indices[[j+(nd+nc)]])],fit$design$Trafo[[j]],ylim=c(-l,l),name=" of age",xlab="Age in years")
j=2
plot_effect_2(fit$design$Xcont[,j],fit$samples$beta[,c(j+nd,fit$design$indices[[j+(nd+nc)]])],fit$design$Trafo[[j]],ylim=c(-l,l),name=" of wbc",xlab="White blood cell count")
j=3
plot_effect_2(fit$design$Xcont[,j],fit$samples$beta[,c(j+nd,fit$design$indices[[j+(nd+nc)]])],fit$design$Trafo[[j]],name=" of tpi",ylim=c(-l,l),xlab="Townsend deprivation index")

### plot baseline hazard rate ###################################################################

par(mfrow=c(1,1))

par(mar=c(4,4,2,1))
plot_hazard_rate(fit,ylim=c(-10,5))
ylab=TeX(r"(log baseline hazard rate   $\widehat{g_0}(t)$)")
title(ylab=ylab, line=2.3, cex.lab=1.5)
title(xlab="Time in days", line=2.3, cex.lab=1.5)

### plot spatial effect ##########################################################################

par(mfrow=c(1,1))
par(mar=c(3.5,5,3,3))
par(cex=1.4)
par(cex.axis=1.5)
par(cex.main=2.2)
par(cex.lab=1.8)
plot_effect_spat(Xspat,fit$samples$beta[,fit$design$indices[[fit$design$J]]],fit$design$Trafos,ncol=30,ngrid=100,xlim=c(-0.1,0.9),ylim=c(-0.07,1.14),cex=0.5)
plotmap(nwengland,add=TRUE)

library(R2BayesX)
nwengland <- read.bnd(system.file("otherdata/nwengland.bnd",
                                  package = "spBayesSurv"))


