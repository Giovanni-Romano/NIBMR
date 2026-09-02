################################################################################
################################################################################
## plot prior effect draws for NBPSS ###########################################
################################################################################
################################################################################

## Note: Here we generate the figures for the prior effect draws
## which are shown in the supplementary material.

neff=30 ## number of effect draws to plot

set.seed(123) ## generate data
x=sort(runif(n=1000))

### draw linear effect NBPSS

ylab="Prior effect draws"

par(mfrow=c(2,2))
par(mar=c(3,4,2,1))
#par(mar=c(4, 4, 2, 1) + 0.1)
par(cex.axis=1.2)
par(cex.main=1.5)
cex.lab=1.3

X=scale(matrix(x))
T=10^3
d=ncol(X)

scale=prior_scaling_sbp(X) ## perform prior scaling for linear effect design matrix

c=scale[1] ## spike

tau2=rbetapr(n=T,shape1=1/2,shape2=5,scale=c)
beta=matrix(rnorm(d*T),T,d) 
beta=Diagonal(x=sqrt(tau2))%*%beta ## the rows of beta contain prior coeff draws

F=beta%*%t(X) ## the rows of F contain prior effect draws

fsup=apply(abs(F),MARGIN=1,FUN=max)
quantile(fsup,0.9)

plot(x,F[1,],ylim=c(-1,1),type="l",main="linear effect, spike, NBPSS",ylab="")
for(j in 1:neff){
  lines(x,F[j,],type="l",ylab="")
}
abline(h=0.1,col="red",lty="dashed",lwd=2)
abline(h=-0.1,col="red",lty="dashed",lwd=2)
title(ylab=ylab, line=2.3, cex.lab=cex.lab)
title(xlab=expression(x[j]), line=2, cex.lab=cex.lab)

c=scale[2] ## slab

tau2=rbetapr(n=T,shape1=1/2,shape2=5,scale=c)
beta=matrix(rnorm(d*T),T,d) 
beta=Diagonal(x=sqrt(tau2))%*%beta ## the rows of beta contain prior coeff draws

F=beta%*%t(X) ## the rows of F contain prior function draws

fsup=apply(abs(F),MARGIN=1,FUN=max)
quantile(fsup,0.9)

plot(x,F[1,],ylim=c(-1,1),type="l",main="linear effect, slab, NBPSS",ylab="")
for(j in 1:neff){
  lines(x,F[j,],type="l",ylab="")
}
abline(h=1,col="red",lty="dashed",lwd=2)
abline(h=-1,col="red",lty="dashed",lwd=2)
title(ylab=ylab, line=2.3, cex.lab=cex.lab)
title(xlab=expression(x[j]), line=2, cex.lab=cex.lab)


### nonlinear effect NBPSS

X=DemmlerReinsch(x,d=9)$X

T=10^3
d=ncol(X)

scale=prior_scaling_sbp(X) ## perform prior scaling for nonlinear effect design matrix 

c=scale[1] ## spike

tau2=rbetapr(n=T,shape1=1/2,shape2=5,scale=c)

beta=matrix(rnorm(d*T),T,d) 
beta=Diagonal(x=sqrt(tau2))%*%beta ## the rows of beta contain prior coeff draws

F=beta%*%t(X) ## the rows of F contain prior function draws

fsup=apply(abs(F),MARGIN=1,FUN=max)
quantile(fsup,0.9)

plot(x,F[1,],ylim=c(-1,1),type="l",main="nonlinear effect, spike, NBPSS",ylab="")
for(j in 1:neff){
  lines(x,F[j,],type="l",ylab="")
}

abline(h=0.1,col="red",lty="dashed",lwd=2)
abline(h=-0.1,col="red",lty="dashed",lwd=2)
title(ylab=ylab, line=2.3, cex.lab=cex.lab)
title(xlab=expression(x[j]), line=2, cex.lab=cex.lab)


c=scale[2] ## slab

tau2=rbetapr(n=T,shape1=1/2,shape2=5,scale=c)
beta=matrix(rnorm(d*T),T,d) 
beta=Diagonal(x=sqrt(tau2))%*%beta ## the rows of beta contain prior coeff draws

F=beta%*%t(X) ## the rows of F contain prior function draws

fsup=apply(abs(F),MARGIN=1,FUN=max)
quantile(fsup,0.9)

plot(x,F[1,],ylim=c(-1,1),type="l",main="nonlinear effect, slab, NBPSS",ylab="")
for(j in 1:neff){
  lines(x,F[j,],ylab="")
}


abline(h=1,col="red",lty="dashed",lwd=2)
abline(h=-1,col="red",lty="dashed",lwd=2)
title(ylab=ylab, line=2.3, cex.lab=cex.lab)
title(xlab=expression(x[j]), line=2, cex.lab=cex.lab)


################################################################################
################################################################################
## plot prior effect draws for SSGL ############################################
################################################################################
################################################################################

### draw linear effect SSGL

X=scale(matrix(x))
T=10^3
d=ncol(X)

scale=prior_scaling_ga(X) ## perform prior scaling for linear effect design matrix

c=scale[1] ## spike

tau2=rgamma(n=T,shape=(d+1)/2,scale=c)
beta=matrix(rnorm(d*T),T,d) 
beta=Diagonal(x=sqrt(tau2))%*%beta ## the rows of beta contain prior coeff draws

F=beta%*%t(X) ## the rows of F contain prior effect draws

fsup=apply(abs(F),MARGIN=1,FUN=max)
quantile(fsup,0.9)

plot(x,F[1,],ylim=c(-1,1),type="l",main="linear effect, spike, SSGL",ylab="")
for(j in 1:neff){
  lines(x,F[j,],type="l",ylab="")
}
abline(h=0.1,col="red",lty="dashed",lwd=2)
abline(h=-0.1,col="red",lty="dashed",lwd=2)
title(ylab=ylab, line=2.3, cex.lab=cex.lab)
title(xlab=expression(x[j]), line=2, cex.lab=cex.lab)

c=scale[2] ## slab

tau2=rgamma(n=T,shape=(d+1)/2,scale=c)
beta=matrix(rnorm(d*T),T,d) 
beta=Diagonal(x=sqrt(tau2))%*%beta ## the rows of beta contain prior coeff draws

F=beta%*%t(X) ## the rows of F contain prior function draws

fsup=apply(abs(F),MARGIN=1,FUN=max)
quantile(fsup,0.9)

plot(x,F[1,],ylim=c(-1,1),type="l",main="linear effect, slab, SSGL",ylab="")
for(j in 1:neff){
  lines(x,F[j,],type="l",ylab="")
}
abline(h=1,col="red",lty="dashed",lwd=2)
abline(h=-1,col="red",lty="dashed",lwd=2)
title(ylab=ylab, line=2.3, cex.lab=cex.lab)
title(xlab=expression(x[j]), line=2, cex.lab=cex.lab)

### nonlinear effect SSGL

X=DemmlerReinsch(x,d=9)$X

T=10^3
d=ncol(X)

scale=prior_scaling_ga(X) ## perform prior scaling for nonlinear effect design matrix 

c=scale[1] ## spike

tau2=rgamma(n=T,shape=(d+1)/2,scale=c)
beta=matrix(rnorm(d*T),T,d) 
beta=Diagonal(x=sqrt(tau2))%*%beta ## the rows of beta contain prior coeff draws

F=beta%*%t(X) ## the rows of F contain prior function draws

fsup=apply(abs(F),MARGIN=1,FUN=max)
quantile(fsup,0.9)

plot(x,F[1,],ylim=c(-1,1),type="l",main="nonlinear effect, spike, SSGL",ylab="")
for(j in 1:neff){
  lines(x,F[j,],type="l",ylab="")
}

abline(h=0.1,col="red",lty="dashed",lwd=2)
abline(h=-0.1,col="red",lty="dashed",lwd=2)
title(ylab=ylab, line=2.3, cex.lab=cex.lab)
title(xlab=expression(x[j]), line=2, cex.lab=cex.lab)

c=scale[2] ## slab

tau2=rgamma(n=T,shape=(d+1)/2,scale=c)
beta=matrix(rnorm(d*T),T,d) 
beta=Diagonal(x=sqrt(tau2))%*%beta ## the rows of beta contain prior coeff draws

F=beta%*%t(X) ## the rows of F contain prior function draws

fsup=apply(abs(F),MARGIN=1,FUN=max)
quantile(fsup,0.9)

plot(x,F[1,],ylim=c(-1,1),type="l",main="nonlinear effect, slab, SSGL",ylab="")
for(j in 1:neff){
  lines(x,F[j,],ylab="")
}


abline(h=1,col="red",lty="dashed",lwd=2)
abline(h=-1,col="red",lty="dashed",lwd=2)
title(ylab=ylab, line=2.3, cex.lab=cex.lab)
title(xlab=expression(x[j]), line=2, cex.lab=cex.lab)


## Crucial difference: For the nonlinear effect under SSGL,
## there is not a single draw from the slab that looks like a draw from the spike

