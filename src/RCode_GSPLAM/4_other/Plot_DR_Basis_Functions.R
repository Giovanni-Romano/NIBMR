
########################################################## plot Demmler-Reinsch basis functions

## Here we generate Figure 2

set.seed(1234)
par(mfrow=c(1,2))
#par(mar=c(2,2,2,2))
par(mar=c(3,4,2,1)+0.1)

n=200

############ uniform design points

x=runif(n)
d=25
#d=17

T=DemmlerReinsch(x,d=d)$Trafo
t=seq(min(x),max(x),length=200)

fit=smoothCon(s(t,bs="ps",k=d+2),data=data.frame(t))
B=fit[[1]]$X
Z=B%*%T

lwd=2
#plot.design.matrix(t,Z,col=TRUE,lwd=lwd,main="DR splines (uniform design)",ylim=c(-2.5,1.2))
plot.design.matrix(t,-Z,col=TRUE,lwd=lwd,main=expression("DR splines for X"[j]*" ~ U[0,1]"))
rug(x)
#points(x,rep(-2.03,n),pch="I")

ylab="linear / nonlinear / overall effect"
ylab="DR spline basis functions"
title(ylab=ylab, line=2.3, cex.lab=1)
title(xlab=expression(x[j]), line=2, cex.lab=1)

############ exponential design points

x=rexp(n)
#x=rbeta(n,2,2)
#x=rbeta(n,2,2)
#fit=smoothCon(s(x,bs="ps",k=d+2),data=data.frame(x),scale.penalty=FALSE)
#X2=fit[[1]]$X

#X2=MixedModel(x,d=9)$X
T=DemmlerReinsch(x,d=d)$Trafo
t=seq(min(x),max(x),length=200)

fit=smoothCon(s(t,bs="ps",k=d+2),data=data.frame(t))
B=fit[[1]]$X
Z=B%*%T

#plot.design.matrix(t,Z,col=TRUE,lwd=lwd,main="DR splines (exponential design)",ylim=c(-2.5,1.2))
plot.design.matrix(t,-Z,col=TRUE,lwd=lwd,main=expression("DR splines for X"[j]*" ~ Exp(1)"))
rug(x)


title(ylab=ylab, line=2.3, cex.lab=1)
title(xlab=expression(x[j]), line=2, cex.lab=1)

