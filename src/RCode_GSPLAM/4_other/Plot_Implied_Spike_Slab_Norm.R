
################################################################################
################################################################################
########## Implied spike and slab on the Euclidean norm ########################
################################################################################
################################################################################

### In this script we generate plots of the marginal spike and slab and 
### the implied spike and slab on the norm for NBPSS. In particular, 
### we generate Figures 4 and 5 of the manuscript.

library(reticulate)

sc<-import("scipy.special")

sc$hyperu(2,3,c(4,5)) ## yep, seems to work fine

dnormNSBP<-function(r,a,b,c,d){
  
  kappa=2*gamma(b+d/2)
  kappa=kappa/((2*c)^(d/2)*gamma(d/2)*beta(a,b)) ## normalizing constant
  
  p=r^(d-1)*sc$hyperu(b+d/2,-a+d/2+1,r^2/(2*c))
  return(kappa*p)
  
} ### density of the Euclidean norm

r=seq(0.001,10,0.01)

pr=dnormNSBP(r,1/2,5,1/2,10)
plot(r,pr,type="l",xlim=c(0,10))

mean(pr)*10

dmargNSBP<-function(x,a,b,c,d){
  
  kappa=gamma(b+d/2)/((2*pi*c)^(d/2)*beta(a,b))  ## normalizing constant
  
  p=sc$hyperu(b+d/2,-a+d/2+1,sum(x^2)/(2*c))
  
  if(is.na(p)){
    p=sc$hyperu(b+d/2,-a+d/2+1+10^-5,sum(x^2)/(2*c))
  }
  
  return(kappa*p)
  
} ### marginal density

################################################################################
################################################################################
################################################################################

## Note: This is Figure 4 in the manuscript

### plot the marginal spike and slab for d=1

par(mfrow=c(1,2))
x=seq(-3,3,0.01)

ylab=TeX("$p(\\beta |  \\gamma=0)$ and $p(\\beta |  \\gamma=1)$")

a=1/2;b=5;c=1;d=1
plot(x,sapply(x,FUN=dmargNSBP,a,b,c,d),type="l",col="red",main="Marginal spike and slab",ylab="")
a=1/2;b=5;c=10;d=1
lines(x,sapply(x,FUN=dmargNSBP,a,b,c,d),type="l")

title(ylab=ylab, line=2.3, cex.lab=1)
title(xlab=expression(beta), line=2, cex.lab=1)
legend(x="topright",legend=c("","spike","","slab"),pch=c(16,16,16,16),col=c("white","red","white","black"),bty="n")

a=1/2;b=5;c=1;d=1
plot(x,log(sapply(x,FUN=dmargNSBP,a,b,c,d)),type="l",col="red",ylab="",main="Marginal spike and slab (log scale)")
a=1/2;b=5;c=10;d=1
lines(x,log(sapply(x,FUN=dmargNSBP,a,b,c,d)))

ylab=TeX("$\\log~p(\\beta|\\gamma=0)$ and $\\log~p(\\beta |  \\gamma=1)$")
title(ylab=ylab, line=2.3, cex.lab=1)
title(xlab=expression(beta), line=2, cex.lab=1)
legend(x="topright",legend=c("","spike","","slab"),pch=c(16,16,16,16),col=c("white","red","white","black"),bty="n")

################################################################################
################################################################################

## Note: This is Figure 5 in the manuscript

## for NBPSS

par(mfrow=c(2,2))
par(mar=c(2,2,2,2))
par(mar=c(3,4,2,1)+0.1)
par(lwd=2)

layout(matrix(c(1,1,2,3,4,5),ncol=2,byrow=TRUE),heights=c(0.5,3,3))
par(mar=c(0,0,0,0))
plot.new()
text(0.5,0.3,"Implied spike and slab of the norm (NBPSS vs. SSGL)",cex=2,font=2)
par(mar=c(3,4,2,1)+0.1)
r=seq(0,2.5,0.0001)
par(cex.axis=1.5)
lwd=3
d=1

a=1/2
b=5

m0=0.1 ## desired expectation under the spike
m1=1 

ylab=TeX("$p(r |  \\gamma=0)$ and $p(r |  \\gamma=1)$")
#ylab="Spike and slab on the norm"

#mean(sqrt(rbetapr(n=10000,shape1=a,shape2=b,scale=c)*rchisq(n=10000,df=d)))
#sqrt(c)*beta(a+1/2,b-1/2)/beta(a,b)*sqrt(2)*gamma((d+1)/2)/gamma(d/2) ## analytical formula for the mean

c0=(m0*beta(a,b)/beta(a+1/2,b-1/2)*gamma(d/2)/gamma((d+1)/2)*1/sqrt(2))^2
c1=(m1*beta(a,b)/beta(a+1/2,b-1/2)*gamma(d/2)/gamma((d+1)/2)*1/sqrt(2))^2

p0=dnormNSBP(r,a,b,c0,d)
p1=dnormNSBP(r,a,b,c1,d)

plot(r,p0,type="l",col="red",main="",ylab="")
lines(r,p1,type="l")
legend(x="topright",legend=c("spike     ","slab"),pch=c(16,16),col=c("red","black"),bty="n",cex=1.5)

title(ylab=ylab, line=2.3, cex.lab=1.5)
title(xlab="r", line=2, cex.lab=1.5)
title(main="NBPSS (d=1)",line=-2,cex.main=1.5)
d=10

#a=1/2;b=10;c0=1;c1=1

c0=(m0*beta(a,b)/beta(a+1/2,b-1/2)*gamma(d/2)/gamma((d+1)/2)*1/sqrt(2))^2
c1=(m1*beta(a,b)/beta(a+1/2,b-1/2)*gamma(d/2)/gamma((d+1)/2)*1/sqrt(2))^2

p0=dnormNSBP(r,a,b,c0,d)
p1=dnormNSBP(r,a,b,c1,d)

plot(r,p0,type="l",col="red",main="",ylab="")
lines(r,p1,type="l")
legend(x="topright",legend=c("spike     ","slab"),pch=c(16,16),col=c("red","black"),bty="n",cex=1.5)

title(ylab=ylab, line=2.3, cex.lab=1.5)
title(xlab="r", line=2, cex.lab=1.5)
title(main="NBPSS (d=10)",line=-2,cex.main=1.5)

## for SSGL #####################################################################

#par(mfrow=c(2,1))
#par(mar=c(2,2,2,2))
#par(mar=c(3,4,2,1)+0.1)

r=seq(0,2.5,0.0001)

d=1
a=(d+1)/2

m1=1 ## expected norm for the slab
m0=0.1 ## expected norm for the spike

plot(r,dgamma(r,shape=d,rate=d/m0),type="l",col="red",main="",ylab="")
lines(r,dgamma(r,shape=d,rate=d/m1),type="l")
legend(x="topright",legend=c("spike     ","slab"),pch=c(16,16),col=c("red","black"),bty="n",cex=1.5)

title(ylab=ylab, line=2.3, cex.lab=1.5)
title(xlab="r", line=2, cex.lab=1.5)
title(main="SSGL (d=1)",line=-2,cex.main=1.5)

d=10
a=(d+1)/2

plot(r,dgamma(r,shape=d,rate=d/m0),type="l",col="red",main="",ylab="")
lines(r,dgamma(r,shape=d,rate=d/m1),type="l")
legend(x="topright",legend=c("spike     ","slab"),pch=c(16,16),col=c("red","black"),bty="n",cex=1.5)

title(ylab=ylab, line=2.3, cex.lab=1.5)
title(xlab="r", line=2, cex.lab=1.5)
title(main="SSGL (d=10)",line=-2,cex.main=1.5)

#title("Implied spike and slab on the Euclidean norm", line = -1.5, outer = TRUE,cex.main=1.5)

#################################################################################
#################################################################################
#################################################################################

