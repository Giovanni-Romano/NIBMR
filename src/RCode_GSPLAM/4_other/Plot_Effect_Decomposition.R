
### Here we create Figure 1 (effect decomposition)

#install.packages("latex2exp")
#library(latex2exp)

par(mfrow=c(1,2))
par(mar=c(2,2,2,2))
par(mar=c(3,3,2,1)+0.1)
par(mar=c(3,4,2,1)+0.1)

x=seq(0,1,length=200)

par(lwd=2)

### uniform design

f=sin(2*pi*x)
flin=-6/pi*(x-1/2)
fnonlin=f-flin
  
plot(x,f,type="l",main=expression("Effect decomposition for X"[j]*" ~ U[0,1]"),ylim=c(-1,1.5),ylab="",xlab="")
rect(0,0,1,1,col="gray",border=0,density=20,lwd=1)
lines(x,f,type="l",main=expression("X"[j]*" ~ U[0,1]"),ylim=c(-1,1.5))
lines(x,flin,lty="dashed")
lines(x,fnonlin,lty="dashed")

ylab="true linear / nonlinear / overall effect"
#ylab="y"
title(ylab=ylab, line=2.3, cex.lab=1)
title(xlab=expression(x[j]), line=2, cex.lab=1)

#lines(x,dbeta(x,1,1),lty="dotted",col="gray")

## beta(2,2) design

f=sin(2*pi*x)
flin=-90/pi^3*(x-1/2)
fnonlin=f-flin

plot(x,f,type="l",main=expression("Effect decomposition for X"[j]*" ~ Beta(2,2)"),ylim=c(-1,1.5),ylab="",xlab="")
polygon(x,dbeta(x,2,2),col="gray",border=0,density=20,lwd=1)
lines(x,f,type="l",main=expression("X"[j]*" ~ Beta(2,2)"),ylim=c(-1,1))
lines(x,flin,lty="dashed")
lines(x,fnonlin,lty="dashed")

#lines(x,dbeta(x,2,2),lty="dotted",col="orange")

title(ylab=ylab, line=2.3, cex.lab=1)
title(xlab=expression(x[j]), line=2, cex.lab=1)


