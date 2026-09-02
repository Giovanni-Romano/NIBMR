##################################################################################
##################################################################################
##################################################################################
##################################################################################

##### In this script we plot the test functions used in the simulation study
##### The resulting figure is shown in the supplementary material. 

n=10^3
p=8
X=matrix(runif(n*p),n,p)

X=matrix(0,n,p)
for(j in 1:p){
  X[,j]=seq(0,1,length=n)
}

################################################### overall effects

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


F=matrix(0,n,p)

for(j in 1:8){
  F[,j]=myfunc(X[,j],nr=j)
}

colMeans(F)

################################################### linear effects

myfuncLin<-function(x,nr=1){
  
  if(nr==1){f=1/sqrt(12)*x}
  if(nr==2){f=-2/sqrt(12)*x}
  if(nr==3){f=0*x}
  if(nr==4){f=0*x}
  if(nr==5){f=-2*(5+exp(3))/(3*sqrt(3)*exp(3))*x}
  if(nr==6){f=-sqrt(3)/pi*x}
  if(nr>6){f=0*x}
  
  return(f)
  
}

FL=matrix(0,n,p)

for(j in 1:8){
  FL[,j]=myfuncLin(sqrt(12)*(X[,j]-1/2),nr=j)
}

################################################### nonlinear effects

myfuncNonLin<-function(x,nr=1){
  
  if(nr==1){f=0*x}
  if(nr==2){f=0*x}
  if(nr==3){f=4*(x-1/2)^2-1/3}
  if(nr==4){f=sin(2*pi*x)+6/pi*(x-1/2)}
  if(nr==5){f=2*exp(-3*x)-(2/3-2/3*exp(-3))--2*(5+exp(3))/(3*sqrt(3)*exp(3))*sqrt(12)*(x-1/2)}
  if(nr==6){f=sin(2*pi*x)+6/pi*(x-1/2)}
  if(nr>6){f=rep(0,length(x))}
  
  return(f)
  
}


FN=matrix(0,n,p)

for(j in 1:8){
  FN[,j]=myfuncNonLin(X[,j],nr=j)
}


################ check that the effect sums agree with the overall effects

max(abs(FN+FL-F))

######################################### plot effects

par(mfrow=c(8,3))
par(mar=c(2,2,2,1))
for(j in 1:8){
  plot(X[,j],FL[,j],ylim=c(-1.2,1.2),main=paste("linear effect of X",j,sep=""),type="l",lwd=2,cex.main=0.9)  
  plot(X[,j],FN[,j],ylim=c(-1.2,1.2),main=paste("nonlinear effect of X",j,sep=""),type="l",lwd=2,cex.main=0.9)  
  plot(X[,j],F[,j],ylim=c(-1.2,1.2),main=paste("overall effect of X",j,sep=""),type="l",lwd=2,cex.main=0.9)  
}


par(mfrow=c(8,3))
par(mar=c(2,2,1,1))
par(cex.lab=0.7)
par(cex.axis=0.7)
par(cex.main=1)


for(j in 1:8){
  
  plot(X[,j],FL[,j],ylim=c(-1.2,1.2),main="",type="l",ylab="",xlab="",yaxt="n",xaxt="n") 
  title(ylab="linear effect",line=1)
  title(xlab=as.expression(bquote(x[.(j)])),line=0.5,cex.lab=1)
  title(main=as.expression(bquote(L[.(j)])),line=0.5)
  axis(1,tck=-0.02,at=c(0,0.25,0.75,1),line=-1,lwd=0)
  axis(1,tck=-0.02,at=c(0,0.25,0.5,0.75,1),lwd=1,label=rep("",5))
  axis(2,tck=-0.02,at=c(-1,0,1),line=-1,lwd=0)
  axis(2,tck=-0.02,at=c(-1,0,1),lwd=1,label=rep("",3))
  
  plot(X[,j],FN[,j],ylim=c(-1.2,1.2),main="",type="l",ylab="",xlab="",xaxt="n",yaxt="n") 
  title(ylab="nonlinear effect",line=1)
  title(xlab=as.expression(bquote(x[.(j)])),line=0.5,cex.lab=1)
  title(main=as.expression(bquote(N[.(j)])),line=0.5)
  axis(1,tck=-0.02,at=c(0,0.25,0.75,1),line=-1,lwd=0)
  axis(1,tck=-0.02,at=c(0,0.25,0.5,0.75,1),lwd=1,label=rep("",5))
  axis(2,tck=-0.02,at=c(-1,0,1),line=-1,lwd=0)
  axis(2,tck=-0.02,at=c(-1,0,1),lwd=1,label=rep("",3))
  
  plot(X[,j],F[,j],ylim=c(-1.2,1.2),main="",type="l",ylab="",xlab="",xaxt="n",yaxt="n")  
  title(ylab="overall effect",line=1)
  title(xlab=as.expression(bquote(x[.(j)])),line=0.5,cex.lab=1)
  title(main=as.expression(bquote(L[.(j)]+N[.(j)])),line=0.5)
  axis(1,tck=-0.02,at=c(0,0.25,0.75,1),line=-1,lwd=0)
  axis(1,tck=-0.02,at=c(0,0.25,0.5,0.75,1),lwd=1,label=rep("",5))
  axis(2,tck=-0.02,at=c(-1,0,1),line=-1,lwd=0)
  axis(2,tck=-0.02,at=c(-1,0,1),lwd=1,label=rep("",3))
}
