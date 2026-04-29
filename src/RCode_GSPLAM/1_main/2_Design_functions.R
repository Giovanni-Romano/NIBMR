
#####################################################################################
#####################################################################################
####################### Design functions ############################################
#####################################################################################
#####################################################################################

## This script contains functions that are relevant for setting up the design matrices
## for the baseline hazard g0 and the geoadditive predictor eta.

## this auxiliary function plots the columns of a design matrix 

plot.design.matrix <- function(x,X,main="Design matrix",ylim=c(min(X),max(X)),col=FALSE,lwd=1,ylab="",xlab="")
{
  if(col==TRUE){
    col=1:dim(X)[2]
  } else {
    col=rep("black",dim(X)[2])
  }
  
  o=order(x)
  x=x[o]
  X=X[o,]
  
  plot(x,X[,1],ylim=ylim,type="l",main=main,col=col[1],lwd=lwd,ylab=ylab,xlab=xlab)
  for(j in 2:dim(X)[2]){
    lines(x,X[,j],col=col[j])
  }
}

## this function returns the nonlinear Demmler-Reinsch spline basis functions

DemmlerReinsch<-function(x,d=20){
  
  ## set up B-spline design
  
  fit=smoothCon(s(x,bs="ps",k=d+2),data=data.frame(x),scale.penalty=FALSE)
  
  B=fit[[1]]$X
  K=fit[[1]]$S[[1]]
  
  n=length(x)
  
  G=crossprod(B)/n
  
  ## realize orthogonality constraints
  
  X0=cbind(rep(1,n),scale(x))
  C=t(X0)%*%B ## constraint matrix for nonlinear effect
  
  #r=as.numeric(rankMatrix(C))
  
  V0=svd(C,nv=d+2)$v[,3:(d+2)]
  
  K2=t(V0)%*%K%*%V0 ## new penalty matrix
  G2=t(V0)%*%G%*%V0 ## new Gramian matrix
  
  ## solve generalized eigenvalue problem
  
  A=geigen(G2,K2)$vectors
  
  Trafo=V0%*%A
  
  X=B%*%Trafo
  
  frob=sqrt(sum(diag(t(X)%*%X/n)))
  
  X=X/frob ### unit Frobenius norm
  
  sum(diag(t(X)%*%X)) ## check
  
  Trafo=Trafo/frob
  
  ret=list()
  ret[[1]]=X
  ret[[2]]=Trafo
  names(ret)=c("X","Trafo")
  
  return(ret)
}

## this function returns the nonlinear mixed model spline basis functions

MixedModel<-function(x,d=20){
  
  ## set up B-spline design
  
  fit=smoothCon(s(x,bs="ps",k=d+2),data=data.frame(x),scale.penalty=FALSE)
  
  B=fit[[1]]$X
  K=fit[[1]]$S[[1]]
  
  E=eigen(K)
  
  U=E$vectors
  l=E$values
  
  Up=U[,-c(d+1,d+2)]
  lp=l[-c(d+1,d+2)]
  
  X=B%*%Up%*%diag(lp^(-1/2))
  
  frob=sqrt(sum(diag(t(X)%*%X/n)))
  
  X=X/frob ### unit Frobenius norm
  
  sum(diag(t(X)%*%X)) ## check
  
  Trafo=Up%*%diag(lp^(-1/2))/frob
  
  ret=list()
  ret[[1]]=X
  ret[[2]]=Trafo
  names(ret)=c("X","Trafo")
  
  return(ret)
} ## Fahrmeir et al. (2004)

## this function returns the nonlinear spline basis functions resulting from Eiler's transformation

Eilers<-function(x,d=20){
  
  ## set up B-spline design
  
  fit=smoothCon(s(x,bs="ps",k=d+2),data=data.frame(x),scale.penalty=FALSE)
  
  B=fit[[1]]$X
  D=fit[[1]]$D
  
  n=length(x)
  
  Trafo=t(solve(D%*%t(D))%*%D)
  
  X=B%*%Trafo
  
  frob=sqrt(sum(diag(t(X)%*%X/n)))
  
  X=X/frob ### unit Frobenius norm
  
  sum(diag(t(X)%*%X)) ## check
  
  Trafo=Trafo/frob
  
  ret=list()
  ret[[1]]=X
  ret[[2]]=Trafo
  names(ret)=c("X","Trafo")
  
  return(ret)
} ## Eilers (1999)

## this function returns the DR basis for continuous spatial effects

DemmlerReinschSpatial<-function(x1,x2,k1=10,k2=10){
  
  fit1=smoothCon(s(x1,bs="ps",k=k1,m=c(3,1)),data=data.frame(x1),scale.penalty=FALSE)
  fit2=smoothCon(s(x2,bs="ps",k=k2,m=c(3,1)),data=data.frame(x2),scale.penalty=FALSE)
  
  B1=fit1[[1]]$X
  B2=fit2[[1]]$X
  
  K1=fit1[[1]]$S[[1]]
  K2=fit2[[1]]$S[[1]]
  
  B=tensor.prod.model.matrix(list(B1,B2))
  K=Reduce("+",tensor.prod.penalties(list(K1,K2)))
  
  G=crossprod(B)/n
  
  C=t(rep(1,n))%*%B
  
  #C=t(cbind(rep(1,n),scale(x1),scale(x2),scale(x1)*scale(x2)))%*%B
  
  V0=svd(C,nv=k1*k2)$v[,2:(k1*k2)]
  #V0=svd(C,nv=100)$v[,5:100]
  
  G2=t(V0)%*%G%*%V0
  K2=t(V0)%*%K%*%V0
  
  fit=geigen(G2,K2)
  
  A=fit$vectors
  
  X=B%*%V0%*%A
  
  frob=sqrt(sum(diag(t(X)%*%X/n)))
  
  X=X/frob ### unit Frobenius norm
  
  sum(diag(t(X)%*%X)) ## check
  
  Trafo=V0%*%A/frob
  
  ret=list()
  ret[[1]]=X
  ret[[2]]=Trafo
  ret[[3]]=B
  names(ret)=c("X","Trafo","B")
  
  return(ret)
}

## integral approximation (for the baseline hazard rate)

## this function augments the observed times y with a fine equidistant grid from (0,ymax)
## we return the corresponding design matrix X0fine, the vector of differences of the time points
## and the indices of the original observations y within the augmented grid ynew

fineGrid<-function(y,length=200,d0=10,bs=bs){
  
  n=length(y)
  ynew=seq(0,max(y),length=length) ## augment observed times with fine equidistant grid for better accuracy 
  ynew=c(y,ynew) ## combine
  ynew=sort(ynew,index.return=TRUE) ## sort
  
  diff=c(0,diff(ynew$x)) ## compute differences for the augmented grid
  
  data=cbind(ynew$ix,1:length(ynew$ix))
  data=data[order(ynew$ix),]
  ind=data[1:n,2] ## these are the indices of y within ynew :-)
  
  ## set up B-spline design matrix on fine grid
  ynew=as.numeric(ynew$x)
  #fitfine=smoothCon(s(ynew,bs="ps",k=d0),data=data.frame(ynew)) 
  fitfine=smoothCon(s(ynew,bs=bs,k=d0),data=data.frame(ynew)) 
  
  #X0fine=Matrix(fitfine[[1]]$X,sparse=TRUE) ## for some reason this makes the sampler much slower
  
  X0fine=fitfine[[1]]$X
  
  ret=list(X0fine,diff,ind,ynew)
}

###################################################################################
###################################################################################
###################################################################################
