
#################################################################################################
#################################################################################################
#################### IWLS functions #############################################################
#################################################################################################
#################################################################################################

### This script contains IWLS functions which are required for efficient MCMC sampling
### in the geoadditive Cox model.

## this function performs a single IWLS step for beta0
perform_iwls_step_beta0<-function(beta0,beta,X0,X,y,delta,K0,tau20,fine,step=1){
  
  ## make proposal for beta0
  
  fit=evaluate_log_FCP_beta0(beta0,beta,X0,X,y,delta,K0,tau20,fine)
  
  Sigma=-solve(fit[[3]])*step
  mu=beta0+Sigma%*%fit[[2]]
  
  beta0.ast=as.vector(rmvnorm(n=1,mean=mu,sigma=Sigma,checkSymmetry=FALSE))
  
  ## compute new moments given the proposal  
  
  fit.ast=evaluate_log_FCP_beta0(beta0.ast,beta,X0,X,y,delta,K0,tau20,fine)
  
  Sigma.ast=-solve(fit.ast[[3]])*step
  mu.ast=beta0.ast+Sigma.ast%*%fit.ast[[2]]
  
  ## compute MH acceptance probability
  
  alpha=fit.ast[[1]]-fit[[1]]
  alpha=alpha+dmvnorm(beta0,mu.ast,Sigma.ast,log=TRUE,checkSymmetry=FALSE)-dmvnorm(beta0.ast,mu,Sigma,log=TRUE,checkSymmetry=FALSE)
  alpha=min(1,exp(alpha))
  
  ## accept or reject proposal
  
  if(rbinom(n=1,size=1,prob=alpha)==1){
    return(beta0.ast)
  } else {
    return(beta0)
  }
  
}

## this function returns the log FCP of beta0 as well as the gradient and Hessian 
evaluate_log_FCP_beta0<-function(beta0,beta,X0,X,y,delta,K0,tau20,fine){
  
  beta0=as.vector(beta0)
  
  ret=list()
  
  ## function value
  
  ret[[1]]=log_lik_cox(beta0,beta,X0,X,y,delta,fine)-t(beta0)%*%K0%*%beta0/(2*tau20)
  
  ## gradient
  
  ret[[2]]=grad_log_lik_cox_beta0(beta0,beta,X0,X,y,delta,fine)-K0%*%beta0/tau20
  
  ## Hessian
  
  H=hess_log_lik_cox_beta0(beta0,beta,X0,X,y,delta,fine)-K0/tau20
  
  E=eigen(H,symmetric=TRUE)
  H=E$vectors%*%diag(pmin(E$values,-0.01))%*%t(E$vectors) ## Hessian modification
  print(max(eigen(H)$values))
  
  ret[[3]]=H
  
  return(ret)
  
}

## this function performs a single IWLS step for beta
perform_iwls_step_beta<-function(beta0,beta,X0,X,y,delta,Omega,fine,step=1){

  ## make proposal for beta
  
  fit=evaluate_log_FCP_beta(beta0,beta,X0,X,y,delta,Omega,fine)
  
  Sigma=-solve(fit[[3]])*step
  mu=beta+Sigma%*%fit[[2]]
  
  beta.ast=as.vector(rmvnorm(n=1,mean=mu,sigma=Sigma,checkSymmetry=FALSE))
  
  ## compute new moments given the proposal  
  
  fit.ast=evaluate_log_FCP_beta(beta0,beta.ast,X0,X,y,delta,Omega,fine)
  
  Sigma.ast=-solve(fit.ast[[3]])*step
  mu.ast=beta.ast+Sigma.ast%*%fit.ast[[2]]
  
  ## compute MH acceptance probability
  
  alpha=fit.ast[[1]]-fit[[1]]
  alpha=alpha+dmvnorm(beta,mu.ast,Sigma.ast,log=TRUE,checkSymmetry=FALSE)-dmvnorm(beta.ast,mu,Sigma,log=TRUE,checkSymmetry=FALSE)
  alpha=min(1,exp(alpha))
  
  ## accept or reject proposal
  
  if(rbinom(n=1,size=1,prob=alpha)==1){
    return(beta.ast)
  } else {
    return(beta)
  }

}

## this function returns the log FCP of beta as well as the gradient and Hessian
## Omega is the precision matrix Omega=Omega(tau21,...tau2J)
evaluate_log_FCP_beta<-function(beta0,beta,X0,X,y,delta,Omega,fine){

  ret=list()
  
  ## function value
  
  ret[[1]]=log_lik_cox(beta0,beta,X0,X,y,delta,fine)-t(beta)%*%Omega%*%beta/2
  
  ## gradient
  
  ret[[2]]=grad_log_lik_cox_beta(beta0,beta,X0,X,y,delta,fine)-Omega%*%beta
  
  ## Hessian
  
  H=hess_log_lik_cox_beta(beta0,beta,X0,X,y,delta,fine)-Omega
  
  E=eigen(H,symmetric=TRUE)
  H=E$vectors%*%diag(pmin(E$values,-1))%*%t(E$vectors) ## Hessian modification
  #print(eigen(H)$values)
  
  ret[[3]]=H
  
  return(ret)

}

############### Cox log-likelihood + gradient + Hessian ########################

## this function returns the Cox log-likelihood
log_lik_cox<-function(beta0,beta,X0,X,y,delta,fine){
  
  beta0=as.vector(beta0)
  
  g0=X0%*%beta0
  eta=X%*%beta
  
  X0fine=fine[[1]]; diff=fine[[2]]; ind=fine[[3]] 
  
  g0fine=X0fine%*%beta0 
  
  I=cumsum(exp(g0fine)*diff)[ind] ## approximate cumulative baseline hazards
  
  ll=t(delta)%*%(g0+eta)-sum(I*exp(eta))
  
  return(ll)
}

## this function returns the gradient of the Cox log-likelihood wrt beta0
grad_log_lik_cox_beta0<-function(beta0,beta,X0,X,y,delta,fine){
  
  eta=X%*%beta
  g0=X0%*%beta0
  
  ## we need the (approximate) Jacobian of I of size n x d0
  
  n=length(y); d0=ncol(X0)
  
  X0fine=fine[[1]]; diff=fine[[2]]; ind=fine[[3]]
  
  g0fine=X0fine%*%beta0 
  
  J=matrix(0,n,d0)
  for(j in 1:d0){
    J[,j]=cumsum(X0fine[,j]*exp(g0fine)*diff)[ind]  
  }
  
  return(t(X0)%*%delta-t(J)%*%exp(eta))
  
}

## this function returns the Hessian of the Cox log-likelihood wrt beta0
hess_log_lik_cox_beta0<-function(beta0,beta,X0,X,y,delta,fine){
  
  eta=X%*%beta
  g0=X0%*%beta0
  
  ## we need the (approximate) Hessian of I(i) of size d0 x d0
  
  d0=ncol(X0)
  
  X0fine=fine[[1]]; diff=fine[[2]]; ind=fine[[3]]
  g0fine=X0fine%*%beta0 
  
  H=matrix(0,d0,d0)
  for(j in 1:d0){
    for(k in 1:j){
      H[j,k]=-sum(cumsum(X0fine[,j]*X0fine[,k]*exp(g0fine)*diff)[ind]*exp(eta))
    }
  }
  
  H[upper.tri(H)]=t(H)[upper.tri(H)]
  
  return(H)
  
}

## this function returns the gradient of the Cox log-likelihood wrt beta
grad_log_lik_cox_beta<-function(beta0,beta,X0,X,y,delta,fine){
  
  eta=X%*%beta
  
  X0fine=fine[[1]]; diff=fine[[2]]; ind=fine[[3]] 
  
  g0fine=X0fine%*%beta0 
  
  I=cumsum(exp(g0fine)*diff)[ind] ## approximate cumulative baseline hazards
  
  return(t(X)%*%(delta-I*exp(eta)))
  
}

## this function returns the Hessian of the Cox log-likelihood wrt beta
hess_log_lik_cox_beta<-function(beta0,beta,X0,X,y,delta,fine){
  
  eta=X%*%beta
  
  X0fine=fine[[1]]; diff=fine[[2]]; ind=fine[[3]] 
  
  g0fine=X0fine%*%beta0
  
  I=cumsum(exp(g0fine)*diff)[ind] ## approximate cumulative baseline hazards
  
  #H=-t(X)%*%Diagonal(x=I*exp(eta))%*%X
  H=-crossprod(Diagonal(x=sqrt(I*exp(eta)))%*%X)
  
  H=as.matrix(H)
  
  return(H)
  
}

#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################

