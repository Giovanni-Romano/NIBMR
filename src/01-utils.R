DemmlerReinsch<-function(x, d=20, rescale = TRUE){
  require(mgcv)
  require(geigen)
  
  ## Set up B-spline design
  fit = smoothCon(s(x, bs="ps", k=d+2), data=data.frame(x), scale.penalty=FALSE)
  B = fit[[1]]$X
  K = fit[[1]]$S[[1]]
  n = length(x)
  G = crossprod(B)/n
  
  ## realize orthogonality constraints
  X0 = cbind(rep(1, n), scale(x))
  C = t(X0)%*%B ## constraint matrix for nonlinear effect
  V0 = svd(C, nv=d+2)$v[,3:(d+2)]
  Ktilde =t (V0) %*% K %*% V0 ## new penalty matrix
  Gtilde = t(V0) %*% G %*% V0 ## new Gramian matrix
  
  ## solve generalized eigenvalue problem
  A = geigen(Gtilde, Ktilde)$vectors
  Tr = V0 %*% A
  Z = B %*% Tr
  
  if (rescale){
    frob = sqrt(sum(Z^2)/n)
    Z = Z/frob ### unit Frobenius norm
    TrM = Tr/frob
  } else {
    TrM = Tr
  }
  
  ret = list()
  ret[[1]] = X
  ret[[2]] = TrM
  names(ret) = c("X", "TrM")
  return(ret)
}

bar_m_n = function(beta, y, B, c, k){
  
  n = length(n)
  
  if (n != nrow(B)){stop("length of y and nrow of B do not coincide.")}
  if (length(beta) != ncol(B)){stop("length of beta and ncol of B do not coincide.")}
  
  r = y - B %*% beta
  s = sqrt(k^2*r^2 + c)
  t = exp(-s) * (k^2*r) / s
  
  out = -(1/n) * t(B) %*% t
  
  return(out)
}

pi = function(beta, lambda, K){
  p = length(lambda)
  
  P = Matrix::bdiag()
  mvnfast::dmvn()
}