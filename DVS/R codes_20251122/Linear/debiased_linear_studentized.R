debiased_linear_studentized=function(Y,Z,U,r,alpha,B){
  
  n=dim(Z)[1];d=dim(Z)[2];q=dim(U)[2]+1;p=d+q
  ###first step:uniform sampling with replacement
  idx=which(runif(n) < r/n)
  z.sel=Z[idx,];u.sel=U[idx,];y.sel=Y[idx]
  ###initial estimator
  x.sel=cbind(z.sel,u.sel)
  tmp=cv.glmnet(x.sel,y.sel,alpha=1)
  mod=glmnet(x.sel,y.sel,lambda=tmp$lambda.min,alpha=1)
  b.hat=as.vector(c(mod$beta[1:d],mod$a0,mod$beta[-(1:d)]))
  b.hat[abs(b.hat)<0.05]=0
  theta.hat=b.hat[1:d]
  gamma.hat=b.hat[-c(1:d)]
  ###approximated projection matrix
  #if (q>1) {
  Hfit=mvr(z.sel,u.sel)
  W.pilot=rbind(as.vector(Hfit$muhat),Hfit$Bhat)
  rreg=cbind(1,U)
  temp1 <- Z%*%theta.hat+rreg%*%gamma.hat
  temp2 <- Z-rreg%*%W.pilot
  
  data.sel <- temp2[idx,]
  G.clime <- clime(x = data.sel, nlambda=10, linsolver = 'simplex', standardize = F, perturb = T)
  G.cv <- cv.clime(G.clime)
  G.hat <- clime(x = data.sel, G.cv$lambdaopt, linsolver = 'simplex', standardize = F, perturb = T)$Omegalist[[1]]
  G.hat = pmin(G.hat, t(G.hat))
  
  residual=as.vector(Y-temp1)
  S=residual*temp2
  theta.final=theta.hat+G.hat%*%apply(S,2,mean)
  
  #sigma.hat = sum((residual)^2)/(n-sum(b.hat!=0))
  #Phi.est=sigma.hat*(t(data.sel)%*%(data.sel))/r
  #std = sqrt(diag(G.hat %*% Phi.est %*% G.hat * sigma.hat))
  std = sqrt(diag(G.hat))
  
  c_n = sapply(1:B,function(xxx,n,std,G.hat,S){
    ei = rnorm(n)
    Q_ne.hat = -diag(std^(-1)) %*% G.hat %*% apply(S*ei, 2, mean)
    boots=max(abs(Q_ne.hat))
    return(boots)
  },n,std,G.hat,S)
  quantiles <- quantile(c_n, probs = 1 - alpha)*std
  
  return(list(est=theta.final,quantiles=quantiles))
  #return(list(est=b.final[(1:d)],quantiles=quantiles))
}