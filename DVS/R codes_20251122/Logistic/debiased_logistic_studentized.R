debiased_logistic_studentized=function(Y,Z,U,r,alpha,B){
  
  n=dim(Z)[1];d=dim(Z)[2];q=dim(U)[2]+1;p=d+q
  #rreg=cbind(1,U)
  
  ###first step:uniform sampling with replacement
  #idx=sample(1:n,r,T)
  idx=which(runif(n) < r/n)
  z.sel=Z[idx,];u.sel=U[idx,];y.sel=Y[idx]
  ###initial estimator
  x.sel=cbind(z.sel,u.sel)
  tmp=cv.glmnet(x.sel,y.sel,family="binomial",alpha=1)
  mod=glmnet(x.sel,y.sel,family="binomial",lambda=tmp$lambda.min,alpha=1)
  b.hat=as.vector(c(mod$beta[1:d],mod$a0,mod$beta[-(1:d)]))
  b.hat[abs(b.hat)<0.05]=0
  theta.hat=b.hat[1:d]
  gamma.hat=b.hat[-c(1:d)]
  ###approximated projection matrix
  rreg=cbind(1,U)
  temp1 <- rreg%*%gamma.hat+Z%*%theta.hat
  th.sel=temp1[idx,]
  b_dev2.sel=as.vector(exp(th.sel)/(1+exp(th.sel))^2)
  z.sel.new=(sqrt(b_dev2.sel))*z.sel
  u.sel.new=(sqrt(b_dev2.sel))*u.sel
  Hfit=mvr(z.sel.new,u.sel.new)
  W.hat=rbind(as.vector(Hfit$muhat),Hfit$Bhat)
  rreg.sel=rreg[idx,]
  temp2 <- Z-rreg%*%W.hat
  
  data.sel <- sqrt(b_dev2.sel)*temp2[idx,]
  G.clime <- clime(x = data.sel, nlambda=10, linsolver = 'simplex', standardize = F, perturb = F)
  G.cv <- cv.clime(G.clime)
  G.hat <- clime(x = data.sel, G.cv$lambdaopt, linsolver = 'simplex', standardize = F, perturb = F)$Omegalist[[1]]
  G.hat = pmin(G.hat, t(G.hat))
  S <- (as.vector(Y-inv.logit(temp1)))*temp2
  theta.final=theta.hat+G.hat%*%apply(S,2,mean)
  std = sqrt(diag(G.hat))
  c_n = sapply(1:B,function(xxx,n,std,G.hat,S){
    ei = rnorm(n)
    Q_ne.hat = - diag(std^(-1))%*% G.hat %*% apply(ei*S, 2, mean)
    quan=max(abs(Q_ne.hat))
    return(quan)
  },n,std,G.hat,S)
  
  quantiles <- quantile(c_n, probs = 1 - alpha)*std
  return(list(est=theta.final,quantiles=quantiles))
}