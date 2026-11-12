###one-step estimation
full_data_estimator_linear=function(Y,Z,U,alpha){
  #t1=proc.time()
  n=dim(Z)[1];d=dim(Z)[2];q=dim(U)[2]+1
  X=cbind(Z,U)
  tmp=cv.glmnet(X,Y,alpha=1)
  mod=glmnet(X,Y,lambda=tmp$lambda.min,alpha=1)
  bhat.step1=as.vector(c(mod$beta[1:d],mod$a0,mod$beta[-(1:d)]))
  bhat.step1[abs(bhat.step1)<0.05]=0
  ###approximated projection matrix
  Hfit=mvr(Z,U)
  W.pilot=rbind(as.vector(Hfit$muhat),Hfit$Bhat)
  rreg=cbind(1,U)
  temp1 <- rreg%*%bhat.step1[-c(1:d)]
  temp2 <- Z-rreg%*%W.pilot
  
  model=full_second_linear(Y,Z,temp1,temp2)
  est.Second=model$est
  sigma.est=sum((Y-Z%*%bhat.step1[1:d]-temp1)^2)/(n-sum(bhat.step1=0))
  #idx <- which(runif(n) < r_est/n)
  #r_real <- length(idx)
  #temp2.sel <- temp2[idx,]
  V <- ginv((t(temp2)%*%(temp2))/(n-1))
  se.est <- sqrt(sigma.est*diag(V)/n)
  return(list(est.ini=bhat.step1[c(1:d)],est=est.Second,se=se.est))
}

full_second_linear=function(Y,Z,temp1,temp2){
  equ=function(par){
    tmp=apply(as.vector((Y-Z%*%par-temp1))*temp2,2,mean)
    return(tmp)
  }
  est.final=nleqslv(rep(0,dim(Z)[2]),equ,control=list(xtol=1e-3,ftol=1e-3,btol=1e-3))$x
  return(list(est=est.final))
}