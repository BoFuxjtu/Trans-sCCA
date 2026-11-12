###full-data estimator
full_data_estimator_Logistic=function(Y,Z,U,alpha){
  t1=proc.time()
  n=dim(Z)[1];d=dim(Z)[2];q=dim(U)[2]+1
  
  X=cbind(Z,U)
  tmp=cv.glmnet(X,Y,family="binomial",alpha=1)
  mod=glmnet(X,Y,family="binomial",alpha=1,lambda=tmp$lambda.min)
  bhat.step1=as.vector(c(mod$beta[1:d],mod$a0,mod$beta[-(1:d)]))
  bhat.step1[abs(bhat.step1)<0.05]=0
  ###approximated projection matrix
  rreg=cbind(1,U)
  temp1 <- rreg%*%bhat.step1[-c(1:d)]
  th1=Z%*%bhat.step1[1:d]+temp1
  b_dev2_step1=as.vector(exp(th1)/(1+exp(th1))^2)
  Z.new=(sqrt(b_dev2_step1))*Z
  U.new=(sqrt(b_dev2_step1))*U
  Hfit=mvr(Z.new,U.new)
  W.pilot=rbind(as.vector(Hfit$muhat),Hfit$Bhat)
  rreg=cbind(1,U)
  temp2 <- Z-rreg%*%W.pilot
  
  model=full_second_Logistic(Y,Z,temp1,temp2)
  est.Second=model$est
  
  V <- ginv((t(temp2)%*%diag(b_dev2_step1)%*%(temp2))/(n-1))
  se.est <- sqrt(diag(V)/n)
  return(list(est.ini=bhat.step1[c(1:d)],est=est.Second,se=se.est))
}

full_second_Logistic=function(Y,Z,temp1,temp2){
  equ=function(par){
    tmp=apply(as.vector((Y-inv.logit(Z%*%par+temp1)))*temp2,2,mean)
    return(tmp)
  }
  est.final=nleqslv(rep(0,dim(Z)[2]),equ,control=list(xtol=1e-3,ftol=1e-3,btol=1e-3))$x
  return(list(est=est.final))
}