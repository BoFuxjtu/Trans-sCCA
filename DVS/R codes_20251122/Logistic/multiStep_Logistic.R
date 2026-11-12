###multi-step estimation
multipleStep_logistic=function(Y,Z,U,r1,alpha,maxiteration,r_est){
  t1=proc.time()
  n=dim(Z)[1];d=dim(Z)[2];q=dim(U)[2]+1
  REG=cbind(Z,1,U)
  ###first sampling
  idx1.step1=which(runif(n) < r1/n)
  z.step1=Z[idx1.step1,];u.step1=U[idx1.step1,];y.step1=Y[idx1.step1]
  ###initial estimator
  x.step1=cbind(z.step1,u.step1)
  tmp=cv.glmnet(x.step1,y.step1,family="binomial",alpha=1)
  mod=glmnet(x.step1,y.step1,family="binomial",lambda=tmp$lambda.min,alpha=1)
  bhat.step1=as.vector(c(mod$beta[1:d],mod$a0,mod$beta[-(1:d)]))
  bhat.step1[abs(bhat.step1)<0.05]=0
  ###approximated projection matrix
  rreg=cbind(1,U)
  temp1 <- rreg%*%bhat.step1[-c(1:d)]
  th1=z.step1%*%bhat.step1[1:d]+temp1[idx1.step1,]
  b_dev2_step1=as.vector(exp(th1)/(1+exp(th1))^2)
  z.step1.new=(sqrt(b_dev2_step1))*z.step1
  u.step1.new=(sqrt(b_dev2_step1))*u.step1
  Hfit=mvr(z.step1.new,u.step1.new)
  W.pilot=rbind(as.vector(Hfit$muhat),Hfit$Bhat)
  rreg.step1=cbind(1,u.step1)
  #Phi.pilot=t(z.step1)%*%(z.step1-rreg.step1%*%W.pilot)
  
  rreg=cbind(1,U)
  temp2 <- Z-rreg%*%W.pilot
  temp2.sel=temp2[idx1.step1,]
  Phi.est=t(temp2.sel)%*%diag(b_dev2_step1)%*%z.step1/r1
  
  Phi_inv.est=ginv(Phi.est)
  iteration <- 0
  error <- rep(NA,maxiteration)
  error[1] <- 1
  
  est_new <- bhat.step1[1:d]
  for (i in 1:maxiteration){
    est_old <- est_new
    est_new <- est_old+Phi_inv.est%*%apply((as.vector(Y-inv.logit(Z%*%est_old+temp1)))*temp2,2,mean)
    error[i+1] <- norm(est_new-est_old,type="2")
    iteration <- iteration+1
    if((error[i+1]-error[i])/error[i]<0.01&&i>1){
      break
    }
  }
  se.est <- sqrt(diag(Phi_inv.est)/n)
  return(list(est1=bhat.step1[1:d],est=est_new,se=se.est,iteration=iteration))
  
}