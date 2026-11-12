###one-step estimation
DVS_Logistic=function(Y,Z,U,r2,r1,alpha,B){
  n=dim(Z)[1];d=dim(Z)[2];q=dim(U)[2]+1
  ###first step:uniform sampling with replacement
  idx.step1=which(runif(n) < r1/n)
  z.step1=Z[idx.step1,];u.step1=U[idx.step1,];y.step1=Y[idx.step1]
  ###initial estimator
  x.step1=cbind(z.step1,u.step1)
  tmp=cv.glmnet(x.step1,y.step1,family="binomial",alpha=1)
  mod=glmnet(x.step1,y.step1,family="binomial",alpha=1,lambda=tmp$lambda.min)
  bhat.step1=as.vector(c(mod$beta[1:d],mod$a0,mod$beta[-(1:d)]))
  bhat.step1[abs(bhat.step1)<0.05]=0
  ###approximated projection matrix
  rreg=cbind(1,U)
  temp1 <- rreg%*%bhat.step1[-c(1:d)]
  th1=z.step1%*%bhat.step1[1:d]+temp1[idx.step1,]
  b_dev2_step1=as.vector(exp(th1)/(1+exp(th1))^2)
  z.step1.new=(sqrt(b_dev2_step1))*z.step1
  u.step1.new=(sqrt(b_dev2_step1))*u.step1
  Hfit=mvr(z.step1.new,u.step1.new)
  W.pilot=rbind(as.vector(Hfit$muhat),Hfit$Bhat)
  rreg.step1=cbind(1,u.step1)
  #Phi.pilot=t(z.step1)%*%(z.step1-rreg.step1%*%W.pilot)
  
  temp2 <- Z-rreg%*%W.pilot
  
  idx=which(runif(n) < r2/n)
  y.sel=Y[idx];z.sel=Z[idx,];u.sel=U[idx,]
  model=second(y.sel,z.sel,temp1[idx,],temp2[idx,])
  est.Second=model$est
  
  
  rreg.sel=cbind(1,u.sel)
  th2=z.sel%*%est.Second+temp1[idx,]
  b_dev2_sel=as.vector(exp(th2)/(1+exp(th2))^2)
  Phi.est=t(temp2[idx,])%*%diag(b_dev2_sel)%*%z.sel/r2
  Phi_inv.sel=ginv(Phi.est)
  
  est.final=est.Second+Phi_inv.sel%*%apply((as.vector(Y-inv.logit(Z%*%est.Second+temp1)))*temp2,2,mean)
  reg.step1=cbind(z.step1,1,u.step1)
  b_dev2.step1=as.vector(exp(th1)/(1+exp(th1))^2)
  Phi.est.step1=t(temp2[idx.step1,])%*%diag(b_dev2.step1)%*%z.step1/r1
  b_dev3.step1=as.vector(exp(th1)*(1-exp(th1))/(1+exp(th1))^3)
  temp2.step1=temp2[idx.step1,]
  r1_real=length(idx.step1)
  T.est.tenser=array(runif(r1_real*d^3), dim = c(r1_real, d, d^2))
  for (i in (1:r1_real)) {
    T.est.tenser[i,,]=b_dev3.step1[i]*as.vector(temp2.step1[i,])%o%as.vector(temp2.step1[i,])
  }
  T.est=apply(T.est.tenser, c(2, 3), mean)
  V11.est=(t(temp2.step1)%*%diag(b_dev2.step1)%*%temp2.step1)/(r1_real-1)
  #se.est=sqrt(diag(ginv(Phi.pilot2)%*%V11.est%*%ginv(Phi.pilot2))/n)
  tensor=array(runif(r1_real*d^4), dim = c(r1_real, d^2, d^2))
  for (i in (1:r1_real)) {
    temp=as.vector((temp2.step1[i,])%o%t(z.step1[i,]))
    tensor[i, , ]=(b_dev2.step1[i])^2*(temp%o%temp)
  }
  V12.est=sqrt(r2/n)*V11.est
  V33.est=apply(tensor, c(2, 3), mean)*(1-r2/n)
  V=rbind(
    cbind(V11.est, V12.est, matrix(0, nrow = d, ncol = d^2)),
    cbind(t(V12.est), V11.est, matrix(0, nrow = d, ncol = d^2)),
    cbind(matrix(0, nrow = d^2, ncol = 2*d), V33.est)
  )
  UC=rmvnorm(B,mean=rep(0, (2*d+d^2)),sigma=V)
  UC1=UC[,(1:d)]
  UC2=UC[,((d+1):(2*d))]
  UC3=UC[,((2*d+1):(2*d+d^2))]
  #UC3=matrix(UC[,((2*d+1):(2*d+d^2))],nrow=d,ncol=d)
  #tensor_UC3 <- array(NA, dim = c(r, d, d))
  #for (i in 1:r2) {
  #  tensor_UC3[i, , ] <- matrix(UC3[i, ], nrow = d, ncol = d)
  #}
  #Phi_inv.est=ginv(V11.est)
  h=array(runif(B*d), dim = c(B, d))
  Phi_inv.est=ginv(Phi.est.step1)
  for (i in (1:B)) {
    Phi_inv_U1=as.vector(Phi_inv.est%*%UC1[i,])
    h[i,]=Phi_inv.est%*%T.est%*%as.vector(Phi_inv_U1%o%Phi_inv_U1)/(2*r2)-Phi_inv.est%*%UC2[i,]/sqrt(n)-Phi_inv.est%*%matrix(UC3[i, ], nrow = d, ncol = d)%*%Phi_inv_U1/r2
  }
  
  lower_quantile=alpha/2
  upper_quantile=1-alpha/2
  quantiles=apply(h, 2, function(x) {
    c(lower = quantile(x, probs = lower_quantile),
      upper = quantile(x, probs = upper_quantile))
  })
  return(list(size=r2,est=est.final,quantiles=quantiles))
}

second=function(y.sel,z.sel,temp1,temp2){
  equ=function(par){
    tmp=apply(as.vector((y.sel-inv.logit(z.sel%*%par+temp1)))*temp2,2,mean)
    return(tmp)
  }
  est.final=nleqslv(rep(0,dim(z.sel)[2]),equ,control=list(xtol=1e-3,ftol=1e-3,btol=1e-3))$x
  return(list(est=est.final))
}