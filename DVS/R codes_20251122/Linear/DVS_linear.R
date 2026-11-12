###one-step estimation
DVS_linear=function(Y,Z,U,r2,r1,alpha,B){
  t1=proc.time()
  n=dim(Z)[1];d=dim(Z)[2];q=dim(U)[2]+1
  
  ###first sampling
  idx.step1=which(runif(n) < r1/n)
  z.step1=Z[idx.step1,];u.step1=U[idx.step1,];y.step1=Y[idx.step1]
  ###initial estimator
  x.step1=cbind(z.step1,u.step1)
  tmp=cv.glmnet(x.step1,y.step1,alpha=1)
  mod=glmnet(x.step1,y.step1,lambda=tmp$lambda.min,alpha=1)
  bhat.step1=as.vector(c(mod$beta[1:d],mod$a0,mod$beta[-(1:d)]))
  bhat.step1[abs(bhat.step1)<0.05]=0
  ###approximated projection matrix
  Hfit=mvr(z.step1,u.step1)
  W.pilot=rbind(as.vector(Hfit$muhat),Hfit$Bhat)
  
  rreg=cbind(1,U)
  temp1 <- rreg%*%bhat.step1[-c(1:d)]
  temp2 <- Z-rreg%*%W.pilot
  temp1.sel=temp1[idx.step1,];temp2.sel=temp2[idx.step1,]
  
  idx=which(runif(n) < r2/n)
  y.sel=Y[idx];z.sel=Z[idx,];
  
  model=second(y.sel,z.sel,temp1[idx,],temp2[idx,])
  est.Second=model$est
  
  #Phi.est=t(temp2.sel)%*%z.step1/r1
  Phi.est=t(temp2[idx,])%*%z.sel/r2
  Phi_inv.est1=ginv(Phi.est)
  
  est.final=est.Second+Phi_inv.est1%*%apply((as.vector(Y-Z%*%est.Second-temp1))*temp2,2,mean)
  r1_real=length(idx.step1)
  sigma.est=sum((Y-Z%*%bhat.step1[1:d]-temp1)^2)/(n-sum(bhat.step1=0))
  V11.est=sigma.est*(t(temp2.sel)%*%(temp2.sel))/(r1_real-1)
  #se.est=sqrt(sigma.est*diag(ginv(Phi.pilot2)%*%V11.est%*%ginv(Phi.pilot2))/n)
  tensor=array(runif(r1_real*d^2), dim = c(r1_real, d^2, d^2))
  for (i in (1:r1_real)) {
    temp=as.vector(temp2.sel[i,]%o%t(z.step1[i,]))
    tensor[i, , ]=temp%o%temp
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
  Phi.est.step1=t(temp2.sel)%*%z.step1/r1
  Phi_inv.est=ginv(Phi.est.step1)
  h=array(runif(B*d), dim = c(B, d))
  for (i in (1:B)) {
    h[i,]=-Phi_inv.est%*%UC2[i,]/sqrt(n)-Phi_inv.est%*%matrix(UC3[i, ], nrow = d, ncol = d)%*%Phi_inv.est%*%UC1[i,]/r2
  }
  
  lower_quantile=alpha/2
  upper_quantile=1-alpha/2
  quantiles=apply(h, 2, function(x) {
    c(lower = quantile(x, probs = lower_quantile),
      upper = quantile(x, probs = upper_quantile))
  })
  t3=proc.time()
  #time.est=(t2-t1)[3][[1]]
  #time.CI=(t3-t1)[3][[1]]
  return(list(size=r2,est=est.final,quantiles=quantiles))
  
}

second=function(y.sel,z.sel,temp1.sel,temp2.sel){
  equ=function(par){
    tmp=apply(as.vector((y.sel-z.sel%*%par-temp1.sel))*temp2.sel,2,mean)
    return(tmp)
  }
  est.final=nleqslv(rep(0,dim(z.sel)[2]),equ,control=list(xtol=1e-3,ftol=1e-3,btol=1e-3))$x
  return(list(est=est.final))
}