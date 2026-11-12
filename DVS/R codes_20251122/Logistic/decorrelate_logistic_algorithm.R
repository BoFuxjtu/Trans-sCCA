w_decorrelate_logistic=function(y.sel,z.sel,u.sel,pi.sel,bhat.step1,W.pilot,r2){
  t1=proc.time()
  d=dim(z.sel)[2]
  rreg.sel=cbind(1,u.sel)
  equ=function(par){
    tmp=apply(as.vector((y.sel-inv.logit(z.sel%*%par+rreg.sel%*%bhat.step1[-c(1:d)]))/pi.sel)*(z.sel-rreg.sel%*%W.pilot),2,mean)
    return(tmp)
  }
  est.final=nleqslv(rep(0,dim(z.sel)[2]),equ,control=list(xtol=1e-3,ftol=1e-3,btol=1e-3))$x
  t2=proc.time()
  
  est1=c(est.final,bhat.step1[-c(1:d)])
  reg.sel=cbind(z.sel,1,u.sel)
  b_dev2=as.vector(exp(reg.sel%*%est1)/(1+exp(reg.sel%*%est1))^2)
  J.est= ginv(t(z.sel)%*%diag(b_dev2/pi.sel)%*%(z.sel-rreg.sel%*%W.pilot)/r2)
  dw=b_dev2/(r2^2*pi.sel^2)-b_dev2/(r2^2*pi.sel)+b_dev2/(r2*pi.sel)
  V.est=t((z.sel-rreg.sel%*%W.pilot))%*%diag(dw)%*%(z.sel-rreg.sel%*%W.pilot)
  se.est=sqrt(diag(J.est%*%V.est%*%t(J.est)))
  t3=proc.time()
  time1=(t2-t1)[3][[1]]
  time2=(t3-t1)[3][[1]]
  return(list(size=r2,est=est.final,se=se.est,time1=time1,time2=time2))
}




uw_decorrelate_logistic=function(y.sel,z.sel,u.sel,pi.sel,bhat.step1,W.pilot,r2){
  t1=proc.time()
  d=dim(z.sel)[2]
  rreg.sel=cbind(1,u.sel)
  equ=function(par){
    tmp=apply(as.vector((y.sel-inv.logit(z.sel%*%par+
              rreg.sel%*%bhat.step1[-c(1:d)])))*(z.sel-rreg.sel%*%W.pilot),2,mean)
    return(tmp)
  }
  est.final=nleqslv(rep(0,dim(z.sel)[2]),equ,control=list(xtol=1e-3,ftol=1e-3,btol=1e-3))$x
  t2=proc.time()
  
  est1=c(est.final,bhat.step1[-c(1:d)])
  reg.sel=cbind(z.sel,1,u.sel)
  b_dev2=as.vector(exp(reg.sel%*%est1)/(1+exp(reg.sel%*%est1))^2)
  duw=pi.sel*b_dev2
  gamma.est=(t(z.sel)%*%diag(b_dev2)%*%(z.sel-rreg.sel%*%W.pilot))
  omega.est=(t(z.sel-rreg.sel%*%W.pilot)%*%diag(duw)%*%(z.sel-rreg.sel%*%W.pilot))
  se.est=sqrt(diag(ginv(gamma.est)+r2*ginv(gamma.est)%*%omega.est%*%t(ginv(gamma.est))))
  t3=proc.time()
  time1=(t2-t1)[3][[1]]
  time2=(t3-t1)[3][[1]]
  return(list(size=r2,est=est.final,se=se.est,time1=time1,time2=time2))
}




two_step_decorrelate_logistic=function(Y,Z,U,r2,r1,optimality,weight){
  t1=proc.time()
  n=dim(Z)[1];d=dim(Z)[2];q=dim(U)[2]+1
  REG=cbind(Z,1,U)
  
  ###first step:uniform sampling with replacement
  idx.step1=sample(1:n,r1,T)
  z.step1=Z[idx.step1,];u.step1=U[idx.step1,];y.step1=Y[idx.step1]
  ###initial estimator
  x.step1=cbind(z.step1,u.step1)
  tmp=cv.glmnet(x.step1,y.step1,family="binomial",alpha=1)
  mod=glmnet(x.step1,y.step1,family="binomial",alpha=1,lambda=tmp$lambda.min)
  bhat.step1=as.vector(c(mod$beta[1:d],mod$a0,mod$beta[-(1:d)]))
  bhat.step1[abs(bhat.step1)<0.05]=0
  ###approximated projection matrix
  b_dev2_full=as.vector(exp(REG%*%bhat.step1)/(1+exp(REG%*%bhat.step1))^2)
  b_dev2_part=b_dev2_full[idx.step1]
  z.step1.new=(sqrt(b_dev2_part))*z.step1
  u.step1.new=(sqrt(b_dev2_part))*u.step1
  Hfit=mvr(z.step1.new,u.step1.new)
  W.pilot=rbind(as.vector(Hfit$muhat),Hfit$Bhat)
  rreg.step1=cbind(1,u.step1)
  Phi.pilot=(t(z.step1)%*%diag(b_dev2_part)%*%(z.step1-rreg.step1%*%W.pilot))
  
  REG=cbind(Z,1,U)
 
  RREG=cbind(1,U)
  if(optimality=="L"){
    PI.Lopt=sqrt(b_dev2_full*rowSums((Z-RREG%*%W.pilot)^2))
    PI.OPT=PI.Lopt/(sum(PI.Lopt))
  }else if(optimality=="A"){
    PI.Aopt=sqrt(b_dev2_full*rowSums(((Z-RREG%*%W.pilot)%*%ginv(Phi.pilot))^2))
    PI.OPT=PI.Aopt/(sum(PI.Aopt))
  }
  
  if (optimality=="U"){PI.OPT=rep(1/n,n)}
  
  idx=sample(1:n,r2,T,PI.OPT)
  pi.sel=PI.OPT[idx];y.sel=Y[idx];z.sel=Z[idx,];u.sel=U[idx,]
  
  t2=proc.time()
  time.est1=(t2-t1)[3][[1]]
  
  if(weight=="yes"){
    model=w_decorrelate_logistic(y.sel,z.sel,u.sel,pi.sel,bhat.step1,W.pilot,r2)
    est.final=model$est;se.est=model$se
  }else if(weight=="no"){
    model=uw_decorrelate_logistic(y.sel,z.sel,u.sel,pi.sel,bhat.step1,W.pilot,r2)
    est.final=model$est;se.est=model$se
  }
  
  time.est=time.est1+model$time1
  time.CI=time.est1+model$time2
  
  return(list(size=r2,crition=optimality,statu=weight,est=est.final,se=se.est,time.est=time.est,time.CI=time.CI))
  
}
