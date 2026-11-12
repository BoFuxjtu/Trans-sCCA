rm(list=ls())
library("rstudioapi")
setwd(dirname(getActiveDocumentContext()$path))
source("decorrelate_linear_algorithm.R")
source("compare_linear_algorithm.R")
source("full_data_estimator_linear.R")
source("debiased_linear_nonstudentized.R")
source("debiased_linear_studentized.R")
source("zhang_debiased_studentized.R")
source("zhang_debiased_nonstudentized.R")

library(MASS)
library(nleqslv)
library(Matrix)
library(glmnet)
library(mvtnorm)
library(mnormt)
library(boot)
library(SparseM)
library(quantreg)
library(pqr)
library(clime)


decorrelate_linear_high=function(n,d,q,r,error,X_type,itr.max,alpha,B){
  
  ###the dimension of x without intercept
  p=d+q
  ###the true coefficient
  #d2=3;bt=c(3,1.5,2,3,1.5,2,rep(0,p-d-d2))
  d2=0;#bt=c(3,1.5,2,3,1.5,2,rep(0,p-6))
  bt=c(sqrt(3),sqrt(3),sqrt(3),rep(0,p-3))
  
  est.duwA=est.duwL=est.dwA=est.dwL=est.uwA=est.uwL=est.wA=est.wL=est.dun=est.un=est.full=matrix(NA,itr.max,d)
  se.duwA=se.duwL=se.dwA=se.dwL=se.uwA=se.uwL=se.wA=se.wL=se.dun=se.un=se.full=matrix(NA,itr.max,d)
  cp.duwA=cp.duwL=cp.dwA=cp.dwL=cp.uwA=cp.uwL=cp.wA=cp.wL=cp.dun=cp.un=cp.full=matrix(NA,itr.max,d)
  
  mse.duwA=mse.duwL=mse.dwA=mse.dwL=mse.uwA=mse.uwL=mse.wA=mse.wL=mse.dun=mse.un=mse.full=rep(NA,itr.max)
  acp.duwA=acp.duwL=acp.dwA=acp.dwL=acp.uwA=acp.uwL=acp.wA=acp.wL=acp.dun=acp.un=acp.full=rep(NA,itr.max)
  al.duwA=al.duwL=al.dwA=al.dwL=al.uwA=al.uwL=al.wA=al.wL=al.dun=al.un=al.full=rep(NA,itr.max)
  #time.duwA=time.duwL=time.dwA=time.dwL=time.uwA=time.uwL=time.wA=time.wL=time.dun=time.un=rep(NA,itr.max)
  
  time.est.dwA=time.CI.dwA=time.est.dwL=time.CI.dwL=time.est.duwA=time.CI.duwA=time.est.duwL=time.CI.duwL=
    time.est.wA=time.CI.wA=time.est.wL=time.CI.wL=time.est.uwA=time.CI.uwA=time.est.uwL=time.CI.uwL=
    time.est.onestep=time.CI.onestep=time.est.dun=time.CI.dun=time.est.un=time.CI.un=time.est.full=time.CI.full=rep(NA,itr.max)
  
  mse.onestep=acp.onestep=acl.onestep=rep(NA,itr.max)
  se.onestep=cp.onestep=est.onestep=matrix(NA,itr.max,d)
  quan.onestep=array(NA, dim = c(itr.max, 2, d))
  
  
  se.debiased=est.debiased=matrix(NA,itr.max,d)
  quan.debiased=mse.debiased=time.debiased=al.debiased=cp.debiased=rep(NA,itr.max)
  
  se.debiased_stu=est.debiased_stu=quan.debiased_stu=matrix(NA,itr.max,d)
  mse.debiased_stu=time.debiased_stu=al.debiased_stu=cp.debiased_stu=rep(NA,itr.max)
  
  se.zhang_nonstu=est.zhang_nonstu=matrix(NA,itr.max,d)
  quan.zhang_nonstu=mse.zhang_nonstu=time.zhang_nonstu=al.zhang_nonstu=cp.zhang_nonstu=rep(NA,itr.max)
  
  se.zhang_stu=est.zhang_stu=quan.zhang_stu=matrix(NA,itr.max,d)
  mse.zhang_stu=time.zhang_stu=al.zhang_stu=cp.zhang_stu=rep(NA,itr.max)
  
  
  meanx=rep(0,p);corr=0.5;sigmax=diag(p)
  for (i in 1:p){
    for (j in 1:p){sigmax[i,j]=corr^(abs(i-j))}}
  
  
  
  for(itr in 1:itr.max){
    
    cat(itr)
    if(X_type=="N"){
      X=rmvnorm(n,mean=meanx,sigma=sigmax)
    }else if(X_type=="t"){ 
      X=rmvt(n,delta=meanx,sigma=sigmax,df=10)
    }
    
    Z=X[,(1:d)];U=X[,-(1:d)]
    ###the distribution of error term 
    if(error=="a"){
      err=rnorm(n,mean=0,sd=1)
    }else if(error=="b"){ 
      err=2*rt(n,df=5)
    }else if(error=="c"){
      err=2*rt(n,df=3)
    }else if(error=="d"){
      err=rnorm(n,mean=-1,sd=4)/2+rnorm(n,mean=1,sd=4)/2
    }else if(error=="e"){ 
      err=rexp(n,rate=0.2)}
    Y=Z%*%bt[1:d]+1+U%*%bt[-(1:d)]+err
    
    t1=proc.time()
    d.debiased=debiased_linear_nonstudentized(Y,Z,U,r,alpha,B)
    t2=proc.time()
    time.debiased=(t2-t1)[3][[1]]
    est.debiased[itr,]=d.debiased$est
    quan.debiased[itr]=d.debiased$quantiles
    cp.debiased[itr] <- (max(abs(bt[1:d]-est.debiased[itr,]))<=quan.debiased[itr])
    mse.debiased[itr] <- c(t(est.debiased[itr,]-bt[1:d])%*%(est.debiased[itr,]-bt[1:d]))
    al.debiased[itr] <- 2*quan.debiased[itr]
    
    t1=proc.time()
    d.debiased_stu=debiased_linear_studentized(Y,Z,U,r,alpha,B)
    t2=proc.time()
    time.debiased_stu=(t2-t1)[3][[1]]
    est.debiased_stu[itr,]=d.debiased_stu$est
    quan.debiased_stu[itr,]=d.debiased_stu$quantiles
    cp.debiased_stu[itr] <- all(abs(bt[1:d]-est.debiased_stu[itr,])<=quan.debiased_stu[itr,])
    mse.debiased_stu[itr] <- c(t(est.debiased_stu[itr,]-bt[1:d])%*%(est.debiased_stu[itr,]-bt[1:d]))
    al.debiased_stu[itr] <- 2*mean(quan.debiased_stu[itr,])
    
    t1=proc.time()
    zhang_nonstu=Zhang_nonstudentized(X,Y,alpha,B,d)
    t2=proc.time()
    time.zhang_nonstu=(t2-t1)[3][[1]]
    est.zhang_nonstu[itr,]=zhang_nonstu$est
    quan.zhang_nonstu[itr]=zhang_nonstu$quantiles
    cp.zhang_nonstu[itr] <- (max(abs(bt[1:d]-est.zhang_nonstu[itr,]))<=quan.zhang_nonstu[itr])
    mse.zhang_nonstu[itr] <- c(t(est.zhang_nonstu[itr,]-bt[1:d])%*%(est.zhang_nonstu[itr,]-bt[1:d]))
    al.zhang_nonstu[itr] <- 2*quan.zhang_nonstu[itr]
    
    t1=proc.time()
    zhang_stu=Zhang_studentized(X,Y,alpha,B,d)
    t2=proc.time()
    time.zhang_stu=(t2-t1)[3][[1]]
    est.zhang_stu[itr,]=zhang_stu$est
    quan.zhang_stu[itr,]=zhang_stu$quantiles
    cp.zhang_stu[itr] <- all(abs(bt[1:d]-est.zhang_stu[itr,])<=quan.zhang_stu[itr,])
    mse.zhang_stu[itr] <- c(t(est.zhang_stu[itr,]-bt[1:d])%*%(est.zhang_stu[itr,]-bt[1:d]))
    al.zhang_stu[itr] <- 2*mean(quan.zhang_stu[itr,])
    
    #d.full=full_data_estimator_linear(Y,Z,U,alpha)
    #est.full[itr,]=d.full$est;se.full[itr,]=d.full$se
    #cp.full[itr,]=abs(est.full[itr,]-bt[1:d])<se.full[itr,]*qnorm(1-alpha/2)
    #cp.full[itr] <- (max(abs(bt[1:d]-est.debiased[itr,]))<=max(se.full[itr,]*qnorm(1-alpha/2)))
    #acp.full[itr]=mean(cp.full[itr,])
    #al.full[itr]=mean(se.full[itr,])*2*1.96
    #time.est.full[itr]=d.full$time.est
    #time.CI.full[itr]=d.full$time.CI
  }
  
  result=list(MSE=c(mean(mse.debiased),mean(mse.debiased_stu),mean(mse.zhang_nonstu),mean(mse.zhang_stu),mean(mse.full)),
              CP=c(mean(cp.debiased),mean(cp.debiased_stu),mean(cp.zhang_nonstu),mean(cp.zhang_stu),mean(acp.full)),
              AL=c(mean(al.debiased),mean(al.debiased_stu),mean(al.zhang_nonstu),mean(al.zhang_stu),mean(al.full)),
              Time=c(mean(time.debiased),mean(time.debiased_stu),mean(time.zhang_nonstu),mean(time.zhang_stu)))
  
  return(result)
  
}

linear_high_result_Na_n50000_d50_q490_r5000<- decorrelate_linear_high(n=5e4,d=50,q=450,r=5000,error="a",X_type = "N",itr.max=500,alpha=0.05,B=500)
cat("linear_high_result_Na_n50000_d50_q490_r5000","\n")

linear_high_result_ta_n50000_d50_q450_r5000<- decorrelate_linear_high(n=5e4,d=50,q=450,r=5000,error="a",X_type = "t",itr.max=500,alpha=0.05,B=500)
cat("linear_high_result_ta_n50000_d50_q450_r5000","\n")

linear_high_result_Nb_n50000_d50_q450_r5000<- decorrelate_linear_high(n=50000,d=50,q=450,r=5000,error="b",X_type = "N",itr.max=500,alpha=0.05,B=500)
cat("linear_high_result_Nb_n50000_d50_q450_r5000","\n")

linear_high_result_tb_n50000_d50_q450_r5000<- decorrelate_linear_high(n=50000,d=50,q=450,r=5000,error="b",X_type = "t",itr.max=500,alpha=0.05,B=500)
cat("linear_high_result_tb_n50000_d50_q450_r5000","\n")


#n:full-sample size 50000
#d:dimension of the parameter z;
#q :dimension of the parameter u, q=50;#p=d+q+1
#r:initial subsample size
#error & X_type:the four different settings
#itr.max: simulation  replications
#alpha:significance level of 0.05