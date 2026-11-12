rm(list=ls())
library("rstudioapi")
setwd(dirname(getActiveDocumentContext()$path))
source("decorrelate_logistic_algorithm.R")
source("compare_logistic_algorithm.R")
source("full_data_estimator_Logistic.R")
source("debiased_logistic_nonstudentized.R")
source("debiased_logistic_studentized.R")
source("van_nonstudentized_logistic.R")
source("van_studentized_logistic.R")


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


decorrelate_logistic_high=function(n,d,q,r,X_type,itr.max,alpha,B){
  
  ###the dimension of x without intercept
  p=d+q
  ###the true coefficient
  #d2=3;bt=c(3,1.5,2,3,1.5,2,rep(0,p-d-d2))
  d2=0;bt=c(1,-1,-1,rep(0,p-3))
  
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
  
  se.van_nonstu=est.van_nonstu=matrix(NA,itr.max,d)
  quan.van_nonstu=mse.van_nonstu=time.van_nonstu=al.van_nonstu=cp.van_nonstu=rep(NA,itr.max)
  
  se.van_stu=est.van_stu=quan.van_stu=matrix(NA,itr.max,d)
  mse.van_stu=time.van_stu=al.van_stu=cp.van_stu=rep(NA,itr.max)
  
  
  for(itr in 1:itr.max){
    
    cat(itr)
    
    meanx=rep(0,p);corr=0.5;sigmax=diag(p)
    for (i in 1:p){
      for (j in 1:p){sigmax[i,j]=corr^(abs(i-j))}} 
    
    
    #block_size <- p / 100
    
    #T_block <- diag(block_size)
    #for (i in 1:block_size) {
    #  for (j in 1:block_size) {
    #    if(i!=j){
    #    k <- abs(i - j)
    #      T_block[i, j] <- 1 * (block_size - k - 1) / (10 * (block_size - 1))}
    #  }
    #}
    #Sigma <- bdiag(replicate(100, T_block, simplify = FALSE))
    #sigmax <- as.matrix(Sigma)
    #browser()
    
    #X=rmvnorm(n,mean=meanx,sigma=sigmax/2)
    if(X_type=="N"){
      X=rmvnorm(n,mean=meanx,sigma=sigmax)/3
    }else if(X_type=="t"){ 
      X=rmvt(n,delta=meanx,sigma=sigmax,df=10)/3
    }
    
    
    Z=X[,(1:d)];U=X[,-(1:d)]
    Y=rbinom(n,1,inv.logit(X%*%bt)) 
    
    #t1=proc.time()
    #d.debiased=debiased_logistic_nonstudentized(Y,Z,U,r,alpha,B)
    #t2=proc.time()
    #time.debiased=(t2-t1)[3][[1]]
    #est.debiased[itr,]=d.debiased$est
    #quan.debiased[itr]=d.debiased$quantiles
    #cp.debiased[itr] <- (max(abs(bt[1:d]-est.debiased[itr,]))<=quan.debiased[itr])
    #mse.debiased[itr] <- c(t(est.debiased[itr,]-bt[1:d])%*%(est.debiased[itr,]-bt[1:d]))
    #al.debiased[itr] <- 2*quan.debiased[itr]
    
    t1=proc.time()
    d.debiased=debiased_logistic_nonstudentized(Y,Z,U,r,alpha,B)
    t2=proc.time()
    time.debiased=(t2-t1)[3][[1]]
    est.debiased[itr,]=d.debiased$est
    quan.debiased[itr]=d.debiased$quantiles
    cp.debiased[itr] <- (max(abs(bt[1:d]-est.debiased[itr,]))<=quan.debiased[itr])
    mse.debiased[itr] <- c(t(est.debiased[itr,]-bt[1:d])%*%(est.debiased[itr,]-bt[1:d]))
    al.debiased[itr] <- 2*quan.debiased[itr]
    
    t1=proc.time()
    d.debiased_stu=debiased_logistic_studentized(Y,Z,U,r,alpha,B)
    t2=proc.time()
    time.debiased_stu=(t2-t1)[3][[1]]
    est.debiased_stu[itr,]=d.debiased_stu$est
    quan.debiased_stu[itr,]=d.debiased_stu$quantiles
    cp.debiased_stu[itr] <- all(abs(bt[1:d]-est.debiased_stu[itr,])<=quan.debiased_stu[itr,])
    mse.debiased_stu[itr] <- c(t(est.debiased_stu[itr,]-bt[1:d])%*%(est.debiased_stu[itr,]-bt[1:d]))
    al.debiased_stu[itr] <- 2*mean(quan.debiased_stu[itr,])
    
    t1=proc.time()
    van_nonstu=van_nonstudentized_logistic(X,Y,alpha,B,d)
    t2=proc.time()
    time.van_nonstu=(t2-t1)[3][[1]]
    est.van_nonstu[itr,]=van_nonstu$est
    quan.van_nonstu[itr]=van_nonstu$quantiles
    cp.van_nonstu[itr] <- (max(abs(bt[1:d]-est.van_nonstu[itr,]))<=quan.van_nonstu[itr])
    mse.van_nonstu[itr] <- c(t(est.van_nonstu[itr,]-bt[1:d])%*%(est.van_nonstu[itr,]-bt[1:d]))
    al.van_nonstu[itr] <- 2*quan.van_nonstu[itr]
    
    t1=proc.time()
    van_stu=van_studentized_logistic(X,Y,alpha,B,d)
    t2=proc.time()
    time.van_stu=(t2-t1)[3][[1]]
    est.van_stu[itr,]=van_stu$est
    quan.van_stu[itr,]=van_stu$quantiles
    cp.van_stu[itr] <- all(abs(bt[1:d]-est.van_stu[itr,])<=quan.van_stu[itr,])
    mse.van_stu[itr] <- c(t(est.van_stu[itr,]-bt[1:d])%*%(est.van_stu[itr,]-bt[1:d]))
    al.van_stu[itr] <- 2*mean(quan.van_stu[itr,])
    
    #d.full=full_data_estimator_Logistic(Y,Z,U,alpha)
    #est.full[itr,]=d.full$est;se.full[itr,]=d.full$se
    #cp.full[itr,]=abs(est.full[itr,]-bt[1:d])<se.full[itr,]*qnorm(1-alpha/2)
    #mse.full[itr]=c(t(est.full[itr,]-bt[1:d])%*%(est.full[itr,]-bt[1:d]))
    #acp.full[itr]=mean(cp.full[itr,])
    #al.full[itr]=mean(se.full[itr,])*2*1.96
    #time.est.full[itr]=d.full$time.est
    #time.CI.full[itr]=d.full$time.CI
  }
  
  result=list(MSE=c(mean(mse.debiased),mean(mse.debiased_stu),mean(mse.van_nonstu),mean(mse.van_stu)),
              CP=c(mean(cp.debiased),mean(cp.debiased_stu),mean(cp.van_nonstu),mean(cp.van_stu)),
              AL=c(mean(al.debiased),mean(al.debiased_stu),mean(al.van_nonstu),mean(al.van_stu)),
              Time=c(mean(time.debiased),mean(time.debiased_stu),mean(time.van_nonstu),mean(time.van_stu)))
  
  return(result)
  
  
}

logistic_high_result_N_n50000_d50_q450_r5000<- decorrelate_logistic_high(n=50000,d=50,q=450,r=5000,X_type = "N",itr.max=500,alpha=0.05,B=500)
cat("logistic_high_result_N_n50000_d10_q490_r5000","\n")

logistic_high_result_t_n50000_d50_q450_r5000<- decorrelate_logistic_high(n=50000,d=50,q=450,r=5000,X_type = "t",itr.max=500,alpha=0.05,B=500)
cat("logistic_high_result_t_n50000_d10_q490_r5000","\n")


#n:full-sample size
#d:dimension of the parameter z;
#q :dimension of the parameter u, q=50,100;#p=d+q+1
#r2:informative subsample size r
#r1:initial subsample size r_p
#error:two settings of x
#itr.max: simulation  replications
#alpha:significance level of 0.05