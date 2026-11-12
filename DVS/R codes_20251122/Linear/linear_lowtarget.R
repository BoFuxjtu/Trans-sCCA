rm(list=ls())
library("rstudioapi")
setwd(dirname(getActiveDocumentContext()$path))
source("decorrelate_linear_algorithm.R")
source("multiStep_linear.R")
source("full_data_estimator_linear.R")
source("onesteplinear.R")
source("DVS_linear.R")

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


decorrelate_subsampling_linear_p400=function(n,d,q,r2,r1,error,X_type,itr.max,alpha){
  
  ###the dimension of x without intercept
  p=d+q
  ###the true coefficient
  #d2=3;bt=c(3,1.5,2,3,1.5,2,rep(0,p-d-d2))
  bt=c(sqrt(3),sqrt(3),sqrt(3),rep(0,p-3))
  
  est.duwA=est.duwL=est.dwA=est.dwL=est.uwA=est.uwL=est.wA=est.wL=est.dun=est.un=matrix(NA,itr.max,d)
  se.duwA=se.duwL=se.dwA=se.dwL=se.uwA=se.uwL=se.wA=se.wL=se.dun=se.un=matrix(NA,itr.max,d)
  cp.duwA=cp.duwL=cp.dwA=cp.dwL=cp.uwA=cp.uwL=cp.wA=cp.wL=cp.dun=cp.un=matrix(NA,itr.max,d)
  
  mse.duwA=mse.duwL=mse.dwA=mse.dwL=mse.uwA=mse.uwL=mse.wA=mse.wL=mse.dun=mse.un=rep(NA,itr.max)
  acp.duwA=acp.duwL=acp.dwA=acp.dwL=acp.uwA=acp.uwL=acp.wA=acp.wL=acp.dun=acp.un=rep(NA,itr.max)
  al.duwA=al.duwL=al.dwA=al.dwL=al.uwA=al.uwL=al.wA=al.wL=al.dun=al.un=rep(NA,itr.max)
  #time.duwA=time.duwL=time.dwA=time.dwL=time.uwA=time.uwL=time.wA=time.wL=time.dun=time.un=rep(NA,itr.max)
  
  time.dwA=time.dwL=time.duwA=time.duwL=
    time.est.wA=time.CI.wA=time.est.wL=time.CI.wL=time.est.uwA=time.CI.uwA=time.est.uwL=time.CI.uwL=
    time.est.onestep=time.CI.onestep=time.dun=time.un=rep(NA,itr.max)
  
  time.duwA=rep(NA,itr.max)
  
  #mse.onestep=acp.onestep=al.onestep=rep(NA,itr.max)
  #se.onestep=cp.onestep=est.onestep=matrix(NA,itr.max,d)
  #quan.onestep=array(NA, dim = c(itr.max, 2, d))
  
  time.full=mse.full=acp.full=al.full=rep(NA,itr.max)
  se.full=cp.full=est.full=matrix(NA,itr.max,d)
  quan.full=array(NA, dim = c(itr.max, 2, d))
  
  time.DVS=mse.DVS=acp.DVS=al.DVS=rep(NA,itr.max)
  se.DVS=cp.DVS=est.DVS=matrix(NA,itr.max,d)
  quan.DVS=array(NA, dim = c(itr.max, 2, d))
  
  est1=est_multiStep=se_multiStep=cp_multiStep=matrix(NA,itr.max,d)
  mse1=mse_multiStep=acp_multiStep=al_multiStep=rep(NA,itr.max)
  iteration_multiStep=time_multiStep=rep(NA,itr.max)
  
  est_multiStep2=se_multiStep2=cp_multiStep2=matrix(NA,itr.max,d)
  mse_multiStep2=acp_multiStep2=al_multiStep2=rep(NA,itr.max)
  iteration_multiStep2=time_multiStep2=rep(NA,itr.max)
  
  mse.onestep=acp.onestep=al.onestep=rep(NA,itr.max)
  se.onestep=cp.onestep=est.onestep=matrix(NA,itr.max,d)
  quan.onestep=time.onestep=rep(NA,itr.max)
  
  
  
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
      err=rt(n,df=5)
    }else if(error=="c"){
      err=2*rt(n,df=3)
    }else if(error=="d"){
      err=rnorm(n,mean=-1,sd=4)/2+rnorm(n,mean=1,sd=4)/2
    }else if(error=="e"){ 
      err=rexp(n,rate=0.2)}
    
    Y=Z%*%bt[1:d]+1+U%*%bt[-(1:d)]+err
    r_est <- 500
    
    
    t1=proc.time()
    d.multiStep=multiStep_linear(Y,Z,U,r1,alpha,100)
    t2=proc.time()
    time_multiStep[itr]=(t2-t1)[3][[1]]
    est1[itr,]=d.multiStep$est1
    mse1[itr]=c(t(est1[itr,]-bt[1:d])%*%(est1[itr,]-bt[1:d]))
    est_multiStep[itr,]=d.multiStep$est
    mse_multiStep[itr]=c(t(est_multiStep[itr,]-bt[1:d])%*%(est_multiStep[itr,]-bt[1:d]))
    se_multiStep[itr,]=d.multiStep$se
    cp_multiStep[itr,]=abs(est_multiStep[itr,]-bt[1:d])<se_multiStep[itr,]*qnorm(1-alpha/2)
    acp_multiStep[itr]=mean(cp_multiStep[itr,])
    al_multiStep[itr]=mean(se_multiStep[itr,])*2*1.96
    iteration_multiStep[itr] <- d.multiStep$iteration
    
    
    t1=proc.time()
    d.DVS=DVS_linear(Y,Z,U,r2,r1,alpha,500)
    t2=proc.time()
    time.DVS[itr]=(t2-t1)[3][[1]]
    est.DVS[itr,]=d.DVS$est
    quan.DVS[itr, , ]=d.DVS$quantiles
    cp.DVS[itr,] <- (bt[1:d] >= est.DVS[itr,]-quan.DVS[itr,2, ]) & (bt[1:d] <= est.DVS[itr,]-quan.DVS[itr,1, ])
    mse.DVS[itr]=c(t(est.DVS[itr,]-bt[1:d])%*%(est.DVS[itr,]-bt[1:d]))
    acp.DVS[itr]=mean(cp.DVS[itr,])
    al.DVS[itr]=mean(quan.DVS[itr,2, ]-quan.DVS[itr,1, ])
    
    
    #Compared decorrelated method
    t1=proc.time()
    duwA=two_step_decorrelate_linear(Y,Z,U,r2,r1,optimality="A",weight="no")
    t2=proc.time()
    time.duwA[itr]=(t2-t1)[3][[1]]
    est.duwA[itr,]=duwA$est;se.duwA[itr,]=duwA$se
    cp.duwA[itr,]=abs(est.duwA[itr,]-bt[1:d])<se.duwA[itr,]*qnorm(1-alpha/2)
    mse.duwA[itr]=c(t(est.duwA[itr,]-bt[1:d])%*%(est.duwA[itr,]-bt[1:d]))
    acp.duwA[itr]=mean(cp.duwA[itr,])
    al.duwA[itr]=mean(se.duwA[itr,])*2*1.96
    
    t1=proc.time()
    duwL=two_step_decorrelate_linear(Y,Z,U,r2,r1,optimality="L",weight="no")
    t2=proc.time()
    time.duwL[itr]=(t2-t1)[3][[1]]
    est.duwL[itr,]=duwL$est;se.duwL[itr,]=duwL$se
    cp.duwL[itr,]=abs(est.duwL[itr,]-bt[1:d])<se.duwL[itr,]*qnorm(1-alpha/2)
    mse.duwL[itr]=c(t(est.duwL[itr,]-bt[1:d])%*%(est.duwL[itr,]-bt[1:d]))
    acp.duwL[itr]=mean(cp.duwL[itr,])
    al.duwL[itr]=mean(se.duwL[itr,])*2*1.96
    
    
    t1=proc.time()
    dwA=two_step_decorrelate_linear(Y,Z,U,r2,r1,optimality="A",weight="yes")
    t2=proc.time()
    time.dwA[itr]=(t2-t1)[3][[1]]
    est.dwA[itr,]=dwA$est;se.dwA[itr,]=dwA$se
    cp.dwA[itr,]=abs(est.dwA[itr,]-bt[1:d])<se.dwA[itr,]*qnorm(1-alpha/2)
    mse.dwA[itr]=c(t(est.dwA[itr,]-bt[1:d])%*%(est.dwA[itr,]-bt[1:d]))
    acp.dwA[itr]=mean(cp.dwA[itr,])
    al.dwA[itr]=mean(se.dwA[itr,])*2*1.96
    
    
    t1=proc.time()
    dwL=two_step_decorrelate_linear(Y,Z,U,r2,r1,optimality="L",weight="yes")
    t2=proc.time()
    time.dwL[itr]=(t2-t1)[3][[1]]
    est.dwL[itr,]=dwL$est;se.dwL[itr,]=dwL$se
    cp.dwL[itr,]=abs(est.dwL[itr,]-bt[1:d])<se.dwL[itr,]*qnorm(1-alpha/2)
    mse.dwL[itr]=c(t(est.dwL[itr,]-bt[1:d])%*%(est.dwL[itr,]-bt[1:d]))
    acp.dwL[itr]=mean(cp.dwL[itr,])
    al.dwL[itr]=mean(se.dwL[itr,])*2*1.96
    
    
    t1=proc.time()
    dun=two_step_decorrelate_linear(Y,Z,U,r2,r1,optimality="U",weight="no")
    t2=proc.time()
    time.dun[itr]=(t2-t1)[3][[1]]
    est.dun[itr,]=dun$est;se.dun[itr,]=dun$se
    cp.dun[itr,]=abs(est.dun[itr,]-bt[1:d])<se.dun[itr,]*qnorm(1-alpha/2)
    mse.dun[itr]=c(t(est.dun[itr,]-bt[1:d])%*%(est.dun[itr,]-bt[1:d]))
    acp.dun[itr]=mean(cp.dun[itr,])
    al.dun[itr]=mean(se.dun[itr,])*2*1.96
    
    
    #t1=proc.time()
    #d.full=onestep_linear(Y,Z,U,alpha,r_est,n)
    #t2=proc.time()
    #time.full[itr]=(t2-t1)[3][[1]]
    #est.full[itr,]=d.full$est;se.full[itr,]=d.full$se
    #mse.full[itr]=c(t(est.full[itr,]-bt[1:d])%*%(est.full[itr,]-bt[1:d]))
    #cp.full[itr,]=abs(est.full[itr,]-bt[1:d])<se.full[itr,]*qnorm(1-alpha/2)
    #acp.full[itr]=mean(cp.full[itr,])
    #al.full[itr]=mean(se.full[itr,])*2*1.96
    
    
    #t1=proc.time()
    #d.full=full_data_estimator_linear(Y,Z,U,alpha)
    #t2=proc.time()
    #time.full[itr]=(t2-t1)[3][[1]]
    #est.full[itr,]=d.full$est;se.full[itr,]=d.full$se
    #mse.full[itr]=c(t(est.full[itr,]-bt[1:d])%*%(est.full[itr,]-bt[1:d]))
    #cp.full[itr,]=abs(est.full[itr,]-bt[1:d])<se.full[itr,]*qnorm(1-alpha/2)
    #acp.full[itr]=mean(cp.full[itr,])
    #al.full[itr]=mean(se.full[itr,])*2*1.96
    
  }
  
  result=list(MSE=c(mean(mse_multiStep),mean(mse.DVS),
                    mean(mse.duwA),mean(mse.duwL),mean(mse.dwA),mean(mse.dwL),mean(mse.dun),
                    mean(mse.full)),
              ACP=c(mean(acp_multiStep),mean(acp.DVS),
                    mean(acp.duwA),mean(acp.duwL),mean(acp.dwA),mean(acp.dwL),mean(acp.dun),
                    mean(acp.full)),
              al=c(mean(al_multiStep),mean(al.DVS),
                   mean(al.duwA),mean(al.duwL),mean(al.dwA),mean(al.dwL),mean(al.dun),
                   mean(al.full)),
              #iteration=c(mean(iteration_multiStep)),
              time=c(mean(time_multiStep),mean(time.DVS),
                     mean(time.duwA),mean(time.duwL), mean(time.dwA),mean(time.dwL),mean(time.dun),
                     mean(time.full)))
  
  
  return(result)
  
  
}

result1_p500_d10_Na<- decorrelate_subsampling_linear_p400(n=1e5,d=10,q=490,r2=1e3,r1=1e3,error="a",X_type = "N",itr.max=200,alpha=0.05)
cat("result1_d","\n")

result1_p500_d5_Na<- decorrelate_subsampling_linear_p400(n=1e5,d=5,q=495,r2=1e3,r1=1e3,error="a",X_type = "N",itr.max=200,alpha=0.05)
cat("result1_d","\n")


result1_p500_d10_Nb<- decorrelate_subsampling_linear_p400(n=1e5,d=10,q=490,r2=1e3,r1=1e3,error="b",X_type = "N",itr.max=200,alpha=0.05)
cat("result1_d","\n")

result1_p500_d5_Nb<- decorrelate_subsampling_linear_p400(n=1e5,d=5,q=495,r2=1e3,r1=1e3,error="b",X_type = "N",itr.max=200,alpha=0.05)
cat("result1_d","\n")



result1_p500_d10_ta<- decorrelate_subsampling_linear_p400(n=1e5,d=10,q=490,r2=1e3,r1=1e3,error="a",X_type = "t",itr.max=200,alpha=0.05)
cat("result1_d","\n")

result1_p500_d5_ta<- decorrelate_subsampling_linear_p400(n=1e5,d=5,q=495,r2=1e3,r1=1e3,error="a",X_type = "t",itr.max=200,alpha=0.05)
cat("result1_d","\n")



result1_p500_d10_tb<- decorrelate_subsampling_linear_p400(n=1e5,d=10,q=490,r2=1e3,r1=1e3,error="b",X_type = "t",itr.max=200,alpha=0.05)
cat("result1_d","\n")

result1_p500_d5_tb<- decorrelate_subsampling_linear_p400(n=1e5,d=5,q=495,r2=1e3,r1=1e3,error="b",X_type = "t",itr.max=200,alpha=0.05)
cat("result1_d","\n")

#n:full-sample size
#d:dimension of the parameter z;
#q :dimension of the parameter u, q+d=500;#p=d+q+1
#r2:informative subsample size r
#r1:initial subsample size r_p
#error & X_type:the four different settings
#itr.max: simulation  replications
#alpha:significance level of 0.05