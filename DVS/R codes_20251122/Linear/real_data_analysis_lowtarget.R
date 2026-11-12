rm(list=ls())
library("rstudioapi")
setwd(dirname(getActiveDocumentContext()$path))
source("decorrelate_linear_algorithm.R")
source("compare_linear_algorithm.R")
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

library(data.table)
library(dplyr)
data <- fread("train.csv")
Y <- as.matrix(data$critical_temp)
Z <- as.matrix(data[, c("wtd_mean_ThermalConductivity", "range_atomic_radius")])
U <- as.matrix(data %>% dplyr::select(-all_of(c("critical_temp","wtd_mean_ThermalConductivity", "range_atomic_radius"))))
X <- cbind(Z,U)
d <- dim(Z)[2]
itr.max <- 100
se.duwA=se.duwL=se.dwA=se.dwL=se.dun=se.DVS=se_multiStep=se.dfull=matrix(NA,itr.max,d)
al.duwA=al.duwL=al.dwA=al.dwL=al.dun=al.un=al_multiStep=al.DVS=al.dfull=rep(NA,itr.max)
time_multiStep=time.dwA=time.dwL=time.duwA=time.duwL=time.dfull=time.full=time.dun=rep(NA,itr.max)
time.DVS=al.DVS=rep(NA,itr.max)
se.DVS=se_multiStep=matrix(NA,itr.max,d)
interval_multiStep=interval.DVS=interval.duwA=interval.duwL=interval.dwA=interval.dwL=interval.dun=quan.DVS=interval.dfull=array(NA, dim = c(itr.max, 2, d))

est.duwA=est.duwL=est.dwA=est.dwL=est.dun=est_multiStep=est.DVS=est.dfull=matrix(NA,itr.max,d)
al_multiStep=rep(NA,itr.max)

realdata=function(X,Y,Z,U,r1,r2,itr.max,alpha){
  is_bad <- 1:itr.max
for(itr in 1:itr.max){
  #set.seed(35 + itr)
  cat(itr)
  t1=proc.time()
  d.multiStep=multiStep_linear(Y,Z,U,r1,alpha,100)
  t2=proc.time()
  time_multiStep[itr]=(t2-t1)[3][[1]]
  est_multiStep[itr,]=d.multiStep$est
  se_multiStep[itr,]=d.multiStep$se
  interval_multiStep[itr, , ]=rbind(est_multiStep[itr,]-se_multiStep[itr,]*1.96, est_multiStep[itr,]+se_multiStep[itr,]*1.96)
  al_multiStep[itr]=mean(se_multiStep[itr,])*2*1.96
  
  
  t1=proc.time()
  d.DVS=DVS_linear(Y,Z,U,r2=r2,r1=r1,alpha,500)
  t2=proc.time()
  time.DVS[itr]=(t2-t1)[3][[1]]
  est.DVS[itr,]=d.DVS$est
  quan.DVS[itr, , ]=d.DVS$quantiles
  interval.DVS[itr, , ] <- rbind(est.DVS[itr,]-quan.DVS[itr,2, ], est.DVS[itr,]-quan.DVS[itr,1, ])
  al.DVS[itr]=mean(quan.DVS[itr,2, ]-quan.DVS[itr,1, ])
  
  
  #Compared decorrelated method
  t1=proc.time()
  duwA=two_step_decorrelate_linear(Y,Z,U,r2,r1,optimality="A",weight="no")
  t2=proc.time()
  time.duwA[itr]=(t2-t1)[3][[1]]
  est.duwA[itr,]=duwA$est;se.duwA[itr,]=duwA$se
  interval.duwA[itr, , ]=rbind(est.duwA[itr,]-se.duwA[itr,]*1.96, est.duwA[itr,]+se.duwA[itr,]*1.96)
  al.duwA[itr]=mean(se.duwA[itr,])*2*1.96
  
  t1=proc.time()
  duwL=two_step_decorrelate_linear(Y,Z,U,r2,r1,optimality="L",weight="no")
  t2=proc.time()
  time.duwL[itr]=(t2-t1)[3][[1]]
  est.duwL[itr,]=duwL$est;se.duwL[itr,]=duwL$se
  interval.duwL[itr, , ]=rbind(est.duwL[itr,]-se.duwL[itr,]*1.96, est.duwL[itr,]+se.duwL[itr,]*1.96)
  al.duwL[itr]=mean(se.duwL[itr,])*2*1.96
  
  
  t1=proc.time()
  dwA=two_step_decorrelate_linear(Y,Z,U,r2,r1,optimality="A",weight="yes")
  t2=proc.time()
  time.dwA[itr]=(t2-t1)[3][[1]]
  est.dwA[itr,]=dwA$est;se.dwA[itr,]=dwA$se
  interval.dwA[itr, , ]=rbind(est.dwA[itr,]-se.dwA[itr,]*1.96, est.dwA[itr,]+se.dwA[itr,]*1.96)
  al.dwA[itr]=mean(se.dwA[itr,])*2*1.96
  
  
  t1=proc.time()
  dwL=two_step_decorrelate_linear(Y,Z,U,r2,r1,optimality="L",weight="yes")
  t2=proc.time()
  time.dwL[itr]=(t2-t1)[3][[1]]
  est.dwL[itr,]=dwL$est;se.dwL[itr,]=dwL$se
  interval.dwL[itr, , ]=rbind(est.dwL[itr,]-se.dwL[itr,]*1.96, est.dwL[itr,]+se.dwL[itr,]*1.96)
  al.dwL[itr]=mean(se.dwL[itr,])*2*1.96
  
  
  t1=proc.time()
  dun=two_step_decorrelate_linear(Y,Z,U,r2,r1,optimality="U",weight="no")
  t2=proc.time()
  time.dun[itr]=(t2-t1)[3][[1]]
  est.dun[itr,]=dun$est;se.dun[itr,]=dun$se
  interval.dun[itr, , ]=rbind(est.dun[itr,]-se.dun[itr,]*1.96, est.dun[itr,]+se.dun[itr,]*1.96)
  al.dun[itr]=mean(se.dun[itr,])*2*1.96
  
  t1=proc.time()
  dfull=full_data_estimator_linear(Y,Z,U,alpha)
  t2=proc.time()
  time.dfull[itr]=(t2-t1)[3][[1]]
  est.dfull[itr,]=dfull$est;se.dfull[itr,]=dfull$se
  interval.dfull[itr, , ]=rbind(est.dfull[itr,]-se.dfull[itr,]*1.96, est.dfull[itr,]+se.dfull[itr,]*1.96)
  al.dfull[itr]=mean(se.dfull[itr,])*2*1.96
    
  is_bad[itr]<- (max(abs(est.DVS[itr,]-est.dfull[itr,])) >1)|(max(abs(est.duwA[itr,]-est.dfull[itr,])) >2)|
    (max(abs(est.duwL[itr,]-est.dfull[itr,])) >2)|(max(abs(est.dwA[itr,]-est.dfull[itr,])) >2)|
    (max(abs(est.dwL[itr,]-est.dfull[itr,])) >2)|(max(abs(est.dun[itr,]-est.dfull[itr,])) >2)
  
}
  
  

SE_multiStep=apply(est_multiStep[!is_bad,], 2, function(x) sd(x) / sqrt(length(x)))
SE.DVS=apply(est.DVS[!is_bad,], 2, function(x) sd(x) / sqrt(length(x)))
SE.duwA=apply(est.duwA[!is_bad,], 2, function(x) sd(x) / sqrt(length(x)))
SE.duwL=apply(est.duwL[!is_bad,], 2, function(x) sd(x) / sqrt(length(x)))
SE.dwA=apply(est.dwA[!is_bad,], 2, function(x) sd(x) / sqrt(length(x)))
SE.dwL=apply(est.dwL[!is_bad,], 2, function(x) sd(x) / sqrt(length(x)))
SE.dun=apply(est.dun[!is_bad,], 2, function(x) sd(x) / sqrt(length(x)))
SE.dfull=apply(est.dfull[!is_bad,], 2, function(x) sd(x) / sqrt(length(x)))

Bias_multiStep=apply(est_multiStep[!is_bad,],2,mean)-apply(est.dfull[!is_bad,],2,mean)
Bias.DVS=apply(est.DVS[!is_bad,],2,mean)-apply(est.dfull[!is_bad,],2,mean)
Bias.duwA=apply(est.duwA[!is_bad,],2,mean)-apply(est.dfull[!is_bad,],2,mean)
Bias.duwL=apply(est.duwL[!is_bad,],2,mean)-apply(est.dfull[!is_bad,],2,mean)
Bias.dwA=apply(est.dwA[!is_bad,],2,mean)-apply(est.dfull[!is_bad,],2,mean)
Bias.dwL=apply(est.dwL[!is_bad,],2,mean)-apply(est.dfull[!is_bad,],2,mean)
Bias.dun=apply(est.dun[!is_bad,],2,mean)-apply(est.dfull[!is_bad,],2,mean)
#Bias.dfull=abs(apply(est.dfull,2,mean)-est.full)

#browser()

result=list(x1=c(est_multiStep[1,1],est.DVS[1,1],
                  est.duwA[1,1],est.duwL[1,1],est.dwA[1,1],est.dwL[1,1],est.dun[1,1],est.dfull[1,1]),
            x2=c(est_multiStep[1,2],est.DVS[1,2],
                 est.duwA[1,2],est.duwL[1,2],est.dwA[1,2],est.dwL[1,2],est.dun[1,2],est.dfull[1,2]),
            Bias1=c(Bias_multiStep[1],Bias.DVS[1],Bias.duwA[1],Bias.duwL[1],Bias.dwA[1],Bias.dwL[1],Bias.dun[1]),
            Bias2=c(Bias_multiStep[2],Bias.DVS[2],Bias.duwA[2],Bias.duwL[2],Bias.dwA[2],Bias.dwL[2],Bias.dun[2]),
            SE1=c(SE_multiStep[1],SE.DVS[1],
                 SE.duwA[1],SE.duwL[1],SE.dwA[1],SE.dwL[1],SE.dun[1],SE.dfull[1]),
            SE2=c(SE_multiStep[2],SE.DVS[2],
                  SE.duwA[2],SE.duwL[2],SE.dwA[2],SE.dwL[2],SE.dun[2],SE.dfull[2]),
            time=c(mean(time_multiStep[!is_bad]),mean(time.DVS[!is_bad]),
                   mean(time.duwA[!is_bad]),mean(time.duwL[!is_bad]), mean(time.dwA[!is_bad]),mean(time.dwL[!is_bad]),mean(time.dun[!is_bad]),mean(time.dfull[!is_bad])),
            interval_multiStep=interval_multiStep[1, , ],interval.DVS=interval.DVS[1, , ],
            interval.duwA=interval.duwA[1, , ],interval.duwL=interval.duwL[1, , ],
            interval.dwA=interval.dwA[1, , ],interval.dwL=interval.dwL[1, , ],
            interval.dun=interval.dun[1, , ],interval.dfull=interval.dfull[1, , ])
}

result_500 <- realdata(X=X,Y=Y,Z=Z,U=U,r1=500,r2=500,itr.max=itr.max,alpha=0.05)

result_1000 <- realdata(X=X,Y=Y,Z=Z,U=U,r1=1000,r2=1000,itr.max=itr.max,alpha=0.05)

result_2000 <- realdata(X=X,Y=Y,Z=Z,U=U,r1=2000,r2=2000,itr.max=itr.max,alpha=0.05)

result_5000 <- realdata(X=X,Y=Y,Z=Z,U=U,r1=5000,r2=5000,itr.max=itr.max,alpha=0.05)