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
library(flare)
library(clime)

library(data.table)
library(dplyr)
data <- fread("train.csv")
Y <- as.matrix(data$critical_temp)
#Z <- as.matrix(data[, c("range_ThermalConductivity", "wtd_std_ThermalConductivity")])
#U <- as.matrix(data %>% dplyr::select(-all_of(c("critical_temp","range_ThermalConductivity", "wtd_std_ThermalConductivity"))))

#Z <- as.matrix(data[, c("wtd_mean_ThermalConductivity", "wtd_std_Valence")])
#U <- as.matrix(data %>% dplyr::select(-all_of(c("critical_temp","wtd_mean_ThermalConductivity", "wtd_std_Valence"))))

#Z <- as.matrix(data[, c("wtd_std_ThermalConductivity", "wtd_entropy_Valence")])
#U <- as.matrix(data %>% dplyr::select(-all_of(c("critical_temp","wtd_std_ThermalConductivity", "wtd_entropy_Valence"))))

Z <- as.matrix(data[, c("range_ThermalConductivity","wtd_std_ThermalConductivity", "range_atomic_radius","wtd_gmean_ThermalConductivity", "std_ThermalConductivity","wtd_entropy_Valence",
                        "wtd_std_ElectronAffinity","wtd_entropy_atomic_mass","wtd_mean_Valence","wtd_gmean_ElectronAffinity",
                        "wtd_range_ElectronAffinity","wtd_mean_ThermalConductivity","wtd_gmean_Valence","std_atomic_mass","std_Density",
                        "wtd_entropy_ThermalConductivity","wtd_range_ThermalConductivity","wtd_mean_atomic_mass","wtd_std_atomic_mass","gmean_Density")])
U <- as.matrix(data %>% dplyr::select(-all_of(c("critical_temp","range_ThermalConductivity","wtd_std_ThermalConductivity", "range_atomic_radius","wtd_gmean_ThermalConductivity", "std_ThermalConductivity","wtd_entropy_Valence",
                                                "wtd_std_ElectronAffinity","wtd_entropy_atomic_mass","wtd_mean_Valence","wtd_gmean_ElectronAffinity",
                                                "wtd_range_ElectronAffinity","wtd_mean_ThermalConductivity","wtd_gmean_Valence","std_atomic_mass","std_Density",
                                                "wtd_entropy_ThermalConductivity","wtd_range_ThermalConductivity","wtd_mean_atomic_mass","wtd_std_atomic_mass","gmean_Density"))))

#Z <- as.matrix(data[, c("range_ThermalConductivity","wtd_std_ThermalConductivity", "range_atomic_radius","wtd_gmean_ThermalConductivity", "std_ThermalConductivity")])
#U <- as.matrix(data %>% dplyr::select(-all_of(c("critical_temp","range_ThermalConductivity","wtd_std_ThermalConductivity", "range_atomic_radius","wtd_gmean_ThermalConductivity", "std_ThermalConductivity"))))

#Z <- as.matrix(data[, c("wtd_mean_ThermalConductivity", "range_atomic_radius")])
#U <- as.matrix(data %>% dplyr::select(-all_of(c("critical_temp","wtd_mean_ThermalConductivity", "range_atomic_radius"))))


X <- cbind(Z,U)

d <- dim(Z)[2]

#X <- scale(X)
#Z <- X[,(1:d)]
#U <- X[,-(1:d)]
#Y <- Y - mean(Y)

itr.max <- 10

#parameter_interest <- which(names(data) %in% c("range_ThermalConductivity", "wtd_std_ThermalConductivity"))

se.debiased=est.debiased=matrix(NA,itr.max,d)
quan.debiased=time.debiased=rep(NA,itr.max)

se.debiased_stu=est.debiased_stu=quan.debiased_stu=matrix(NA,itr.max,d)
time.debiased_stu=rep(NA,itr.max)

se.zhang_nonstu=est.zhang_nonstu=matrix(NA,itr.max,d)
quan.zhang_nonstu=time.zhang_nonstu=rep(NA,itr.max)

se.zhang_stu=est.zhang_stu=quan.zhang_stu=matrix(NA,itr.max,d)
time.zhang_stu=rep(NA,itr.max)

interval.debiased=interval.debiased_stu=interval.zhang_nonstu=interval.zhang_stu=array(NA, dim = c(itr.max, 2, d))

realdata=function(X,Y,Z,U,r,itr.max,alpha,B){
  
  for(itr in 1:itr.max){
    
    cat(itr)
    
    t1=proc.time()
    d.debiased=debiased_linear_nonstudentized(Y,Z,U,r,alpha,B)
    t2=proc.time()
    time.debiased[itr]=(t2-t1)[3][[1]]
    est.debiased[itr,]=d.debiased$est
    quan.debiased[itr]=d.debiased$quantiles
    interval.debiased[itr, , ] <- rbind(est.debiased[itr,]-quan.debiased[itr], est.debiased[itr,]+quan.debiased[itr])
    
    t1=proc.time()
    d.debiased_stu=debiased_linear_studentized(Y,Z,U,r,alpha,B)
    t2=proc.time()
    time.debiased_stu[itr]=(t2-t1)[3][[1]]
    est.debiased_stu[itr,]=d.debiased_stu$est
    quan.debiased_stu[itr,]=d.debiased_stu$quantiles
    interval.debiased_stu[itr, , ] <- rbind(est.debiased_stu[itr,]-quan.debiased_stu[itr,], est.debiased_stu[itr,]+quan.debiased_stu[itr,])
    
    t1=proc.time()
    zhang_nonstu=Zhang_nonstudentized(X,Y,alpha,B,d)
    t2=proc.time()
    time.zhang_nonstu[itr]=(t2-t1)[3][[1]]
    est.zhang_nonstu[itr,]=zhang_nonstu$est
    quan.zhang_nonstu[itr]=zhang_nonstu$quantiles
    interval.zhang_nonstu[itr, , ] <- rbind(est.zhang_nonstu[itr,]-quan.zhang_nonstu[itr], est.zhang_nonstu[itr,]+quan.zhang_nonstu[itr])
    
    t1=proc.time()
    zhang_stu=Zhang_studentized(X,Y,alpha,B,d)
    t2=proc.time()
    time.zhang_stu[itr]=(t2-t1)[3][[1]]
    est.zhang_stu[itr,]=zhang_stu$est
    quan.zhang_stu[itr,]=zhang_stu$quantiles
    interval.zhang_stu[itr, , ] <- rbind(est.zhang_stu[itr,]-quan.zhang_stu[itr,], est.zhang_stu[itr,]+quan.zhang_stu[itr,])
  }
  
  result=list(time=c(mean(time.debiased),mean(time.debiased_stu),mean(time.zhang_nonstu),mean(time.zhang_stu)),
              interval.debiased=interval.debiased[1, , ],interval.debiased_stu=interval.debiased_stu[1, , ],
              interval.zhang_nonstu=interval.zhang_nonstu[1, , ],interval.zhang_stu=interval.zhang_stu[1, , ])
}

result <- realdata(X=X,Y=Y,Z=Z,U=U,r=5000,itr.max=itr.max,alpha=0.05,B=500)