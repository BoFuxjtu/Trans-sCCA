
This is the readme file for the R scripts of the paper “Unifiedly Efficient Inference on All-Dimensional Targets for Large-Scale GLMs” by Bo Fu, and Dandan Jiang. 
Please make sure the packages needed for each R script are installed. 

Install the following R packages:   "MASS", "nleqslv", "Matrix", "glmnet", "mvtnorm","mnormt", "boot", "SparseM", "quantreg","pqr","clime".
Please make sure the packages needed for each R script are installed. 

###Some steps for installing the package "pqr" in R 4.3.2
###install Rtools 4.3 (http://cran.r-project.org/bin/windows/Rtools)
###please make sure that the R tools 4.3 is installed in the same directory as the R 4.3.2 location
###install.packages("devtools")
###library(devtools)
###install_github("xliusufe/pqr")
###More details about the installation of the package "pqr" can be found at "https://github.com/xliusufe/pqr".


Explanation of files:

1. Linear: 
The "linear_lowtarget.R" is the main R code for simulations of linear regression with single coefficient inference in Section 5.1.
The "linear_hightarget.R" is the main R code for simulations of linear regression with high-dimensional simultaneous inference in Section 5.2.
The "DVS_linear.R" and "multiStep_linear.R" are used for our propsoed algorithms in the R functions "linear_lowtarget.R", while
 the "debiased_linear_nonstudentized.R" and "debiased_linear_nonstudentized.R" correspond to the non-studentized and studentized form of our methods implemented in
 "linear_hightarget.R" .
The "decorrelate_linear_algorithm.R" is used for the baseline algorithm in Section 5.1; The "full data estimator_linea.R" is used for the full-data-based decorrelated algorithm in 
 Section 5.1; the "zhang_debiased_nonstudentized.R" and "zhang_debiased_studentized.R" correspond to the non-studentized and studentized form of the baseline
 algorithm in Section 5.2.
The "train.csv" is the superconductivity dataset for real data analysis. The "real data analysis.R" is aplied to the inference for low parameters of interest while "real data analysis2.R" is aplied to simultaneous inference. 
The "Linear Figure 1.R" is to reproduce Figure 1 in Section 5.1.

2. Logistic: 
The "Logistic_lowtargets.R" is the main R code for simulations of logistic regression with single coefficient inference in Section A1.
The "Logistic_hightargets.R" is the main R code for simulations of logistic regression with high-dimensional simultaneous inference in Section A2.
The "DVS_Logistic.R" and "multiStep_Logistic.R" are used for our propsoed algorithms in the R functions "Logistic_lowtarget.R", while
 the "debiased_logistic_nonstudentized.R" and "debiased_logistic_nonstudentized.R" correspond to the non-studentized and studentized form of our methods implemented in
 "Logistic_hightarget.R" .
The "decorrelate_logistic_algorithm.R" is used for the baseline algorithm in Section A1; The "full data estimator_linea.R" is used for the full-data-based decorrelated algorithm in Section A1; the "van_nonstudentized logisitic.R" and "van_studentized logisitic.R" correspond to the non-studentized and studentized form of the baseline
 algorithm in Section A2.
The "Logistic Figure 1.R" is to reproduce Figure A1 in Section A1.
