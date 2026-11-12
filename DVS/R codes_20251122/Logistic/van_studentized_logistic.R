van_studentized_logistic = function(X,Y,alpha=0.05,B = 500,d){
  n = nrow(X)
  p = ncol(X)
  REG <- cbind(X,1)
  tmp <- cv.glmnet(X,Y,family="binomial",alpha = 1)
  mod <- glmnet(X,Y,family="binomial",lambda=tmp$lambda.min,alpha=1)
  beta_hat <- as.vector(c(mod$beta[1:p],mod$a0))
  beta_hat[abs(beta_hat)<0.05]=0
  REG <- cbind(X,1)
  temp1 <- REG%*%beta_hat
  b_dev2=as.vector(exp(temp1)/(1+exp(temp1))^2)
  residual <- as.vector(Y-inv.logit(REG%*%beta_hat))*X
  X.new <- sqrt(b_dev2) * X
  W_hat <- matrix(0, nrow = d, ncol = p)
  for (j in 1:d){
    tmp_j <- cv.glmnet(X.new[,-j], X.new[,j], alpha = 1)
    gamma_j <- coef(tmp_j, s = "lambda.min")[-1]
    T_j_square <- sum((X.new[,j]-X.new[,-j]%*%gamma_j)^2)/n + tmp_j$lambda.min * sum(abs(gamma_j))
    W_hat[j, -j] <- -gamma_j / T_j_square
    W_hat[j, j] <- 1 / T_j_square
  }
  theta_tilde = beta_hat[1:d] + W_hat%*%apply(residual,2,mean)
  std=sqrt(diag(W_hat))
  c_n = sapply(1:B,function(xxx,n,std,W_hat,residual){
    ei = rnorm(n)
    boots = diag(std^(-1)) %*% W_hat %*% apply(ei*residual,2,mean)
    quan=max(abs(boots))
    return(quan)
  },n,std,W_hat,residual)
  quantiles = quantile(c_n,1-alpha)*std[1:d]
  return(list(est=theta_tilde,quantiles=quantiles))
}