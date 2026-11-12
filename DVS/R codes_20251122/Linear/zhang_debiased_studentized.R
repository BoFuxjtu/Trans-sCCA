Zhang_studentized = function(X,Y,alpha=0.05,B = 500,d){
  n = nrow(X)
  p = ncol(X)
  tmp <- cv.glmnet(X, Y, alpha = 1)
  mod <- glmnet(X,Y,lambda=tmp$lambda.min,alpha=1)
  beta_hat <- as.vector(c(mod$beta[1:p],mod$a0))
  beta_hat[abs(beta_hat)<0.05]=0
  REG <- cbind(X,1)
  residual <- Y-REG%*%beta_hat
  sigma_hat = sum((residual)^2)/(n-sum(beta_hat[-1]!=0))
  W_hat <- matrix(0, nrow = d, ncol = p)
  for (j in 1:d){
    tmp_j <- cv.glmnet(X[,-j], X[,j], alpha = 1)
    gamma_j <- coef(tmp_j, s = "lambda.min")[-1]
    T_j_square <- sum((X[,j]-X[,-j]%*%gamma_j)^2)/n + tmp_j$lambda.min * sum(abs(gamma_j))
    W_hat[j, -j] <- -gamma_j / T_j_square
    W_hat[j, j] <- 1 / T_j_square
  }
  theta_tilde = beta_hat[1:d] + W_hat%*%t(X)%*%(residual)/n
  std = sqrt(diag(W_hat %*% (t(X)%*%X/n) %*% t(W_hat) * sigma_hat))
  c_n = sapply(1:B,function(xxx,n,std,W_hat,X,sigma_hat){
    ei = rnorm(n)
    boots = diag(std^(-1)) %*% W_hat %*% t(X) %*% ei*sqrt(sigma_hat)/n
    quan=max(abs(boots))
    return(quan)
  },n,std,W_hat,X,sigma_hat)
  Q = quantile(c_n,1-alpha)
  quantiles=Q*std
  return(list(est=theta_tilde,quantiles=quantiles))
}