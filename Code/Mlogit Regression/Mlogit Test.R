# covariate matrix
X = matrix(rnorm(75000), 5000, 15)

# coefficients for each choice
vCoef1 = rep(0, 15)
vCoef2 = c(2.5, 1.6, -2, -4, 2.8, 0, 0 , 0, 0, 0 , 0 , 0, 0, 0 ,0)
vCoef3 = c(1.3, -4, 3.6, 4, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

# vector of probabilities
sum_probs = exp(X%*%vCoef1) + exp(X%*%vCoef2) + exp(X%*%vCoef3)
vProb = cbind(exp(X%*%vCoef1)/sum_probs, exp(X%*%vCoef2)/sum_probs, exp(X%*%vCoef3)/sum_probs)

# multinomial draws
mChoices = t(apply(vProb, 1, rmultinom, n = 1, size = 1))
df = cbind.data.frame(y = as.factor(apply(mChoices, 1, function(x) which(x==1))), X)


X = model.matrix(~.,data=df[,2:16])
Y = model.matrix(~0+y,data=df)

beta = solve(t(X) %*% X) %*% t(X) %*% Y

logit.nll = function (beta, X, Y) {
  
  beta = matrix(beta,ncol=3)    
  P = as.matrix(rowSums(exp(X %*% beta))); #Sum_(h=1)^3 exp(X * Beta_(h))    
  Pr_1 = exp(X %*% beta[,2])/(1 + P); #P(y = 2 | X)
  Pr_2 = exp(X %*% beta[,3])/(1 + P); #P(y = 3 | X)    
  Pr_0 = 1/(1+P);#P(y = 1 | X)
  
  LL = (colSums(Y[,1] * log(Pr_0)) + colSums(Y[,2] * log(Pr_1)) + colSums(Y[,3] * log(Pr_2))) - 5*sum(abs(beta))#log-likelihood
  print(LL)
  return(-LL)
  
}

res = optim(beta, logit.nll, X = X, Y = Y, method = "BFGS")
res
