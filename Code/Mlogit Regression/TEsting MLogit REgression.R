iris

rDiscreta<-function(p){
  u<-runif(1)
  P<-cumsum(p)
  val<-sum(P<u)+1
  return(val)
}

N=500 #Tamanho da amostra Binomial
T=1200 #Cumprimento da cadeia simulada
K=3   #Numero de estados ocultos
D=10   #Quantidade de Covariaveis
P0=rep(1/K,K) #Inicializamos vetor de probabilidades inciais para o HMM
theta=c(0.10,0.5,0.9) # vetor com a probabilidade de sucesso das 2 distribuiçoes Binomiais
Nt=rep(N,T) # número de ensaios de Bernoulli associado a dada uma das T variáveis Binomiais. Cada coloquei tudo igual mas eles podem diferentes.
Betas=array(0, dim=c(K,D,K)) # valores de beta utilizados na geração dos valores (consideranda intercepto e duas covariáveis)

S<-NULL #Inicializamos a sequência de estados não observaveis
Y<-NULL #Inicializamos a sequência de valores observaveis
X<-NULL #Inicializamos o vetor de covariaveis
S_treino<-NULL #Inicializamos a sequência S de treinamento
X <- matrix(rnorm(D*T), nrow = T, ncol = D)
X[,1] = 1



S[1]<-rDiscreta(P0) #O valor para o primeiro estado oculto
Y[1]<-rbinom(1,Nt[1],theta[S[1]])# O valor para o primeiro valor observavel
#Geramos a Sequencia S de cumprimento T e um valor observavel Y para cada S.

for (t in 2:T){
  prob<-NULL
  for (i in 1:K) prob[i]<-exp(X[t,]%*%matrix(Betas[i,,S[t-1]],ncol=1))
  prob<-prob/sum(prob)
  S[t]<-rDiscreta(prob)
  Y[t]<-rbinom(1,Nt[t],theta[S[t]])
}



data <- cbind(S,X)
data<-as.data.frame(data)
X = model.matrix(~.,data=data[,2:11])
S = model.matrix(~0+S,data=data)





beta = solve(t(X) %*% X) %*% t(X) %*% Y

logit.nll = function (beta, X, Y) {
  
  beta = matrix(beta,ncol=3)    
  P = as.matrix(rowSums(exp(X %*% beta))); #Sum_(h=1)^3 exp(X * Beta_(h))    
  Pr_1 = exp(X %*% beta[,2])/(1 + P); #P(y = 2 | X)
  Pr_2 = exp(X %*% beta[,3])/(1 + P); #P(y = 3 | X)    
  Pr_0 = 1/(1+P);#P(y = 1 | X)
  
  LL = (colSums(Y[,1] * log(Pr_0)) + colSums(Y[,2] * log(Pr_1)) + colSums(Y[,3] * log(Pr_2))) #log-likelihood
  print(LL)
  return(-LL)
  
}

res = optim(beta, logit.nll, X = X, Y = Y, method = "BFGS")
res
