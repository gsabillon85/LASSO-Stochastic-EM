
#####Função para gerar valores uniformes discretos
rDiscreta<-function(p){
  u<-runif(1)
  P<-cumsum(p)
  val<-sum(P<u)+1
  return(val)}
#####


zero_threshold = 0.05
R <- 100 # Numero de Replicas
T=1000 #Cumprimento da cadeia simulada
K=3   #Numero de estados ocultos
D=8   #Quantidade de Covariaveis
tol<-0.0000001 #Nivel de tolerancia que estabelecemos como criterio de parada do EM Est
tolval=NULL
tolval[1]=1
optim_algo = "BFGS" #Algorithm to use in the optimization process
n_max_iter_EM = 50
Tempo <- NULL

set.seed(2)
seeds <-sample(110000000,R) # Seed number para conseguer os mesmos valores simulados

P0=rep(1/K,K) #Inicializamos vetor de probabilidades inciais para o HMM
Betas=array(0, dim=c(K,D,K)) # valores de beta utilizados na geração dos valores (consideranda intercepto e duas covariáveis)
Real=NULL


#Outros 
Betas[2,1,1]=-1
Betas[2,3,1]=1.5
Betas[2,8,1]=2.3

Betas[2,1,2]=0.6
Betas[2,4,2]=1.2
Betas[2,7,2]=-2.4

Betas[2,1,3]=-0.4
Betas[2,5,3]=-1.4
Betas[2,6,3]=3.4

# Betas transção 3
Betas[3,1,1]=-1
Betas[3,2,1]=0
Betas[3,4,1]=0

Betas[3,1,2]=-0.6
Betas[3,6,2]=1.2
Betas[3,3,2]=5.4

Betas[3,1,3]=0.5
Betas[3,5,3]=2
Betas[3,3,3]=-1.5

# -- Caractersticas Especificas para a distribuição das VA observaveis
mu = c(50,70,90) # vetor com media para as duas Normais
sigma = c(,1.5,1) #Vetor com os desvios padrões para as duas normais
RealTransCount <- array(rep(0,K*K*R), dim=c(K,K,R))

for (p in 1:R){
  set.seed(seeds[p])
  #   INICIO DE SIMULAÇÃO DOS DADOS ##
  #########################################
  print("Simulando Dados...")
  S<-NULL #Inicializamos a sequência de estados não observaveis
  Y<-NULL #Inicializamos a sequência de valores observaveis
  X<-NULL #Inicializamos o vetor de covariaveis
  S_treino<-NULL #Inicializamos a sequência S de treinamento
  X <- matrix(rnorm(D*T), nrow = T, ncol = D)
  X[,1] = 1
  
  S[1]<-rDiscreta(P0) #O valor para o primeiro estado oculto
  Y[1]<-rnorm(1,mu[S[1]],sigma[S[1]])# O valor para o primeiro valor observavel
  #Geramos a Sequencia S de cumprimento T e um valor observavel Y para cada S.
  for (t in 2:T){
    prob<-NULL
    for (i in 1:K) prob[i]<-exp(X[t,]%*%matrix(Betas[i,,S[t-1]],ncol=1))
    prob<-prob/sum(prob)
    S[t]<-rDiscreta(prob)
    Y[t]<-rnorm(1,mu[S[t]],sigma[S[t]])
  }
  
  for (i in 2:T) {
    for (j in 1:K) {
      for (k in 1:K) {
        if (S[i]==j && S[i-1]==k)
          RealTransCount[k,j,p]=RealTransCount[k,j,p]+1
      }
    }
  }
}

Final_TransCount <- matrix(nrow = K, ncol = K)
Final_TransCount_sd <- matrix(nrow = K, ncol = K)
for (j in 1:K) {
  for (k in 1:K) {
    Final_TransCount[k,j] <- round(mean(RealTransCount[k,j,])/length(S)*100,2)
    Final_TransCount_sd[k,j] <- round(sd(RealTransCount[k,j,])/length(S)*100,2)
  }
}
table(S)
Final_TransCount
Final_TransCount_sd