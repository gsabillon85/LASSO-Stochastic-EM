library('label.switching')
#library("r2excel")
#library("xlsx")
library(e1071)

#COnseguimos Zerar para 3 covariaveis e 3 estados ocultos, porme o codigo nao esta muito generalizado


options(digits=8)
options(scipen=999)
num<-sample(110000000,1) # Seed number para conseguer os mesmos valores simulados
#num = 70431463
set.seed(num)
num
tempo_inicial<-proc.time()#Calcularemos quanto demoro todo o proceso

N=500 #Tamanho da amostra Binomial
T=1200 #Cumprimento da cadeia simulada
K=3   #Numero de estados ocultos
D=10   #Quantidade de Covariaveis
tol<-0.0000001 #Nivel de tolerancia que estabelecemos como criterio de parada do EM Est
tolval=NULL
tolval[1]=1

lambda = 0.5

path = "/home/gustavo/Documentos/Programas PosMestrado/NHMM PosMestrado/K3/EM Estocastico/K3/Sem_Replicas/Output"
Output_Number = 4
path_completo = paste(path, toString(Output_Number),sep = "")
ifelse(!dir.exists(path_completo), dir.create(path_completo), FALSE)
VerProx<-NULL
VerAct<-NULL
arqname<-c(path_completo,"/Output_Estocastico_",toString(T),".xlsx")
nome_arquivo<-paste(arqname,collapse = "") 
nomefolha<-paste("Saida T=",toString(T))
header<-paste("Resultados para K=",toString(K), "e T=",toString(T), collapse = "")

####### Simulação #######
TransCount <- matrix(data = c(rep(0,K^2)), nrow = K, ncol = K)
P0=rep(1/K,K) #Inicializamos vetor de probabilidades inciais para o HMM
theta=c(0.10,0.5,0.9) # vetor com a probabilidade de sucesso das 2 distribuiçoes Binomiais
Nt=rep(N,T) # número de ensaios de Bernoulli associado a dada uma das T variáveis Binomiais. Cada coloquei tudo igual mas eles podem diferentes.
theta_hat = NULL #Variavel para estimar os thetas em cada iteração do EM Estocastico
BetaArray = array(0, dim=c(K,D,K))
Betas=array(0, dim=c(K,D,K)) # valores de beta utilizados na geração dos valores (consideranda intercepto e duas covariáveis)

#Fazemos o valor iniciais dos Betas especificos igual a 0, para ter a função
#de Ligação mlogit



Betas[1,1,1]=0
Betas[1,2,1]=0
Betas[1,3,1]=0
Betas[2,1,1]=0
Betas[2,2,1]=-1.9
Betas[2,3,1]=0
Betas[3,1,1]=-2.3
Betas[3,2,1]=0
Betas[3,3,1]=0

Betas[1,1,2]=0
Betas[1,2,2]=0
Betas[1,3,2]=0
Betas[2,1,2]=2.4
Betas[2,2,2]=0
Betas[2,3,2]=5
Betas[3,1,2]=0
Betas[3,2,2]=1
Betas[3,3,2]=2.2

Betas[1,1,3]=0
Betas[1,2,3]=0
Betas[1,3,3]=0
Betas[2,1,3]=0
Betas[2,2,3]=0
Betas[2,3,3]=-1.3
Betas[3,1,3]=0
Betas[3,2,3]=1.9
Betas[3,3,3]=1.7

#####Função para gerar valores uniformes discretos
rDiscreta<-function(p){
  u<-runif(1)
  P<-cumsum(p)
  val<-sum(P<u)+1
  return(val)}
#####

#######   Escrevemos as funções que serão o objetivo da optimização   ######
# Com o temos um array de Betas, utilizaremos tres funções para achar os valores otimos
# Uma para a matriz Betas[,,1] uma para a matriz Betas[,,2] e uma para 
# a matriz Betas[,,3]

# FSM1 <-function(params){#função a maximizar para achar os Betas_1
#   resp <- sum(1 - log(1 + exp(params[1]*Xtemp11[,1] + params[2]*Xtemp11[,2] + params[3]*Xtemp11[,3])+ exp(params[4]*Xtemp11[,1]+params[5]*Xtemp11[,2] + params[6]*Xtemp11[,3]))) + sum(params[1]*Xtemp12[,1] + params[2]*Xtemp12[,2] + params[3]*Xtemp12[,3] - log(1 + exp(params[1]*Xtemp12[,1] + params[2]*Xtemp12[,2] + params[3]*Xtemp12[,3]) + exp(params[4]*Xtemp12[,1]+params[5]*Xtemp12[,2] + params[6]*Xtemp12[,3] ))) + sum(params[4]*Xtemp13[,1] + params[5]*Xtemp13[,2] + params[6]*Xtemp13[,3] - log(1 + exp(params[1]*Xtemp13[,1] + params[2]*Xtemp13[,2] + params[3]*Xtemp13[,3]) + exp(params[4]*Xtemp13[,1]+params[5]*Xtemp13[,2] + params[6]*Xtemp13[,3]))) - 20*sum(abs(params[2:D])) 
# }

FSM1 <-function(params){#função a maximizar para achar os Betas_1
  resp <- (sum(1 - log(1 + exp(Xtemp11%*%params[1:D])+ exp(Xtemp11%*%params[(D+1):(2*D)]))) + sum((Xtemp12%*%params[1:D]) - log( 1 + exp(Xtemp12%*%params[1:D])+ exp(Xtemp12%*%params[(D+1):(2*D)]) )) + sum((Xtemp13%*%params[(D+1):(2*D)]) - log( 1 + exp(Xtemp13%*%params[1:D])+ exp(Xtemp13%*%params[(D+1):(2*D)]) ))) + 0 + ( sum(1 - log(1 + exp(Xtemp21%*%params[(2*D+1):(2*D+D)])+ exp(Xtemp21%*%params[((2*D)+(D+1)):((2*D)+2*D)]))) + sum((Xtemp22%*%params[((2*D)+1):((2*D)+D)]) - log( 1 + exp(Xtemp22%*%params[((2*D)+1):((2*D)+D)])+ exp(Xtemp22%*%params[((2*D)+(D+1)):((2*D)+2*D)]) )) + sum((Xtemp23%*%params[((2*D)+D+1):((2*D)+2*D)]) - log( 1 + exp(Xtemp23%*%params[((2*D)+1):((2*D)+D)])+ exp(Xtemp23%*%params[((2*D)+D+1):((2*D)+2*D)]))))  + 0 +  (   sum(1 - log(1 + exp(Xtemp31%*%params[((4*D)+1):((4*D)+D)])+ exp(Xtemp31%*%params[((4*D)+D+1):((4*D)+2*D)]))) + sum((Xtemp32%*%params[((4*D)+1):((4*D)+D)]) - log( 1 + exp(Xtemp32%*%params[((4*D)+1):((4*D)+D)])+ exp(Xtemp32%*%params[((4*D)+D+1):((4*D)+2*D)]) )) + sum((Xtemp33%*%params[((4*D)+D+1):((4*D)+2*D)]) - log( 1 + exp(Xtemp33%*%params[((4*D)+1):((4*D)+D)])+ exp(Xtemp33%*%params[((4*D)+D+1):((4*D)+2*D)]) )) ) - lambda*(sum(abs(params[2:6*D])))
  print(resp)
}

################################
###        SIMULAÇÃO         ###
################################

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

#########################################
#####   Procedimento de Estimação   #####
#########################################

#######Primeiro geramos uma sequência não observavel de treinamento######
P_Treino=rep(1/K,K) #Vetor de probabilidade utilizadas para gerar a sequência de treino
S_treino<-NULL # Inicializamos a sequência oculta de treinamento
init1 = c(rnorm(D*K*(K-1)))#Valores iniciais para os Betas_1

#Geramos uma sequência de treinamento
for (i in 1:T) {
  S_treino[i] = rDiscreta(P_Treino)
}

## Escrevemos a função para recalcular a matriz de transição em cada iteração
## do algoritmo EM Estocástico.
Mat_trans <-function(covar){
  B = matrix(nrow=K, ncol=K)
  for (j in 1:K) {
    for (i in 1:K){
      numerator = exp(covar%*%BetaArray[i,,j])
      denom = 0
      for (l in 1:K){
        denom = denom + exp(covar%*%BetaArray[l,,j])
      }
      B[i,j] = numerator/denom
    }  
  }
  return(B)
}

val=1
tolval[1]=1
D
#Agora executamos o Algoritmo EM Estocástico
while ( (abs(tolval[val]))>tol ){
  #VeroSimActual=VeroSimProxima
  #Aqui devemos calcular a diferença entre a L.V. em na iteração atual e na anterior  
  acVS1=0
  acVS2=0
  acVS3=0
  VeroSimActual=0
  
  for (k in 1:K) {id = S_treino == k
  theta_hat[k] = sum(id*Y)/sum(id*Nt)
  }
  
  for (i in 1:T) {#Calculo do primeiro segmento da LL
     acVS1 = acVS1 + Y[i]*log(theta_hat[S_treino[i]])
  }
  for (i in 1:T) {#Calculo do segundo segmento da LL
    acVS2 = acVS2 + (Nt[i]-Y[i])*(log(1-theta_hat[S_treino[i]]))
  }
  temp=NULL
  for (i in 2:T) {#Calculo do terceiro segmento da LL
    for (g in 1:K) {
      temp[g]<-exp(X[i,]%*%matrix(BetaArray[g,,S_treino[t-1]],ncol=1))
  }
    acVS3 = acVS3 + (X[i,]%*%matrix(BetaArray[S_treino[i],,S_treino[i-1]]) - log(sum(temp), base = exp(1)))
  }
  VeroSimActual <- sum(log(choose(Nt,Y), base = exp(1))) + log(P0[S_treino[1]]) + acVS1 + acVS2 + acVS3 #calculo da LogVerosim
  
  #Calculamos a sequência S_treino utilizando os Betas
  #Atualizados na iteração passada e os valores observados Y
  S_treino[1]=which.max(dbinom(Y[1], Nt[1], theta_hat))
  for (i in 2:T) {
    A_hat_t = Mat_trans(X[i,])
    prob<-(A_hat_t[S_treino[i], ]*dbinom(Y[i], Nt[i], theta_hat))/sum(A_hat_t[S_treino[i], ]*dbinom(Y[i], Nt[i], theta_hat))
    S_treino[i]=rDiscreta(prob)
    #S_treino[i]=which.max(A_hat_t[S_treino[i-1], ]*dbinom(Y[i], Nt[i], theta_hat))
  }
  
  #####################################
  #Este segmento de codigo testa se aconteceram todas as transições possiveis
  #No caso que elas não tinham acontecido, as que
  #não aconteceram são forçadas a acontecer
  TransCount <- matrix(data = c(0,0,0,0,0,0,0,0,0), nrow = K, ncol = K)
  for (i in 2:T) {
    for (j in 1:K) {
      for (k in 1:K) {
        if (S_treino[i]==j && S_treino[i-1]==k)
          TransCount[k,j]=TransCount[k,j]+1
      }
    }
  }
  
  for (j in 1:K) {
    for (k in 1:K) {
      if (TransCount[k,j]==0){
        positions = sample(2:T, 4)
        for (d in 1:4) {
          S_treino[positions[d]]=j
          S_treino[positions[d]-1]=k
        }
      }
    }
  }
  
  #### Aqui inicia a filtragem dos dados para cada iteração
  Xtemp11<-NULL
  Xtemp12<-NULL
  Xtemp13<-NULL
  Xtemp21<-NULL
  Xtemp22<-NULL
  Xtemp23<-NULL
  Xtemp31<-NULL
  Xtemp32<-NULL
  Xtemp33<-NULL
  
  for (t in 2:T) {
    #filtros indo para o Estado # 1
    if(S_treino[t]%in%1 && S_treino[t-1]%in%1)
      Xtemp11<-rbind(Xtemp11, X[t,])
    
    if(S_treino[t]%in%1 && S_treino[t-1]%in%2)
      Xtemp21<-rbind(Xtemp21, X[t,])
    
    if(S_treino[t]%in%1 && S_treino[t-1]%in%3)
      Xtemp31<-rbind(Xtemp31, X[t,])
    
    #Filtros indo para o Estado # 2
    if(S_treino[t]%in%2 && S_treino[t-1]%in%1)
      Xtemp12<-rbind(Xtemp12, X[t,])
    
    if(S_treino[t]%in%2 && S_treino[t-1]%in%2)
      Xtemp22<-rbind(Xtemp22, X[t,])
    
    if(S_treino[t]%in%2 && S_treino[t-1]%in%3)
      Xtemp32<-rbind(Xtemp32, X[t,])
    
    #Filtros indo para o Estado # 3
    if(S_treino[t]%in%3 && S_treino[t-1]%in%1)
      Xtemp13<-rbind(Xtemp13, X[t,])
    
    if(S_treino[t]%in%3 && S_treino[t-1]%in%2)
      Xtemp23<-rbind(Xtemp23, X[t,])
    
    if(S_treino[t]%in%3 && S_treino[t-1]%in%3)
      Xtemp33<-rbind(Xtemp33, X[t,])
  }

  
  ##O ajuste para estimar os parâmetros de transição é
  ##feito aqui usando a função optim e os valores das
  #covariaveis filtradas
  fit1 <- optim(par = init1, fn = FSM1, control = list(fnscale=-1), method = "BFGS", hessian = FALSE)

  # Aqui atribuimos os valores estimados dos parâmetros de 
  # transição a um array que sera utilizado para recalcular 
  # a sequência S_treino na seguinte iteração do EM Est. 
  # Em outras palavras, aqui acontece a ATUALIZAÇÃO dos parâmetros de transição.
  
  for (i in 1:K){
    for (d in 1:D){
      if (i == 1){
        BetaArray[i,d,1]=0
      } else if (i == 2){
        BetaArray[i,d,1]=fit1$par[d]
      } else if (i == 3){
        BetaArray[i,d,1]=fit1$par[D+d]
      }
      
    }
  }
  
  for (i in 1:K){
    for (d in 1:D){
      if (i == 1){
        BetaArray[i,d,2]=0
      } else if (i == 2){
        BetaArray[i,d,2]=fit1$par[(2*D + d)]
      } else if (i == 3){
        BetaArray[i,d,2]=fit1$par[((2*D) + D+d)]
      }
      
    }
  }
  
  for (i in 1:K){
    for (d in 1:D){
      if (i == 1){
        BetaArray[i,d,3]=0
      } else if (i == 2){
        BetaArray[i,d,3]=fit1$par[(4*D + d)]
      } else if (i == 3){
        BetaArray[i,d,3]=fit1$par[(4*D + D+d)]
      }
      
    }
  }
  
  
  # BetaArray[1,1,1]=0
  # BetaArray[1,2,1]=0
  # BetaArray[1,2,1]=0
  # BetaArray[2,1,1]=fit1$par[1]
  # BetaArray[2,2,1]=fit1$par[2]
  # BetaArray[2,3,1]=fit1$par[3]
  # BetaArray[3,1,1]=fit1$par[4]
  # BetaArray[3,2,1]=fit1$par[5]
  # BetaArray[3,3,1]=fit1$par[6]
  # 
  # BetaArray[1,1,2]=0
  # BetaArray[1,2,2]=0
  # BetaArray[1,3,2]=0
  # BetaArray[2,1,2]=fit2$par[1]
  # BetaArray[2,2,2]=fit2$par[2]
  # BetaArray[2,3,2]=fit2$par[3]
  # BetaArray[3,1,2]=fit2$par[4]
  # BetaArray[3,2,2]=fit2$par[5]
  # BetaArray[3,3,2]=fit2$par[6]
  # 
  # Betas
  # 
  # BetaArray[1,1,3]=0
  # BetaArray[1,2,3]=0
  # BetaArray[1,3,3]=0
  # BetaArray[2,1,3]=fit3$par[1]
  # BetaArray[2,2,3]=fit3$par[2]
  # BetaArray[2,3,3]=fit3$par[3]
  # BetaArray[3,1,3]=fit3$par[4]
  # BetaArray[3,2,3]=fit3$par[5]
  # BetaArray[3,3,3]=fit3$par[6]
  
  ac2VS1=0
  ac2VS2=0
  ac2VS3=0
  VeroSimProxima=0
  
  for (i in 1:T) {#Calculo do primeiro segmento da LL
    ac2VS1 = ac2VS1 + Y[i]*log(theta_hat[S_treino[i]])
  }
  for (i in 1:T) {#Calculo do segundo segmento da LL
    ac2VS2 = ac2VS2 + (Nt[i]-Y[i])*(log(1-theta_hat[S_treino[i]]))
  }
  temp=NULL
  for (i in 2:T) {#Calculo do terceiro segmento da LL
    for (g in 1:K) {
      temp[g]<-exp(X[i,]%*%matrix(BetaArray[g,,S_treino[t-1]],ncol=1))
    }
    ac2VS3 = ac2VS3 + (X[i,]%*%matrix(BetaArray[S_treino[i],,S_treino[i-1]]) - log(sum(temp), base = exp(1)))
  }
  VeroSimProxima <- sum(log(choose(Nt,Y), base = exp(1))) + log(P0[S_treino[1]]) + ac2VS1 + ac2VS2 + ac2VS3 #calculo da LogVerosim
  
  val=val+1
  VerAct[val]<-VeroSimActual
  VerProx[val]<-VeroSimProxima
  tolval[val]<-VeroSimProxima - VeroSimActual
  print(tolval[val])
  
}#######Fim da primeira rodada do EM Estocastico#######

#Criar algumas matrizes para fazer calculos e manipular a saida MCMC
#nestas matrizes, as estimativas serão reordenadas usando o metodo ECR
mat_thetar<-matrix(nrow = 1, ncol = K)
reorder_S<-matrix(nrow = 1, ncol = T)
mat_S<-matrix(nrow = 1, ncol = T)
mat_S[1,]<-S_treino
zpvt_S = S #Como pivot para o metodo ECR usamos o S original
perms_S = ecr(zpivot = zpvt_S, z = mat_S, K = 3)# aplicamos o metodo ECR que retornara as permutações das dos estados ocultos que devem ser utilizadas para reordenar a saida do algoritmo bayesiano

# Reordenamos a saido do algoritmo EMEst usando as 
# permutações fornecidas pelo ECR para K=3 
# só rerotulamos a Sequência S_treino, e reordenamos os Thetas
# Os Betas serão estimados usando a sequência S_Treino rerotulada
# e os Thetas, na segunda etapa do EMEst

for (i in 1:1) {
  for (j in 1:T) {
    if(S_treino[j]!=S[j] && ((perms_S$permutations[i,1]==2 && perms_S$permutations[i,2]==3 && perms_S$permutations[i,3]==1) | (perms_S$permutations[i,1]==3 && perms_S$permutations[i,2]==1 && perms_S$permutations[i,3]==2))){
      S_treino[j]=perms_S$permutations[i,perms_S$permutations[i,S_treino[j]]]
    }
    
    else {
      S_treino[j]=perms_S$permutations[i,perms_S$permutations[i,S[j]]]
    }
  }
  theta_hat<-theta_hat[perms_S$permutations[i,]]
}

# repetimos o EM Estocastico, porque para K=3
# Entraremos com a sequência S_treino estimada na rodada 
# anterior E ja rotulada corretamente usando o ECR Para resolver
# o problema dos Parametros Fantasmas e a troca de rotulos

VeroSimProxima=1
VeroSimActual=0
val=1
tolval[1]=1
while ( abs(tolval[val])>tol ) {
  #Aqui devemos calcular a diferença entre a L.V. em na iteração atual e na anterior  
  acVS1=0
  acVS2=0
  acVS3=0
  VeroSimActual=0
  
  for (i in 1:T) {#Calculo do primeiro segmento da LL
    acVS1 = acVS1 + Y[i]*log(theta_hat[S_treino[i]])
  }
  for (i in 1:T) {#Calculo do segundo segmento da LL
    acVS2 = acVS2 + (Nt[i]-Y[i])*(log(1-theta_hat[S_treino[i]]))
  }
  temp=NULL
  for (i in 2:T) {#Calculo do terceiro segmento da LL
    for (g in 1:K) {
      temp[g]<-exp(X[i,]%*%matrix(BetaArray[g,,S_treino[t-1]],ncol=1))
    }
    acVS3 = acVS3 + (X[i,]%*%matrix(BetaArray[S_treino[i],,S_treino[i-1]]) - log(sum(temp), base = exp(1)))
  }
  VeroSimActual <- sum(log(choose(Nt,Y), base = exp(1))) + log(P0[S_treino[1]]) + acVS1 + acVS2 + acVS3 #calculo da LogVerosim
  
  #filtragem dos dados
  Xtemp11<-NULL
  Xtemp12<-NULL
  Xtemp13<-NULL
  Xtemp21<-NULL
  Xtemp22<-NULL
  Xtemp23<-NULL
  Xtemp31<-NULL
  Xtemp32<-NULL
  Xtemp33<-NULL
  
  for (t in 2:T) {
    #filtros indo para o Estado # 1
    if(S_treino[t]%in%1 && S_treino[t-1]%in%1)
      Xtemp11<-rbind(Xtemp11, X[t,])
    
    if(S_treino[t]%in%1 && S_treino[t-1]%in%2)
      Xtemp21<-rbind(Xtemp21, X[t,])
    
    if(S_treino[t]%in%1 && S_treino[t-1]%in%3)
      Xtemp31<-rbind(Xtemp31, X[t,])
    
    #Filtros indo para o Estado # 2
    if(S_treino[t]%in%2 && S_treino[t-1]%in%1)
      Xtemp12<-rbind(Xtemp12, X[t,])
    
    if(S_treino[t]%in%2 && S_treino[t-1]%in%2)
      Xtemp22<-rbind(Xtemp22, X[t,])
    
    if(S_treino[t]%in%2 && S_treino[t-1]%in%3)
      Xtemp32<-rbind(Xtemp32, X[t,])
    
    #Filtros indo para o Estado # 3
    if(S_treino[t]%in%3 && S_treino[t-1]%in%1)
      Xtemp13<-rbind(Xtemp13, X[t,])
    
    if(S_treino[t]%in%3 && S_treino[t-1]%in%2)
      Xtemp23<-rbind(Xtemp23, X[t,])
    
    if(S_treino[t]%in%3 && S_treino[t-1]%in%3)
      Xtemp33<-rbind(Xtemp33, X[t,])
  }
  
  fit1 <- optim(par = init1, fn = FSM1, control = list(fnscale=-1), method = "Nelder-Mead", hessian = FALSE)
 
  for (i in 1:K){
    for (d in 1:D){
      if (i == 1){
        BetaArray[i,d,1]=0
      } else if (i == 2){
        BetaArray[i,d,1]=fit1$par[d]
      } else if (i == 3){
        BetaArray[i,d,1]=fit1$par[D+d]
      }
      
    }
  }
  
  for (i in 1:K){
    for (d in 1:D){
      if (i == 1){
        BetaArray[i,d,2]=0
      } else if (i == 2){
        BetaArray[i,d,2]=fit1$par[(2*D + d)]
      } else if (i == 3){
        BetaArray[i,d,2]=fit1$par[((2*D) + D+d)]
      }
      
    }
  }
  
  for (i in 1:K){
    for (d in 1:D){
      if (i == 1){
        BetaArray[i,d,3]=0
      } else if (i == 2){
        BetaArray[i,d,3]=fit1$par[(4*D + d)]
      } else if (i == 3){
        BetaArray[i,d,3]=fit1$par[(4*D + D+d)]
      }
      
    }
  }
  
  
  # BetaArray[1,1,1]=0
  # BetaArray[1,2,1]=0
  # BetaArray[1,2,1]=0
  # BetaArray[2,1,1]=fit1$par[1]
  # BetaArray[2,2,1]=fit1$par[2]
  # BetaArray[2,3,1]=fit1$par[3]
  # BetaArray[3,1,1]=fit1$par[4]
  # BetaArray[3,2,1]=fit1$par[5]
  # BetaArray[3,3,1]=fit1$par[6]
  # 
  # BetaArray[1,1,2]=0
  # BetaArray[1,2,2]=0
  # BetaArray[1,3,2]=0
  # BetaArray[2,1,2]=fit2$par[1]
  # BetaArray[2,2,2]=fit2$par[2]
  # BetaArray[2,3,2]=fit2$par[3]
  # BetaArray[3,1,2]=fit2$par[4]
  # BetaArray[3,2,2]=fit2$par[5]
  # BetaArray[3,3,2]=fit2$par[6]
  # 
  # BetaArray[1,1,3]=0
  # BetaArray[1,2,3]=0
  # BetaArray[1,3,3]=0
  # BetaArray[2,1,3]=fit3$par[1]
  # BetaArray[2,2,3]=fit3$par[2]
  # BetaArray[2,3,3]=fit3$par[3]
  # BetaArray[3,1,3]=fit3$par[4]
  # BetaArray[3,2,3]=fit3$par[5]
  # BetaArray[3,3,3]=fit3$par[6] 
  
  ac2VS1=0
  ac2VS2=0
  ac2VS3=0
  VeroSimProxima=0
  
  for (i in 1:T) {#Calculo do primeiro segmento da LL
    ac2VS1 = ac2VS1 + Y[i]*log(theta_hat[S_treino[i]])
  }
  for (i in 1:T) {#Calculo do segundo segmento da LL
    ac2VS2 = ac2VS2 + (Nt[i]-Y[i])*(log(1-theta_hat[S_treino[i]]))
  }
  temp=NULL
  for (i in 2:T) {#Calculo do terceiro segmento da LL
    for (g in 1:K) {
      temp[g]<-exp(X[i,]%*%matrix(BetaArray[g,,S_treino[t-1]],ncol=1))
    }
    ac2VS3 = ac2VS3 + (X[i,]%*%matrix(BetaArray[S_treino[i],,S_treino[i-1]]) - log(sum(temp), base = exp(1)))
  }
  VeroSimProxima <- sum(log(choose(Nt,Y), base = exp(1))) + P0[S_treino[1]] + ac2VS1 + ac2VS2 + ac2VS3 #calculo da LogVerosim
  VerAct[val]<-VeroSimActual
  VerProx[val]<-VeroSimProxima
  tolval[val]<-VeroSimProxima-VeroSimActual
  print(tolval[val])
  val=val+1
}###fim da segunda rodada do EM Estocastico###

Beta_Post_Array=array(round(runif((K*D*K),0,1),1),c(K,D,K))#Utilizamos este arreglo para reorganizar os Betas numa estrutura conveniente na hora de calcular a LL e o BIC
#Inicializamos algumas variaveis para
#almacenar calculos finais importantes
Porcentagem_Acertos=NULL
Vies_Betas=NULL
Vies_Thetas=NULL
Thetas_Finais=NULL
Betas_Finais=NULL

#Comparamos S com S_treino para saber a porcentagem de acertos
acertos=0
for (i in 1:T) {
  if (S[i]==S_treino[i]){
    acertos = acertos + 1
  }
}

Porcentagem_Acertos<-acertos/T

for (i in 1:K){
  for (d in 1:D){
    if (i == 1){
      Beta_Post_Array[i,d,1]=0
    } else if (i == 2){
      Beta_Post_Array[i,d,1]=fit1$par[d]
    } else if (i == 3){
      Beta_Post_Array[i,d,1]=fit1$par[D+d]
    }
    
  }
}

for (i in 1:K){
  for (d in 1:D){
    if (i == 1){
      Beta_Post_Array[i,d,2]=0
    } else if (i == 2){
      Beta_Post_Array[i,d,2]=fit1$par[(2*D + d)]
    } else if (i == 3){
      Beta_Post_Array[i,d,2]=fit1$par[((2*D) + D+d)]
    }
    
  }
}

for (i in 1:K){
  for (d in 1:D){
    if (i == 1){
      Beta_Post_Array[i,d,3]=0
    } else if (i == 2){
      Beta_Post_Array[i,d,3]=fit1$par[(4*D + d)]
    } else if (i == 3){
      Beta_Post_Array[i,d,3]=fit1$par[(4*D + D+d)]
    }
    
  }
}


Thetas_Finais<-theta_hat



#Calculo do AICc e BIC

acumLL1=0 #usaremos para calcular 1era parte da LL
acumLL2=0 #usaremos para calcular segunda parte da LL
acumLL3=0 #usaremos para calcular 3era parte da LL
LL=0 #usaremos para calcular a logverosim
AICc=0
BIC=0

for (i in 1:T) {#Calculo do primeiro segmento da LL
  acumLL1 = acumLL1 + Y[i]*log(Thetas_Finais[S_treino[i]])
}
for (i in 1:T) {#Calculo do segundo segmento da LL
  acumLL2 = acumLL2 + (Nt[i]-Y[i])*(log(1-Thetas_Finais[S_treino[i]], base = exp(1)))
}
temp=NULL
for (i in 2:T) {#Calculo do terceiro segmento da LL
  for (g in 1:K) {
    temp[g]<-exp(X[i,]%*%matrix(Beta_Post_Array[g,,S_treino[t-1]],ncol=1))
  }
  acumLL3 = acumLL3 + (X[i,]%*%matrix(Beta_Post_Array[S_treino[i],,S_treino[i-1]]) - log(sum(temp), base = exp(1)))
}
LL <- sum(log(choose(Nt,Y), base = exp(1))) + sum(log(P0[])) + acumLL1 + acumLL2 + acumLL3 #calculo da LogVerosim

BIC<--2*LL + 8*log(T) #Calculo do BIC para cada Replica. O 8 aqui representa a quantidade de parametros do modelo
AICc<--2*(LL) + 2*8 + (2*8^2 + 2*8)/(T - 8 - 1) #Calculo do AIC corregido
#####################################################################
#####################################################################
theta
theta_hat


Real<-NULL
Estimado <- NULL
for (i in 2:K) {
  for (j in 1:K){
    for (d in 1:D){
      Real <- c(Real, Betas[i,d,j]) 
      Estimado <- c(Estimado, Beta_Post_Array[i,d,j]) 
    }
  }
}

TP = NULL
TN = NULL
FP = NULL
FN = NULL

for (i in 1:(D*K*(K-1))){
  if ((Real[i]==0) && (abs(Estimado[i])<=0.05)){
    TP[i] = 1
  }
  if ((Real[i]!=0) && (abs(Estimado[i])>0.05)){
    TN[i]=1
  }
  if ((Real[i]!=0) && (abs(Estimado[i])<=0.05)){
    FP[i]=1
  }
  if ((Real[i]==0) && (abs(Estimado[i])>0.05)){
    FN[i]=1
  }
}

Precision <- NULL
Sensitivity <- NULL
Specificity <- NULL
Accuracy <- NULL

Precision = sum(TP, na.rm = TRUE) / ( sum(TP, na.rm = TRUE) +  sum(FP, na.rm = TRUE) )  
Sensitivity = sum(TP, na.rm = TRUE) / ( sum(TP, na.rm = TRUE) +  sum(FN, na.rm = TRUE) )
Specificity = sum(TN, na.rm = TRUE) / ( sum(TN, na.rm = TRUE) + sum(FP, na.rm = TRUE) )
Accuracy = (sum(TN, na.rm = TRUE) + sum(TP, na.rm = TRUE)) / ( sum(TP, na.rm = TRUE) + sum(TN, na.rm = TRUE) + sum(FP, na.rm = TRUE) + sum(FN, na.rm = TRUE))

Precision
Specificity
Sensitivity
Accuracy
  

df1<-data.frame(Real, Estimado)

df1



tempo_final<-proc.time()
tempoTOTAL<-(tempo_final - tempo_inicial)/60


# wb<-createWorkbook(type = "xlsx")#Criamos um livro de trabalho de excel
# sheet<-createSheet(wb, sheetName = nomefolha) #criamos uma folha no livro de trabalh
# 
# xlsx.addHeader(wb, sheet, value=header,level=1, color="black", underline=1) # escrevemos um titulo na folha 
# xlsx.addLineBreak(sheet, 1) #Um linha em branco embaixo do titulo
# xlsx.addTable(wb, sheet, df1, startCol=2)#Criamos a primeira tabela usando o primeiro dataframe
# xlsx.addLineBreak(sheet, 2)#insertamos 2 linhas em branco
# xlsx.addTable(wb, sheet, df2, startCol=2)#criamos a segunda Tabela usando o segundo dataframe
# saveWorkbook(wb, nome_arquivo)#guardamos o arquivo de excel.

