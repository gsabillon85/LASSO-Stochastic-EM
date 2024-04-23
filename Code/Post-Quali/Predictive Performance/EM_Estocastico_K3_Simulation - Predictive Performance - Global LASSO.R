library('label.switching')
#library("r2excel")
#library("xlsx")
library(e1071)

options(digits=8)
options(scipen=999)


#####Função para gerar valores uniformes discretos
rDiscreta<-function(p){
  u<-runif(1)
  P<-cumsum(p)
  val<-sum(P<u)+1
  return(val)}
#####


train_size = 0.70
validation_size = 0.20
test_size = 0.10


zero_threshold = 0.05
R <- 2 # Numero de Replicas
N=500 #Tamanho da amostra Binomial
T=400 #Cumprimento da cadeia simulada
K=3   #Numero de estados ocultos
D=10   #Quantidade de Covariaveis
tol<-0.0000001 #Nivel de tolerancia que estabelecemos como criterio de parada do EM Est
tolval=NULL
tolval[1]=1
optim_algo = "BFGS" #Algorithm to use in the optimization process


mainDir = "/home/gustavo/Projects/LASSO-NHMM/Code/Testing - Regularization Accuracy/Shrinkage Accuracy - Results"
subDir = paste("Resultados_T",toString(T),"_D",toString(D),"_zero-threshold",toString(zero_threshold),"optimethod",toString(optim_algo))
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
setwd(file.path(mainDir, subDir))

# Generate 20 random seeds to generate 20 random samples
set.seed(10)
seeds <-sample(110000000,R) # Seed number para conseguer os mesmos valores simulados
lambdas <- c(0, 0.0001, 0.001, 0.01, 0.1, 0.25, 0.5, 1, 2.5, 5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 30, 50, 75, 100, 200)
#lambdas = lambdas_large

# 

#Metricas de Performance Preditiva 

EQM_Preditivo_Validação <- NULL
EQM_Preditivo_Teste <- NULL

#Acertos na Sequência Oculta
Acertos_S_Teste <- NULL


#Metricas de Zerado dos Coeficientes
Sensitivity = NULL
Specificity = NULL
Accuracy = NULL
Pctgm_Zerado= NULL
RMSE_Parameters = NULL
Parameters = matrix(nrow = R, ncol = D*K*(K-1) )


## SEÇÃO DE DEFINICAÇÃO DOS PARAMETROS PARA SIMULAÇÃO DE DADOS ##
################################################################
P0=rep(1/K,K) #Inicializamos vetor de probabilidades inciais para o HMM
Betas=array(0, dim=c(K,D,K)) # valores de beta utilizados na geração dos valores (consideranda intercepto e duas covariáveis)
Real=NULL

#Fazemos o valor iniciais dos Betas da transição 1 igual a 0, para ter a função de Ligação mlogit
#Além disso setamos alguns Betas manualmente para controla a simulação
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
Betas[2,3,2]=-1.5
Betas[3,1,2]=0
Betas[3,2,2]=0
Betas[3,3,2]=2.2

Betas[1,1,3]=0
Betas[1,2,3]=0
Betas[1,3,3]=0
Betas[2,1,3]=0
Betas[2,2,3]=0
Betas[2,3,3]=-1.3
Betas[3,1,3]=0
Betas[3,2,3]=-2.3
Betas[3,3,3]=1.7

# -- Caracteristicas Especificas para a distribuição das VA observaveis
theta=c(0.25,0.5,0.75) # vetor com a probabilidade de sucesso das K distribuiçoes Binomiais
Nt=rep(N,T) # número de ensaios de Bernoulli associado a dada uma das T variáveis Binomiais. Cada coloquei tudo igual mas eles podem diferentes.

#Criar vetor de valores reais de Betas para comparar com o vetor de betas estimados
for (i in 2:K) {
  for (j in 1:K){
    for (d in 1:D){
      Real <- c(Real, Betas[i,d,j]) 
    }
  }
}
################################################################
## FIM DA SEÇÃO DE DEFINICAÇÃO DOS PARAMETROS PARA SIMULAÇÃO DE DADOS ##


tempo_inicial<-proc.time()
for (p in 1:R){
  set.seed(seeds[p])
  #   INICIO DE SIMULAÇÃO DOS DADOS ##
  #########################################
  
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
  ###########################################
  #   FIM DA SIMULAÇÂO DOS DADOS ##
  
  #   SEPARAR BASES EM TREINO< VALIDAÇÃO E TESTE
  ##############################################
  
  #Calcula os indices de corte
  cutoff_treino = length(Y)*train_size
  cutoff_validation = length(Y)*(train_size+validation_size)
  
  #Cria as bases 
  S_training = S[1:cutoff_treino] 
  Y_training = Y[1:cutoff_treino]
  X_training = X[1:cutoff_treino, ]
  
  S_validation = S[(cutoff_treino+1):cutoff_validation] 
  Y_validation = Y[(cutoff_treino+1):cutoff_validation]
  X_validation = X[(cutoff_treino+1):cutoff_validation, ]
  
  S_test = S[(cutoff_validation+1):T] 
  Y_test = Y[(cutoff_validation+1):T]
  X_test = X[(cutoff_validation+1):T, ]
  
  Nt_training = rep(N,cutoff_treino)
  Nt_validation = rep(N,length(S_validation))
  Nt_validation = rep(N,length(S_test))
  ##############################################
  # FIM DE SEPARAÇÃO DAS BASES EM TREINO, VALIDATION E TESTE
  
  
  
  #######Primeiro geramos uma sequência não observavel de treinamento######
  P_Treino=rep(1/K,K) #Vetor de probabilidade utilizadas para gerar a sequência de treino
  S_treino<-NULL # Inicializamos a sequência oculta de treinamento
  
  init1 = c(rnorm(D*(K-1), 0, 5))#Valores iniciais para os Betas_1
  init2 = c(rnorm(D*(K-1), 0, 5))#Valores iniciais para os Betas_2
  init3 = c(rnorm(D*(K-1), 0, 5))#Valores iniciais para os Betas_3
  
  lasso_iterator = 1
  
  lasso_RMSE <- NULL
  lasso_theta_hat_estimates <- matrix(nrow = length(lambdas), ncol = K)
  lasso_Beta_estimates <- matrix(nrow = length(lambdas), ncol = D*K*(K-1))
  lasso_Beta_arrays <- array(rep(0,K*D*K*length(lambdas)), dim=c(K,D,K,length(lambdas)))
  for (h in 1:length(lambdas)){
    tempo_inicial<-proc.time()#Calcularemos quanto demoro todo o proceso
    #Estruturas necessarias no processo de estimação
    theta_hat = NULL #Variavel para estimar os thetas em cada iteração do EM Estocastico
    BetaArray = array(0, dim=c(K,D,K)) #Estrutura para guardar as estimativas dos Betas em cada iteração do EM
    
    
    VerProx<-NULL
    VerAct<-NULL
    
    ####### Simulação #######
    
    lambda = lambdas[lasso_iterator]
    
    #######   Escrevemos as funções que serão o objetivo da optimização   ######
    # Com o temos um array de Betas, utilizaremos tres funções para achar os valores otimos
    # Uma para a matriz Betas[,,1] uma para a matriz Betas[,,2] e uma para 
    # a matriz Betas[,,3]
    lambda <- lambdas[h]
    FSM1 <-function(params){#função a maximizar para achar os Betas_1
      resp <- sum(1 - log(1 + exp(Xtemp11%*%params[1:D])+ exp(Xtemp11%*%params[(D+1):(2*D)]))) + sum((Xtemp12%*%params[1:D]) - log( 1 + exp(Xtemp12%*%params[1:D])+ exp(Xtemp12%*%params[(D+1):(2*D)]) )) + sum((Xtemp13%*%params[(D+1):(2*D)]) - log( 1 + exp(Xtemp13%*%params[1:D])+ exp(Xtemp13%*%params[(D+1):(2*D)]) )) - lambda*(sum(abs(params[2:D])) + sum(abs(params[(D+2):(2*D)])))
    }
    
    FSM2 <-function(params){#função a maximizar para achar os Betas_2
      resp <- sum(1 - log(1 + exp(Xtemp21%*%params[1:D])+ exp(Xtemp21%*%params[(D+1):(2*D)]))) + sum((Xtemp22%*%params[1:D]) - log( 1 + exp(Xtemp22%*%params[1:D])+ exp(Xtemp22%*%params[(D+1):(2*D)]) )) + sum((Xtemp23%*%params[(D+1):(2*D)]) - log( 1 + exp(Xtemp23%*%params[1:D])+ exp(Xtemp23%*%params[(D+1):(2*D)]) )) - lambda*(sum(abs(params[2:D])) + sum(abs(params[(D+2):(2*D)])))
    }
    
    FSM3 <-function(params){#função a maximizar para achar os Betas_3
      resp <- sum(1 - log(1 + exp(Xtemp31%*%params[1:D])+ exp(Xtemp31%*%params[(D+1):(2*D)]))) + sum((Xtemp32%*%params[1:D]) - log( 1 + exp(Xtemp32%*%params[1:D])+ exp(Xtemp32%*%params[(D+1):(2*D)]) )) + sum((Xtemp33%*%params[(D+1):(2*D)]) - log( 1 + exp(Xtemp33%*%params[1:D])+ exp(Xtemp33%*%params[(D+1):(2*D)]) )) - lambda*(sum(abs(params[2:D])) + sum(abs(params[(D+2):(2*D)])))
    }
    
    
    
    #########################################
    #####   Procedimento de Estimação   #####
    #########################################
    
    
  
    #Geramos uma sequência de treinamento
    for (i in 1:length(S_training)) {
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
    #Agora executamos o Algoritmo EM Estocástico
    while ( tolval[val]>tol ){
      #print(tolval[val])
      #VeroSimActual=VeroSimProxima
      #Aqui devemos calcular a diferença entre a L.V. em na iteração atual e na anterior  
      acVS1=0
      acVS2=0
      acVS3=0
      VeroSimActual=0
      
      for (k in 1:K) {
        id <- S_treino == k
        theta_hat[k] <- sum(id*Y_training)/sum(id*Nt_training)
      }
      
      for (i in 1:length(S_training)) {#Calculo do primeiro segmento da LL
        acVS1 = acVS1 + Y_training[i]*log(theta_hat[S_treino[i]])
      }
      for (i in 1:length(S_training)) {#Calculo do segundo segmento da LL
        acVS2 = acVS2 + (Nt_training[i]-Y_training[i])*(log(1-theta_hat[S_treino[i]]))
      }
      temp=NULL
      for (i in 2:length(S_training)) {#Calculo do terceiro segmento da LL
        for (g in 1:K) {
          temp[g]<-exp(X[i,]%*%matrix(BetaArray[g,,S_treino[i-1]],ncol=1))
        }
        acVS3 = acVS3 + (X[i,]%*%matrix(BetaArray[S_treino[i],,S_treino[i-1]]) - log(sum(temp), base = exp(1)))
      }
      VeroSimActual <- sum(log(choose(Nt_training,Y_training), base = exp(1))) + log(P0[S_treino[1]]) + acVS1 + acVS2 + acVS3 #calculo da LogVerosim
      
      #Calculamos a sequência S_treino utilizando os Betas
      #Atualizados na iteração passada e os valores observados Y
      S_treino[1]=which.max(dbinom(Y[1], Nt_training[1], theta_hat))
      for (i in 2:length(S_training)) {
        A_hat_t = Mat_trans(X[i,])
        prob<-(A_hat_t[S_treino[i], ]*dbinom(Y[i], Nt_training[i], theta_hat))/sum(A_hat_t[S_treino[i], ]*dbinom(Y[i], Nt_training[i], theta_hat))
        S_treino[i]=rDiscreta(prob)
        #S_treino[i]=which.max(A_hat_t[S_treino[i-1], ]*dbinom(Y[i], Nt[i], theta_hat))
      }
      
      S_treino[is.na(S_treino)] <- 1
      
      if (length(S_treino[is.na(S_treino)]) > 0){
        print(length(S_treino[is.na(S_treino)]))
      }
        
      #####################################
      #Este segmento de codigo testa se aconteceram todas as transições possiveis
      #No caso que elas não tinham acontecido, as que
      #não aconteceram são forçadas a acontecer
      TransCount <- matrix(data = c(rep(0,K^2)), nrow = K, ncol = K)
      for (i in 2:length(S_treino)) {
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
            positions = sample(2:length(S_treino), 4)
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
      
      for (t in 2:length(S_training)) {
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
      fit1 <- optim(par = init1, fn = FSM1, control = list(fnscale=-1), method = optim_algo, hessian = FALSE)
      fit2 <- optim(par = init2, fn = FSM2, control = list(fnscale=-1), method = optim_algo, hessian = FALSE)
      fit3 <- optim(par = init3, fn = FSM3, control = list(fnscale=-1), method = optim_algo, hessian = FALSE)
      
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
            BetaArray[i,d,2]=fit2$par[d]
          } else if (i == 3){
            BetaArray[i,d,2]=fit2$par[D+d]
          }
          
        }
      }
      
      for (i in 1:K){
        for (d in 1:D){
          if (i == 1){
            BetaArray[i,d,3]=0
          } else if (i == 2){
            BetaArray[i,d,3]=fit3$par[d]
          } else if (i == 3){
            BetaArray[i,d,3]=fit3$par[D+d]
          }
          
        }
      }
      
      
      
      ac2VS1=0
      ac2VS2=0
      ac2VS3=0
      VeroSimProxima=0
      
      for (i in 1:length(S_training)) {#Calculo do primeiro segmento da LL
        ac2VS1 = ac2VS1 + Y_training[i]*log(theta_hat[S_treino[i]])
      }
      for (i in 1:length(S_training)) {#Calculo do segundo segmento da LL
        ac2VS2 = ac2VS2 + (Nt_training[i]-Y_training[i])*(log(1-theta_hat[S_treino[i]]))
      }
      temp=NULL
      for (i in 2:length(S_training)) {#Calculo do terceiro segmento da LL
        for (g in 1:K) {
          temp[g]<-exp(X[i,]%*%matrix(BetaArray[g,,S_treino[i-1]],ncol=1))
        }
        ac2VS3 = ac2VS3 + (X[i,]%*%matrix(BetaArray[S_treino[i],,S_treino[i-1]]) - log(sum(temp), base = exp(1)))
      }
      VeroSimProxima <- sum(log(choose(Nt_training,Y_training), base = exp(1))) + log(P0[S_treino[1]]) + ac2VS1 + ac2VS2 + ac2VS3 #calculo da LogVerosim
      
      val=val+1
      VerAct[val]<-VeroSimActual
      VerProx[val]<-VeroSimProxima
      tolval[val]<-VeroSimProxima - VeroSimActual
      #print(tolval[val])
      
    }#######Fim da primeira rodada do EM Estocastico#######
    
    #Criar algumas matrizes para fazer calculos e manipular a saida MCMC
    #nestas matrizes, as estimativas serão reordenadas usando o metodo ECR
    mat_thetar<-matrix(nrow = 1, ncol = K)
    reorder_S<-matrix(nrow = 1, ncol = length(S_training))
    mat_S<-matrix(nrow = 1, ncol = length(S_training))
    mat_S[1,]<-S_treino
    zpvt_S = S #Como pivot para o metodo ECR usamos o S original
    perms_S = ecr(zpivot = zpvt_S, z = mat_S, K = 3)# aplicamos o metodo ECR que retornara as permutações das dos estados ocultos que devem ser utilizadas para reordenar a saida do algoritmo bayesiano
    
    # Reordenamos a saido do algoritmo EMEst usando as 
    # permutações fornecidas pelo ECR para K=3 
    # só rerotulamos a Sequência S_treino, e reordenamos os Thetas
    # Os Betas serão estimados usando a sequência S_Treino rerotulada
    # e os Thetas, na segunda etapa do EMEst
    
    for (i in 1:1) {
      for (j in 1:length(S_training)) {
        if(S_treino[j]!=S_training[j] && ((perms_S$permutations[i,1]==2 && perms_S$permutations[i,2]==3 && perms_S$permutations[i,3]==1) | (perms_S$permutations[i,1]==3 && perms_S$permutations[i,2]==1 && perms_S$permutations[i,3]==2))){
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
    tolval=NULL
    tolval[1]=3
    tol2 = 2
    while ( abs(tolval[val])>tol2 ) {
      #print(tolval[val])
      #Aqui devemos calcular a diferença entre a L.V. em na iteração atual e na anterior  
      acVS1=0
      acVS2=0
      acVS3=0
      VeroSimActual=0
      
      for (i in 1:length(S_training)) {#Calculo do primeiro segmento da LL
        acVS1 = acVS1 + Y_training[i]*log(theta_hat[S_treino[i]])
      }
      for (i in 1:length(S_training)) {#Calculo do segundo segmento da LL
        acVS2 = acVS2 + (Nt_training[i]-Y_training[i])*(log(1-theta_hat[S_treino[i]]))
      }
      temp=NULL
      for (i in 2:length(S_training)) {#Calculo do terceiro segmento da LL
        for (g in 1:K) {
          temp[g]<-exp(X[i,]%*%matrix(BetaArray[g,,S_treino[t-1]],ncol=1))
        }
        acVS3 = acVS3 + (X[i,]%*%matrix(BetaArray[S_treino[i],,S_treino[i-1]]) - log(sum(temp), base = exp(1)))
      }
      VeroSimActual <- sum(log(choose(Nt_training,Y_training), base = exp(1))) + log(P0[S_treino[1]]) + acVS1 + acVS2 + acVS3 #calculo da LogVerosim
      
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
      
      for (t in 2:length(S_training)) {
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
      
      fit1 <- optim(par = init1, fn = FSM1, control = list(fnscale=-1), method = optim_algo, hessian = FALSE)
      fit2 <- optim(par = init2, fn = FSM2, control = list(fnscale=-1), method = optim_algo, hessian = FALSE)
      fit3 <- optim(par = init3, fn = FSM3, control = list(fnscale=-1), method = optim_algo, hessian = FALSE)
      
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
            BetaArray[i,d,2]=fit2$par[d]
          } else if (i == 3){
            BetaArray[i,d,2]=fit2$par[D+d]
          }
          
        }
      }
      
      for (i in 1:K){
        for (d in 1:D){
          if (i == 1){
            BetaArray[i,d,3]=0
          } else if (i == 2){
            BetaArray[i,d,3]=fit3$par[d]
          } else if (i == 3){
            BetaArray[i,d,3]=fit3$par[D+d]
          }
          
        }
      }
      
      
      ac2VS1=0
      ac2VS2=0
      ac2VS3=0
      VeroSimProxima=0
      
      for (i in 1:length(S_training)) {#Calculo do primeiro segmento da LL
        ac2VS1 = ac2VS1 + Y_training[i]*log(theta_hat[S_treino[i]])
      }
      for (i in 1:length(S_training)) {#Calculo do segundo segmento da LL
        ac2VS2 = ac2VS2 + (Nt_training[i]-Y_training[i])*(log(1-theta_hat[S_treino[i]]))
      }
      temp=NULL
      for (i in 2:length(S_training)) {#Calculo do terceiro segmento da LL
        for (g in 1:K) {
          temp[g]<-exp(X[i,]%*%matrix(BetaArray[g,,S_treino[t-1]],ncol=1))
        }
        ac2VS3 = ac2VS3 + (X[i,]%*%matrix(BetaArray[S_treino[i],,S_treino[i-1]]) - log(sum(temp), base = exp(1)))
      }
      VeroSimProxima <- sum(log(choose(Nt_training,Y_training), base = exp(1))) + P0[S_treino[1]] + ac2VS1 + ac2VS2 + ac2VS3 #calculo da LogVerosim
      val=val+1
      VerAct[val]<-VeroSimActual
      VerProx[val]<-VeroSimProxima
      tolval[val]<-VeroSimProxima-VeroSimActual
      #print(tolval[val])
    }###fim da segunda rodada do EM Estocastico###
    
    Y_hat_validation = NULL
    S_hat_validation = NULL
    S_hat_validation[1]<-S_validation[1] #O valor para o primeiro estado oculto
    Y_hat_validation[1]<-rbinom(1,Nt_training[1],theta_hat[S_hat_validation[1]]) # O valor para o primeiro valor observavel
    for (t in 2:length(Y_validation)){
      prob<-NULL
      for (i in 1:K) prob[i]<-exp(X_validation[t,]%*%matrix(BetaArray[i,,S_hat_validation[t-1]],ncol=1))
      prob<-prob/sum(prob)
      S_hat_validation[t]<-rDiscreta(prob)
      Y_hat_validation[t]<-rbinom(1,Nt[t],theta_hat[S_hat_validation[t]])
    }
    
    Beta_Estimates <- NULL
    for (i in 2:K) {
      for (j in 1:K){
        for (d in 1:D){
          Beta_Estimates <- c(Beta_Estimates, BetaArray[i,d,j]) 
        }
      }
    }
    lasso_Beta_estimates[lasso_iterator,] <- Beta_Estimates
    lasso_theta_hat_estimates[lasso_iterator,] <- theta_hat
    lasso_RMSE[lasso_iterator] <- (sum((Y_hat_validation - Y_validation)^2))/length(Y_validation)
    lasso_Beta_arrays[,,,lasso_iterator] <- BetaArray
    lasso_iterator = lasso_iterator+1
  } # FIM DO LASSO
  
  # CAPTURA DE METRICAS PARA CADA REPLICA (Usando conjunto de Teste)
  ##################################################################
  
  
  
  
  
} #FIM DAS REPLICAS 
  
lasso_Beta_estimates[which.min(lasso_RMSE),]
Real

df_avgs<-data.frame(cbind(lambdas, avg_Sensitivity, avg_Specificity, avg_Accuracy, avg_pctgm_zerado, avg_RMSE))
df_sd <- data.frame(cbind(lambdas, sd_Sensitivity, sd_Specificity, sd_Accuracy, sd_pctgm_zerado, sd_RMSE))
df_full <- data.frame(cbind(lambdas, avg_Sensitivity, sd_Sensitivity, avg_Specificity, sd_Specificity, avg_Accuracy, sd_Accuracy, avg_pctgm_zerado, sd_pctgm_zerado, avg_RMSE, sd_RMSE))


tempo_final<-proc.time()
tempoTOTAL<-(tempo_final - tempo_inicial)/60

df_full  




######################################
#Sensitividade
sample.n <- 30 
lower.bound = NULL
upper.bound = NULL
for (z in 1:length(lambdas)){
  sample.mean <- df_avgs$avg_Sensitivity[z]
  sample.se <- df_sd$sd_Sensitivity[z]
  
  alpha = 0.05
  degrees.freedom = sample.n - 1
  t.score = qt(p=alpha/2, df=degrees.freedom,lower.tail=F)
  
  margin.error <- t.score * sample.se
  print(margin.error)
  lower.bound[z] <- sample.mean - margin.error
  upper.bound[z] <- sample.mean + margin.error
}

graph_df <- data.frame(x =1:length(lambdas),
                       F = df_avgs$avg_Sensitivity*100,
                       L =lower.bound*100,
                       U =upper.bound*100)


jpeg(paste("Sensitivity_T-",toString(T),"_D-",toString(D),"_zerothreshold-",toString(zero_threshold),"_optimethod-",toString(optim_algo),".jpg", sep = ""), width = 1400, height = 700)

plot(graph_df$x, xaxt="n", graph_df$F, ylim = c(0,100),  type = "l", xlab = "Valor de Lambda", ylab = "Porcentagem")
axis(1, at=1:22, labels=lambdas, cex.axis=0.75)
title(paste("Sensitivity T=",toString(T)," D=",toString(D)," zero-threshold=",toString(zero_threshold)," optim-method=",toString(optim_algo),".jpg", sep = ""))
#make polygon where coordinates start with lower limit and 
# then upper limit in reverse order
polygon(c(graph_df$x,rev(graph_df$x)),c(graph_df$L,rev(graph_df$U)),col = "grey75", border = FALSE)
lines(graph_df$x, graph_df$F, lwd = 1)
#add red lines on borders of polygon
lines(graph_df$x, graph_df$U, col="red",lty=2)
lines(graph_df$x, graph_df$L, col="red",lty=2) 
dev.off() 


#Specificidade
lower.bound = NULL
upper.bound = NULL
for (z in 1:length(lambdas)){
  sample.mean <- df_avgs$avg_Specificity[z]
  sample.se <- df_sd$sd_Specificity[z]
  
  alpha = 0.05
  degrees.freedom = sample.n - 1
  t.score = qt(p=alpha/2, df=degrees.freedom,lower.tail=F)
  
  margin.error <- t.score * sample.se
  print(margin.error)
  lower.bound[z] <- sample.mean - margin.error
  upper.bound[z] <- sample.mean + margin.error
}

graph_df <- data.frame(x =1:length(lambdas),
                       F = df_avgs$avg_Specificity*100,
                       L =lower.bound*100,
                       U =upper.bound*100)


jpeg(paste("Specificity_T-",toString(T),"_D-",toString(D),"_zerothreshold-",toString(zero_threshold),"_optimethod-",toString(optim_algo),".jpg", sep = ""), width = 1400, height = 700)

plot(graph_df$x, xaxt="n", graph_df$F, ylim = c(0,100),  type = "l", xlab = "Valor de Lambda", ylab = "Porcentagem")
axis(1, at=1:22, labels=lambdas, cex.axis=0.75)
title(paste("Specificity T=",toString(T)," D=",toString(D)," zero-threshold=",toString(zero_threshold)," optim-method=",toString(optim_algo),".jpg", sep = ""))
#make polygon where coordinates start with lower limit and 
# then upper limit in reverse order
polygon(c(graph_df$x,rev(graph_df$x)),c(graph_df$L,rev(graph_df$U)),col = "grey75", border = FALSE)
lines(graph_df$x, graph_df$F, lwd = 1)
#add red lines on borders of polygon
lines(graph_df$x, graph_df$U, col="red",lty=2)
lines(graph_df$x, graph_df$L, col="red",lty=2) 
dev.off() 


#Acuracia
lower.bound = NULL
upper.bound = NULL
for (z in 1:length(lambdas)){
  sample.mean <- df_avgs$avg_Accuracy[z]
  sample.se <- df_sd$sd_Accuracy[z]
  
  alpha = 0.05
  degrees.freedom = sample.n - 1
  t.score = qt(p=alpha/2, df=degrees.freedom,lower.tail=F)
  
  margin.error <- t.score * sample.se
  print(margin.error)
  lower.bound[z] <- sample.mean - margin.error
  upper.bound[z] <- sample.mean + margin.error
}

graph_df <- data.frame(x =1:length(lambdas),
                       F = df_avgs$avg_Accuracy*100,
                       L =lower.bound*100,
                       U =upper.bound*100)


jpeg(paste("Accuracy_T-",toString(T),"_D-",toString(D),"_zerothreshold-",toString(zero_threshold),"_optimethod-",toString(optim_algo),".jpg", sep = ""), width = 1400, height = 700)

plot(graph_df$x, xaxt="n", graph_df$F, ylim = c(0,100),  type = "l", xlab = "Valor de Lambda", ylab = "Porcentagem")
axis(1, at=1:22, labels=lambdas, cex.axis=0.75)
title(paste("Accuracy T=",toString(T)," D=",toString(D)," zero-threshold=",toString(zero_threshold)," optim-method=",toString(optim_algo),".jpg", sep = ""))
#make polygon where coordinates start with lower limit and 
# then upper limit in reverse order
polygon(c(graph_df$x,rev(graph_df$x)),c(graph_df$L,rev(graph_df$U)),col = "grey75", border = FALSE)
lines(graph_df$x, graph_df$F, lwd = 1)
#add red lines on borders of polygon
lines(graph_df$x, graph_df$U, col="red",lty=2)
lines(graph_df$x, graph_df$L, col="red",lty=2) 
dev.off() 


#Taxa de Zerado
lower.bound = NULL
upper.bound = NULL
for (z in 1:length(lambdas)){
  sample.mean <- df_avgs$avg_pctgm_zerado[z]
  sample.se <- df_sd$sd_pctgm_zerado[z]
  
  alpha = 0.05
  degrees.freedom = sample.n - 1
  t.score = qt(p=alpha/2, df=degrees.freedom,lower.tail=F)
  
  margin.error <- t.score * sample.se
  print(margin.error)
  lower.bound[z] <- sample.mean - margin.error
  upper.bound[z] <- sample.mean + margin.error
}

graph_df <- data.frame(x =1:length(lambdas),
                       F = df_avgs$avg_pctgm_zerado*100,
                       L =lower.bound*100,
                       U =upper.bound*100)


jpeg(paste("PctgmZerado_T-",toString(T),"_D-",toString(D),"_zerothreshold-",toString(zero_threshold),"_optimethod-",toString(optim_algo),".jpg", sep = ""), width = 1400, height = 700)

plot(graph_df$x, xaxt="n", graph_df$F, ylim = c(0,100),  type = "l", xlab = "Valor de Lambda", ylab = "Porcentagem")
axis(1, at=1:22, labels=lambdas, cex.axis=0.75)
title(paste("PctgmZerado T=",toString(T)," D=",toString(D)," zero-threshold=",toString(zero_threshold)," optim-method=",toString(optim_algo),".jpg", sep = ""))
#make polygon where coordinates start with lower limit and 
# then upper limit in reverse order
polygon(c(graph_df$x,rev(graph_df$x)),c(graph_df$L,rev(graph_df$U)),col = "grey75", border = FALSE)
lines(graph_df$x, graph_df$F, lwd = 1)
#add red lines on borders of polygon
lines(graph_df$x, graph_df$U, col="red",lty=2)
lines(graph_df$x, graph_df$L, col="red",lty=2) 
dev.off() 


#RMSE
lower.bound = NULL
upper.bound = NULL
for (z in 1:length(lambdas)){
  sample.mean <- df_avgs$avg_RMSE[z]
  sample.se <- df_sd$sd_RMSE[z]
  
  alpha = 0.05
  degrees.freedom = sample.n - 1
  t.score = qt(p=alpha/2, df=degrees.freedom,lower.tail=F)
  
  margin.error <- t.score * sample.se
  print(margin.error)
  lower.bound[z] <- sample.mean - margin.error
  upper.bound[z] <- sample.mean + margin.error
}

graph_df <- data.frame(x =1:length(lambdas),
                       F = df_avgs$avg_RMSE,
                       L =lower.bound,
                       U =upper.bound)


jpeg(paste("RMSE_T-",toString(T),"_D-",toString(D),"_zerothreshold-",toString(zero_threshold),"_optimethod-",toString(optim_algo),".jpg", sep = ""), width = 1400, height = 700)

plot(graph_df$x, xaxt="n", graph_df$F,  type = "l", xlab = "Valor de Lambda", ylab = "RMSE")
axis(1, at=1:22, labels=lambdas, cex.axis=0.75)
title(paste("RMSE T=",toString(T)," D=",toString(D)," zero-threshold=",toString(zero_threshold)," optim-method=",toString(optim_algo),".jpg", sep = ""))
#make polygon where coordinates start with lower limit and 
# then upper limit in reverse order
polygon(c(graph_df$x,rev(graph_df$x)),c(graph_df$L,rev(graph_df$U)),col = "grey75", border = FALSE)
lines(graph_df$x, graph_df$F, lwd = 1)
#add red lines on borders of polygon
lines(graph_df$x, graph_df$U, col="red",lty=2)
lines(graph_df$x, graph_df$L, col="red",lty=2) 
dev.off() 


df_params <- data.frame(avg_Parameters)
write.csv(df_full, paste("ResultsData_T-",toString(T),"_D-",toString(D),"_zerothreshold-",toString(zero_threshold),"_optimethod-",toString(optim_algo),".csv", sep = ""), row.names=FALSE)
write.csv(df_params, paste("ParameterEstimates_T-",toString(T),"_D-",toString(D),"_zerothreshold-",toString(zero_threshold),"_optimethod-",toString(optim_algo),".csv", sep = ""), row.names=FALSE)


# Create an empty plot
jpeg(paste("Coeficientes_T-",toString(T),"_D-",toString(D),"_zerothreshold-",toString(zero_threshold),"_optimethod-",toString(optim_algo),".jpg", sep = ""), width = 1400, height = 1000)
plot(lambdas, xaxt="n", avg_Parameters[,1],  type = "l", xlim = c(0, 25), ylim = c(-3, 3), 
     xlab = "Lambda", ylab = "Valor do Coeficiente", main = "Coeficientes de Transição")
axis(1, at=1:22, labels=lambdas, cex.axis=0.75)
for (v in 1:(D*K*(K-1))){
  if (Real[v] == 0)
  color = "red"
  else {
    color = "blue"
  }
  lines(lambdas, avg_Parameters[,v], type = "l", col = color)
}
dev.off() 

