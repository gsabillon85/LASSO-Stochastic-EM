library('label.switching')
library(glmnet) # Regressão Linear Penalizada
library(forecast) # Auto Regressive Integrated Moving Average
#library("r2excel")
#library("xlsx")
library(e1071)

options(digits=4)
options(scipen=999)



#####Função para gerar valores uniformes discretos
rDiscreta<-function(p){
  u<-runif(1)
  P<-cumsum(p)
  val<-sum(P<u)+1
  return(val)}
#####


train_size = 0.80
validation_size = 0.15
test_size = 0.05


zero_threshold = 0.05
R <- 30 # Numero de Replicas
T=2000 #Cumprimento da cadeia simulada
K=3   #Numero de estados ocultos
D=6   #Quantidade de Covariaveis
tol<-0.0000001 #Nivel de tolerancia que estabelecemos como criterio de parada do EM Est
tolval=NULL
tolval[1]=1
optim_algo = "BFGS" #Algorithm to use in the optimization process
n_max_iter_EM = 100
Tempo <- NULL

cenario <- "TESTANDO_CODIGO"
mainDir = paste("/home/gustavo/Projects/LASSO-NHMM/Code/Post-Quali/Predictive Performance - Normal/Global LASSO Results/K3/",cenario,sep = "")
subDir = paste("Resultados_T",toString(T),"_D",toString(D),"_zero-threshold",toString(zero_threshold),"optimethod",toString(optim_algo))
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
setwd(file.path(mainDir, subDir))

set.seed(2)
seeds <-sample(110000000,R) # Seed number para conseguer os mesmos valores simulados
lambdas <- c(0, 0.0001, 0.001, 0.01, 0.1, 0.25, 0.5, 1, 2.5, 5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 30, 50, 75)


# seeds[2] = 305
# seeds[5] = 312
# seeds[6] = 308
# seeds[8] = 300
# seeds[9] = 3005
# seeds[10] = 30106
# seeds[11] = 3008
# seeds[14] = 30031
# seeds[15] = 860
# seeds[16] = 301065
# seeds[17] = 3010878
# seeds[17] = 3010873
# seeds[19] = 2310
# seeds[21] = 353005
# seeds[22] = 47
# seeds[23] = 443008
# seeds[26] = 30090
# seeds[27] = 30051
# seeds[28] = 317
# seeds[29] = 43611


# 
#Metricas de Performance de Estimação dos ParÂmetros das VA observáveis
Mu_Rep=matrix(nrow = R, ncol = K )
Sigma_Rep=matrix(nrow = R, ncol = K )
Best_Beta_Arrays <- array(rep(0,K*D*K*R), dim=c(K,D,K,R))

#Metricas de Performance Preditiva 
MSPE_Validação <- NULL
MSPE_Teste <- NULL
Y_test_DF <- matrix(nrow = R, ncol = (test_size*T))
Y_hat_test_NHMM_DF <- matrix(nrow = R, ncol = (test_size*T))
Y_hat_test_glmnet_DF <- matrix(nrow = R, ncol = (test_size*T))
Y_hat_test_arima_DF <- matrix(nrow = R, ncol = (test_size*T))

#Acertos na Sequência Oculta
Acertos_S_Validação <- rep(0, R)
Acertos_S_Teste <- rep(0, R)


#Metricas de Zerado dos Coeficientes
Sensitivity = NULL
Specificity = NULL
Accuracy = NULL
Pctgm_Zerado= NULL
RMSE_Parameters = NULL
Best_Beta_Estimates = matrix(nrow = R, ncol = D*K*(K-1) )

#Outras informações de possivel interesse
Best_S = matrix(nrow = R, ncol = (validation_size*T))
Best_Y = matrix(nrow = R, ncol = (validation_size*T))
True_Y_Test = matrix(nrow = R, ncol = (T))
True_S_Test = matrix(nrow = R, ncol = (T))
RealTransCount <- array(rep(0,K*K*R), dim=c(K,K,R))


## SEÇÃO DE DEFINICAÇÃO DOS PARAMETROS PARA SIMULAÇÃO DE DADOS ##
################################################################
P0=rep(1/K,K) #Inicializamos vetor de probabilidades inciais para o HMM
Betas=array(0, dim=c(K,D,K)) # valores de beta utilizados na geração dos valores (consideranda intercepto e duas covariáveis)
Real=NULL

Betas[2,1,1]=-1.5
Betas[2,2,1]=-1.5
Betas[2,3,1]=-2.6

Betas[2,1,2]=-2.0
Betas[2,2,2]=2.6
Betas[2,3,2]=1.4

Betas[2,1,3]=-2.4
Betas[2,2,3]=2.1
Betas[2,3,3]=-1.5

# Betas transição 3
Betas[3,1,1]=-1.3
Betas[3,2,1]=-3.2
Betas[3,3,1]=-2.4

Betas[3,1,2]=-1.3
Betas[3,2,2]=1.7
Betas[3,3,2]=1.3

Betas[3,1,3]=-1.3
Betas[3,2,3]=-2.7
Betas[3,3,3]=-2.5


# -- Caracteristicas Especificas para a distribuição das VA observaveis
mu = c(20,150,250) # vetor com media para as duas Normais
sigma = c(1,1.5,0.5) #Vetor com os desvios padrões para as duas normais


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


##########################
# glmnet
MSPE_Teste_glmnet <- NULL
##########################

##########################
# arima 

MSPE_Teste_arima <- NULL
##########################

for (p in 1:R){
  tempo_inicial<-Sys.time()
  set.seed(seeds[p])
  cat(paste("\nNumero de Replica:",toString(p),"\n", collapse = "")) #Messagem indicando o numero da replica atual
  
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
  
  
  ###########################################
  #   FIM DA SIMULAÇÂO DOS DADOS ##
  
  print("Segmentando Base de Dados em Treino, Validação e Teste...")
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
  
  ##############################################
  # FIM DE SEPARAÇÃO DAS BASES EM TREINO, VALIDATION E TESTE
  
  
  # INICIO DO PROCESSO DE ESTIMAÇÃO
  ##########################################
  # Primeiro geramos uma sequência não observavel de treinamento
  P_Treino=rep(1/K,K) #Vetor de probabilidade utilizadas para gerar a sequência de treino
  S_treino<-NULL # Inicializamos a sequência oculta de treinamento
  
  init1 = c(rnorm(D*(K-1), 0, 5))#Valores iniciais para os Betas_1
  init2 = c(rnorm(D*(K-1), 0, 5))#Valores iniciais para os Betas_2
  init3 = c(rnorm(D*(K-1), 0, 5))#Valores iniciais para os Betas_3
  
  lasso_iterator = 1 #Criamos um contador para iterar a traves dos valores de lambda
  
  # Algumas estruturas para almacenar valores gerados pelo LASSO
  lasso_RMSE <- NULL
  lasso_S <- matrix(nrow = length(lambdas), ncol = length(S_validation))
  lasso_Y <- matrix(nrow = length(lambdas), ncol = length(S_validation))
  lasso_mu_hat_estimates <- matrix(nrow = length(lambdas), ncol = K)
  lasso_sigma_hat_estimates <- matrix(nrow = length(lambdas), ncol = K)
  lasso_Beta_estimates <- matrix(nrow = length(lambdas), ncol = D*K*(K-1))
  lasso_Beta_arrays <- array(rep(0,K*D*K*length(lambdas)), dim=c(K,D,K,length(lambdas)))
  
  # INICIO DO LASSO
  ###############################################
  pb <- txtProgressBar(min = 1, max = length(lambdas), style = 3)
  print('Executando LASSO...')
  for (h in 1:length(lambdas)){
    #setTxtProgressBar(pb, lasso_iterator)
    
    #Estruturas necessarias no processo de estimação
    mu_hat = NULL #Variavel para estimar os mus em cada iteração do EM Estocástico
    sigma_hat = NULL #Variavel para estimar os sigmas em cada iteração do EM Estocástico
    BetaArray = array(0, dim=c(K,D,K)) #Estrutura para guardar as estimativas dos Betas em cada iteração do EM
    
    
    VerProx<-NULL
    VerAct<-NULL
  
    lambda = lambdas[lasso_iterator]
    
    #######   Escrevemos as funções que serão o objetivo da optimização   ######
    # Com o temos um array de Betas, utilizaremos tres funções para achar os valores otimos
    # Uma para a matriz Betas[,,1] uma para a matriz Betas[,,2] e uma para 
    # a matriz Betas[,,3]
    FSM1 <-function(params){#função a maximizar para achar os Betas_1
      resp <- sum(1 - log(1 + exp(Xtemp11%*%params[1:D])+ exp(Xtemp11%*%params[(D+1):(2*D)]))) + sum((Xtemp12%*%params[1:D]) - log( 1 + exp(Xtemp12%*%params[1:D])+ exp(Xtemp12%*%params[(D+1):(2*D)]) )) + sum((Xtemp13%*%params[(D+1):(2*D)]) - log( 1 + exp(Xtemp13%*%params[1:D])+ exp(Xtemp13%*%params[(D+1):(2*D)]) )) - lambda*(sum(abs(params[2:D])) + sum(abs(params[(D+2):(2*D)])))
    }
    
    FSM2 <-function(params){#função a maximizar para achar os Betas_2
      resp <- sum(1 - log(1 + exp(Xtemp21%*%params[1:D])+ exp(Xtemp21%*%params[(D+1):(2*D)]))) + sum((Xtemp22%*%params[1:D]) - log( 1 + exp(Xtemp22%*%params[1:D])+ exp(Xtemp22%*%params[(D+1):(2*D)]) )) + sum((Xtemp23%*%params[(D+1):(2*D)]) - log( 1 + exp(Xtemp23%*%params[1:D])+ exp(Xtemp23%*%params[(D+1):(2*D)]) )) - lambda*(sum(abs(params[2:D])) + sum(abs(params[(D+2):(2*D)])))
    }
    
    FSM3 <-function(params){#função a maximizar para achar os Betas_3
      resp <- sum(1 - log(1 + exp(Xtemp31%*%params[1:D])+ exp(Xtemp31%*%params[(D+1):(2*D)]))) + sum((Xtemp32%*%params[1:D]) - log( 1 + exp(Xtemp32%*%params[1:D])+ exp(Xtemp32%*%params[(D+1):(2*D)]) )) + sum((Xtemp33%*%params[(D+1):(2*D)]) - log( 1 + exp(Xtemp33%*%params[1:D])+ exp(Xtemp33%*%params[(D+1):(2*D)]) )) - lambda*(sum(abs(params[2:D])) + sum(abs(params[(D+2):(2*D)])))
    }
    
    
    
    
    #   Procedimento de Estimação   
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
    while ( abs(tolval[val])>tol && val < n_max_iter_EM){
      #print(val)
      #print(tolval[val])
      #VeroSimActual=VeroSimProxima
      #Aqui devemos calcular a diferença entre a L.V. em na iteração atual e na anterior  
      LL_parte1 = 0
      LL_parte2 = 0
      LL_parte3 = 0
      LL_parte4 = 0
      VeroSimActual=0
      
      for (k in 1:K){
        id = S_treino == k
        mu_hat[k] = sum(id*Y_training)/sum(id)
        Y_id_list = split(Y_training,id)
        Y_id = unlist(Y_id_list[2], use.names = FALSE)
        sigma_hat[k] = sqrt((sum((Y_id - mu_hat[k])^2)) / (sum(id) - 1)) #DECIDIR SOBRE O ESTIMADOR DA VARIANCIA (VICIADO OU NṼICIADO)
      }
      
      #print(mu_hat)
      #print(sigma_hat)
      #Calculo da Verosimilhança como valor de tolerança
      LL_parte1 = -.5*length(S_training)*log(2*pi)
      
      for (i in 1:length(S_training)) {#Calculo do primeiro segmento da LL
        LL_parte2 = LL_parte2 -.5*log(sigma_hat[S_treino[i]]) 
      }
      for (i in 1:length(S_training)) {#Calculo do segundo segmento da LL
        LL_parte3 = LL_parte3 -(1/(2*sigma_hat[S_treino[i]]))*((Y_training[i]-mu_hat[S_treino[i]])^2)
      }
      temp=NULL
      for (i in 2:length(S_training)) {#Calculo do terceiro segmento da LL
        for (g in 1:K) {
          temp[g]<-exp(X[i,]%*%matrix(BetaArray[g,,S_treino[i-1]],ncol=1))
        }
        LL_parte4 = LL_parte4 + (X[i,]%*%matrix(BetaArray[S_treino[i],,S_treino[i-1]]) - log(sum(temp), base = exp(1)))
      }
      VeroSimActual <- log(P0[S_treino[1]]) + LL_parte1 + (LL_parte2 + LL_parte3) + LL_parte4 #calculo da LogVerosim
      
      
      #Calculamos a sequência S_treino utilizando os Betas
      #Atualizados na iteração passada e os valores observados Y
      S_treino[1]=which.max(dnorm(Y[1], mu_hat, sigma_hat))
      for (i in 2:length(S_training)) {
        A_hat_t = Mat_trans(X[i,])
        if (any(is.na(A_hat_t))){
          print("NaN encountered in Transition Matrix Calculation")
          A_hat_t[is.nan(A_hat_t)] = 1 
        }
        prob<-(A_hat_t[S_treino[i], ]*dnorm(Y_training[i], mu_hat, sigma_hat))/sum(A_hat_t[S_treino[i], ]*dnorm(Y_training[i], mu_hat, sigma_hat))
        #S_treino[i]=rDiscreta(prob)
        if (any(is.na(prob))){
          print("NaN encountered in S_treino update")
          S_treino[i]=which.max(A_hat_t[S_treino[i], ])
        } else {
          #S_treino[i]=which.max(prob)  
          S_treino[i]=which.max(prob)  
        }
      }
      
      S_treino[is.na(S_treino)] <- 1
      
      if (length(S_treino[is.na(S_treino)]) > 0){
        print(length(S_treino[is.na(S_treino)]))
      }
        
      
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
      
      fit1 <- optim(par = init1, fn = FSM1, control = list(fnscale=-1), method = "Nelder-Mead", hessian = FALSE)
      fit2 <- optim(par = init2, fn = FSM2, control = list(fnscale=-1), method = "Nelder-Mead", hessian = FALSE)
      fit3 <- optim(par = init3, fn = FSM3, control = list(fnscale=-1), method = "Nelder-Mead", hessian = FALSE)
      
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
      
      
      LL2_parte1 = 0
      LL2_parte2 = 0
      LL2_parte3 = 0
      LL2_parte4 = 0
      VeroSimProxima=0
      
      #Calculo da Verosimilhança como valor de tolerança
      LL2_parte1 = -.5*length(S_training)*log(2*pi)
      
      for (i in 1:length(S_training)) {#Calculo do primeiro segmento da LL
        LL2_parte2 = LL2_parte2 -.5*log(sigma_hat[S_treino[i]]) 
      }
      for (i in 1:length(S_training)) {#Calculo do segundo segmento da LL
        LL2_parte3 = LL2_parte3 -(1/(2*sigma_hat[S_treino[i]]))*((Y_training[i]-mu_hat[S_treino[i]])^2)
      }
      temp=NULL
      for (i in 2:length(S_training)) {#Calculo do terceiro segmento da LL
        for (g in 1:K) {
          temp[g]<-exp(X[i,]%*%matrix(BetaArray[g,,S_treino[i-1]],ncol=1))
        }
        LL2_parte4 = LL2_parte4 + (X[i,]%*%matrix(BetaArray[S_treino[i],,S_treino[i-1]]) - log(sum(temp), base = exp(1)))
      }
      VeroSimProxima <- log(P0[S_treino[1]]) + LL2_parte1 + (LL2_parte2 + LL2_parte3) + LL2_parte4 #calculo da LogVerosim
      
      val=val+1
      VerAct[val]<-VeroSimActual
      VerProx[val]<-VeroSimProxima
      tolval[val]<-VeroSimProxima - VeroSimActual
     # print(tolval[val])
      
      message(paste('\r',"Lasso iteration # ",toString(lasso_iterator),"; Valor de Lambda = ",toString(lambdas[lasso_iterator]),"; Mu_hat:",toString(round(mu_hat,3)),". Sigma_hat:",toString(round(sigma_hat,3)),"                  ", collapse = ""), appendLF = FALSE) #Messagem indicando o numero da replica atual
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
      mu_hat<-mu_hat[perms_S$permutations[i,]]
      sigma_hat<-sigma_hat[perms_S$permutations[i,]]
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
    
    while (tolval[val]>tol2) {
      #print(tolval[val])
      #Aqui devemos calcular a diferença entre a L.V. em na iteração atual e na anterior  
      LL_parte1 = 0
      LL_parte2 = 0
      LL_parte3 = 0
      LL_parte4 = 0
      VeroSimActual=0
      
      LL_parte1 = -.5*T*log(2*pi)
      
      for (i in 1:length(S_training)) {#Calculo do primeiro segmento da LL
        LL_parte2 = LL_parte1 +.5*log(sigma_hat[S_treino[i]]) 
      }
      for (i in 1:length(S_training)) {#Calculo do segundo segmento da LL
        LL_parte3 = LL_parte3 +(1/(2*sigma_hat[S_treino[i]]))*((Y_training[i]-mu_hat[S_treino[i]])^2)
      }
      temp=NULL
      for (i in 2:length(S_training)) {#Calculo do terceiro segmento da LL
        for (g in 1:K) {
          temp[g]<-exp(X[i,]%*%matrix(BetaArray[g,,S_treino[i-1]],ncol=1))
        }
        LL_parte4 = LL_parte4 + (X[i,]%*%matrix(BetaArray[S_treino[i],,S_treino[i-1]]) - log(sum(temp), base = exp(1)))
      }
      VeroSimActual <- log(P0[S_treino[1]]) + LL_parte1 - (LL_parte2 + LL_parte3) + LL_parte4 #calculo da LogVerosim
      
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
      
      
      LL2_parte1 = 0
      LL2_parte2 = 0
      LL2_parte3 = 0
      LL2_parte4 = 0
      VeroSimProxima=0
      
      #Calculo da Verosimilhança como valor de tolerança
      LL2_parte1 = -.5*length(S_training)*log(2*pi)
      
      for (i in 1:length(S_training)) {#Calculo do primeiro segmento da LL
        LL2_parte2 = LL2_parte2 +.5*log(sigma_hat[S_treino[i]]) 
      }
      for (i in 1:length(S_training)) {#Calculo do segundo segmento da LL
        LL2_parte3 = LL2_parte3 +(1/(2*sigma_hat[S_treino[i]]))*((Y_training[i]-mu_hat[S_treino[i]])^2)
      }
      temp=NULL
      for (i in 2:length(S_training)) {#Calculo do terceiro segmento da LL
        for (g in 1:K) {
          temp[g]<-exp(X[i,]%*%matrix(BetaArray[g,,S_treino[i-1]],ncol=1))
        }
        LL2_parte4 = LL2_parte4 + (X[i,]%*%matrix(BetaArray[S_treino[i],,S_treino[i-1]]) - log(sum(temp), base = exp(1)))
      }
      VeroSimProxima <- log(P0[S_treino[1]]) + LL2_parte1 - (LL2_parte2 + LL2_parte3) + LL2_parte4 #calculo da LogVerosim
      
      val=val+1
      tolval[val]<-VeroSimProxima-VeroSimActual
      
     # print(tolval[val])
    }###fim da segunda rodada do EM Estocastico###
    
    Y_hat_validation = NULL
    S_hat_validation = NULL
    S_hat_validation[1]<-S_validation[1] #O valor para o primeiro estado oculto
    Y_hat_validation[1]<-rnorm(1,mu_hat[S_hat_validation[1]],sigma_hat[S_hat_validation[1]])# O valor para o primeiro valor observavel
    for (t in 2:length(Y_validation)){
      prob<-NULL
      for (i in 1:K) prob[i]<-exp(X_validation[t,]%*%matrix(BetaArray[i,,S_hat_validation[t-1]],ncol=1))
      prob<-prob/sum(prob)
      S_hat_validation[t]<-which.max(prob)
      Y_hat_validation[t]<-sum(prob * mu_hat)
    }
    
    Beta_Estimates <- NULL
    for (i in 2:K) {
      for (j in 1:K){
        for (d in 1:D){
          Beta_Estimates <- c(Beta_Estimates, BetaArray[i,d,j]) 
        }
      }
    }
    lasso_Y[lasso_iterator,] <- Y_hat_validation
    lasso_S[lasso_iterator,] <- S_hat_validation
    lasso_Beta_estimates[lasso_iterator,] <- Beta_Estimates
    lasso_mu_hat_estimates[lasso_iterator,] <- mu_hat
    lasso_sigma_hat_estimates[lasso_iterator,] <- sigma_hat
    lasso_RMSE[lasso_iterator] <- (sum((Y_hat_validation - Y_validation)^2))/length(Y_validation)
    lasso_Beta_arrays[,,,lasso_iterator] <- BetaArray
    lasso_iterator = lasso_iterator+1
  } ##################################################
  # FIM DO PROCESSO DE ESTIMAÇÃO (LASSO)
  
  # CAPTURA INDICE DO LAMBDA COM MELHORES RESULTADOS
  min_index = which.min(lasso_RMSE)
 
  # CAPTURA DE METRICAS PARA CADA REPLICA
  ##################################################################

  
  # COLETANDO VALORES NO CONJUNTO DE VALIDAÇÃO
  
  # Valor de Lambda optimo
  
  
  # Coletar valores estimados dos parâmetros das VA observaveis
  Mu_Rep[p,] <- lasso_mu_hat_estimates[min_index,]
  Sigma_Rep[p,] <- lasso_sigma_hat_estimates[min_index,]
  Best_Beta_Estimates[p,] <- lasso_Beta_estimates[min_index,]
  Best_Beta_Arrays[,,,p] <- lasso_Beta_arrays[,,,min_index]
  
  # Coletar o valor da melhor sequência S e Y no conjunto de Validação
  Best_S[p, ] <- lasso_S[min_index, ]
  Best_Y[p, ] <- lasso_Y[min_index, ]
  
  #Metricas de Performance Preditiva 
  MSPE_Validação[p] <- lasso_RMSE[min_index] #Mean Square Predictive Error para o melhor lambda
  
  #Metricas de Predição da sequência S
  for (i in 1:length(S_validation)) {
    if (S_validation[i] == Best_S[p,i]){
      Acertos_S_Validação[p] =  Acertos_S_Validação[p] + 1
    }
  }
  
  #Metricas de Zerado dos Coeficientes
  TP = NULL
  TN = NULL
  FP = NULL
  FN = NULL
  
  for (j in 1:(D*K*(K-1))){
    if ((Real[j]==0) && (abs(Best_Beta_Estimates[p,j])<=zero_threshold)){
      TP[j] = 1
    }
    if ((Real[j]!=0) && (abs(Best_Beta_Estimates[p,j])>zero_threshold)){
      TN[j]=1
    }
    if ((Real[j]!=0) && (abs(Best_Beta_Estimates[p,j])<=zero_threshold)){
      FP[j]=1
    }
    if ((Real[j]==0) && (abs(Best_Beta_Estimates[p,j])>zero_threshold)){
      FN[j]=1
    }
  }
  
  Pctgm_Zerado[p] = sum(abs(Best_Beta_Estimates[p,])<=zero_threshold)/(D*K*(K-1))
  Sensitivity[p] = sum(TP, na.rm = TRUE) / ( sum(TP, na.rm = TRUE) +  sum(FN, na.rm = TRUE) )
  Specificity[p] = sum(TN, na.rm = TRUE) / ( sum(TN, na.rm = TRUE) + sum(FP, na.rm = TRUE) )
  Accuracy[p] = (sum(TN, na.rm = TRUE) + sum(TP, na.rm = TRUE)) / ( sum(TP, na.rm = TRUE) + sum(TN, na.rm = TRUE) + sum(FP, na.rm = TRUE) + sum(FN, na.rm = TRUE))
  RMSE_Parameters[p] = sum((Real - Best_Beta_Estimates[p,])^2)
  
  
  # Train dataset for other models (Concatenation of train and validation)
  X_tr <- NULL
  Y_tr <- NULL
  
  X_tr <- rbind(X_training, X_validation)
  Y_tr <- c(Y_training, Y_validation)
  
  ############################################
  # GLMNET
  
 
  glmnet_mod <- cv.glmnet(X_tr, Y_tr)
  Y_hat_test_glmnet <- predict(glmnet_mod, newx = X_test, s = "lambda.min")
  MSPE_Teste_glmnet[p] <- sum((Y_test - Y_hat_test_glmnet)^2)/length(Y_test) 
  Y_hat_test_glmnet_DF[p,] <- Y_hat_test_glmnet
  #############################################
  
  ############################################
  # ARIMA
  
  arima_mod <- try(auto.arima(y=Y_tr, xreg = data.matrix(X_tr)))
  if (inherits(arima_mod,"try-error")){
    Y_hat_test_arima <- Y_hat_test_arima_DF[p-1,]
    Y_hat_test_arima_DF[p,] <- Y_hat_test_arima
  } else {
    Y_hat_test_arima <- forecast(arima_mod,xreg=data.matrix(X_test))  
    Y_hat_test_arima_DF[p,] <- Y_hat_test_arima$mean
  }
  
  MSPE_Teste_arima[p] <- sum((Y_test - Y_hat_test_arima_DF[p,])^2)/length(Y_test)
  ############################################
  
  
  ##########################################################
  #           AVALIAÇÃO NO CONJUNTO DE TESTE
  #--------------------------------------------------------#
  Y_hat_test <- NULL
  S_hat_test <- NULL
  
  S_hat_test[1]<-S_test[1] #O valor para o primeiro estado oculto
  Y_hat_test[1]<-rnorm(1,mu_hat[S_hat_test[1]],sigma_hat[S_hat_test[1]])# O valor para o primeiro valor observavel
  for (t in 2:length(Y_test)){
    prob<-NULL
    for (i in 1:K) prob[i]<-exp(X_test[t,]%*%matrix(Best_Beta_Arrays[i,,S_hat_test[t-1],p],ncol=1))
    prob<-prob/sum(prob)
    S_hat_test[t]<-which.max(prob)
    Y_hat_test[t]<-sum(prob * mu_hat)
  }
  
  MSPE_Teste[p] <- (sum((Y_hat_test - Y_test)^2))/length(Y_test)
  Y_test_DF[p,] <- Y_test
  Y_hat_test_NHMM_DF[p,] <- Y_hat_test
  for (i in 1:length(S_test)) {
    if (S_test[i] == S_hat_test[i]){
      Acertos_S_Teste[p] =  Acertos_S_Teste[p] + 1
    }
  }
  tempo_final<-Sys.time()
  Tempo[p] <- difftime(tempo_final, tempo_inicial, units = "secs")[[1]]/60
} #FIM DAS REPLICAS 



#########################################################
#  METRICAS DOS PARAMETROS DAS DISTRIBUÇÕES OBSERVAVEIS #
#########################################################
EQM_Mu=NULL
EQM_Sigma=NULL
Vies_Mu=NULL
Vies_Sigma=NULL
Mu_Finais_Mean=NULL
Sigma_Finais_Mean=NULL
Mu_Finais_Median=NULL
Sigma_Finais_Median=NULL
SD_Mu=NULL
SD_Sigma=NULL
AICc_Final=NULL
BIC_Final=NULL


#calculo do Valor EStimado Medio Final
for (i in 1:K) {
  Mu_Finais_Mean[i]=mean(Mu_Rep[,i])
  Sigma_Finais_Mean[i]=mean(Sigma_Rep[,i])
  Mu_Finais_Median[i]=median(Mu_Rep[,i])
  Sigma_Finais_Median[i]=median(Sigma_Rep[,i])
}

#Calculo do Desvio Padrão ao longo das R replicas
for (i in 1:K) {
  SD_Mu[i]=sd(Mu_Rep[,i])
  SD_Sigma[i]=sd(Sigma_Rep[,i])
}

#Calculo do Vies 
for (i in 1:K) {
  Vies_Mu[i]=mu[i]-Mu_Finais_Mean[i]
  Vies_Sigma[i]=sigma[i]-Sigma_Finais_Mean[i]
}

#Calculo do Erro Quadratico Medio dos parametros
for (i in 1:K) {
  EQM_Mu[i]=Vies_Mu[i]^2 + SD_Mu[i]^2
  EQM_Sigma[i]=Vies_Sigma[i]^2 + SD_Sigma[i]^2
}

#calculo dos quanties para os parametros
quantiles<-matrix(nrow = 6, ncol = 2)
quantiles[1,]<-quantile(Mu_Rep[,1],probs = c(0.025, 0.975))
quantiles[2,]<-quantile(Mu_Rep[,2],probs = c(0.025, 0.975))
quantiles[3,]<-quantile(Mu_Rep[,3],probs = c(0.025, 0.975))
quantiles[4,]<-quantile(Sigma_Rep[,1],probs = c(0.025, 0.975))
quantiles[5,]<-quantile(Sigma_Rep[,2],probs = c(0.025, 0.975))
quantiles[6,]<-quantile(Sigma_Rep[,3],probs = c(0.025, 0.975))

#Juntamos as Metricas dos parametros das distribuições observaveis num dataframe
Parametro<-c("Mu_1","Mu_2","Mu_3", "Sigma_1", "Sigma_2", "Sigma_3")
Real<-c(mu[1],mu[2],mu[3], sigma[1],sigma[2], sigma[3])
Estimado_Mean<-c(Mu_Finais_Mean[1],Mu_Finais_Mean[2],Mu_Finais_Mean[3], Sigma_Finais_Mean[1], Sigma_Finais_Mean[2], Sigma_Finais_Mean[3])
Estimado_Median<-c(Mu_Finais_Median[1],Mu_Finais_Median[2],Mu_Finais_Median[3], Sigma_Finais_Median[1], Sigma_Finais_Median[2], Sigma_Finais_Median[3])
SD<-c(SD_Mu[1],SD_Mu[2], SD_Mu[3], SD_Sigma[1],SD_Sigma[2],SD_Sigma[3])
Vies<-c(Vies_Mu[1],Vies_Mu[2],Vies_Mu[3], Vies_Sigma[1], Vies_Sigma[2],Vies_Sigma[3])
EQM<-c(EQM_Mu[1],EQM_Mu[2],EQM_Mu[3], EQM_Sigma[1], EQM_Sigma[2],EQM_Sigma[3])
Assimetria<-c(skewness(Mu_Rep[,1]), skewness(Mu_Rep[,2]),skewness(Mu_Rep[,3]), skewness(Sigma_Rep[,1]), skewness(Sigma_Rep[,2]), skewness(Sigma_Rep[,3]))
IC_95<-c(toString(quantiles[1,]),toString(quantiles[2,]),toString(quantiles[3,]),toString(quantiles[4,]),toString(quantiles[5,]),toString(quantiles[6,]))

#Criamos o DAtaframe
df1<-data.frame(Parametro,Real,Estimado_Mean, Estimado_Median,SD,Vies,EQM,Assimetria,IC_95)



#############################################################
## METRICAS PARA COMPARAR/AVALIAR A QUALIDADE DA ESTIMAÇÃO ##
#############################################################
EQM_Validação_Median_Final = NULL
EQM_Validação_Mean_Final = NULL
EQM_Validação_SD = NULL

EQM_Test_Median_Final = NULL
EQM_Test_Mean_Final = NULL
EQM_Test_SD = NULL

EQM_Test_Median_glmnet = NULL
EQM_Test_Mean_glmnet = NULL
EQM_Test_SD_glmnet = NULL

EQM_Test_Median_arima = NULL
EQM_Test_Mean_arima = NULL
EQM_Test_SD_arima = NULL

EQM_Validação_Median_Final <- median(MSPE_Validação)
EQM_Validação_Mean_Final <- mean(MSPE_Validação)
EQM_Validação_SD <- sd(MSPE_Validação)

EQM_Test_Median_Final <- median(MSPE_Teste)
EQM_Test_Mean_Final = mean(MSPE_Teste)
EQM_Test_SD = sd(MSPE_Teste)

EQM_Test_Median_arima = median(MSPE_Teste_arima)
EQM_Test_Mean_arima = mean(MSPE_Teste_arima)
EQM_Test_SD_arima = sd(MSPE_Teste_arima)

EQM_Test_Median_glmnet = median(MSPE_Teste_glmnet)
EQM_Test_Mean_glmnet = mean(MSPE_Teste_glmnet)
EQM_Test_SD_glmnet = sd(MSPE_Teste_glmnet)


#Criamos o Dataframe
Chain_Length <- paste("T =",toString(T))
df2<-data.frame(Chain_Length, EQM_Validação_Median_Final, EQM_Validação_Mean_Final, EQM_Validação_SD, EQM_Test_Median_Final, EQM_Test_Mean_Final, EQM_Test_SD, EQM_Test_Median_arima, EQM_Test_Mean_arima, EQM_Test_SD_arima, EQM_Test_Median_glmnet, EQM_Test_Mean_glmnet, EQM_Test_SD_glmnet)


#####################################################################
## METRICAS PARA COMPARAR/AVALIAR A QUANTO ACERTAMOS NA SEQUÊNCIA S ##
#####################################################################
Acertos_Validação_Mean_Final <- mean(Acertos_S_Validação/length(S_validation))
Acertos_Validação_Median_Final <- median(Acertos_S_Validação/length(S_validation))
Acertos_Validação_SD <- sd(Acertos_S_Validação/length(S_validation))

Acertos_Teste_Mean_Final <- mean(Acertos_S_Teste/length(S_test))
Acertos_Teste_Median_Final <- median(Acertos_S_Teste/length(S_test))
Acertos_Teste_SD <- sd(Acertos_S_Teste/length(S_test))

Chain_Length <- paste("T =",toString(T))
df3<-data.frame(Chain_Length, Acertos_Validação_Mean_Final, Acertos_Validação_Median_Final, Acertos_Validação_SD,  Acertos_Teste_Mean_Final, Acertos_Teste_Median_Final, Acertos_Teste_SD)

#############################################################
### Calculo das metricas de Shrinkage (Zerou ou Não Zerou?) #
#############################################################
Specificity_Final_Median = NULL
Specificity_Final = NULL
Specificity_SD = NULL

Sensitivity_Final_Median = NULL
Sensitivity_Final = NULL
Sensitivity_SD = NULL

Accuracy_Final_Median = NULL
Accuracy_Final = NULL
Accuracy_SD = NULL

#Especificidade
Specificity_Final = mean(Specificity)
Specificity_Final_Median = median(Specificity)
Specificity_SD = sd(Specificity)

#Sensitividade
Sensitivity_Final = mean(Sensitivity)
Sensitivity_Final_Median = median(Sensitivity)
Sensitivity_SD = sd(Sensitivity)

#Acuracia
Accuracy_Final = mean(Accuracy)
Accuracy_Final_Median = median(Accuracy)
Accuracy_SD = sd(Accuracy)

Chain_Length <- paste("T =",toString(T))
df4<-data.frame(Chain_Length, Specificity_Final_Median, Specificity_Final,Specificity_SD, Sensitivity_Final_Median, Sensitivity_Final,Sensitivity_SD, Accuracy_Final_Median, Accuracy_Final,Accuracy_SD)

###################################################
#     Tempos de Execução  
###################################################
Tempo_total <- sum(Tempo)
Tempo_medio <- mean(Tempo)
tempo_median <- median(Tempo)
Sd_tempo <- sd(Tempo)

df5 <- data.frame(Tempo_total, Tempo_medio, tempo_median, Sd_tempo)

#####################################################
# Quantidade de Transições por Replica
#####################################################
Final_TransCount <- matrix(nrow = K, ncol = K)
Final_TransCount_sd <- matrix(nrow = K, ncol = K)
for (j in 1:K) {
  for (k in 1:K) {
    Final_TransCount[k,j] <- round(mean(RealTransCount[k,j,])/length(S)*100,2)
    Final_TransCount_sd[k,j] <- round(sd(RealTransCount[k,j,])/length(S)*100,2)
  }
}

write.csv(df1, paste("1_ObservableParams_T-",toString(T),"_D-",toString(D),"_zerothreshold-",toString(zero_threshold),"_optimethod-",toString(optim_algo),".csv", sep = ""), row.names=FALSE)
write.csv(df2, paste("2_MSPE_T-",toString(T),"_D-",toString(D),"_zerothreshold-",toString(zero_threshold),"_optimethod-",toString(optim_algo),".csv", sep = ""), row.names=FALSE)
write.csv(df3, paste("3_AcertosSequênciaOculta_T-",toString(T),"_D-",toString(D),"_zerothreshold-",toString(zero_threshold),"_optimethod-",toString(optim_algo),".csv", sep = ""), row.names=FALSE)
write.csv(df4, paste("4_Shrinkage_T-",toString(T),"_D-",toString(D),"_zerothreshold-",toString(zero_threshold),"_optimethod-",toString(optim_algo),".csv", sep = ""), row.names=FALSE)
write.csv(df5, paste("5_ExecutionTime_T-",toString(T),"_D-",toString(D),"_zerothreshold-",toString(zero_threshold),"_optimethod-",toString(optim_algo),".csv", sep = ""), row.names=FALSE)


write.csv(data.frame(Mu_Rep) , paste("Mu_hat_Replicas_T-",toString(T),"_D-",toString(D),"_zerothreshold-",toString(zero_threshold),"_optimethod-",toString(optim_algo),".csv", sep = ""), row.names=FALSE)
write.csv(data.frame(Sigma_Rep), paste("Sigma_hat_Replicas_T-",toString(T),"_D-",toString(D),"_zerothreshold-",toString(zero_threshold),"_optimethod-",toString(optim_algo),".csv", sep = ""), row.names=FALSE)
write.csv(data.frame(Best_Beta_Estimates), paste("Best_Betas_T-",toString(T),"_D-",toString(D),"_zerothreshold-",toString(zero_threshold),"_optimethod-",toString(optim_algo),".csv", sep = ""), row.names=FALSE)
write.csv(data.frame(Best_Y), paste("Y_BestValidation_T-",toString(T),"_D-",toString(D),"_zerothreshold-",toString(zero_threshold),"_optimethod-",toString(optim_algo),".csv", sep = ""), row.names=FALSE)
write.csv(data.frame(Best_S), paste("S_BestValidation_T-",toString(T),"_D-",toString(D),"_zerothreshold-",toString(zero_threshold),"_optimethod-",toString(optim_algo),".csv", sep = ""), row.names=FALSE)

write.csv(data.frame(round(Y_hat_test_arima_DF,2)), paste("Y_Test_hat_ARIMA_T-",toString(T),"_D-",toString(D),"_zerothreshold-",toString(zero_threshold),"_optimethod-",toString(optim_algo),".csv", sep = ""), row.names=FALSE)
write.csv(data.frame(round(Y_hat_test_glmnet_DF,2)), paste("Y_Test_hat_GLMNET_T-",toString(T),"_D-",toString(D),"_zerothreshold-",toString(zero_threshold),"_optimethod-",toString(optim_algo),".csv", sep = ""), row.names=FALSE)
write.csv(data.frame(round(Y_hat_test_NHMM_DF,2)), paste("Y_Test_hat_NHMM_T-",toString(T),"_D-",toString(D),"_zerothreshold-",toString(zero_threshold),"_optimethod-",toString(optim_algo),".csv", sep = ""), row.names=FALSE)
write.csv(data.frame(round(Y_test_DF,2)), paste("Real_Y_Test_T-",toString(T),"_D-",toString(D),"_zerothreshold-",toString(zero_threshold),"_optimethod-",toString(optim_algo),".csv", sep = ""), row.names=FALSE)

write.csv(data.frame(Final_TransCount), paste("Real_TransitionsCount_T-",toString(T),"_D-",toString(D),"_zerothreshold-",toString(zero_threshold),"_optimethod-",toString(optim_algo),".csv", sep = ""), row.names=FALSE)
write.csv(data.frame(Final_TransCount_sd), paste("Real_Transitions_SD_T-",toString(T),"_D-",toString(D),"_zerothreshold-",toString(zero_threshold),"_optimethod-",toString(optim_algo),".csv", sep = ""), row.names=FALSE)


df1
df2
df3
df4
df5

# plot(Y_hat_test_NHMM_DF[14,],type = 'l')
# plot(Y_hat_test_arima_DF[14,],type = 'l')
# plot(Y_test_DF[14,], type = 'l')
# MSPE_Teste
# Tempo
# Sigma_Rep
