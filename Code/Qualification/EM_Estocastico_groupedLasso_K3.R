library('label.switching')
library("ggplot2")
#library("r2excel")
#library("xlsx")
library(e1071)
library(EnvStats)
library("Matrix")

get.lambda.l1  <-  function(x,y,nlambda) {
  xs     <-  as(x, "dgCMatrix")
  ## currently not robust to missing values in xs or y
  ybar <- mean(y,na.rm=TRUE)
  xbar <- colMeans(xs,na.rm=TRUE)
  x      <-  list(Xi = xs@i, Xj = xs@p, Xnp = diff(xs@p), Xx = xs@x)
  xty    <-  drop(crossprod(y-ybar,scale(xs,xbar,FALSE)))
  lmax  <-  max(abs(xty))
  return(10^seq(log10(lmax), log(1), len=nlambda))
} 



################################################
########           FUNÇÕES            ########## 
################################################

#####Função para gerar valores uniformes discretos
rDiscreta<-function(v){
  u<-runif(1)
  P<-cumsum(v)
  val<-sum(P<u)+1
  return(val)}


#função para calcular Matriz de Transição em cada instante de tempo
# Mat_trans <-function(covar){
#   B = matrix(nrow=K, ncol=K)
#   for (j in 1:K) {
#     for (i in 1:K){
#       B[i,j] = exp(covar%*%BetaArray[i,,j])/(exp(covar%*%BetaArray[1,,j])+exp(covar%*%BetaArray[2,,j]))+exp(covar%*%BetaArray[3,,j]))
#     }  
#   }
#   return(B)
# }

#função para calcular Matriz de Transição em cada instante de tempo
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


T=600 #Cumprimento da cadeia simulada
K=3   #Numero de estados ocultos
D=6   #Quantidade de Covariaveis
R<-10 #Numero de Replicas que serão executadas
tol<-0.0001 #Nivel de tolerancia que estabelecemos como criterio de parada do EM Est

##Tamanhos das bases
train_size = 0.8
test_size = 0.2

Y_hat_teste_minEQM_Rep <- matrix(nrow = R, ncol = T*test_size)
Y_teste_Rep <- matrix(nrow = R, ncol = T*test_size)

#####################################################
## Criação de Directorio para salvar os Resultados ##
#####################################################
direc = dirname(rstudioapi::getActiveDocumentContext()$path)
folder_name = paste("/",toString(D),"Bcovs_Output_EM-Estocastico_GroupedLASSO_","K",toString(K),"_","T",toString(T), sep = "")
path_completo = paste(direc, folder_name,sep = "")
ifelse(!dir.exists(path_completo), dir.create(path_completo), FALSE)

arqname<-c(path_completo,folder_name,".xlsx")
nome_arquivo<-paste(arqname,collapse = "") 
nomefolha<-paste("Saida T=",toString(T))
header<-paste("Resultados para K=", toString(K), "e T=",toString(T), collapse = "")




set.seed(3)
seeds <-sample(110000000,R) # Seed number para conseguer os mesmos valores simulados
options(digits=6)
options(scipen=999)


#Inicializamos algumas estruturas para
#almacenar valores importantes ao longo das
#"R" replicas que serão feitas
Tempo = NULL
AICc_Acum<-NULL
BIC_Acum<-NULL
Betas_Rep<-matrix(nrow = R, ncol = D)
mu_Rep<-matrix(nrow = R, ncol = K)
sigma_Rep<-matrix(nrow = R, ncol = K)
acertos<-rep.int(0,R)

EQM_mean_LASSO_Rep = NULL
SD_EQM_LASSO_Rep = NULL
EQM_min_LASSO_Rep = NULL
EQM_No_LASSO_Rep = NULL

Acertos_teste_min_Rep = NULL
Acertos_teste_mean_Rep = NULL
SD_Acertos_teste_Rep = NULL
Acertos_teste_NoLasso_Rep = NULL

specificity <-rep.int(0,R)
sensitivity <-rep.int(0,R)
accuracy <-rep.int(0,R)

S_Table<-matrix(nrow = T, ncol = R)
S_Post_Table<-matrix(nrow = T, ncol = R)
CritInfo_Table<-matrix(nrow = 2, ncol = R)

P0 = rep(1/K,K) #Inicializamos vetor de probabilidades inciais para o HMM
mu = c(20,150,250) # vetor com media para as duas Normais
sigma = c(1,2,3) #Vetor com os desvios padrões para as duas normais
Betas=array(dim=c(K,D,K)) # valores de beta utilizados na geração dos valores (consideranda intercepto e duas covariáveis)


#Fazemos o valor iniciais dos Betas especificos igual a 0, para ter a função
#de Ligação mlogit
for (i in 1:K){
  for (j in 1:D){
    Betas[1,j,i] = 0
  }
}


# # Betas transição 2
Betas[2,1,1]=-1.5
Betas[2,2,1]=1.5
Betas[2,3,1]=-2.6

Betas[2,1,2]=-2.0
Betas[2,2,2]=4.6
Betas[2,3,2]=1.4

Betas[2,1,3]=2.4
Betas[2,2,3]=1.1
Betas[2,3,3]=-1.5

# Betas transição 3
Betas[3,1,1]=1.3
Betas[3,2,1]=3.2
Betas[3,3,1]=2.4

Betas[3,1,2]=-3.3
Betas[3,2,2]=1.7
Betas[3,3,2]=1.3

Betas[3,1,3]=-2.3
Betas[3,2,3]=-2.7
Betas[3,3,3]=-3.5



qtd_nonzero_coefs = 3
for (i in 1:K){
  for (j in 2:K){
    for (d in (qtd_nonzero_coefs+1):D){
      Betas[j,d,i] = 0
    }
  }
}

########Inicializamos estruturas para calcular specificidade e sensitividade
is_non_zero_real = array(dim=c(K,D,K))

#Marcamos quais coeficientes realmente são diferentes de 0
for (i in 1:K){
  for (j in 2:K)
    for (d in 1:D){
      if(Betas[j,d,i] != 0){
        is_non_zero_real[j,d,i] = 1
      } else {
        is_non_zero_real[j,d,i] = 0
      }
    }
}


T1 = proc.time()
#Inicio das Replicas
for (p in 1:R){
  tempo_inicial<-proc.time()#Calcularemos quanto demoro todo o proceso
  set.seed(seeds[p]) 
  cat(paste("\nNumero de Replica:",toString(p),"\n", collapse = "")) #Messagem indicando o numero da replica atual
  
  TransCount <- matrix(data = c(rep(0,K^2)), nrow = K, ncol = K) #matriz que almacena a quantidade de vezes que cada transição ocorre
  mu_hat = NULL #Variavel para estimar os mus em cada iteração do EM Estocástico
  sigma_hat = NULL #Variavel para estimar os sigmas em cada iteração do EM Estocástico
  BetaArray = array(c(rep(0,(K*D*K))), dim = c(K,D,K))
  S_treino<-NULL #Inicializamos a sequência S de treinamento
  tolval=NULL #Inicializamos a variavel que almacenara a diferença nas verosim no paso i e i+1
  tolval[1]=1 #Para garantir que a primeira iteração do EMEst rode, atribuimos 1 a esse valor
  
  
  print("Simulando Dados...")
  ########################################
  ###        SIMULAÇÃO DE DADOS        ###
  ########################################
  S<-NULL #Inicializamos a sequência de estados não observaveis
  Y<-NULL #Inicializamos a sequência de valores observaveis
  X<-NULL #Inicializamos o vetor de covariaveis
  
  #Geração de covariáveis
  X <- matrix(rnorm(D*T), nrow = T, ncol = D)
  X[,1] = 1
  
  ### Finalmente, Aqui simulamos um NHMM com os 
  ### parâmetros de simulação anteriormente definidos
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

  ###########################
  #So pra enxergar as caracteristicas da sequencia S
  #E avaliar se ela é facil ou dificil de estimar
  S_TransCount <- matrix(data = c(rep(0,K^2)), nrow = K, ncol = K)
  for (i in 2:length(S)) {
    for (j in 1:K) {
      for (k in 1:K) {
        if (S[i]==j && S[i-1]==k)
          S_TransCount[k,j]=S_TransCount[k,j]+1
      }
    }
  }
  S_TransCount

  #Calcula os indices de corte
  cutoff_treino = length(Y)*train_size
  
  
  
  #Cria as bases 
  S_training = S[0:cutoff_treino] 
  Y_training = Y[0:cutoff_treino]
  X_training = X[0:cutoff_treino, ]
  
  S_teste = S[(cutoff_treino+1):T] 
  Y_teste = Y[(cutoff_treino+1):T]
  X_teste = X[(cutoff_treino+1):T, ]
  
  
  #########################################
  #####   Procedimento de Estimação   #####
  #########################################
  
  # geramos uma sequência não observavel de treinamento
  P_Treino=rep(1/K,K) #Vetor de probabilidade utilizadas para gerar a sequência de treino
  S_treino<-NULL # Inicializamos a sequência oculta de treinamento
  
  #Geramos uma sequência de treinamento
  for (i in 1:length(S_training)) {
    S_treino[i] = rDiscreta(P_Treino)
  }
  
  #Valores iniciais para a estimação dos Betas
  init1 = c(c(rnorm(D*2, 3, 3)))#Valores iniciais para os Betas_1 
  init2 = c(c(rnorm(D*2, 3, 3)))#Valores iniciais para os Betas_2
  init3 = c(c(rnorm(D*2, 3, 3)))#Valores iniciais para os Betas_3
  
  
  #Antes de Rodar o EM Estocástico, precisaremos de algumas estruturas para
  #almacenar os valores de mu_hat, sigma_hat, S_treino
  #lambdas = seq(log(exp(0.001)), log(10000000), len=19)
  #lambdas = append(0,lambdas)
  lambdas <- c(0, 0.0001, 0.001, 0.01, 0.1, 0.25, 0.5, 1, 2.5, 5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 30, 50, 75, 100, 200)
  
  #lambdas = get.lambda.l1(X_training,Y_training,50)
  #lambdas <- c(0.01, 0.05, 0.1, 0.15, 0.5, 0.75, 1, 1.5, 2, 4, 5, 10)
  mu_teste_vector = matrix(nrow = length(lambdas)^K, ncol = K)
  sigma_teste_vector = matrix(nrow = length(lambdas)^K, ncol = K)
  Beta21_teste_vector = matrix(nrow = length(lambdas)^K, ncol = D)
  Beta22_teste_vector = matrix(nrow = length(lambdas)^K, ncol = D)
  Beta23_teste_vector = matrix(nrow = length(lambdas)^K, ncol = D)
  Beta31_teste_vector = matrix(nrow = length(lambdas)^K, ncol = D)
  Beta32_teste_vector = matrix(nrow = length(lambdas)^K, ncol = D)
  Beta33_teste_vector = matrix(nrow = length(lambdas)^K, ncol = D)
  Y_hat_teste_vector = matrix(nrow = length(lambdas)^K, ncol = T*test_size)
  acertos_teste <- rep.int(0,length(lambdas)^K)
  EQM_Teste = NULL
  # O EMEst sera rodado em duas etapas. Uma para estimar 
  # S e theta, e outra  para estimar os Betas
  loop_count = 0
  print('Running LASSO')
  for (lasso_iter3 in 1:length(lambdas)){
    for (lasso_iter2 in 1:length(lambdas)){
      for (lasso_iter1 in 1:length(lambdas)){
        message('\r', paste("Lasso iteration # = ",toString(loop_count)), appendLF = FALSE)
        FSM1 <-function(params){#função a maximizar para achar os Betas_1
          resp <- (sum(1 - log(1 + exp(Xtemp11%*%params[1:D]) + exp(Xtemp11%*%params[(D+1):(2*D)]))) + sum(Xtemp12%*%params[1:D] - log(1 + exp(Xtemp12%*%params[1:D]) + exp(Xtemp12%*%params[(D+1):(2*D)]))) + sum(Xtemp13%*%params[(D+1):(2*D)] - log(1 + exp(Xtemp13%*%params[1:D]) + exp(Xtemp13%*%params[(D+1):(2*D)])))) - (lambdas[lasso_iter3]*sum(abs(params[2:D])) + lambdas[lasso_iter3]*sum(abs(params[(D+2):(2*D)]))) 
        }
        
        FSM2 <-function(params){#função a maximizar para achar os Betas_2
          resp <- (sum(1 - log(1 + exp(Xtemp21%*%params[1:D]) + exp(Xtemp21%*%params[(D+1):(2*D)]))) + sum(Xtemp22%*%params[1:D] - log(1 + exp(Xtemp22%*%params[1:D]) + exp(Xtemp22%*%params[(D+1):(2*D)]))) + sum(Xtemp23%*%params[(D+1):(2*D)] - log(1 + exp(Xtemp23%*%params[1:D]) + exp(Xtemp23%*%params[(D+1):(2*D)])))) - (lambdas[lasso_iter2]*sum(abs(params[2:D])) + lambdas[lasso_iter2]*sum(abs(params[(D+2):(2*D)]))) 
        }
        
        FSM3 <-function(params){#função a maximizar para achar os Betas_2
          resp <- (sum(1 - log(1 + exp(Xtemp31%*%params[1:D]) + exp(Xtemp31%*%params[(D+1):(2*D)]))) + sum(Xtemp32%*%params[1:D] - log(1 + exp(Xtemp32%*%params[1:D]) + exp(Xtemp32%*%params[(D+1):(2*D)]))) + sum(Xtemp33%*%params[(D+1):(2*D)] - log(1 + exp(Xtemp33%*%params[1:D]) + exp(Xtemp33%*%params[(D+1):(2*D)])))) - (lambdas[lasso_iter1]*sum(abs(params[2:D])) + lambdas[lasso_iter1]*sum(abs(params[(D+2):(2*D)]))) 
        }
        
        val=1 # Esse é um contador para controlar o calculo do tolval em cada iteração
        #Agora executamos a primeira etapa do EM Estocástico
        while (tolval[val]>tol){
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
          S_treino[1]=which.max(dnorm(Y_training[1], mu_hat, sigma_hat))
          for (i in 2:length(S_training)) {
            A_hat_t = Mat_trans(X[i,])
            prob<-(A_hat_t[S_treino[i], ]*dnorm(Y_training[i], mu_hat, sigma_hat))/sum(A_hat_t[S_treino[i], ]*dnorm(Y_training[i], mu_hat, sigma_hat))
            S_treino[i]=which.max(prob)
          }
          
          
          #####################################
          #Este segmento de codigo testa se aconteceram todas as transições possiveis
          #No caso que elas não tinham acontecido, as que
          #não aconteceram são forçadas a acontecer
          TransCount <- matrix(data = c(rep(0,K^2)), nrow = K, ncol = K)
          for (i in 2:length(S_training)) {
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
                positions = sample(2:length(S_training), 2)
                for (d in 1:2) {
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
          fit2 <- optim(par = init2, fn = FSM2, control = list(fnscale=-1), method = "BFGS", hessian = FALSE)
          fit3 <- optim(par = init3, fn = FSM3, control = list(fnscale=-1), method = "BFGS", hessian = FALSE)
          
          for (i in 1:K){
            for (d in 1:D){
              if (i == 1){
                BetaArray[i,d,1]=0
              } else {
                BetaArray[i,d,1]=fit1$par[d]
              }
              
            }
          }
          
          for (i in 1:K){
            for (d in 1:D){
              if (i == 1){
                BetaArray[i,d,2]=0
              } else {
                BetaArray[i,d,2]=fit2$par[d]
              }
              
            }
          }
          
          for (i in 1:K){
            for (d in 1:D){
              if (i == 1){
                BetaArray[i,d,3]=0
              } else {
                BetaArray[i,d,3]=fit2$par[d]
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
          tolval[val]<-VeroSimProxima - VeroSimActual
        }#########Fim da primeira rodada do EM Estocastico#########
        
        # Criar algumas matrizes para fazer calculos e manipular 
        # a saida MCMC nestas matrizes, as estimativas serão 
        # reordenadas usando o metodo ECR
        
        mat_thetar<-matrix(nrow = 1, ncol = K)
        reorder_S<-matrix(nrow = 1, ncol = length(S_training))
        mat_S<-matrix(nrow = 1, ncol = length(S_training))
        mat_S[1,]<-S_treino
        zpvt_S = S #Como pivot para o metodo ECR usamos o S original
        perms_S = ecr(zpivot = zpvt_S, z = mat_S, K = 3)# aplicamos o metodo ECR que retornara as permutações das dos estados ocultos que devem ser utilizadas para reordenar a saida do algoritmo bayesiano
        
        for (i in 1:1) {
          for (j in 1:length(S_training)) {
            if(S_treino[j]!=S[j] && ((perms_S$permutations[i,1]==2 && perms_S$permutations[i,2]==3 && perms_S$permutations[i,3]==1) | (perms_S$permutations[i,1]==3 && perms_S$permutations[i,2]==1 && perms_S$permutations[i,3]==2))){
              S_treino[j]=perms_S$permutations[i,perms_S$permutations[i,S_treino[j]]]
            }
            
            else {
              S_treino[j]=perms_S$permutations[i,perms_S$permutations[i,S[j]]]
            }
          }
          mu_hat<-mu_hat[perms_S$permutations[i,]]
          sigma_hat<-sigma_hat[perms_S$permutations[i,]]
        }
        
        VeroSimProxima=1
        VeroSimActual=0
        val=1
        tolval[1]=1
        
        ###############################################
        #É possivel que na rerotulagem tinha ficado alguma transição com 0 ocorrencias então caso
        #Seja necessario, forçamos essa transição
        TransCount <- matrix(data = c(rep(0,K^2)), nrow = K, ncol = K)
        for (i in 2:length(S_training)) {
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
              positions = sample(2:length(S_training), 2)
              for (d in 1:2) {
                S_treino[positions[d]]=j
                S_treino[positions[d]-1]=k
              }
            }
          }
        }
        
        
        ##############################################
        
        
        ########Inicia segunda rodada do EM Estocastico#########
        while (tolval[val]>tol) {
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
          
          ###############################################
          #É possivel que na rerotulagem tinha ficado alguma transição com 0 ocorrencias então caso
          #Seja necessario, forçamos essa transição
          TransCount <- matrix(data = c(rep(0,K^2)), nrow = K, ncol = K)
          for (i in 2:length(S_training)) {
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
                positions = sample(2:length(S_training), 3)
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
          fit2 <- optim(par = init2, fn = FSM2, control = list(fnscale=-1), method = "BFGS", hessian = FALSE)
          fit3 <- optim(par = init3, fn = FSM3, control = list(fnscale=-1), method = "BFGS", hessian = FALSE)
          
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
          
          tolval[val]<-VeroSimProxima-VeroSimActual
          val=val+1
          
        }###fim da segunda rodada do EM Estocastico##
        
        
        mu_teste_vector[loop_count,] = mu_hat
        sigma_teste_vector[loop_count, ] = sigma_hat
        Beta21_teste_vector[loop_count, ] = BetaArray[2,,1]
        Beta22_teste_vector[loop_count, ] = BetaArray[2,,2]
        Beta23_teste_vector[loop_count, ] = BetaArray[2,,3]
        Beta31_teste_vector[loop_count, ] = BetaArray[3,,1]
        Beta32_teste_vector[loop_count, ] = BetaArray[3,,2]
        Beta33_teste_vector[loop_count, ] = BetaArray[3,,3]
        
        T_teste = length(S_teste)
        Y_hat_teste = NULL
        S_hat_teste = NULL
        S_hat_teste[1]<-S_teste[1] #O valor para o primeiro estado oculto
        Y_hat_teste[1]<-rnorm(1,mu_hat[S_hat_teste[1]],sigma_hat[S_hat_teste[1]])# O valor para o primeiro valor observavel
        for (t in 2:T_teste){
          prob<-NULL
          for (i in 1:K) prob[i]<-exp(X_teste[t,]%*%matrix(BetaArray[i,,S_hat_teste[t-1]],ncol=1))
          prob<-prob/sum(prob)
          S_hat_teste[t]<-which.max(prob)
          #S_hat_teste[t]<-rDiscreta(prob)
          Y_hat_teste[t] = sum(prob * mu_hat)
        }
        Y_hat_teste_vector[loop_count,] = Y_hat_teste
        for (i in 1:length(S_teste)) {
          if (S_hat_teste[i]==S_teste[i]){
            acertos_teste[loop_count] = acertos_teste[loop_count] + 1
          }
        }
        EQM_Teste[loop_count] = sum((Y_hat_teste - Y_teste)^2)/length(S_teste)
        loop_count = loop_count + 1
      }
    }
  }
  acertos_teste[1:length(lambdas)^3]/length(S_teste)
  ############################################################
  #Usamos a Base de Teste para escolher calcular or erro quadratico
  #medio na base de teste e desta forma determinaremos qual é o melhor
  #valor de lambda
  ############################################################
  Betas_teste21 = NULL
  Betas_teste22 = NULL
  Betas_teste23 = NULL
  Betas_teste31 = NULL
  Betas_teste32 = NULL
  Betas_teste33 = NULL
  acertos_S_teste_min = NULL
  is_non_zero_teste = array(dim=c(K,D,K))
  index_min = which.min(EQM_Teste)
  acertos_S_min_teste = acertos_teste[index_min]/length(Y_teste)
  Betas_teste21 = Beta21_teste_vector[index_min,]
  Betas_teste22 = Beta22_teste_vector[index_min,]
  Betas_teste23 = Beta23_teste_vector[index_min,]
  Betas_teste31 = Beta31_teste_vector[index_min,]
  Betas_teste32 = Beta32_teste_vector[index_min,]
  Betas_teste33 = Beta33_teste_vector[index_min,]
  
  
  #Marcamos quais coeficientes foram zerados pelo LASSO
  for (d in 1:D){
    if(abs(Betas_teste21[d]) >= 0.02){
      is_non_zero_teste[2,d,1] = 1
    } else {
      is_non_zero_teste[2,d,1] = 0
    }
  }
  
  for (d in 1:D){
    if(abs(Betas_teste22[d]) >= 0.02){
      is_non_zero_teste[2,d,2] = 1
    } else {
      is_non_zero_teste[2,d,2] = 0
    }
  }
  
  for (d in 1:D){
    if(abs(Betas_teste23[d]) >= 0.02){
      is_non_zero_teste[2,d,3] = 1
    } else {
      is_non_zero_teste[2,d,3] = 0
    }
  }
  
  for (d in 1:D){
    if(abs(Betas_teste31[d]) >= 0.02){
      is_non_zero_teste[3,d,1] = 1
    } else {
      is_non_zero_teste[3,d,1] = 0
    }
  }
  
  for (d in 1:D){
    if(abs(Betas_teste32[d]) >= 0.02){
      is_non_zero_teste[3,d,2] = 1
    } else {
      is_non_zero_teste[3,d,2] = 0
    }
  }
  
  for (d in 1:D){
    if(abs(Betas_teste33[d]) >= 0.02){
      is_non_zero_teste[3,d,3] = 1
    } else {
      is_non_zero_teste[3,d,3] = 0
    }
  }
  
  
  TP = 0
  TN = 0
  FN = 0
  FP = 0
  
  #Calculamos as quantidade de TP, TN, FP e FN
  for (i in 1:K){
    for (j in 2:K)
      for (d in 1:D){
        if((is_non_zero_teste[j,d,i] == 1) && (is_non_zero_real[j,d,i]==1)){
          TP = TP + 1
        } else if ((is_non_zero_teste[j,d,i] == 0) && (is_non_zero_real[j,d,i]==0)) {
          TN = TN + 1
        } else if ((is_non_zero_teste[j,d,i] == 0) && (is_non_zero_real[j,d,i]==1)) {
          FN = FN + 1
        } else if ((is_non_zero_teste[j,d,i] == 1) && (is_non_zero_real[j,d,i]==0)) {
          FP = FP + 1
        }
      }
  }
  
  #Qual foi a sequencia Y que teve o menor EQM
  Y_hat_teste_minEQM_Rep[p,] = Y_hat_teste_vector[index_min,]
  sum((Y_hat_teste_minEQM_Rep[1,] - Y_teste)^2)/length(Y_teste)
  #Calculamos as metricas de classificação
  accuracy[p] = (TN + TP)/(TP + FP +TN +FN)
  specificity[p] = TN / (TN + FP)
  sensitivity[p] = TP / (TP + FN)
  
  
  #############
  # Agora almacenaremos todas as quantidades de interesse 
  # em cada Replica, produtos do processo de estimação
  
  #Calculando valores de EQM do Lasso
  EQM_min_LASSO_Rep[p] = EQM_Teste[index_min] #Menor EQM obtido no Lasso
  EQM_No_LASSO_Rep[p] = EQM_Teste[1] #EQM quando Lambda = 0 (Equivalente a não ter penalização)
  EQM_mean_LASSO_Rep[p] = mean(EQM_Teste) #EQM medio de todos os Lassos (Para comparar)
  SD_EQM_LASSO_Rep[p] = sd(EQM_Teste) #Desvio padrão do EQM do Lasso (Para comparar)
  
  Acertos_teste_min_Rep[p] = acertos_teste[index_min]/length(S_teste)
  Acertos_teste_NoLasso_Rep[p] = acertos_teste[1]/length(S_teste)
  Acertos_teste_mean_Rep[p] = mean(acertos_teste/length(S_teste))
  SD_Acertos_teste_Rep[p] = sd(acertos_teste/length(S_teste))
 
  mu_Rep[p,1] = mu_hat[1]
  mu_Rep[p,2] = mu_hat[2]
  mu_Rep[p,3] = mu_hat[3]
  
  sigma_Rep[p,1] = sigma_hat[1]
  sigma_Rep[p,2] = sigma_hat[2]
  sigma_Rep[p,3] = sigma_hat[3]
  
  
  LL2_parte1 = 0
  LL2_parte2 = 0
  LL2_parte3 = 0
  LL2_parte4 = 0
  LL=0
  
  
  #Calculo do AICc e BIC
  LL2_parte1 = -.5*length(S_training)*log(2*pi)
  
  for (i in 1:T) {#Calculo do primeiro segmento da LL
    LL2_parte2 = LL2_parte2 -.5*log(sigma_hat[S_treino[i]]) 
  }
  for (i in 1:T) {#Calculo do segundo segmento da LL
    LL2_parte3 = LL2_parte3 -(1/(2*sigma_hat[S_treino[i]]))*((Y_training[i]-mu_hat[S_treino[i]])^2)
  }
  temp=NULL
  for (i in 2:length(S_training)) {#Calculo do terceiro segmento da LL
    for (g in 1:K) {
      temp[g]<-exp(X[i,]%*%matrix(BetaArray[g,,S_treino[i-1]],ncol=1))
    }
    LL2_parte4 = LL2_parte4 + (X[i,]%*%matrix(BetaArray[S_treino[i],,S_treino[i-1]]) - log(sum(temp), base = exp(1)))
  }
  LL <- log(P0[S_treino[1]]) + LL2_parte1 + LL2_parte2 + LL2_parte3 + LL2_parte4 #calculo da LogVerosim
  
  #Almacenamos os valores calculados dos criterios de info
  #Nas CritInfo_Tables, as almacenamos para poder ter acesso
  #A cada valor em cada replica, e não so à media.
  qtd_params = D*K*(K-1) + K*2
  CritInfo_Table[1,p]<-BIC_Acum[p]<--2*LL + qtd_params*log(T) #Calculo do BIC para cada Replica. O 8 aqui representa a quantidade de parametros do modelo
  CritInfo_Table[2,p]<-AICc_Acum[p]<--2*(LL) + 2*qtd_params + (2*qtd_params^2 + 2*qtd_params)/(T - qtd_params - 1) #Calculo do AIC corregido
  
  tempo_final<-proc.time()
  Tempo[p]<-(tempo_final[1] - tempo_inicial[1])/60
}# Fim das REplicas
T2 = proc.time()
T2 - T1
#Inicializamos algumas variaveis para
#almacenar calculos finais importantes



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
  Mu_Finais_Mean[i]=mean(mu_Rep[,i])
  Sigma_Finais_Mean[i]=mean(sigma_Rep[,i])
  Mu_Finais_Median[i]=median(mu_Rep[,i])
  Sigma_Finais_Median[i]=median(sigma_Rep[,i])
}

#Calculo do Desvio Padrão ao longo das R replicas
for (i in 1:K) {
  SD_Mu[i]=sd(mu_Rep[,i])
  SD_Sigma[i]=sd(sigma_Rep[,i])
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
quantiles[1,]<-quantile(mu_Rep[,1],probs = c(0.025, 0.975))
quantiles[2,]<-quantile(mu_Rep[,2],probs = c(0.025, 0.975))
quantiles[3,]<-quantile(mu_Rep[,3],probs = c(0.025, 0.975))
quantiles[4,]<-quantile(sigma_Rep[,1],probs = c(0.025, 0.975))
quantiles[5,]<-quantile(sigma_Rep[,2],probs = c(0.025, 0.975))
quantiles[6,]<-quantile(sigma_Rep[,3],probs = c(0.025, 0.975))

#Juntamos as Metricas dos parametros das distribuições observaveis num dataframe
Parametro<-c("Mu_1","Mu_2","Mu_3", "Sigma_1", "Sigma_2", "Sigma_3")
Real<-c(mu[1],mu[2],mu[3], sigma[1],sigma[2], sigma[3])
Estimado_Mean<-c(Mu_Finais_Mean[1],Mu_Finais_Mean[2],Mu_Finais_Mean[3], Sigma_Finais_Mean[1], Sigma_Finais_Mean[2], Sigma_Finais_Mean[3])
Estimado_Median<-c(Mu_Finais_Median[1],Mu_Finais_Median[2],Mu_Finais_Median[3], Sigma_Finais_Median[1], Sigma_Finais_Median[2], Sigma_Finais_Median[3])
SD<-c(SD_Mu[1],SD_Mu[2], SD_Mu[3], SD_Sigma[1],SD_Sigma[2],SD_Sigma[3])
Vies<-c(Vies_Mu[1],Vies_Mu[2],Vies_Mu[3], Vies_Sigma[1], Vies_Sigma[2],Vies_Sigma[3])
EQM<-c(EQM_Mu[1],EQM_Mu[2],EQM_Mu[3], EQM_Sigma[1], EQM_Sigma[2],EQM_Sigma[3])
Assimetria<-c(skewness(mu_Rep[,1]), skewness(mu_Rep[,2]),skewness(mu_Rep[,3]), skewness(sigma_Rep[,1]), skewness(sigma_Rep[,2]), skewness(sigma_Rep[,3]))
IC_95<-c(toString(quantiles[1,]),toString(quantiles[2,]),toString(quantiles[3,]),toString(quantiles[4,]),toString(quantiles[5,]),toString(quantiles[6,]))

#Criamos o DAtaframe
df1<-data.frame(Parametro,Real,Estimado_Mean, Estimado_Median,SD,Vies,EQM,Assimetria,IC_95)
df1


#############################################################
## METRICAS PARA COMPARAR/AVALIAR A QUALIDADE DA ESTIMAÇÃO ##
#############################################################
EQM_Min_Final_Median = NULL
EQM_Min_Final = NULL
EQM_Min_SD = NULL

EQM_No_LASSO_Final = NULL
EQM_No_LASSO_Final_Median = NULL
EQM_No_LASSO_SD = NULL

EQM_Mean_LASSO_Final = NULL
EQM_Mean_LASSO_Final_Median = NULL
EQM_Mean_LASSO_SD = NULL


#Media e DP do menor EQM obtido ao longo do Lasso em cada Replica
EQM_Min_Final_Median = median(EQM_min_LASSO_Rep)
EQM_Min_Final = mean(EQM_min_LASSO_Rep)
EQM_Min_SD = sd(EQM_min_LASSO_Rep)

#Media e DP do EQM medio obtido ao longo do Lasso em cada Replica
EQM_Mean_LASSO_Final = mean(EQM_mean_LASSO_Rep)
EQM_Mean_LASSO_Final_Median = median(EQM_mean_LASSO_Rep)
EQM_Mean_LASSO_SD = sd(EQM_mean_LASSO_Rep)

#Media e DP do EQM para modelo sem Lasso 
EQM_No_LASSO_Final = mean(EQM_No_LASSO_Rep)
EQM_No_LASSO_Final_Median = median(EQM_No_LASSO_Rep)
EQM_No_LASSO_SD = sd(EQM_No_LASSO_Rep)

#Criamos p DAtaframe
Chain_Length <- paste("T =",toString(T))
df2<-data.frame(Chain_Length, EQM_Min_Final_Median, EQM_Min_Final, EQM_Min_SD, EQM_No_LASSO_Final_Median, EQM_No_LASSO_Final, EQM_No_LASSO_SD, EQM_Mean_LASSO_Final_Median, EQM_Mean_LASSO_Final, EQM_Mean_LASSO_SD)
df2


#####################################################################
## METRICAS PARA COMPARAR/AVALIAR A QUANTO ACERTAMOS NA SEQUÊNCIA S ##
#####################################################################
Acertos_teste_min_Final = NULL
Acertos_teste_min_Final_Median = NULL
Acertos_teste_min_SD = NULL

Acertos_teste_NoLasso_Final = NULL
Acertos_teste_NoLasso_Final_Median = NULL
Acertos_teste_NoLasso_SD = NULL

Acertos_teste_mean_Final = NULL
Acertos_teste_mean_Final_Median = NULL
Acertos_teste_mean_SD = NULL


#Acertos em S_teste com o modelo com o Menor EQM
Acertos_teste_min_Final_Median = median(Acertos_teste_min_Rep)
Acertos_teste_min_Final = mean(Acertos_teste_min_Rep)
Acertos_teste_min_SD = sd(Acertos_teste_min_Rep)


#Acertos em S_teste com o modelo sem LASSO 
Acertos_teste_NoLasso_Final = mean(Acertos_teste_NoLasso_Rep)
Acertos_teste_NoLasso_Final_Median = median(Acertos_teste_NoLasso_Rep)
Acertos_teste_NoLasso_SD = sd(Acertos_teste_NoLasso_Rep)

#Acertos medios em S_teste ao longo do LASSO em cada Replica
Acertos_teste_mean_Final = mean(Acertos_teste_mean_Rep)
Acertos_teste_mean_Final_Median = median(Acertos_teste_mean_Rep)
Acertos_teste_mean_SD = sd(Acertos_teste_mean_Rep)


#Criamos o Dataframe
Chain_Length <- paste("T =",toString(T))
df3<-data.frame(Chain_Length, Acertos_teste_min_Final_Median, Acertos_teste_min_Final, Acertos_teste_min_SD,  Acertos_teste_NoLasso_Final_Median, Acertos_teste_NoLasso_Final, Acertos_teste_NoLasso_SD, Acertos_teste_mean_Final_Median, Acertos_teste_mean_Final, Acertos_teste_mean_SD)
df3

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
Specificity_Final = mean(specificity)
Specificity_Final_Median = median(specificity)
Specificity_SD = sd(specificity)

#Sensitividade
Sensitivity_Final = mean(sensitivity)
Sensitivity_Final_Median = median(sensitivity)
Sensitivity_SD = sd(sensitivity)

#Acuracia
Accuracy_Final = mean(accuracy)
Accuracy_Final_Median = median(accuracy)
Accuracy_SD = sd(accuracy)

Chain_Length <- paste("T =",toString(T))
df4<-data.frame(Chain_Length, Specificity_Final_Median, Specificity_Final,Specificity_SD, Sensitivity_Final_Median, Sensitivity_Final,Sensitivity_SD, Accuracy_Final_Median, Accuracy_Final,Accuracy_SD)
df4




#######################
## Outras Metricas ####
#######################


#Calculo da media dos AICc para o modelo
AICc_Final=mean(AICc_Acum)

#Calculo da media dos AICc para o modelo
BIC_Final=mean(BIC_Acum)


Tempo_medio = mean(Tempo)
Sd_Tempo = sd(Tempo)
Tempo_Total = sum(Tempo)

Indicadores<-c("Tempo Medio por Replica", "Desvio Padrão do Tempo Medio", "Tempo Total de Procesamento")
Valor<-c(Tempo_medio, Sd_Tempo, Tempo_Total)
df5<-data.frame(Indicadores,Valor)

df1
df2
df3
df4
df5

library("r2excel")
wb<-createWorkbook(type = "xlsx")   #Criamos um livro de trabalho de excel
sheet<-createSheet(wb, sheetName = nomefolha) #criamos uma folha no livro de trabalho que tinhamos criado
xlsx.addHeader(wb, sheet, value=header,level=1, color="black", underline=1) # escrevemos um titulo na folha 
xlsx.addLineBreak(sheet, 1) #Um linha em ranco embaixo do titulo
xlsx.addParagraph(wb, sheet, "Metricas dos Parametros das Distribuições Observaveis", fontColor = "#FFFFFF", fontSize = 12, 
                  backGroundColor = "#FFFFFF", isBold = TRUE, isItalic = TRUE,
                  startRow = 3, startCol = 3, colSpan = 10, rowSpan = 2)
xlsx.addTable(wb, sheet, df1, startCol=2)#Criamos a primeira tabela usando o primeiro dataframe
xlsx.addLineBreak(sheet, 4)#insertamos 2 linhas em branco

xlsx.addParagraph(wb, sheet, "Metricas do Erro Quadratico Medio (Qualidade da Previsão)", fontColor = "#FFFFFF", fontSize = 12, 
                  backGroundColor = "#FFFFFF", isBold = TRUE, isItalic = TRUE,
                  startRow = 14, startCol = 3, colSpan = 10, rowSpan = 2)
xlsx.addTable(wb, sheet, df2, startCol=2)#criamos a segunda Tabela usando o segundo dataframe
xlsx.addLineBreak(sheet, 4)#insertamos 2 linhas em branco

xlsx.addParagraph(wb, sheet, "Metricas de Acertos em S Teste", fontColor = "#FFFFFF", fontSize = 12, 
                  backGroundColor = "#FFFFFF", isBold = TRUE, isItalic = TRUE,
                  startRow = 22, startCol = 3, colSpan = 10, rowSpan = 2)
xlsx.addTable(wb, sheet, df3, startCol=2)#criamos a segunda Tabela usando o segundo dataframe
xlsx.addLineBreak(sheet, 4)#insertamos 2 linhas em branco
# 

xlsx.addParagraph(wb, sheet, "Metricas de Acuracia do LASSO (Zerou ou Não Zerou?)", fontColor = "#FFFFFF", fontSize = 12, 
                  backGroundColor = "#FFFFFF", isBold = TRUE, isItalic = TRUE,
                  startRow = 30, startCol = 3, colSpan = 10, rowSpan = 2)
xlsx.addLineBreak(sheet, 1)#insertamos 2 linhas em branco
xlsx.addTable(wb, sheet, df4, startCol=2)#criamos a segunda Tabela usando o segundo dataframe
xlsx.addLineBreak(sheet, 7)#insertamos 2 linhas em branco


xlsx.addTable(wb, sheet, df5, startCol=2)#criamos a segunda Tabela usando o segundo dataframe
xlsx.addLineBreak(sheet, 3)#insertamos 2 linhas em branco
saveWorkbook(wb, nome_arquivo)#guardamos o arquivo de excel.

sigma_Rep
