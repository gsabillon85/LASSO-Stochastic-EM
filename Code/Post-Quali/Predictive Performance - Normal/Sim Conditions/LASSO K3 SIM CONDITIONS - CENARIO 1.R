###############################################################
#                                                             #
#            GLOBAL LASSO - 3 ESTADOS OCULTOS                 #
#                                                             #
###############################################################

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
R <- 6 # Numero de Replicas
T=1200 #Cumprimento da cadeia simulada
K=3   #Numero de estados ocultos
D=8   #Quantidade de Covariaveis
tol<-0.0000001 #Nivel de tolerancia que estabelecemos como criterio de parada do EM Est
tolval=NULL
tolval[1]=1
optim_algo = "BFGS" #Algorithm to use in the optimization process
n_max_iter_EM = 50
Tempo <- NULL

cenario <- "TESTANDO_CODIGO"
mainDir = paste("/home/gustavo/Projects/LASSO-NHMM/Code/Post-Quali/Predictive Performance - Normal/Global LASSO Results/K3/",cenario,sep = "")
subDir = paste("Resultados_T",toString(T),"_D",toString(D),"_zero-threshold",toString(zero_threshold),"optimethod",toString(optim_algo))
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
setwd(file.path(mainDir, subDir))

set.seed(12)
seeds <-sample(110000000,R) # Seed number para conseguer os mesmos valores simulados
#lambdas <- c(0, 0.0001, 0.001, 0.01, 0.1, 0.25, 0.5, 1, 2.5, 5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 30)
lambdas <- seq(0, 20, by=0.5)

# seeds[1] = 307
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
Best_Lambdas <- NULL

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
sigma = c(2,1.5,3) #Vetor com os desvios padrões para as duas normais


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





###############################################################
#                                                             #
#            INDIVUAL LASSO - 3 ESTADOS OCULTOS              #
#                                                             #
###############################################################

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
R <- 6 # Numero de Replicas
T=1200 #Cumprimento da cadeia simulada
K=3   #Numero de estados ocultos
D=8   #Quantidade de Covariaveis
tol<-0.0000001 #Nivel de tolerancia que estabelecemos como criterio de parada do EM Est
tolval=NULL
tolval[1]=1
optim_algo = "BFGS" #Algorithm to use in the optimization process
n_max_iter_EM = 50
Tempo <- NULL

cenario <- "TESTANDO_CODIGO"
mainDir = paste("/home/gustavo/Projects/LASSO-NHMM/Code/Post-Quali/Predictive Performance - Normal/Global LASSO Results/K3/",cenario,sep = "")
subDir = paste("Resultados_T",toString(T),"_D",toString(D),"_zero-threshold",toString(zero_threshold),"optimethod",toString(optim_algo))
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
setwd(file.path(mainDir, subDir))

set.seed(2)
seeds <-sample(110000000,R) # Seed number para conseguer os mesmos valores simulados
lambdas <- c(0, 0.25, 0.5, 1, 2.5, 5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 30)
lambdas <- seq(2.5, 12.5, by=0.75)

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
Best_Lambdas <- matrix(nrow = R, ncol = K)

## SEÇÃO DE DEFINICAÇÃO DOS PARAMETROS PARA SIMULAÇÃO DE DADOS ##
################################################################
P0=rep(1/K,K) #Inicializamos vetor de probabilidades inciais para o HMM
Betas=array(0, dim=c(K,D,K)) # valores de beta utilizados na geração dos valores (consideranda intercepto e duas covariáveis)
Real=NULL



Betas[2,1,1]=-1.4
Betas[2,3,1]=2.5
Betas[2,8,1]=1.7

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
sigma = c(2,1.5,3) #Vetor com os desvios padrões para as duas normais


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