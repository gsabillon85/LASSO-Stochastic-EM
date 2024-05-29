library('label.switching')
#library("r2excel")
#library("xlsx")
library(e1071)
zero_threshold = 0.05
T=300 #Cumprimento da cadeia simulada
K=2   #Numero de estados ocultos
D=5   #Quantidade de Covariaveis
tol<-0.0000001 #Nivel de tolerancia que estabelecemos como criterio de parada do EM Est
tolval=NULL
tolval[1]=1

optim_algo = "BFGS"

# Generate 30 random seeds to generate 20 random samples
n_datasets = 30
set.seed(10)
seeds <-sample(110000000,n_datasets) # Seed number para conseguer os mesmos valores simulados
lambdas <- c(0, 0.0001, 0.001, 0.01, 0.1, 0.25, 0.5, 1, 2.5, 5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 30, 50, 75, 100, 200)
#lambdas = lambdas_large

RealTransCount = array(rep(0,K*K*length(lambdas)), dim=c(K,K,length(lambdas)))
avg_Sensitivity = NULL
sd_Sensitivity = NULL
avg_Specificity = NULL
sd_Specificity = NULL
avg_Accuracy = NULL
sd_Accuracy= NULL
avg_pctgm_zerado= NULL
sd_pctgm_zerado = NULL
avg_RMSE = NULL
sd_RMSE = NULL
avg_Parameters = matrix(nrow = length(lambdas), ncol = D*K*(K-1) )

mainDir = "/home/gustavo/Projects/LASSO-NHMM/Code/Testing - Regularization Accuracy/Shrinkage Accuracy - Results"
subDir = paste("Resultados_T",toString(T),"_D",toString(D),"_zero-threshold",toString(zero_threshold),"optimethod",toString(optim_algo))
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
setwd(file.path(mainDir, subDir))


tempo_inicial<-proc.time()
for (p in 1:length(lambdas)){
  Parameters <- matrix(nrow = length(seeds), ncol = D*K*(K-1) )
  pctgm_zerado <- NULL
  Precision <- NULL
  Sensitivity <- NULL
  Specificity <- NULL
  Accuracy <- NULL
  RMSE<-NULL
  # Create Lambda vector
  
  options(digits=8)
  options(scipen=999)
  cat(paste("\nValor de Lambda: ",toString(lambdas[p]), collapse = "")) #Messagem indicando o numero da replica atual
  for (h in 1:length(seeds)){
    cat(paste("\nReplica de Dados Numero: ",toString(h), collapse = "")) #Messagem indicando o numero da replica atual
    set.seed(seeds[h])
    tempo_inicial<-proc.time()#Calcularemos quanto demoro todo o proceso
    
    lambda = lambdas[p]
    
    VerProx<-NULL
    VerAct<-NULL
    
    ####### Simulação #######
    TransCount <- matrix(data = c(rep(0,K^2)), nrow = K, ncol = K)
    P0=rep(1/K,K) #Inicializamos vetor de probabilidades inciais para o HMM
    mu = c(50,150) # vetor com media para as duas Normais
    sigma = c(0.75,1.25) #Vetor com os desvios padrões para as duas normais
    mu_hat <- NULL
    sigma_hat <- NULL
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
    
    Betas[1,1,2]=0
    Betas[1,2,2]=0
    Betas[1,3,2]=0
    Betas[2,1,2]=2.4
    Betas[2,2,2]=0
    Betas[2,3,2]=-1.5
    
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
      resp <- (sum(1 - log(1 + exp(Xtemp11%*%params))) + sum(Xtemp12%*%params - log(1 + exp(Xtemp12%*%params)))) - lambda*sum(abs(params[2:D])) 
    }
    
    FSM2 <-function(params){#função a maximizar para achar os Betas_2
      resp <- (sum(1 - log(1 + exp(Xtemp21%*%params))) + sum(Xtemp22%*%params - log(1 + exp(Xtemp22%*%params))))  - lambda*sum(abs(params[2:D]))
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
    
    #########################################
    #####   Procedimento de Estimação   #####
    #########################################
    
    #######Primeiro geramos uma sequência não observavel de treinamento######
    P_Treino=rep(1/K,K) #Vetor de probabilidade utilizadas para gerar a sequência de treino
    S_treino<-NULL # Inicializamos a sequência oculta de treinamento
    init1 = c(rnorm(D*(K-1)))#Valores iniciais para os Betas_1
    init2 = c(rnorm(D*(K-1)))#Valores iniciais para os Betas_2

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
    #Agora executamos o Algoritmo EM Estocástico
    while ( tolval[val]>tol ){
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
        mu_hat[k] = sum(id*Y)/sum(id)
        Y_id_list = split(Y,id)
        Y_id = unlist(Y_id_list[2], use.names = FALSE)
        sigma_hat[k] = sqrt((sum((Y_id - mu_hat[k])^2)) / (sum(id) - 1)) #DECIDIR SOBRE O ESTIMADOR DA VARIANCIA (VICIADO OU NṼICIADO)
      }
      
      #print(mu_hat)
      #print(sigma_hat)
      #Calculo da Verosimilhança como valor de tolerança
      LL_parte1 = -.5*length(S)*log(2*pi)
      
      for (i in 1:length(S)) {#Calculo do primeiro segmento da LL
        LL_parte2 = LL_parte2 -.5*log(sigma_hat[S_treino[i]]) 
      }
      for (i in 1:length(S)) {#Calculo do segundo segmento da LL
        LL_parte3 = LL_parte3 -(1/(2*sigma_hat[S_treino[i]]))*((Y[i]-mu_hat[S_treino[i]])^2)
      }
      temp=NULL
      for (i in 2:length(S)) {#Calculo do terceiro segmento da LL
        for (g in 1:K) {
          temp[g]<-exp(X[i,]%*%matrix(BetaArray[g,,S_treino[i-1]],ncol=1))
        }
        LL_parte4 = LL_parte4 + (X[i,]%*%matrix(BetaArray[S_treino[i],,S_treino[i-1]]) - log(sum(temp), base = exp(1)))
      }
      VeroSimActual <- log(P0[S_treino[1]]) + LL_parte1 + (LL_parte2 + LL_parte3) + LL_parte4 #calculo da LogVerosim
      
      
      #Calculamos a sequência S_treino utilizando os Betas
      #Atualizados na iteração passada e os valores observados Y
      S_treino[1]=which.max(dnorm(Y[1], mu_hat, sigma_hat))
      for (i in 2:length(S)) {
        A_hat_t = Mat_trans(X[i,])
        if (any(is.na(A_hat_t))){
          print("NaN encountered in Transition Matrix Calculation")
          A_hat_t[is.nan(A_hat_t)] = 1 
        }
        prob<-(A_hat_t[S_treino[i], ]*dnorm(Y[i], mu_hat, sigma_hat))/sum(A_hat_t[S_treino[i], ]*dnorm(Y[i], mu_hat, sigma_hat))
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
        
      #####################################
      #Este segmento de codigo testa se aconteceram todas as transições possiveis
      #No caso que elas não tinham acontecido, as que
      #não aconteceram são forçadas a acontecer
      TransCount <- matrix(data = c(rep(0, K*K)), nrow = K, ncol = K)
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
      Xtemp21<-NULL
      Xtemp22<-NULL
  
     
      
      for (t in 2:T) {
        #filtros indo para o Estado # 1
        if(S_treino[t]%in%1 && S_treino[t-1]%in%1)
          Xtemp11<-rbind(Xtemp11, X[t,])
        
        if(S_treino[t]%in%1 && S_treino[t-1]%in%2)
          Xtemp21<-rbind(Xtemp21, X[t,])
        
        
        #Filtros indo para o Estado # 2
        if(S_treino[t]%in%2 && S_treino[t-1]%in%1)
          Xtemp12<-rbind(Xtemp12, X[t,])
        
        if(S_treino[t]%in%2 && S_treino[t-1]%in%2)
          Xtemp22<-rbind(Xtemp22, X[t,])
      }
      
      
      ##O ajuste para estimar os parâmetros de transição é
      ##feito aqui usando a função optim e os valores das
      #covariaveis filtradas
      fit1 <- optim(par = init1, fn = FSM1, control = list(fnscale=-1), method = optim_algo, hessian = FALSE)
      fit2 <- optim(par = init2, fn = FSM2, control = list(fnscale=-1), method = optim_algo, hessian = FALSE)
          
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
      
      
      
      LL2_parte1 = 0
      LL2_parte2 = 0
      LL2_parte3 = 0
      LL2_parte4 = 0
      VeroSimProxima=0
      
      #Calculo da Verosimilhança como valor de tolerança
      LL2_parte1 = -.5*length(S)*log(2*pi)
      
      for (i in 1:length(S)) {#Calculo do primeiro segmento da LL
        LL2_parte2 = LL2_parte2 -.5*log(sigma_hat[S_treino[i]]) 
      }
      for (i in 1:length(S)) {#Calculo do segundo segmento da LL
        LL2_parte3 = LL2_parte3 -(1/(2*sigma_hat[S_treino[i]]))*((Y[i]-mu_hat[S_treino[i]])^2)
      }
      temp=NULL
      for (i in 2:length(S)) {#Calculo do terceiro segmento da LL
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
      #print(tolval[val])
      
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
    while (tolval[val]>tol2 ) {
      #print(tolval[val])
      #Aqui devemos calcular a diferença entre a L.V. em na iteração atual e na anterior  
      LL_parte1 = 0
      LL_parte2 = 0
      LL_parte3 = 0
      LL_parte4 = 0
      VeroSimActual=0
      
      LL_parte1 = -.5*T*log(2*pi)
      
      for (i in 1:length(S)) {#Calculo do primeiro segmento da LL
        LL_parte2 = LL_parte1 +.5*log(sigma_hat[S_treino[i]]) 
      }
      for (i in 1:length(S)) {#Calculo do segundo segmento da LL
        LL_parte3 = LL_parte3 +(1/(2*sigma_hat[S_treino[i]]))*((Y[i]-mu_hat[S_treino[i]])^2)
      }
      temp=NULL
      for (i in 2:length(S)) {#Calculo do terceiro segmento da LL
        for (g in 1:K) {
          temp[g]<-exp(X[i,]%*%matrix(BetaArray[g,,S_treino[i-1]],ncol=1))
        }
        LL_parte4 = LL_parte4 + (X[i,]%*%matrix(BetaArray[S_treino[i],,S_treino[i-1]]) - log(sum(temp), base = exp(1)))
      }
      VeroSimActual <- log(P0[S_treino[1]]) + LL_parte1 - (LL_parte2 + LL_parte3) + LL_parte4 #calculo da LogVerosim
      
      #filtragem dos dados
      Xtemp11<-NULL
      Xtemp12<-NULL
      Xtemp21<-NULL
      Xtemp22<-NULL
      
      for (t in 2:T) {
        #filtros indo para o Estado # 1
        if(S_treino[t]%in%1 && S_treino[t-1]%in%1)
          Xtemp11<-rbind(Xtemp11, X[t,])
        
        if(S_treino[t]%in%1 && S_treino[t-1]%in%2)
          Xtemp21<-rbind(Xtemp21, X[t,])
        
        #Filtros indo para o Estado # 2
        if(S_treino[t]%in%2 && S_treino[t-1]%in%1)
          Xtemp12<-rbind(Xtemp12, X[t,])
        
        if(S_treino[t]%in%2 && S_treino[t-1]%in%2)
          Xtemp22<-rbind(Xtemp22, X[t,])
      }
      
      fit1 <- optim(par = init1, fn = FSM1, control = list(fnscale=-1), method = optim_algo, hessian = FALSE)
      fit2 <- optim(par = init2, fn = FSM2, control = list(fnscale=-1), method = optim_algo, hessian = FALSE)
 
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
      
      LL2_parte1 = 0
      LL2_parte2 = 0
      LL2_parte3 = 0
      LL2_parte4 = 0
      VeroSimProxima=0
      
      #Calculo da Verosimilhança como valor de tolerança
      LL2_parte1 = -.5*length(S)*log(2*pi)
      
      for (i in 1:length(S)) {#Calculo do primeiro segmento da LL
        LL2_parte2 = LL2_parte2 +.5*log(sigma_hat[S_treino[i]]) 
      }
      for (i in 1:length(S)) {#Calculo do segundo segmento da LL
        LL2_parte3 = LL2_parte3 +(1/(2*sigma_hat[S_treino[i]]))*((Y[i]-mu_hat[S_treino[i]])^2)
      }
      temp=NULL
      for (i in 2:length(S)) {#Calculo do terceiro segmento da LL
        for (g in 1:K) {
          temp[g]<-exp(X[i,]%*%matrix(BetaArray[g,,S_treino[i-1]],ncol=1))
        }
        LL2_parte4 = LL2_parte4 + (X[i,]%*%matrix(BetaArray[S_treino[i],,S_treino[i-1]]) - log(sum(temp), base = exp(1)))
      }
      VeroSimProxima <- log(P0[S_treino[1]]) + LL2_parte1 - (LL2_parte2 + LL2_parte3) + LL2_parte4 #calculo da LogVerosim
      
      val=val+1
      tolval[val]<-VeroSimProxima-VeroSimActual
      #print(tolval[val])
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
          Beta_Post_Array[i,d,2]=fit2$par[d]
        } else if (i == 3){
          Beta_Post_Array[i,d,2]=fit2$par[D+d]
        }
        
      }
    }

    #Calculo do AICc e BIC
    
    #####################################################################
    #####################################################################
    
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
    
    for (j in 1:(D*K*(K-1))){
      if ((Real[j]==0) && (abs(Estimado[j])<=zero_threshold)){
        TP[j] = 1
      }
      if ((Real[j]!=0) && (abs(Estimado[j])>zero_threshold)){
        TN[j]=1
      }
      if ((Real[j]!=0) && (abs(Estimado[j])<=zero_threshold)){
        FP[j]=1
      }
      if ((Real[j]==0) && (abs(Estimado[j])>zero_threshold)){
        FN[j]=1
      }
    }
    
    for (w in 1:length(Estimado)){
      Parameters[h,w] = Estimado[w]
    }
    
    pctgm_zerado[h] = sum(abs(Estimado<=zero_threshold))/(D*K*(K-1))
    Precision[h] = sum(TP, na.rm = TRUE) / ( sum(TP, na.rm = TRUE) +  sum(FP, na.rm = TRUE) )  
    Sensitivity[h] = sum(TP, na.rm = TRUE) / ( sum(TP, na.rm = TRUE) +  sum(FN, na.rm = TRUE) )
    Specificity[h] = sum(TN, na.rm = TRUE) / ( sum(TN, na.rm = TRUE) + sum(FP, na.rm = TRUE) )
    Accuracy[h] = (sum(TN, na.rm = TRUE) + sum(TP, na.rm = TRUE)) / ( sum(TP, na.rm = TRUE) + sum(TN, na.rm = TRUE) + sum(FP, na.rm = TRUE) + sum(FN, na.rm = TRUE))
    RMSE[h] = sum( (Real - Estimado)^2  )
  }
  avg_RMSE[p] = mean(RMSE)
  sd_RMSE[p] = sd(RMSE)
  avg_Sensitivity[p] = mean(Sensitivity)
  sd_Sensitivity[p] = sd(Sensitivity)
  avg_Specificity[p] = mean(Specificity)
  sd_Specificity[p] = sd(Specificity)
  avg_Accuracy[p] = mean(Accuracy)
  sd_Accuracy[p] = sd(Accuracy)
  avg_pctgm_zerado[p] = mean(pctgm_zerado)
  sd_pctgm_zerado[p] = sd(pctgm_zerado)
  for (o in 1:(D*K*(K-1))){
    avg_Parameters[p,o] <- mean(Parameters[,o])
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
write.csv(data.frame(Final_TransCount), paste("Real_TransitionsCount_T-",toString(T),"_D-",toString(D),"_zerothreshold-",toString(zero_threshold),"_optimethod-",toString(optim_algo),".csv", sep = ""), row.names=FALSE)
write.csv(data.frame(Final_TransCount_sd), paste("Real_Transitions_SD_T-",toString(T),"_D-",toString(D),"_zerothreshold-",toString(zero_threshold),"_optimethod-",toString(optim_algo),".csv", sep = ""), row.names=FALSE)


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

