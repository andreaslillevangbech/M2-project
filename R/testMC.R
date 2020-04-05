
execute = function(DGP, Pc=50,  N=200, T=180){
  
  library(reshape2)
  library(tidyverse)
  library(abind)
  library(glmnet)
  library(pls)
  library(lars)
  library(mvtnorm)
  library(gglasso)
  library(splines)
  library(leaps)
  library(tree)
  library(randomForest)
  library(gbm)
  
  #Creating c_bar for firm characteristics.
  rho_c=runif(n = Pc, min = 0.9, max = 1)
  eps_c = rnorm(T*N*Pc)
  dim(eps_c)=c(N,Pc,T)
  c_bar=array(dim = c(N, Pc, T))
  c=array(dim = c(N, Pc, T))
  c_bar[,,1]=eps_c[,,1]
  for (i in 1:N){
    for (j in 1:Pc){
      for (t in 2:T){
        c_bar[i,j,t]=rho_c[j]*c_bar[i,j,t-1]+eps_c[i,j,t]
        
      }
    }
  }
  
  
  
  rank=rank(c(c_bar))
  dim(rank)=c(N,Pc,T)
  n=N*T*Pc
  
  for (i in 1:N){
    for (j in 1:Pc){
      for (t in 1:T){
        c[i,j,t] = 2*rank[i,j,t]/(n+1) -1
        
      }
    }
  }
  
  beta=c[,1:3,]
  
  #Creating x, which is just a common timeseries. AR(1) no constant.
  rho=0.95
  u=rnorm(T+200, sd = 1-rho^2)
  x=numeric(T+200)
  x[1]=0
  for (t in 2:(T+200)){
    x[t]=rho*x[t-1]+u[t]
  }
  x=x[(length(x)-T+1):length(x)]
  
  #Create vector z which is interaction of all x with c
  c2=array(dim=dim(c))
  for (t in 1:T){
    c2[,,t]=x[t]*c[,,t]
  }
  z=abind::abind(c,c2,along= 2)
  
  #Finally create the returns from the characteristics
  #v=mvrnorm(n = T,mu= c(0,0,0),Sigma = diag(0.05^2,3,3))
  v=mvtnorm::rmvnorm(n=T,mean = c(0,0,0),sigma = diag(0.05^2,3,3))
  eps=rt(n = T*N,df = 5)*0.05
  dim(eps)=c(N,T)
  
  #There are two cases for the function g(z), case (a) and (b).
  g=matrix(nrow = N,ncol = T)
  
  #for DGP (a)
  if(DGP=="a"){
    theta=c(rep(0.02,3))
    for (i in 1:N){
      for (t in 1:T){
        g[i,t]=z[i,1,t]*theta[1]+z[i,2,t]*theta[2]+z[i,(Pc+3),t]*theta[3]
      }
    }
  }
  
  #for DGP (b)
  if(DGP=="b"){
    theta=c(0.04,0.03,0.012)
    for (i in 1:N){
      for (t in 1:T){
        g[i,t]=(z[i,1,t]^2)*theta[1]+z[i,1,t]*z[i,2,t]*theta[2]+sign(z[i,(Pc+3),t])*theta[3]
      }
    }
  }
  
  
  #Finally we calculate the returns, one for dgp (a) and one for (b)
  e=matrix(nrow = N,ncol = T)
  for (i in 1:N){
    for (t in 1:T){
      e[i,t] = beta[i,,t]%*%v[t,]+eps[i,t]
    }
  }
  
  #The response is thus
  r=g+e
  #Note that inputs range from period 1 to T and response ranges from 2 to 
  #T+1
  
  #make inputs to a panel of data instead of 3d array
  #likewise make response into vector instead of 2d array
  #the rows are ordered st. all firms in one period are together
  Y=c(r)
  
  input <- lapply(seq_len(T), function(i) data.frame(z[,,i])) %>% 
    bind_rows(.) 
  
  df <- cbind(Y,input)
  df$time = sapply(1:T, function(i) rep(i,N)) %>% c(.)
  df$id=rep(1:N,T)
  
  ########################################################################
  ################## ------- SIMULATION OVER ----------- ###################
  ########################################################################
  
  initial_T=T/2
  val_T = T/4
  sample.size = N*T
  test.sample = which(df$time>T-val_T)
  test.sample.size = length(df$time[df$time > (T-val_T)])
  
  train.sample = which(df$time<=T-val_T)
  historical.mean = mean(df$Y[train.sample])
  
  tot.test = mean((df$Y[test.sample] - historical.mean)^2)
  
  tot.test.no.demean = mean(df$Y[test.sample]^2)
  
  lambda.grid = 10^seq(2,-5, length = 100)
  
  # For Neural Network
  build_model = function(layers=2, lr = 0.0001, lambda=lambda.grid.NN[i]){
    
    ## One layer
    
    if (layers == 1){
      model <- keras_model_sequential() %>%
        layer_dense(units = 32, 
                    kernel_regularizer = regularizer_l1(lambda),
                    input_shape = dim(x_train)[2]) %>%
        layer_activation("relu") %>% 
        layer_batch_normalization() %>%
        layer_dense(units = 1, kernel_regularizer = regularizer_l1(lambda))
    }
    
    if(layers == 2){
      model <- keras_model_sequential() %>%
        layer_dense(units = 32, 
                    kernel_regularizer = regularizer_l1(lambda),
                    input_shape = dim(x_train)[2]) %>%
        layer_activation("relu") %>% 
        layer_batch_normalization() %>%
        
        layer_dense(units = 16,
                    kernel_regularizer = regularizer_l1(lambda)) %>%
        layer_activation("relu") %>% 
        layer_batch_normalization() %>% 
        layer_dense(units = 1, kernel_regularizer = regularizer_l1(lambda))
    }
    
    if(layers == 3){
      model <- keras_model_sequential() %>%
        layer_dense(units = 32, 
                    kernel_regularizer = regularizer_l1(lambda),
                    input_shape = dim(x_train)[2]) %>%
        layer_activation("relu") %>% 
        layer_batch_normalization() %>%
        
        layer_dense(units = 16,
                    kernel_regularizer = regularizer_l1(lambda)) %>%
        layer_activation("relu") %>% 
        layer_batch_normalization() %>% 
        layer_dense(units = 8, 
                    kernel_regularizer = regularizer_l1(lambda)) %>% 
        layer_activation("relu") %>% 
        layer_batch_normalization() %>% 
        layer_dense(units = 1, kernel_regularizer = regularizer_l1(lambda))
      
    }
    if(layers >= 4){
      model <- keras_model_sequential() %>%
        layer_dense(units = 32, 
                    kernel_regularizer = regularizer_l1(lambda),
                    input_shape = dim(x_train)[2]) %>%
        layer_activation("relu") %>% 
        layer_batch_normalization() %>%
        
        layer_dense(units = 16,
                    kernel_regularizer = regularizer_l1(lambda)) %>%
        layer_activation("relu") %>% 
        layer_batch_normalization() %>% 
        layer_dense(units = 8, 
                    kernel_regularizer = regularizer_l1(lambda)) %>% 
        layer_activation("relu") %>% 
        layer_batch_normalization() %>% 
        layer_dense(units = 4, 
                    kernel_regularizer = regularizer_l1(lambda)) %>% 
        layer_activation("relu") %>% 
        layer_batch_normalization() %>% 
        layer_dense(units = 1, kernel_regularizer = regularizer_l1(lambda))
      
    }
    
    model %>% compile(
      loss = "mean_squared_error",
      optimizer = optimizer_adam(lr = lr),
      metrics = list("mean_absolute_error")
    )
    
    model
    
    #End function  
  }
  
  
  ##Vectors for predictions
  lm.pred = numeric(N*T)
  oracle.pred = numeric(N*T)
  ridge.pred = numeric(N*T)
  lasso.pred = numeric(N*T)
  enet.pred = numeric(N*T)
  pcr.pred = numeric(N*T)
  plsr.pred = numeric(N*T)
  glm.pred = numeric(N*T)
  randomF.pred = numeric(N*T)
  boost.pred = numeric(N*T)
  NN1.pred = numeric(N*T)
  NN2.pred = numeric(N*T)
  NN3.pred = numeric(N*T)
  
  ########################################################################
  ################## ------- STARTING LOOP ----------- ###################
  ########################################################################
  
  
  ### Loop it over
  seq = seq(initial_T, T - val_T -1, by = 12)
  
  #Creating a loop over the periods to reestimate the model. 
  #Using expanding window like in the paper
  for (t in seq){
    
    #The current window in terms of time index as well as the corresponding observations
    test.window = seq(t+val_T+1, ifelse(t+val_T+12<= T, t+val_T+12, T) )
    test.window.obs = which(df$time %in% test.window)
    
    #Create the three samples for training, tuning and testing
    train= df %>% filter(time %in% 1:t) %>% select(-c(id, time))
    val= df %>% filter(time %in% (t+1):(t+val_T)) %>% select(-c(id, time))
    test= df %>% 
      filter( time %in% test.window ) %>% 
      select(-c(id,time))
    
    full.train = rbind(train,val)
    
    ################################################################################
    ## ------ ESTIMATING MODELS FOR THE GIVEN TRAINING AND VAL SAMPLE --------
    ################################################################################
    
    ### ------------------------------------------------------------------
    # Standard Linear Model
    ###--------------------------------------------------------------------
    fit = lm(Y ~ . , data=rbind(train,val))
    lm.pred[test.window.obs] = predict(fit, newdata = test[,-1])
    
    ## IS
    lm.pred.IS = predict(fit, newdata = full.train[,-1])
    
    ## What variables are in the oracle
    
    fml = as.formula(str_c("Y ~ " , "X1 + ", "X2 + ", "X", Pc+3 ))
    fit.oracle = lm(fml, data=rbind(train,val) )
    oracle.pred[test.window.obs] = predict(fit.oracle, newdata = test[,-1])
    oracle.pred.IS = predict(fit.oracle, newdata = full.train[,-1])
    
    
    
    ### ------------------------------------------------------------------
    # Shrinking / selection parameter framework
    ###--------------------------------------------------------------------
    
    ### ------------------------------------------------------------------
    #Ridge regression
    ##--------------------------------------------------------------------
    ridge = glmnet(as.matrix(train[,-1]), as.matrix(train[,1]),
                   family = "gaussian", alpha=0 , lambda = lambda.grid)
    
    pre.ridge.val = predict(ridge , newx = as.matrix(val[,-1]))
    mse.ridge.val = numeric(length(pre.ridge.val[1,]))
    
    for (i  in 1:length(pre.ridge.val[1,])){
      mse.ridge.val[i] = mean((pre.ridge.val[,i] - val[,1])^2)
    }
    min.lambda.ridge = which.min(mse.ridge.val) %>% ridge$lambda[.]
    ridge.pred[test.window.obs] = predict(ridge, newx=as.matrix(test[,-1]), 
                                          s = min.lambda.ridge)
    ridge.pred.IS = predict(ridge, newx = as.matrix(full.train[,-1]),
                            s = min.lambda.ridge)
    
    
    
    lm.mse.IS = mean((lm.pred.IS - full.train$Y)^2)
    lm.R2.IS = (1 - lm.mse.IS/(tot.test))*100
    
    oracle.mse.IS = mean((oracle.pred.IS - full.train$Y)^2)
    oracle.R2.IS = (1 - oracle.mse.IS/(tot.test))*100
    
    ridge.mse.IS = mean((ridge.pred.IS - full.train$Y)^2)
    ridge.R2.IS = (1 - ridge.mse.IS/(tot.test))*100
    
    ##############################################################################
    ###  -------------- END OF LOOP CURLY BRACKET ----------------------
    ##############################################################################
  }
  
  
  #######
  # Calculation of MSE and R^2 for all model based on predictions from loop
  
  #######
  
  lm.mse = mean((lm.pred[test.sample] - df$Y[test.sample])^2)
  lm.R2 = (1 - lm.mse/(tot.test))*100
  
  oracle.mse = mean((oracle.pred[test.sample] - df$Y[test.sample])^2)
  oracle.R2 = (1 - oracle.mse/(tot.test))*100
  
  ridge.mse = mean((ridge.pred[test.sample] - df$Y[test.sample])^2)
  ridge.R2 = (1 - ridge.mse/(tot.test) )*100
  
  
  ##################################################################
  ########### ------- END OF SIMULATION ---------------- ###########
  ##################################################################
  return(list("OLS"=c("R2-OOS"=lm.R2, "R2-IS"=lm.R2.IS, "MSE-OOS"=lm.mse, "MSE-IS"=lm.mse.IS), 
              "Ridge"= c("R2-OOS"=ridge.R2, "R2-IS"=ridge.R2.IS, "MSE-OOS"=ridge.mse, "MSE-IS"=ridge.mse.IS), 
              "Oracle" = c("R2-OOS"=oracle.R2, "R2-IS"=oracle.R2.IS, "MSE-OOS"=oracle.mse, "MSE-IS"=oracle.mse.IS )
  ))
  
}

mytest = replicate(3, execute(DGP = "b"))
