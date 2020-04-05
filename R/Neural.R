### Neural Networks
library(parallel)
#library(snowfall)

set.seed(2019)
rm(list=ls())
execute = function(DGP, Pc=50,  N=100, T=90){
  
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
  initial_T=floor(T/2)
  val_T = floor(T/4)
  sample.size = N*T
  test.sample = which(df$time>initial_T + val_T)
  test.sample.size = length(test.sample)
  
  
  train.sample = which(df$time<=initial_T+val_T)
  historical.mean = mean(df$Y[train.sample])
  
  tot.test = mean((df$Y[test.sample] - historical.mean)^2)
  
  tot.test.no.demean = mean(df$Y[test.sample]^2)
  
  lambda.grid = 10^seq(2,-5, length = 100)
  
  
  
  
  #Create the three samples for training, tuning and testing
  train= df %>% filter(time %in% 1:initial_T) %>% select(-c(id, time))
  val= df %>% filter(time %in% (initial_T+1):(initial_T+val_T)) %>% select(-c(id, time))
  test= df %>% 
    filter( time > initial_T+val_T ) %>% 
    select(-c(id,time))
  
  full.train = rbind(train,val)
  
  
  
  ################################################################################
  ## ------ ESTIMATING MODELS FOR THE GIVEN TRAINING AND VAL SAMPLE --------
  ################################################################################
  
  ### ------------------------------------------------------------------
  # Neural Networks
  ##--------------------------------------------------------------------
  library(keras)
  
  x_train = as.matrix(train[,-1])
  y_train = as.matrix(train[,1])
  x_val = as.matrix(val[,-1])
  y_val = as.matrix(val[,1])
  x_test = as.matrix(test[,-1])
  y_test = as.matrix(test[,1])
  
  # Normalize training data
  x_train <- scale(x_train) 
  
  # Use means and standard deviations from training set to normalize test set
  col_means_train <- attr(x_train, "scaled:center") 
  col_stddevs_train <- attr(x_train, "scaled:scale")
  x_val <- scale(x_val, center = col_means_train, scale = col_stddevs_train)
  x_test <- scale(x_test, center = col_means_train, scale = col_stddevs_train)
  
  #Grid for penalty param
  lambda.grid.NN = c(0, 0.0001,  0.01)
  
  
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
  
  
  #######
  ## One layer
  #######
  
  R2.NN1 = numeric(length(lambda.grid.NN))
  pred.NN1.test = list(0)
  NN1.mse = numeric(length(lambda.grid.NN))
  NN1.mse.IS = numeric(length(lambda.grid.NN))
  
  for (i in 1:length(lambda.grid.NN)){
    
    model = build_model(layers = 1)
    model %>% summary()
    
    print_dot_callback <- callback_lambda(
      on_epoch_end = function(epoch, logs) {
        if (epoch %% 80 == 0) cat("\n")
        cat(".")
      }
    )   
    
    # The patience parameter is the amount of epochs to check for improvement.
    early_stop <- callback_early_stopping(monitor = "val_loss", 
                                          patience = 10, min_delta = 0.001, 
                                          mode = "min",
                                          restore_best_weights = TRUE)
    
    history <- model %>% fit(
      x_train,
      y_train,
      epochs = 500,
      batch_size = 100,
      shuffle = FALSE,
      validation_data = list(x_val, y_val),
      verbose = 0,
      callbacks = list(early_stop, print_dot_callback)
    )
    
    plot(history, metrics = "loss", smooth = FALSE) 
    #coord_cartesian(xlim = c(0, 50), ylim = c(0, 3))
    
    pred.NN1.test[[i]] <- model %>% predict(x_test)
    NN1.mse[i] = mean((pred.NN1.test[[i]] - test[,1])^2)
    R2.NN1[i] = (1 - NN1.mse[i]/tot.test)*100
    NN1.mse.IS[i] = history$metrics$loss %>% tail(1)
  }
  
  NN1.pred = pred.NN1.test[[which.min(NN1.mse)]]
  NN1.mse.IS = NN1.mse.IS[which.min(NN1.mse)]
  
  #######
  ## two layerss
  #######
  
  R2.NN2 = numeric(length(lambda.grid.NN))
  pred.NN2.test = list(0)
  NN2.mse = numeric(length(lambda.grid.NN))
  NN2.mse.IS = numeric(length(lambda.grid.NN))
  
  for (i in 1:length(lambda.grid.NN)){
    
    model = build_model(layers = 2)
    model %>% summary()
    
    print_dot_callback <- callback_lambda(
      on_epoch_end = function(epoch, logs) {
        if (epoch %% 80 == 0) cat("\n")
        cat(".")
      }
    )   
    
    # The patience parameter is the amount of epochs to check for improvement.
    early_stop <- callback_early_stopping(monitor = "val_loss", 
                                          patience = 10, min_delta = 0.001, 
                                          mode = "min",
                                          restore_best_weights = TRUE)
    
    history <- model %>% fit(
      x_train,
      y_train,
      epochs = 500,
      batch_size = 100,
      shuffle = FALSE,
      validation_data = list(x_val, y_val),
      verbose = 0,
      callbacks = list(early_stop, print_dot_callback)
    )
    
    plot(history, metrics = "loss", smooth = FALSE) 
    #coord_cartesian(xlim = c(0, 50), ylim = c(0, 3))
    
    pred.NN2.test[[i]] <- model %>% predict(x_test)
    NN2.mse[i] = mean((pred.NN2.test[[i]] - test[,1])^2)
    R2.NN2[i] = (1 - NN2.mse[i]/tot.test)*100
    NN2.mse.IS[i] = history$metrics$loss %>% tail(1)
  }
  
  NN2.pred = pred.NN2.test[[which.min(NN2.mse)]]
  NN2.mse.IS = NN2.mse.IS[which.min(NN2.mse)]
  
  
  
  #######
  ## Three layers
  #######
  
  R2.NN3 = numeric(length(lambda.grid.NN))
  pred.NN3.test = list(0)
  NN3.mse = numeric(length(lambda.grid.NN))
  NN3.mse.IS = numeric(length(lambda.grid.NN))
  
  for (i in 1:length(lambda.grid.NN)){
    
    model = build_model(layers = 3)
    model %>% summary()
    
    print_dot_callback <- callback_lambda(
      on_epoch_end = function(epoch, logs) {
        if (epoch %% 80 == 0) cat("\n")
        cat(".")
      }
    )   
    
    # The patience parameter is the amount of epochs to check for improvement.
    early_stop <- callback_early_stopping(monitor = "val_loss", 
                                          patience = 10, min_delta = 0.001, 
                                          mode = "min",
                                          restore_best_weights = TRUE)
    
    history <- model %>% fit(
      x_train,
      y_train,
      epochs = 500,
      batch_size = 100,
      shuffle = FALSE,
      validation_data = list(x_val, y_val),
      verbose = 0,
      callbacks = list(early_stop, print_dot_callback)
    )
    
    pred.NN3.test[[i]] <- model %>% predict(x_test)
    NN3.mse[i] = mean((pred.NN3.test[[i]] - test[,1])^2)
    R2.NN3[i] = (1 - NN3.mse[i]/tot.test)*100
    NN3.mse.IS[i] = history$metrics$loss %>% tail(1)
  }
  
  NN3.pred = pred.NN3.test[[which.min(NN3.mse)]]
  NN3.mse.IS = NN3.mse.IS[which.min(NN3.mse)]

  
  ## performance
  
  NN1.R2.IS = (1 - NN1.mse.IS/(tot.test))*100
  
  NN2.R2.IS = (1 - NN2.mse.IS/(tot.test))*100
  
  NN3.R2.IS = (1 - NN3.mse.IS/(tot.test))*100
  

  
  NN1.mse = mean((NN1.pred - df$Y[test.sample])^2)
  NN1.R2 = (1 - NN1.mse/(tot.test) )*100
  
  NN2.mse = mean((NN2.pred - df$Y[test.sample])^2)
  NN2.R2 = (1 - NN2.mse/(tot.test) )*100
  
  NN3.mse = mean((NN3.pred - df$Y[test.sample])^2)
  NN3.R2 = (1 - NN3.mse/(tot.test) )*100
  
  

  
  
  ##################################################################
  ########### ------- END OF EXCUTE ---------------- ###########
  ##################################################################
  return(list("NN1" = c("R2-OOS"=NN1.R2, "R2-IS"=NN1.R2.IS, "MSE-OOS"=NN1.mse, "MSE-IS"=NN1.mse.IS ),
              "NN2" = c("R2-OOS"=NN2.R2, "R2-IS"=NN2.R2.IS, "MSE-OOS"=NN2.mse, "MSE-IS"=NN2.mse.IS ),
              "NN3" = c("R2-OOS"=NN3.R2, "R2-IS"=NN3.R2.IS, "MSE-OOS"=NN3.mse, "MSE-IS"=NN3.mse.IS )
  ))
  
}



library(pbapply)



clust = makeCluster(no_cores)

#DGP (a) and Pc = 50
MC.a.50 = pblapply(rep("a", 50), execute , cl = clust) 

#DGP (b) and Pc = 50
MC.b.50 = pblapply(rep("b", 50), execute , cl = clust)

#DGP (a) and Pc = 100
MC.a.100 = pblapply(rep("a", 50), execute , Pc = 100, cl = clust)

#DGP (b) and Pc = 100
MC.b.100 = pblapply(rep("b", 50), execute , Pc = 100, cl = clust )


stopCluster(clust)

setwd("C:/Users/bech/Documents/Dropbox/Project")
save.image(file = "Neural.RData")

















