## MC Quarterly

library(parallel)
#library(snowfall)

set.seed(2019)

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
  
  ####################################################################  
  ### Now I will change the data frame to do quarterly predictions..
  ####################################################################
  
  df = df %>% mutate(Y = Y + lead(Y, N) + lead(Y,2*N))  
  
  df = df[complete.cases(df),]
  
  T = max(df$time)
  
  
  initial_T=floor(T/2)
  val_T = floor(T/4)
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
    lm.mse.IS = mean((lm.pred.IS - full.train$Y)^2)
    lm.R2.IS = (1 - lm.mse.IS/(tot.test))*100
    
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
    
    
    ### ------------------------------------------------------------------
    #Lasso regression
    ##--------------------------------------------------------------------
    lasso = glmnet(as.matrix(train[,-1]), as.matrix(train[,1]),
                   family = "gaussian", alpha=1, lambda = lambda.grid)
    #x
    #plot(x)
    #coef(x)[,45]
    pre.lasso.val = predict(lasso , newx = as.matrix(val[,-1]), type = "link" )
    mse.lasso.val = numeric(length(pre.lasso.val[1,]))
    for (i  in 1:length(pre.lasso.val[1,])){
      mse.lasso.val[i] = mean((pre.lasso.val[,i] - val[,1])^2)
    }
    min.lambda.lasso = which.min(mse.lasso.val) %>% lasso$lambda[.]
    lasso.pred[test.window.obs] = predict(lasso, newx = as.matrix(test[,-1]),
                                          s = min.lambda.lasso)
    lasso.pred.IS = predict(lasso, newx = as.matrix(full.train[,-1]),
                            s = min.lambda.lasso)
    ### ------------------------------------------------------------------
    #Parameter search for alpha and lambda in elastic net.
    ##--------------------------------------------------------------------
    
    al <- seq(0,1,length.out = 5)
    lol = lapply(al, 
                 function(i) glmnet(as.matrix(train[,-1]), as.matrix(train[,1]), 
                                    alpha=i, lambda = lambda.grid))
    what = lapply(lol  ,
                  function(i) predict(i, newx=as.matrix(val[,-1])))
    what2 = sapply(what, function(i) colMeans((i - val[,1])^2) %>% 
                     which.min(.) )
    
    mse.enet.vect = numeric(length(al))
    for (i in 1:length(al)){
      fit = lol[[i]]$lambda[what2[[i]]]
      mse.enet.vect[i] = mean((predict(lol[[i]], s=fit, newx=as.matrix(val[,-1])) - val[,1])^2)
    }
    which.model =which.min(mse.enet.vect)
    enet.alpha = which.model %>% al[.]
    enet.lambda = lol[[which.model]]$lambda[ what2[[which.model]] ] 
    
    ##calculate the mse on the test sample ... 
    glmnet.model = glmnet(as.matrix(train[,-1]), as.matrix(train[,1]),
                          family = "gaussian", 
                          alpha=enet.alpha , lambda = enet.lambda)
    
    enet.pred[test.window.obs] = predict(glmnet.model, newx = as.matrix(test[,-1]))
    enet.pred.IS = predict(glmnet.model, newx = as.matrix(full.train[,-1]))
    
    ### ------------------------------------------------------------------
    #Parameter search PLS and PCR. Predictor averaging techniques. 
    ##--------------------------------------------------------------------
    pcr.fit = pcr(Y ~ . , data = train ,  scale = TRUE, 
                  validation = "none")
    mse.pcr = numeric(pcr.fit$ncomp)
    for (i in 1:pcr.fit$ncomp){
      pcr.pred.val = predict(pcr.fit, val[,-1], ncomp = i)
      mse.pcr[i] = mean((pcr.pred.val - val[,1])^2)
    }
    pcr.comp = which.min(mse.pcr)
    pcr.pred[test.window.obs] = predict(pcr.fit, test[,-1], ncomp = pcr.comp)
    pcr.pred.IS = predict(pcr.fit, full.train[,-1], ncomp = pcr.comp)
    
    
    plsr.fit = plsr(Y ~ . , data = train, scale = TRUE,
                    validation = "none")
    mse.plsr.val = numeric(plsr.fit$ncomp)
    for (i in 1:plsr.fit$ncomp){
      plsr.pred.val = predict(plsr.fit, val[,-1], ncomp = i)
      mse.plsr.val[i] = mean((plsr.pred.val - val[,1])^2)
    }
    plsr.comp = which.min(mse.plsr.val)
    plsr.pred[test.window.obs] = predict(plsr.fit , test[,-1], ncomp = plsr.comp)
    plsr.pred.IS = predict(plsr.fit, full.train[,-1], ncomp = pcr.comp)
    
    
    ### ------------------------------------------------------------------
    #Generalized linear model (GAM) penalized by group lasso.
    ##--------------------------------------------------------------------
    num.pred = length(train[1,-1])
    lambda.grid.glm = 10^seq(3,-7, length = 100)
    
    predictor.matrix.train = lapply(2:5, function(x) lapply(seq(2,num.pred+1), 
                                                            function(i) bs(train[,i], degree = 2, df = x)) %>% 
                                      do.call(cbind, .) )
    
    predictor.matrix.val = lapply(2:5, function(x) lapply(seq(2,num.pred+1), 
                                                          function(i) bs(val[,i], degree = 2, df = x)) %>% 
                                    do.call(cbind, .) )
    predictor.matrix.test = lapply(2:5, function(x) lapply(seq(2,num.pred+1), 
                                                           function(i) bs(test[,i], degree = 2, df = x)) %>% 
                                     do.call(cbind, .) )
    predictor.matrix.fulltrain = lapply(2:5, function(x) lapply(seq(2,num.pred+1), 
                                                                function(i) bs(full.train[,i], degree = 2, df = x)) %>% 
                                          do.call(cbind, .) )
    
    groups = lapply(2:5, function(i) rep(1:num.pred, each = i
    ) )
    
    ## TRain Models
    
    glm.models = lapply(1:4, function(i) gglasso(predictor.matrix.train[[i]], 
                                                 as.matrix(train[,1]),
                                                 group = groups[[i]], 
                                                 loss = "ls", lambda = lambda.grid.glm) )
    
    #plot(glm.models[[9]])
    #coef(glm.models[[9]], s = lambda.grid.glm[4])
    
    ## Which lambda minimizes the mse for each of the models with different df (no. of knots)
    which.lambda.glm =   sapply(1:length(glm.models)  , 
                                function(i) predict(glm.models[[i]], newx=predictor.matrix.val[[i]],type="link") %>% 
                                  colMeans((. - val[,1])^2) %>%
                                  which.min() ) 
    
    #glm.models[[i]]$lambda[which.lambda.glm[[i]]]
    
    ## Evaluate on  validation set
    glm.mse.val = numeric(length(glm.models))
    for (i in 1:length(glm.models)){
      glm.pred.val = predict(glm.models[[i]],
                             newx=predictor.matrix.val[[i]], type="link",
                             s = glm.models[[i]]$lambda[which.lambda.glm[[i]]])
      glm.mse.val[i] = mean((glm.pred.val - val[,1])^2)
    }
    #plot(glm.mse.val)
    best.model = which.min(glm.mse.val)
    
    glm.pred[test.window.obs] = predict(glm.models[[best.model]], newx=predictor.matrix.test[[best.model]], type="link",
                                        s = glm.models[[best.model]]$lambda[which.lambda.glm[[best.model]]] )
    
    glm.pred.IS = predict(glm.models[[best.model]], 
                          newx=predictor.matrix.fulltrain[[best.model]], 
                          type="link",
                          s = glm.models[[best.model]]$lambda[which.lambda.glm[[best.model]]] )
    
    ### ------------------------------------------------------------------
    #Tree based models (boosting and bagging).
    ##--------------------------------------------------------------------
    
    ### Random Forest
    library(ranger)
    no.pred = dim(train)[2] - 1
    
    RF_grid <- expand.grid(
      mtry       = c(floor(sqrt(no.pred)), floor(no.pred/3)),
      node_size  = c(10, 1000 ),
      OOB_RMSE  = 0
    )
    
    ## Use all data as we can use OOB error to evaluate tuning parameters
    
    train.rf = rbind(train,val)
    
    for(i in 1:nrow(RF_grid)) {
      
      # train model
      model <- ranger(
        formula         = Y ~ ., 
        data            = train.rf, 
        num.trees       = 500,
        mtry            = RF_grid$mtry[i],
        min.node.size   = RF_grid$node_size[i]
      )
      
      # add OOB error to grid
      RF_grid$OOB_RMSE[i] <- sqrt(model$prediction.error)
    }
    
    best.model.RF = which.min(RF_grid$OOB_RMSE) %>% RF_grid[.,]
    
    
    model.RF = ranger(
      formula         = Y ~ ., 
      data            = train.rf, 
      num.trees       = 500,
      mtry            = best.model.RF$mtry,
      min.node.size   = best.model.RF$node_size
    )
    
    randomF.pred[test.window.obs] = predict(model.RF, test[,-1])$predictions
    randomF.pred.IS = predict(model.RF, full.train[,-1])$predictions
    
    ## GRadient Boosted regression trees
    
    library(gbm)
    
    boost_grid <- expand.grid(
      depth       = c(1,4),
      shrink  = c(0.01, 0.2)
    )
    boost_grid_list = apply(boost_grid, MARGIN = 1, as.list)
    
    
    model.boosting = lapply(boost_grid_list, function(x)
      gbm(Y~. , distribution = "gaussian", data = train,
          interaction.depth = x$depth,
          shrinkage = x$shrink,
          n.trees =500 ) )
    
    boost.pred.val = lapply(model.boosting , 
                            function(x) predict(x, newdata=val[,-1], n.trees=500))
    boost.pred.mse = sapply(boost.pred.val, function(x) mean(( x - val$Y)^2) )   
    
    opt.boost = which.min(boost.pred.mse)
    
    boost.pred[test.window.obs] = predict(model.boosting[[opt.boost]], newdata = test[,-1], n.trees =500)
    boost.pred.IS = predict(model.boosting[[opt.boost]], newdata = full.train[,-1], n.trees =500)
    
    
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
                                            patience = 2, 
                                            mode = "min",
                                            restore_best_weights = TRUE)
      
      history <- model %>% fit(
        x_train,
        y_train,
        epochs = 500,
        batch_size = 64,
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
    
    NN1.pred[test.window.obs] = pred.NN1.test[[which.min(NN1.mse)]]
    NN1.mse.IS = NN1.mse.IS[which.min(NN1.mse)]
    
    
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
                                            patience = 2, 
                                            mode = "min",
                                            restore_best_weights = TRUE)
      
      history <- model %>% fit(
        x_train,
        y_train,
        epochs = 500,
        batch_size = 64,
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
    
    NN3.pred[test.window.obs] = pred.NN3.test[[which.min(NN3.mse)]]
    NN3.mse.IS = NN3.mse.IS[which.min(NN3.mse)]
    
    
    ## In sample numbers
    
    lm.mse.IS = mean((lm.pred.IS - full.train$Y)^2)
    lm.R2.IS = (1 - lm.mse.IS/(tot.test))*100
    
    oracle.mse.IS = mean((oracle.pred.IS - full.train$Y)^2)
    oracle.R2.IS = (1 - oracle.mse.IS/(tot.test))*100
    
    ridge.mse.IS = mean((ridge.pred.IS - full.train$Y)^2)
    ridge.R2.IS = (1 - ridge.mse.IS/(tot.test))*100
    
    lasso.mse.IS = mean((lasso.pred.IS - full.train$Y)^2)
    lasso.R2.IS = (1 - lasso.mse.IS/(tot.test))*100
    
    enet.mse.IS = mean((enet.pred.IS - full.train$Y)^2)
    enet.R2.IS = (1 - enet.mse.IS/(tot.test))*100
    
    plsr.mse.IS = mean((plsr.pred.IS - full.train$Y)^2)
    plsr.R2.IS = (1 - plsr.mse.IS/(tot.test))*100
    
    pcr.mse.IS = mean((pcr.pred.IS - full.train$Y)^2)
    pcr.R2.IS = (1 - pcr.mse.IS/(tot.test))*100
    
    glm.mse.IS = mean((glm.pred.IS - full.train$Y)^2)
    glm.R2.IS = (1 - glm.mse.IS/(tot.test))*100
    
    randomF.mse.IS = mean((randomF.pred.IS - full.train$Y)^2)
    randomF.R2.IS = (1 - randomF.mse.IS/(tot.test))*100
    
    boost.mse.IS = mean((boost.pred.IS - full.train$Y)^2)
    boost.R2.IS = (1 - boost.mse.IS/(tot.test))*100
    
    NN1.R2.IS = (1 - NN1.mse.IS/(tot.test))*100
    
    
    NN3.R2.IS = (1 - NN3.mse.IS/(tot.test))*100
    
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
  
  lasso.mse = mean((lasso.pred[test.sample] - df$Y[test.sample])^2)
  lasso.R2 = (1 - lasso.mse/(tot.test) )*100
  
  enet.mse = mean((enet.pred[test.sample] - df$Y[test.sample])^2)
  enet.R2 = (1 - enet.mse/(tot.test) )*100
  
  pcr.mse = mean((pcr.pred[test.sample] - df$Y[test.sample])^2)
  pcr.R2 = (1 - pcr.mse/(tot.test) )*100
  
  plsr.mse = mean((plsr.pred[test.sample] - df$Y[test.sample])^2)
  plsr.R2 = (1 - plsr.mse/(tot.test) )*100
  
  glm.mse = mean((glm.pred[test.sample] - df$Y[test.sample])^2)
  glm.R2 = (1 - glm.mse/(tot.test) )*100
  
  randomF.mse = mean((randomF.pred[test.sample] - df$Y[test.sample])^2)
  randomF.R2 = (1 - randomF.mse/(tot.test) )*100
  
  boost.mse = mean((boost.pred[test.sample] - df$Y[test.sample])^2)
  boost.R2 = (1 - boost.mse/(tot.test) )*100
  
  NN1.mse = mean((NN1.pred[test.sample] - df$Y[test.sample])^2)
  NN1.R2 = (1 - NN1.mse/(tot.test) )*100
  
  NN2.mse = mean((NN2.pred[test.sample] - df$Y[test.sample])^2)
  NN2.R2 = (1 - NN2.mse/(tot.test) )*100
  
  NN3.mse = mean((NN3.pred[test.sample] - df$Y[test.sample])^2)
  NN3.R2 = (1 - NN3.mse/(tot.test) )*100
  
  ##################################################################
  ########### ------- END OF SIMULATION ---------------- ###########
  ##################################################################
  return(list("OLS"=c("R2-OOS"=lm.R2, "R2-IS"=lm.R2.IS, "MSE-OOS"=lm.mse, "MSE-IS"=lm.mse.IS), 
              "Ridge"= c("R2-OOS"=ridge.R2, "R2-IS"=ridge.R2.IS, "MSE-OOS"=ridge.mse, "MSE-IS"=ridge.mse.IS), 
              "Lasso"= c("R2-OOS"=lasso.R2, "R2-IS"=lasso.R2.IS, "MSE-OOS"=lasso.mse, "MSE-IS"=lasso.mse.IS ),
              "ENet" = c("R2-OOS"=enet.R2, "R2-IS"=enet.R2.IS, "MSE-OOS"=enet.mse, "MSE-IS"=enet.mse.IS ),
              "PCR" = c("R2-OOS"=pcr.R2, "R2-IS"=pcr.R2.IS, "MSE-OOS"=pcr.mse, "MSE-IS"=pcr.mse.IS ),
              "PLSR" = c("R2-OOS"=plsr.R2, "R2-IS"=plsr.R2.IS, "MSE-OOS"=plsr.mse, "MSE-IS"=plsr.mse.IS ),
              "GLM" = c("R2-OOS"=glm.R2, "R2-IS"=glm.R2.IS, "MSE-OOS"=glm.mse, "MSE-IS"=glm.mse.IS ),
              "RandomF"= c("R2-OOS"=randomF.R2, "R2-IS"=randomF.R2.IS, "MSE-OOS"=randomF.mse, "MSE-IS"=randomF.mse.IS ),
              "GBRT" = c("R2-OOS"=boost.R2, "R2-IS"=boost.R2.IS, "MSE-OOS"=boost.mse, "MSE-IS"=boost.mse.IS ), 
              "NN1" = c("R2-OOS"=NN1.R2, "R2-IS"=NN1.R2.IS, "MSE-OOS"=NN1.mse, "MSE-IS"=NN1.mse.IS ),
              "NN2" = c("R2-OOS"=NN2.R2, "R2-IS"=NN2.R2.IS, "MSE-OOS"=NN2.mse, "MSE-IS"=NN2.mse.IS ),
              "NN3" = c("R2-OOS"=NN3.R2, "R2-IS"=NN3.R2.IS, "MSE-OOS"=NN3.mse, "MSE-IS"=NN3.mse.IS ),
              "Oracle" = c("R2-OOS"=oracle.R2, "R2-IS"=oracle.R2.IS, "MSE-OOS"=oracle.mse, "MSE-IS"=oracle.mse.IS )
  ))
  
}

### USing the Parallel package to utilize all cores on the machine

#TEST RUN
#no_cores = detectCores()
#clust = makeCluster(no_cores)
#MC.b.50 = parLapply(clust, rep("b", 10), execute , Pc = 50 )
#stopCluster(clust)





#########  
## For Markdown
#########
no_cores = detectCores()
clust = makeCluster(no_cores)

#DGP (a) and Pc = 50
MC.a.50.quart = parLapply(clust, rep("a", 100), execute )

#DGP (b) and Pc = 50
MC.b.50.quart = parLapply(clust, rep("b", 100), execute )

#DGP (a) and Pc = 100
MC.a.100.quart = parLapply(clust, rep("a", 100), execute , Pc = 100)

#DGP (b) and Pc = 100
MC.b.100.quart = parLapply(clust, rep("b", 100), execute , Pc = 100)

stopCluster(clust)


