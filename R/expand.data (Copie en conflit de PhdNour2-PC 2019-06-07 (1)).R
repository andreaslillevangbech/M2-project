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
  library(Hmisc)
  
  rm(list=ls())
  
  library(readxl)
  ## Using PredictorData2017 from Amit Goyals' website (University of lausanne)
  ## data is from paper ("http://www.hec.unil.ch/agoyal/")
  ## "A Comprehensive Look at The Empirical Performance of Equity Premium Prediction"
  ## Data is updated up to 2017
  ## 
  monthly <- 
    read_excel("Dropbox/Project/PredictorData2017.xlsx" , sheet = "Monthly",
               na = "NaN")
  quart <- 
    read_excel("Dropbox/Project/PredictorData2017.xlsx" , 
               sheet = "Quarterly",
               na = "NaN")
  annual <- 
    read_excel("Dropbox/Project/PredictorData2017.xlsx" , sheet = "Annual",
               na = "NaN")
  
  monthly$yyyymm = monthly$yyyymm %>% as.character()
  quart$yyyyq = quart$yyyyq %>% as.character()
  annual$yyyy = annual$yyyy %>% as.character()
  
  # Adding quarters to monthly
  df = select(quart, -c(names(monthly)[-1], yyyyq)) %>% 
    slice(rep(1:n(), each=3)) %>% 
    apply(., MARGIN = 2, dplyr::lag, n=2) %>% 
    cbind(monthly, . ) %>% as_tibble
  
  # Adding Annual (only the eqis variable is missing)
  df$eqis = rep(annual$eqis, each=12) %>% lag(n=11)
  
  # Gonna Remove cay, E3, D3, cause they dont really belong
  df = select(df,-c(cay, D3, E3,csp))
  
  #Now I have
  names(df)
  # Index"       "D12"        "E12"        "b/m"
  # "tbl"        "AAA"        "BAA"       
  # "lty"        "ntis"       "Rfree"      "infl"      
  # "ltr"        "corpr"      "svar"       
  # "CRSP_SPvw"  "CRSP_SPvwx" "ik"         "eqis"
  
  # Create a date variable
  df$date = ymd(paste0(df$yyyymm, "01"))
  
  ## Creating new variables
  df = mutate(df, dp = log(D12) - log(Index))
  df = mutate(df, dy = log(D12) - lag(log(Index)))
  df = mutate(df, ep = log(E12) - log(Index))
  df = mutate(df, dfy = BAA - AAA)
  df = mutate(df, dfr = corpr - ltr)
  ## The returns, notice timing is such that r_t+1 is next to p_t
  df = mutate(df, r = c(diff(log(Index)), NA))
  ## Creating momentum variables
  df$mom3 = c(rep(NA,3), embed(df$r, 4)[,2:4] %>% rowSums())
  df$mom6 = c(rep(NA,6), embed(df$r, 7)[,2:7] %>% rowSums())
  df$mom612 = c(rep(NA,12), embed(df$r, 13)[,2:13] %>% rowSums())
  
  # full availability from 1947-03-01 , i/k not available before.
  plot(df$date, complete.cases(df),type = "l") 
  
  # Saw off dataframe to omit all NA rows.
  df = df[complete.cases(df),]
  #So the sample runs from 1947-03-01 to 2017-06-01
  
  ## Now remove superflous variables so not to have multicoll
  df = select(df, -c(date, Index, yyyymm, D12,
                        E12, AAA, BAA, corpr, CRSP_SPvwx,
                        CRSP_SPvw, Rfree))
  
  df = select(df, r, everything())
  names(df)[names(df)=="b/m"] = "bm"
  T= nrow(df)
  df$time=1:T
  
  df = as.data.frame(df)
  
  
  
#expanding_window = function(){
  
    
  
  initial_T=floor(T/2)
  val_T=floor(T/4)
  
  sample.size = T
  test.sample = which(df$time>T-val_T)
  test.sample.size = length(df$time[df$time > (T-val_T)])
  
  train.sample = which(df$time<=T-val_T)
  historical.mean = mean(df$r[train.sample])
  
  tot.test = mean((df$r[test.sample] - historical.mean)^2)
  
  tot.test.no.demean = mean(df$r[test.sample]^2)
  
  
  
  
  ########################################################################
  ## FOr use in loop
  ########################################################################
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
  lm.pred = numeric(T)
  ridge.pred = numeric(T)
  lasso.pred = numeric(T)
  enet.pred = numeric(T)
  pcr.pred = numeric(T)
  plsr.pred = numeric(T)
  glm.pred = numeric(T)
  randomF.pred = numeric(T)
  boost.pred = numeric(T)
  NN1.pred = numeric(T)
  NN2.pred = numeric(T)
  NN3.pred = numeric(T)
  
  ## Lists for model fits for each model
  ridge.list = list(0)
  lasso.list = list(0)
  enet.list = list(0)
  pcr.list = list(0)
  plsr.list = list(0)
  glm.list = list(0)
  RF.list = list(0)
  GBRT.list = list(0)
  NN1.list = list(0)
  NN2.list = list(0)
  NN3.list = list(0)
  
  ########################################################################
  ################## ------- STARTING LOOP ----------- ###################
  ########################################################################
  
  
  ### Loop it over
  seq = seq(initial_T, T - val_T -1, by = 12)
  no.of.models = length(seq)
  
  
  #Creating a loop over the periods to reestimate the model. 
  #Using expanding window like in the paper
  for (i in 1:length(seq)){
    
    t = seq[i]
    #The current window in terms of time index as well as the corresponding observations
    test.window = seq(t+val_T+1, ifelse(t+val_T+12<= T, t+val_T+12, T) )
    test.window.obs = which(df$time %in% test.window)
    
    #Create the three samples for training, tuning and testing
    train= df %>% filter(time %in% 1:t) %>% select(-time)
    val= df %>% filter(time %in% (t+1):(t+val_T)) %>% select(-time)
    test= df %>% 
      filter( time %in% test.window ) %>% 
      select(-time)
    
    full.train = rbind(train,val)
    
    ################################################################################
    ## ------ ESTIMATING MODELS FOR THE GIVEN TRAINING AND VAL SAMPLE --------
    ################################################################################
    
    ### ------------------------------------------------------------------
    # Standard Linear Model
    ###--------------------------------------------------------------------
    fit = lm(r ~ . , data=rbind(train,val))
    lm.pred[test.window.obs] = predict(fit, newdata = test[,-1])
    
    ## IS
    lm.pred.IS = predict(fit, newdata = full.train[,-1])
    lm.mse.IS = mean((lm.pred.IS - full.train$r)^2)
    lm.R2.IS = (1 - lm.mse.IS/(tot.test))*100
    
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
    ridge.list[[i]] = list("model" = ridge, "lambda" = min.lambda.ridge)
    
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
    
    lasso.list[[i]] = list("model" = lasso, "lambda" = min.lambda.lasso)
    
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
    
    enet.list[[i]] = list("model" = glmnet.model,
                          "lambda" = enet.lambda, "alpha"=enet.alpha)
    ### ------------------------------------------------------------------
    #Parameter search PLS and PCR. Predictor averaging techniques. 
    ##--------------------------------------------------------------------
    pcr.fit = pcr(r ~ . , data = train ,  scale = TRUE, 
                  validation = "none")
    mse.pcr = numeric(pcr.fit$ncomp)
    for (i in 1:pcr.fit$ncomp){
      pcr.pred.val = predict(pcr.fit, val[,-1], ncomp = i)
      mse.pcr[i] = mean((pcr.pred.val - val[,1])^2)
    }
    pcr.comp = which.min(mse.pcr)
    pcr.pred[test.window.obs] = predict(pcr.fit, test[,-1], ncomp = pcr.comp)
    pcr.pred.IS = predict(pcr.fit, full.train[,-1], ncomp = pcr.comp)
    
    
    plsr.fit = plsr(r ~ . , data = train, scale = TRUE,
                    validation = "none")
    mse.plsr.val = numeric(plsr.fit$ncomp)
    for (i in 1:plsr.fit$ncomp){
      plsr.pred.val = predict(plsr.fit, val[,-1], ncomp = i)
      mse.plsr.val[i] = mean((plsr.pred.val - val[,1])^2)
    }
    plsr.comp = which.min(mse.plsr.val)
    plsr.pred[test.window.obs] = predict(plsr.fit , test[,-1], ncomp = plsr.comp)
    plsr.pred.IS = predict(plsr.fit, full.train[,-1], ncomp = pcr.comp)
    
    pcr.list[[i]] = list("model" = pcr.fit, "ncomp" = pcr.comp)
    plsr.list[[i]] = list("model" = plsr.fit, "ncomp" = plsr.comp)
    
    
    ### ------------------------------------------------------------------
    #Generalized linear model (GAM) penalized by group lasso.
    ##--------------------------------------------------------------------
    num.pred = length(train[1,-1])
    lambda.grid.glm = 10^seq(3,-7, length = 100)
    
    predictor.matrix.train = lapply(2:10, function(x) lapply(seq(2,num.pred+1), 
                                                             function(i) bs(train[,i], degree = 2, df = x)) %>% 
                                      do.call(cbind, .) )
    
    predictor.matrix.val = lapply(2:10, function(x) lapply(seq(2,num.pred+1), 
                                                           function(i) bs(val[,i], degree = 2, df = x)) %>% 
                                    do.call(cbind, .) )
    predictor.matrix.test = lapply(2:10, function(x) lapply(seq(2,num.pred+1), 
                                                            function(i) bs(test[,i], degree = 2, df = x)) %>% 
                                     do.call(cbind, .) )
    predictor.matrix.fulltrain = lapply(2:10, function(x) lapply(seq(2,num.pred+1), 
                                                                 function(i) bs(full.train[,i], degree = 2, df = x)) %>% 
                                          do.call(cbind, .) )
    
    groups = lapply(2:10, function(i) rep(1:num.pred, each = i
    ) )
    
    ## TRain Models
    
    glm.models = lapply(1:9, function(i) gglasso(predictor.matrix.train[[i]], 
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
      node_size  = c(10,50,100,200,300,400,500,600,700),
      OOB_RMSE  = 0
    )
    
    ## Use all data as we can use OOB error to evaluate tuning parameters
    
    train.rf = rbind(train,val)
    
    for(i in 1:nrow(RF_grid)) {
      
      # train model
      model <- ranger(
        formula         = r ~ ., 
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
      formula         = r ~ ., 
      data            = train.rf, 
      num.trees       = 500,
      mtry            = best.model.RF$mtry,
      min.node.size   = best.model.RF$node_size,
      importance      = 'impurity'
    )
    
    randomF.pred[test.window.obs] = predict(model.RF, test[,-1])$predictions
    randomF.pred.IS = predict(model.RF, full.train[,-1])$predictions
    
    RF.list[[i]] = list("model" = model.RF)
    
    ## GRadient Boosted regression trees
    
    library(gbm)
    
    boost_grid <- expand.grid(
      depth       = seq(1,4),
      shrink  = c(0.01, 0.2)
    )
    boost_grid_list = apply(boost_grid, MARGIN = 1, as.list)
    
    
    model.boosting = lapply(boost_grid_list, function(x)
      gbm(r~. , distribution = "gaussian", data = train,
          interaction.depth = x$depth,
          shrinkage = x$shrink,
          n.trees =500 ) )
    
    boost.pred.val = lapply(model.boosting , 
                            function(x) predict(x, newdata=val[,-1], n.trees=500))
    boost.pred.mse = sapply(boost.pred.val, function(x) mean(( x - val$r)^2) )   
    
    opt.boost = which.min(boost.pred.mse)
    
    boost.pred[test.window.obs] = predict(model.boosting[[opt.boost]], newdata = test[,-1], n.trees =500)
    boost.pred.IS = predict(model.boosting[[opt.boost]], newdata = full.train[,-1], n.trees =500)
    
    var.select.boost = sapply(1:500, function(x) pretty.gbm.tree(model.boosting[[opt.boost]],
                                              i.tree = x)[1,1])
    var.select.boost = table(var.select.boost) %>% as_tibble
    
    GBRT.list[[i]] = list("model" = model.boosting[[opt.boost]],
                          "pretty.tree" = var.select.boost)
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
    ## Two layers
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
      
      pred.NN2.test[[i]] <- model %>% predict(x_test)
      NN2.mse[i] = mean((pred.NN2.test[[i]] - test[,1])^2)
      R2.NN2[i] = (1 - NN2.mse[i]/tot.test)*100
      NN2.mse.IS[i] = history$metrics$loss %>% tail(1)
    }
    
    NN2.pred[test.window.obs] = pred.NN2.test[[which.min(NN2.mse)]]
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
    
    lm.mse.IS = mean((lm.pred.IS - full.train$r)^2)
    lm.R2.IS = (1 - lm.mse.IS/(tot.test))*100
    
    ridge.mse.IS = mean((ridge.pred.IS - full.train$r)^2)
    ridge.R2.IS = (1 - ridge.mse.IS/(tot.test))*100
    
    lasso.mse.IS = mean((lasso.pred.IS - full.train$r)^2)
    lasso.R2.IS = (1 - lasso.mse.IS/(tot.test))*100
    
    enet.mse.IS = mean((enet.pred.IS - full.train$r)^2)
    enet.R2.IS = (1 - enet.mse.IS/(tot.test))*100
    
    plsr.mse.IS = mean((plsr.pred.IS - full.train$r)^2)
    plsr.R2.IS = (1 - plsr.mse.IS/(tot.test))*100
    
    pcr.mse.IS = mean((pcr.pred.IS - full.train$r)^2)
    pcr.R2.IS = (1 - pcr.mse.IS/(tot.test))*100
    
    glm.mse.IS = mean((glm.pred.IS - full.train$r)^2)
    glm.R2.IS = (1 - glm.mse.IS/(tot.test))*100
    
    randomF.mse.IS = mean((randomF.pred.IS - full.train$r)^2)
    randomF.R2.IS = (1 - randomF.mse.IS/(tot.test))*100
    
    boost.mse.IS = mean((boost.pred.IS - full.train$r)^2)
    boost.R2.IS = (1 - boost.mse.IS/(tot.test))*100
    
    NN1.R2.IS = (1 - NN1.mse.IS/(tot.test))*100
    
    NN2.R2.IS = (1 - NN2.mse.IS/(tot.test))*100
    
    NN3.R2.IS = (1 - NN3.mse.IS/(tot.test))*100
    
    ##############################################################################
    ###  -------------- END OF LOOP CURLY BRACKET ----------------------
    ##############################################################################
  }
  
  #######
  # Calculation of MSE and R^2 for all model based on predictions from loop
  #######
  
  lm.mse = mean((lm.pred[test.sample] - df$r[test.sample])^2)
  lm.R2 = (1 - lm.mse/(tot.test))*100

  ridge.mse = mean((ridge.pred[test.sample] - df$r[test.sample])^2)
  ridge.R2 = (1 - ridge.mse/(tot.test) )*100
  
  lasso.mse = mean((lasso.pred[test.sample] - df$r[test.sample])^2)
  lasso.R2 = (1 - lasso.mse/(tot.test) )*100
  
  enet.mse = mean((enet.pred[test.sample] - df$r[test.sample])^2)
  enet.R2 = (1 - enet.mse/(tot.test) )*100
  
  pcr.mse = mean((pcr.pred[test.sample] - df$r[test.sample])^2)
  pcr.R2 = (1 - pcr.mse/(tot.test) )*100
  
  plsr.mse = mean((plsr.pred[test.sample] - df$r[test.sample])^2)
  plsr.R2 = (1 - plsr.mse/(tot.test) )*100
  
  glm.mse = mean((glm.pred[test.sample] - df$r[test.sample])^2)
  glm.R2 = (1 - glm.mse/(tot.test) )*100
  
  randomF.mse = mean((randomF.pred[test.sample] - df$r[test.sample])^2)
  randomF.R2 = (1 - randomF.mse/(tot.test) )*100
  
  boost.mse = mean((boost.pred[test.sample] - df$r[test.sample])^2)
  boost.R2 = (1 - boost.mse/(tot.test) )*100
  
  NN1.mse = mean((NN1.pred[test.sample] - df$r[test.sample])^2)
  NN1.R2 = (1 - NN1.mse/(tot.test) )*100
  
  NN2.mse = mean((NN2.pred[test.sample] - df$r[test.sample])^2)
  NN2.R2 = (1 - NN2.mse/(tot.test) )*100
  
  NN3.mse = mean((NN3.pred[test.sample] - df$r[test.sample])^2)
  NN3.R2 = (1 - NN3.mse/(tot.test) )*100
  
  ##################################################################
  ########### ------- END OF SIMULATION ---------------- ###########
  ##################################################################
#  return(list("OLS"=c("R2-OOS"=lm.R2, "R2-IS"=lm.R2.IS, "MSE-OOS"=lm.mse, "MSE-IS"=lm.mse.IS), 
#              "Ridge"= c("R2-OOS"=ridge.R2, "R2-IS"=ridge.R2.IS, "MSE-OOS"=ridge.mse, "MSE-IS"=ridge.mse.IS), 
#              "Lasso"= c("R2-OOS"=lasso.R2, "R2-IS"=lasso.R2.IS, "MSE-OOS"=lasso.mse, "MSE-IS"=lasso.mse.IS ),
#              "ENet" = c("R2-OOS"=enet.R2, "R2-IS"=enet.R2.IS, "MSE-OOS"=enet.mse, "MSE-IS"=enet.mse.IS ),
#              "PCR" = c("R2-OOS"=pcr.R2, "R2-IS"=pcr.R2.IS, "MSE-OOS"=pcr.mse, "MSE-IS"=pcr.mse.IS ),
#              "PLSR" = c("R2-OOS"=plsr.R2, "R2-IS"=plsr.R2.IS, "MSE-OOS"=plsr.mse, "MSE-IS"=plsr.mse.IS ),
#              "GLM" = c("R2-OOS"=glm.R2, "R2-IS"=glm.R2.IS, "MSE-OOS"=glm.mse, "MSE-IS"=glm.mse.IS ),
#              "RandomF"= c("R2-OOS"=randomF.R2, "R2-IS"=randomF.R2.IS, "MSE-OOS"=randomF.mse, "MSE-IS"=randomF.mse.IS ),
#              "GBRT" = c("R2-OOS"=boost.R2, "R2-IS"=boost.R2.IS, "MSE-OOS"=boost.mse, "MSE-IS"=boost.mse.IS ), 
#              "NN1" = c("R2-OOS"=NN1.R2, "R2-IS"=NN1.R2.IS, "MSE-OOS"=NN1.mse, "MSE-IS"=NN1.mse.IS ),
#              "NN2" = c("R2-OOS"=NN2.R2, "R2-IS"=NN2.R2.IS, "MSE-OOS"=NN2.mse, "MSE-IS"=NN2.mse.IS ),
#              "NN3" = c("R2-OOS"=NN3.R2, "R2-IS"=NN3.R2.IS, "MSE-OOS"=NN3.mse, "MSE-IS"=NN3.mse.IS ),
#  ))
  
#}

output = expanding_window()

  
  
  
  
### ------------------------------------------------------------------
#Plots of model complexity over time
##--------------------------------------------------------------------

sum(var.select.boost$n == 0)











### ------------------------------------------------------------------
#Estimating Model's and giving a nice plot from the Book Solutions.
##--------------------------------------------------------------------
mod.ls <- lm(response~ . -1 ,data = df)
mod.ridge <- lm.ridge(response~ .,data = df)
mod.pcr <- pcr(formula=response ~ ., data=df, validation="CV") 
mod.plsr <- plsr(formula=response~ ., data=df, validation="CV") 
mod.lars <- lars(as.matrix(df[,2:ncol(df)]),
                 df[,1], type="lar")
mod.lasso <- lars(as.matrix(df[,2:ncol(df)]),
                  df[,1], type="lasso")

mods.coeffs <- data.frame(ls=mod.ls$coef , ridge=mod.ridge$coef ,
                          lasso=mod.lasso$beta[10,])

mods.coeffs$xs = row.names(mods.coeffs)
mods.coeffs <- as_tibble(mods.coeffs)

#mods.coeffs %>% gather(variable, value, -xs)

plot.data <- melt(mods.coeffs , id="xs")
ggplot(data=plot.data, aes(x=factor(xs), y=value, group=variable , colour=variable)) +
  geom_line() +
  geom_point() +
  xlab("Factor") +
  ylab("Regression Coefficient") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())











