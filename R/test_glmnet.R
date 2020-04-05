library(reshape2)
library(tidyverse)
library(abind)
library(mvtnorm)
library(leaps)
library(glmnet)
library(pls)
library(lars)
library(elasticnet)
library(gglasso)
library(splines)

rm(list = ls())
simulation <- function(DGP, Pc = 50, N=200, T=180){
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
  z=abind(c,c2,along= 2)
  
  #Finally create the returns from the characteristics
  #v=mvrnorm(n = T,mu= c(0,0,0),Sigma = diag(0.05^2,3,3))
  v=rmvnorm(n=T,mean = c(0,0,0),sigma = diag(0.05^2,3,3))
  eps=rt(n = T*N,df = 5)*0.05
  dim(eps)=c(N,T)
  
  #There are two cases for the function g(z), case (a) and (b).
  g=matrix(nrow = N,ncol = T)
  
  if(DGP=="a"){
    theta.a=c(rep(0.02,3))
    for (i in 1:N){
      for (t in 1:T){
        g[i,t]=z[i,1,t]*theta.a[1]+z[i,2,t]*theta.a[2]+z[i,(Pc+3),t]*theta.a[3]
      }
    }
  }
    
  if(DGP=="b"){
    theta=c(0.04,0.03,0.012)
    for (i in 1:N){
      for (t in 1:T){
        g[i,t]=(z[i,1,t]^2)*theta[1]+
          z[i,1,t]*z[i,2,t]*theta[2]+sign(z[i,(Pc+3),t])*theta[3]
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
#  df = df[ , c(53,52,1:51)]

train_T=T/2
val_T=train_T+T/4

train= df %>% filter(time %in% 1:train_T ) %>% select(-c(id, time))
val= df %>% filter(time %in% (train_T+1):val_T) %>% select(-c(id, time))
test= df %>% filter(time %in% (val_T+1):T) %>% select(-c(id,time))


return(list("train"=train, "val"=val,"test"=test, "df"=df))
###############################
##        SIMULATION OVER
###############################
}
sim = simulation(DGP = "b")
df = sim$df
train = sim$train
val = sim$val
test = sim$test

### Test sample total sum of squares (using historical mean)
(tot.test = mean((test[,1] - mean(rbind(train,val)[,1]))^2))

### --------------------------------------------------------------------
# OLS with backward selection
### --------------------------------------------------------------------
fit.OLS = regsubsets(Y ~ . , data=train, method = "forward", nvmax = 80)

val.mat = model.matrix(Y ~ . , data = val)

OLS.mse.val = numeric(length(fit.OLS))
for (i in seq(1,length(fit.OLS))){
  coefi = coef(fit.OLS, id = i)
  pred.OLS.val = val.mat[,names(coefi)]%*%coefi
  OLS.mse.val[i] = mean((pred.OLS.val - val[,1])^2) 
}
plot(OLS.mse.val, type = "l")
subset.OLS = which.min(OLS.mse.val)

## On the test set
test.mat = model.matrix(Y ~ . , data = test)

coef.OLS = coef(fit.OLS, id = subset.OLS)
pred.OLS.test = test.mat[,names(coef.OLS)]%*%coef.OLS
mse.lm = mean((pred.OLS.test - test[,1])^2)
(R2.lm =( 1 - mse.lm/(tot.test) )*100     )



reg.summary = summary(fit.OLS)
#Adjusted Rsquares on the training sample, plot.
par(mfrow=c(1,2))
plot(reg.summary$rss ,xlab="Number of Variables ",ylab="RSS",
       type="l")
plot(reg.summary$adjr2 ,xlab="Number of Variables ",
       ylab="Adjusted RSq",type="l")
red.dot = which.max(reg.summary$adjr2)
points(red.dot,reg.summary$adjr2[red.dot], col="red",cex=2,pch=20)



# Shrinking / selection parameter framework
lambda.grid = 10^seq(2,-5, length = 100)
### ------------------------------------------------------------------
#Ridge regression
##--------------------------------------------------------------------
ridge = glmnet(as.matrix(train[,-1]), as.matrix(train[,1]),
           family = "gaussian", alpha=0 , lambda = lambda.grid)
#ridge
#plot(x)
#coef(x)[,45]
pre.ridge = predict(ridge , newx = as.matrix(val[,-1]) ) # , type = "link"
mse.ridge = numeric(length(pre.ridge[1,]))
#Calculate mse for each lambda
for (i  in 1:length(pre.ridge[1,])){
  mse.ridge[i] = mean((pre.ridge[,i] - val[,1])^2)
}
plot(ridge$lambda, mse.ridge)
plot(mse.ridge)
min_lambda = which.min(mse.ridge) %>% ridge$lambda[.]

#Get test error for the ridge with the lambda just found
pre.ridge.test = predict(ridge, newx = as.matrix(test[,-1]), s = min_lambda)
mse.ridge.test = mean((pre.ridge.test - test[,1])^2)
(R2.ridge.test = 1 - mse.ridge.test/tot.test)

### ------------------------------------------------------------------
#Lasso regression
##--------------------------------------------------------------------
lasso = glmnet(as.matrix(train[,-1]), as.matrix(train[,1]),
               family = "gaussian", alpha=1, lambda = lambda.grid)
#lasso
#plot(lasso)
#coef(lasso)[,45]
pre.lasso = predict(lasso , newx = as.matrix(val[,-1]))
mse.lasso = numeric(length(pre.lasso[1,]))
for (i  in 1:length(pre.lasso[1,])){
  mse.lasso[i] = mean((pre.lasso[,i] - val[,1])^2)
}
plot(lasso$lambda, mse.lasso)
plot(mse.lasso)
min_lambda.lasso = which.min(mse.lasso) %>% lasso$lambda[.]

pre.lasso.test = predict(lasso, newx = as.matrix(test[,-1]), s = min_lambda.lasso)
mse.lasso.test = mean((pre.lasso.test - test[,1])^2)
(R2.lasso.test = 1 - mse.lasso.test/tot.test)

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
plot(al, mse.enet.vect)
which.model =which.min(mse.enet.vect)
enet.alpha = which.model %>% al[.]
enet.lambda = lol[[which.model]]$lambda[ what2[[which.model]] ] 
  
##calculate the mse on the test sample ... 
glmnet.model = glmnet(as.matrix(train[,-1]), as.matrix(train[,1]),
                  family = "gaussian", 
                    alpha=enet.alpha , lambda = enet.lambda)

enet.pred = predict(glmnet.model, newx = as.matrix(test[,-1]))


### ------------------------------------------------------------------
#Parameter search PLS and PCR. Predictor averaging techniques. 
##--------------------------------------------------------------------
pcr.fit = pcr(Y ~ . , data = train ,  scale = TRUE, 
              validation = "none")
mse.pcr = numeric(pcr.fit$ncomp)
for (i in 1:pcr.fit$ncomp){
  pcr.pred = predict(pcr.fit, val[,-1], ncomp = i)
  mse.pcr[i] = mean((pcr.pred - val[,1])^2)
}
plot(1:pcr.fit$ncomp ,mse.pcr)
(pcr.comp = which.min(mse.pcr))
predict.pcr = predict(pcr.fit , test[,-1], ncomp = pcr.comp)
(mse.pcr.test=mean((predict.pcr - test[,1])^2) )
 
plsr.fit = plsr(Y ~ . , data = train, scale = TRUE,
               validation = "none")
mse.plsr = numeric(plsr.fit$ncomp)
for (i in 1:plsr.fit$ncomp){
  plsr.pred = predict(plsr.fit, val[,-1], ncomp = i)
  mse.plsr[i] = mean((plsr.pred - val[,1])^2)
}
plot(1:plsr.fit$ncomp ,mse.plsr)
(plsr.comp = which.min(mse.plsr))
plsr.pred = predict(plsr.fit , test[,-1], ncomp = plsr.comp)
(mse.plsr.test = mean((plsr.pred - test[,1])^2) )


### ------------------------------------------------------------------
# Generalized Linear Model
##--------------------------------------------------------------------
#using gglasso pacakage
#Create a matrix of transformed predicters, a list with each element a matrix for given df for model
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

groups = lapply(2:10, function(i) rep(1:num.pred, each = i) )

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

## Evaluate on validation set
glm.mse.val = numeric(length(glm.models))
for (i in 1:length(glm.models)){
  glm.pred.val = predict(glm.models[[i]],
                         newx=predictor.matrix.val[[i]], type="link",
                         s = glm.models[[i]]$lambda[which.lambda.glm[[i]]])
  glm.mse.val[i] = mean((glm.pred.val - val[,1])^2)
}
plot(glm.mse.val)
best.model = which.min(glm.mse.val)

pred.glm.test = predict(glm.models[[best.model]], newx=predictor.matrix.test[[best.model]], type="link",
          s = glm.models[[best.model]]$lambda[which.lambda.glm[[best.model]]] )





### ------------------------------------------------------------------
#Tree based Models (Using tree package)
##--------------------------------------------------------------------


library(reshape2)
library(tidyverse)
library(abind)
library(mvtnorm)
rm(list = ls())
simulation <- function(DGP, Pc = 50, N=200, T=180){
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
  z=abind(c,c2,along= 2)
  
  #Finally create the returns from the characteristics
  #v=mvrnorm(n = T,mu= c(0,0,0),Sigma = diag(0.05^2,3,3))
  v=rmvnorm(n=T,mean = c(0,0,0),sigma = diag(0.05^2,3,3))
  eps=rt(n = T*N,df = 5)*0.05
  dim(eps)=c(N,T)
  
  #There are two cases for the function g(z), case (a) and (b).
  g=matrix(nrow = N,ncol = T)
  
  if(DGP=="a"){
    theta.a=c(rep(0.02,3))
    for (i in 1:N){
      for (t in 1:T){
        g[i,t]=z[i,1,t]*theta.a[1]+z[i,2,t]*theta.a[2]+z[i,(Pc+3),t]*theta.a[3]
      }
    }
  }
  
  if(DGP=="b"){
    theta=c(0.04,0.03,0.012)
    for (i in 1:N){
      for (t in 1:T){
        g[i,t]=(z[i,1,t]^2)*theta[1]+
          z[i,1,t]*z[i,2,t]*theta[2]+sign(z[i,(Pc+3),t])*theta[3]
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
  #  df = df[ , c(53,52,1:51)]
  
  train_T=T/2
  val_T=train_T+T/4
  
  train= df %>% filter(time %in% 1:train_T ) %>% select(-c(id, time))
  val= df %>% filter(time %in% (train_T+1):val_T) %>% select(-c(id, time))
  test= df %>% filter(time %in% (val_T+1):T) %>% select(-c(id,time))
  
  
  return(list("train"=train, "val"=val,"test"=test, "df"=df))
  ###############################
  ##        SIMULATION OVER
  ###############################
}
sim = simulation(DGP = "b")
df = sim$df
train = sim$train
val = sim$val
test = sim$test

### Test sample total sum of squares
(tot.test = mean((test[,1] - mean(test[,1]))^2))

library(ranger)
no.pred = dim(train)[2] - 1

RF_grid <- expand.grid(
  mtry       = c(floor(sqrt(no.pred)), floor(no.pred/3)),
  node_size  = c(10, 100, 1000),
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

RF_grid %>% 
  dplyr::arrange(OOB_RMSE) %>%
  head(10)

best.model.RF = which.min(RF_grid$OOB_RMSE) %>% RF_grid[.,]


model.RF = ranger(
  formula         = Y ~ ., 
  data            = train.rf, 
  num.trees       = 500,
  mtry            = best.model.RF$mtry,
  min.node.size   = best.model.RF$node_size
)

pred.RF.test = predict(model.RF, test[,-1])
mse.RF.test = mean((pred.RF.test$predictions - test$Y)^2)
(R2.RF = (1 -  mse.RF.test/tot.test)*100)


#### Boosting
library(gbm)

boost_grid <- expand.grid(
  depth       = seq(1,4),
  shrink  = c(0.01, 0.2)
)
boost_grid_list = apply(boost_grid, MARGIN = 1, as.list)


#boost = gbm(Y~. , distribution = "gaussian", data = train,
#            interaction.depth = 4,
#            n.trees = 5000)


model.boosting = lapply(boost_grid_list, function(x)
                        gbm(Y~. , distribution = "gaussian", data = train,
                            interaction.depth = x$depth,
                            shrinkage = x$shrink,
                            n.trees =500 ) )

boost.pred.val = lapply(model.boosting , 
                       function(x) predict(x, newdata=val[,-1], n.trees=500))
boost.pred.mse = sapply(boost.pred.val, function(x) mean(( x - val$Y)^2) )   

opt.boost = which.min(boost.pred.mse)

pred.boost = predict(model.boosting[[opt.boost]], newdata = test[,-1], n.trees =500)
mse.boost = mean((pred.boost - test$Y)^2)
R2.boost = ( 1 - mse.boost/tot.test )*100 

### ------------------------------------------------------------------
#Neural Network
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

R2.NN = numeric(length(lambda.grid.NN))
pred.NN.test = list(0)
NN.mse = numeric(length(lambda.grid.NN))
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

for (i in 1:length(lambda.grid.NN)){

model = build_model()
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

pred.NN.test[[i]] <- model %>% predict(x_test)
NN.mse[i] = mean((pred.NN.test[[i]] - test[,1])^2)
R2.NN[i] = (1 - NN.mse[i]/tot.test)*100
}

pred.NN = pred.NN.test[[which.min(NN.mse)]]












softplus <- function(x) log(1+exp(x))
relu <- function(x) sapply(x, function(z) max(0,z))
x <- seq(from=-5, to=5, by=0.1)

fits <- data.frame(x=x, softplus = softplus(x), relu = relu(x))
long <- melt(fits, id.vars="x")
ggplot(data=long, aes(x=x, y=value, group=variable, colour=variable))+
  geom_line(size=1) +
  ggtitle("ReLU & Softplus") +
  theme(plot.title = element_text(size = 14)) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size = 8))


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
                          lasso=mod.lasso$beta[10,], enet= glmnet.model.i.have)

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

grid=10^seq(10,-2,length=100)
ridge.mod=glmnet(input,response_a,alpha=0,lambda=grid)
coef(ridge.mod)[,50]











