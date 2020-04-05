###TEsts not needed

###Elasticnet package
enet.model =enet(as.matrix(train[,-1]), as.matrix(train[,1]), lambda = 0)
attributes(enet.model)
enet.model$beta.pure
plot(enet.model)

elastic_net.ifelse <- function(v) {
  if (between(v[1],0,1) & v[2]>0){
    model = glmnet(as.matrix(train[,-1]), as.matrix(train[,1]), 
                   family = "gaussian", alpha=v[1], lambda = v[2])
    yhat = predict(model, newx=as.matrix(val[,-1]),  type = "link" )
    mse = mean((as.numeric(yhat) - val[,1])^2)
  } else {mse = Inf}
  
  return(mse)
}

#####################
#### THis part is not really important
####################

alpha.grid2 = seq(optimal_par$alpha - 0.05, optimal_par$alpha + 0.05,length = 11)
glmnet.grid2 = expand.grid(alpha = alpha.grid2,  lambda = lambda.grid)
elastic_overgrid = sapply(1:length(glmnet.grid2[,1]), 
                          function(i) elastic_net( as.matrix(glmnet.grid2[i,]) ) )
plot(1:length(glmnet.grid2[,1]), elastic_overgrid)
optimal_par = which.min(elastic_overgrid) %>% glmnet.grid2[., ]

glmnet.model = glmnet(as.matrix(train[,-1]), as.matrix(train[,1]),
                      family = "gaussian", 
                      alpha=optimal_par$alpha , lambda = optimal_par$lambda)
#coef(glmnet.model)
glmnet.predict = predict(glmnet.model, newx = as.matrix(test[,-1]))


#colmeans = numeric(length(al))
#for (i in seq(1, length(al))){
#  colmeans[i] = colMeans((what[[i]] - val[,1])^2) %>% 
#  which.min(.)
#}

#mse.colmeans = numeric(length(al))
#for (i in 1:length(al)){
#  fit = lol[[i]]$lambda[colmeans[i]]
#  mse.colmeans[i] = mean((predict(lol[[i]], s=fit, newx=as.matrix(val[,-1])) - val[,1])^2)
#}



elastic_net <- function(v) {
  model = glmnet(as.matrix(train[,-1]), as.matrix(train[,1]), 
                 family = "gaussian", alpha=v[2], lambda = v[1])
  yhat = predict(model, newx=as.matrix(val[,-1]))
  mse = mean((as.numeric(yhat) - val[,1])^2)
  return(mse)
}

elastic_overgrid = sapply(1:length(glmnet.grid[,1]), 
                          function(i) elastic_net( as.matrix(glmnet.grid[i,]) ) )
(optimal_par = which.min(elastic_overgrid) %>% glmnet.grid[., ])
plot(1:length(glmnet.grid[,1]), elastic_overgrid)

glmnet.model = glmnet(as.matrix(train[,-1]), as.matrix(train[,1]),
                      family = "gaussian", 
                      alpha=optimal_par$alpha , lambda = optimal_par$lambda)

enet.pred = predict(glmnet.model, newx = as.matrix(test[,-1]))

##One dimensional functio for glmnet
afternoon = function(l){
  if (l>0){
    x = glmnet(as.matrix(train[,-1]), as.matrix(train[,1]),
               family = "gaussian", alpha=0, lambda = l)
    pre = predict(x , newx = as.matrix(val[,-1]), type = "link" )
    mse = mean((as.numeric(pre) - val[,1])^2)
  } else {mse=Inf}
  return(mse)
}

##Testing the caret package ..
alpha.grid = seq(0,1,length = 11)
lambda.grid = 10^seq(2,-2, length = 100)
glmnet.grid = expand.grid(alpha = alpha.grid,  lambda = lambda.grid)

elastic_overgrid = sapply(1:1100, function(i) elastic_net(  as.matrix(grid[i,])   )   )
plot(1:1100, elastic_overgrid)
(optimal_par = which.min(elastic_overgrid) %>% grid[., ] )

glmnet.model = glmnet(as.matrix(train[,-1]), as.matrix(train[,1]),
                      family = "gaussian", 
                      alpha=optimal_par$alpha , lambda = optimal_par$lambda)
coef(glmnet.model)
glmnet.predict = predict(glmnet.model, newx = as.matrix(test[,-1]))


my.coef <- glmnet(as.matrix(rbind(train, val)[,-1]), as.matrix(rbind(train, val)[,1]),
                  family = "gaussian", alpha=my.train$bestTune$alpha,
                  lambda = my.train$bestTune$lambda) %>%
  coef(.)

##Something a bit different
al <- seq(0,1,length.out = 3)
gnet = function(i){
  glmnet(as.matrix(train[,-1]),
         as.matrix(train[,1]), alpha=i, lambda = lambda.grid)
}
x = lapply(al, gnet)
what = lapply(x  ,
              function(i) predict(i, newx=as.matrix(val[,-1])) %>% 
                colMeans((. - val[,1])^2) %>% 
                which.min(.))
mse = numeric(11)
for (i in 1:11){
  fit = x[[i]]$lambda[what[[i]]]
  mse[i] = mean((predict(x[[i]], s=fit, newx=as.matrix(val[,-1])) - val[,1])^2)
}
plot(al,mse)
##Yields the same result
##If I try with a finder grid for alpha I get a sligtly different result.
##Whether it has any relevance in terms of MSE I don't know.
al <- seq(0,1,length.out = 100)
x = lapply(al, function(i) glmnet(as.matrix(train[,-1]), as.matrix(train[,1]), alpha=i, lambda = lambda.grid))
what = lapply(x  ,
              function(i) predict(i, newx=as.matrix(val[,-1])) %>% colMeans((. - val[,1])^2) %>% 
                which.min(.))
mse = numeric(100)
for (i in 1:100){
  fit = x[[i]]$lambda[what[[i]]]
  mse[i] = mean((predict(x[[i]], s=fit, newx=as.matrix(val[,-1])) - val[,1])^2)
}



#with caret package, does not give me the same.
trnCtrl <- trainControl(method = "timeslice", 
                        initialWindow = 18000,  horizon = 9000, skip = 9000)
control <- trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 10,
  ## repeated ten times
  repeats = 2)
set.seed(2019)
my.train = train(Y ~. - id - time, data=df[1:27000,] , method = "glmnet",
                 tuneGrid = glmnet.grid , trControl = trnCtrl , standardize = TRUE, 
                 maxit = 100000)
plot(my.train)
attributes(my.train)
my.train$bestTune
my.glmnet <- my.train$finalModel
coef.glmnet.caret <-  coef(my.glmnet, s= my.train$bestTune$lambda) # Selects only the first two predicters
