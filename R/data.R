library(reshape2)
library(tidyverse)
library(lubridate)
library(abind)
library(leaps)
library(glmnet)
library(pls)
library(lars)
library(mvtnorm)
library(elasticnet)
library(gglasso)
library(splines)
library(Hmisc)

rm(list=ls())

library(readxl)
## Using PredictorData2017 from Amit Goyals' website (University of lausanne)
## data is from paper ("http://www.hec.unil.ch/agoyal/")
## "A Comprehensive Look at The Empirical Performance of Equity Premium Prediction"
## Data is updated up to 2017
## 
monthly <- 
  read_excel("~/Dropbox/Project/PredictorData2017.xlsx" , sheet = "Monthly",
             na = "NaN")
quart <- 
  read_excel("~/Dropbox/Project/PredictorData2017.xlsx" , 
             sheet = "Quarterly",
             na = "NaN")
annual <- 
  read_excel("~/Dropbox/Project/PredictorData2017.xlsx" , sheet = "Annual",
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

# full availability from 1947-03-01 , i/k not available before.
plot(df$date, complete.cases(df),type = "l") 

# Saw off dataframe to omit all NA rows.
df = df[complete.cases(df),]
#So the sample runs from 1947-03-01 to 2017-06-01

## Now remove superflous variables so not to have multicoll
final = select(df, -c(date, Index, yyyymm, D12,
                     E12, AAA, BAA, corpr, CRSP_SPvwx,
                     CRSP_SPvw, Rfree))

final = select(final, r, everything())
names(final)[names(final)=="b/m"] = "bm"

pairs(select(final, r, bm, lty,ntis, dfy, bm), col = "slateblue3", cex = sample(6)/15)


## Do sample splitting

T= floor(nrow(final))
train_T=floor(T/2)
val_T=floor(train_T+T/4)

train = final %>% slice(1:train_T ) %>% as.data.frame()
val = final %>% slice((train_T+1):val_T) %>% as.data.frame()
test = final %>% slice((val_T+1):(val_T+T)) %>% as.data.frame()




### Test sample total sum of squares (using historical mean)
tot.test = mean((test[,1] - mean(rbind(train,val)[,1]))^2)

### --------------------------------------------------------------------
# OLS
### --------------------------------------------------------------------
fit.OLS = lm(r ~ . , data=rbind(train,val))

summary(fit.OLS)
coef.OLS = coef(fit.OLS)
pred.OLS.test = predict(fit.OLS, test[,-1])
mse.lm = mean((pred.OLS.test - test[,1])^2)
(R2.lm =( 1 - mse.lm/(tot.test) )*100 )




# Shrinking / selection parameter framework
lambda.grid = 10^seq(5,-5, length = 100)
### ------------------------------------------------------------------
#Ridge regression
##--------------------------------------------------------------------
ridge = glmnet(as.matrix(train[,-1]), as.matrix(train[,1]),
               family = "gaussian", alpha=0 , lambda = lambda.grid)


plot(ridge)

pre.ridge = predict(ridge , newx = as.matrix(val[,-1]) ) 
mse.ridge = numeric(length(pre.ridge[1,]))
#Calculate mse for each lambda
for (i  in 1:length(pre.ridge[1,])){
  mse.ridge[i] = mean((pre.ridge[,i] - val[,1])^2)
}
plot(ridge$lambda, mse.ridge)
plot(mse.ridge)
min_lambda = which.min(mse.ridge) %>% ridge$lambda[.]

coef.ridge = coef(ridge)

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
pcr.fit = pcr(r ~ . , data = train ,  scale = TRUE, 
              validation = "none")
mse.pcr = numeric(pcr.fit$ncomp)
for (i in 1:pcr.fit$ncomp){
  pcr.pred = predict(pcr.fit, val[,-1], ncomp = i)
  mse.pcr[i] = mean((pcr.pred - val[,1])^2)
}
plot(1:pcr.fit$ncomp ,mse.pcr)
(pcr.comp = which.min(mse.pcr))
predict.pcr = predict(pcr.fit , test[,-1], ncomp = pcr.comp)
mse.pcr.test=mean((predict.pcr - test[,1])^2) 
(R2.pcr = (1 - mse.pcr.test/tot.test)*100)


plsr.fit = plsr(r ~ . , data = train, scale = TRUE,
                validation = "none")
mse.plsr = numeric(plsr.fit$ncomp)
for (i in 1:plsr.fit$ncomp){
  plsr.pred = predict(plsr.fit, val[,-1], ncomp = i)
  mse.plsr[i] = mean((plsr.pred - val[,1])^2)
}
plot(1:plsr.fit$ncomp ,mse.plsr)
(plsr.comp = which.min(mse.plsr))
plsr.pred = predict(plsr.fit , test[,-1], ncomp = plsr.comp)
mse.plsr.test = mean((plsr.pred - test[,1])^2) 
(R2.plsr = (1 - mse.plsr.test/tot.test)*100)


### ------------------------------------------------------------------
# Generalized Linear Model
##--------------------------------------------------------------------
#using gglasso pacakage
#Create a matrix of transformed predicters, a list with each element a matrix for given df for model
num.pred = length(train[1,-1])
lambda.grid.glm = 10^seq(5,-7, length = 100)

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

plot(glm.models[[1]])
coef(glm.models[[1]], s = lambda.grid.glm[4])

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

mse.glm.test = mean((pred.glm.test - test$r)^2)
(R2.glm.test = (1 - mse.glm.test/tot.test)*100   )

###############################################################
## Tree based models
################################################################
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
    formula         = r ~ ., 
    data            = train.rf, 
    num.trees       = 5000,
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
  formula         = r ~ ., 
  data            = train.rf, 
  num.trees       = 5000,
  mtry            = best.model.RF$mtry,
  min.node.size   = best.model.RF$node_size
)

pred.RF.test = predict(model.RF, test[,-1])
mse.RF.test = mean((pred.RF.test$predictions - test$r)^2)
(R2.RF = (1 -  mse.RF.test/tot.test)*100)


#### Boosting
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

pred.boost = predict(model.boosting[[opt.boost]], newdata = test[,-1], n.trees =500)
mse.boost = mean((pred.boost - test$r)^2)
R2.boost = ( 1 - mse.boost/tot.test )*100 



###############################################################
## Neural Network
################################################################
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

NN3 <- keras_model_sequential() %>%
  layer_dense(units = 32, activation = "relu", 
              kernel_regularizer = regularizer_l1(0.01),
              input_shape = dim(x_train)[2]) %>%
  layer_dense(units = 16, activation = "relu",
              kernel_regularizer = regularizer_l1(0.01)) %>%
  layer_dense(units = 8, activation = "relu",
              kernel_regularizer = regularizer_l1(0.01)) %>% 
  layer_dense(units = 1, kernel_regularizer = regularizer_l1(0.01))

NN2 <- keras_model_sequential() %>%
  layer_dense(units = 16, activation = "relu", 
              kernel_regularizer = regularizer_l1(0.01),
              input_shape = dim(x_train)[2]) %>%
  layer_dense(units = 8, activation = "relu",
              kernel_regularizer = regularizer_l1(0.01)) %>%
  layer_dense(units = 1, kernel_regularizer = regularizer_l1(0.01))

NN3 %>% compile(
  loss = "mean_squared_error",
  optimizer = optimizer_rmsprop(),
  metrics = list("mean_absolute_error")
)

NN3 %>% summary()

print_dot_callback <- callback_lambda(
  on_epoch_end = function(epoch, logs) {
    if (epoch %% 80 == 0) cat("\n")
    cat(".")
  }
)   

# The patience parameter is the amount of epochs to check for improvement.
early_stop <- callback_early_stopping(monitor = "val_loss", 
                                      patience = 2)

history <- NN3 %>% fit(
  x_train,
  y_train,
  epochs = 500,
  batch_size = 16,
  shuffle = FALSE,
  validation_data = list(x_val, y_val),
  verbose = 0,
  callbacks = list(early_stop, print_dot_callback)
)

plot(history, metrics = "mean_absolute_error", smooth = FALSE) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 0.4))

pred.NN3.test <- NN3 %>% predict(x_test)
NN3.mse = mean((pred.NN3.test - test[,1])^2)
(1 - NN3.mse/tot.test)*100






