library(parallel)
#library(snowfall)

set.seed(2019)

execution = function(DGP, Pc,  N=200, T=180){
  
  
  library(reshape2)
  library(tidyverse)
  library(abind)
  library(MonteCarlo)
  library(glmnet)
  library(pls)
  library(lars)
  library(mvtnorm)
  library(caret)
  library(elasticnet)
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
  

    
    #Create the three samples for training, tuning and testing
    train= df %>% filter(time %in% 1:initial_T) %>% select(-c(id, time))
    val= df %>% filter(time %in% (initial_T+1):(initial_T+val_T)) %>% select(-c(id, time))
    test= df %>% 
      filter( time > initial_T+val_T ) %>% 
      select(-c(id,time))
    
    ################################################################################
    ## ------ ESTIMATING MODELS FOR THE GIVEN TRAINING AND VAL SAMPLE --------
    ################################################################################
    tot.test = mean(test[,1]^2)
    
    lambda.grid = 10^seq(1,-4, length = 80)
    
    ### ------------------------------------------------------------------
    # Standard Linear Model
    ###--------------------------------------------------------------------
    fit = lm(Y ~ . , data=train)
    lm.pred = predict(fit, newdata = test[,-1])
    
    
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
    ridge.pred = predict(ridge, newx=as.matrix(test[,-1]), s = min.lambda.ridge)
    
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
    lasso.pred = predict(lasso, newx = as.matrix(test[,-1]), s = min.lambda.lasso)
    
    ### ------------------------------------------------------------------
    #Parameter search for alpha and lambda in elastic net.
    ##--------------------------------------------------------------------
    
    
    #elastic_overgrid = sapply(1:length(glmnet.grid[,1]), 
    #                          function(i) elastic_net( as.matrix(glmnet.grid[i,]) ) )
    #optimal_par = which.min(elastic_overgrid) %>% glmnet.grid[., ]
    
    #glmnet.model = glmnet(as.matrix(train[,-1]), as.matrix(train[,1]),
    #                      family = "gaussian", 
    #                      alpha=optimal_par$alpha , lambda = optimal_par$lambda)
    
    #enet.pred[test.window.obs] = predict(glmnet.model, newx = as.matrix(test[,-1]))
    
    ##### WITH VECOTOIZATION
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
      pcr.pred.val = predict(pcr.fit, val[,-1], ncomp = i)
      mse.pcr[i] = mean((pcr.pred.val - val[,1])^2)
    }
    pcr.comp = which.min(mse.pcr)
    pcr.pred = predict(pcr.fit, test[,-1], ncomp = pcr.comp)
    
    
    plsr.fit = plsr(Y ~ . , data = train, scale = TRUE,
                    validation = "none")
    mse.plsr.val = numeric(plsr.fit$ncomp)
    for (i in 1:plsr.fit$ncomp){
      plsr.pred.val = predict(plsr.fit, val[,-1], ncomp = i)
      mse.plsr.val[i] = mean((plsr.pred.val - val[,1])^2)
    }
    plsr.comp = which.min(mse.plsr.val)
    plsr.pred = predict(plsr.fit , test[,-1], ncomp = plsr.comp)
    
    
  #######
  # Calculation of MSE and R^2 for all model based on predictions from loop
  #######
  
  lm.mse = mean((lm.pred - test[,1])^2)
  lm.R2 = (1 - lm.mse/(tot.test))*100
  
  ridge.mse = mean((ridge.pred - test[,1])^2)
  ridge.R2 = (1 - ridge.mse/(tot.test) )*100
  
  lasso.mse = mean((lasso.pred - test[,1])^2)
  lasso.R2 = (1 - lasso.mse/(tot.test) )*100
  
  enet.mse = mean((enet.pred - test[,1])^2)
  enet.R2 = (1 - enet.mse/(tot.test) )*100
  
  pcr.mse = mean((pcr.pred - test[,1])^2)
  pcr.R2 = (1 - pcr.mse/(tot.test) )*100
  
  plsr.mse = mean((plsr.pred - test[,1])^2)
  plsr.R2 = (1 - plsr.mse/(tot.test) )*100
  
  
  ##################################################################
  ########### ------- END OF SIMULATION ---------------- ###########
  ##################################################################
  return(list("OLS"=lm.R2, "Ridge"= ridge.R2, "Lasso"= lasso.R2, "ENet" = enet.R2, 
              "PCR" = pcr.R2, "PLSR" = plsr.R2))
  
}

######################################
### Using replicate from base R
######################################
lapply(1, execution, DGP = "a",  Pc = 50)

no_cores = detectCores() - 1
clust = makeCluster(no_cores)
MC.a.50 = parLapply(clust, rep("b", 100), execution, Pc = 50 )
stopCluster(clust)

clust = makeCluster(no_cores)
MC.b.50 = parLapply(clust, rep("b", 30), execution , Pc = 100)
stopCluster(clust)

output.a.50 = sapply(MC.a.50, unlist) %>% t() #%>% as_tibble
output.b.50 = sapply(MC.b.50, unlist) %>% t() #%>% as_tibble

colMeans(output.b.50)
