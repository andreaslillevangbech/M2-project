## Diebold Mariano simulations

library(tidyverse)
library(xtable)

load("/Users/alexanderbech/Dropbox/Project/diebold_sim.RData")

DB.a.50 = sapply(1:50, function(i) MC.a.50[[i]]$`Diebold Mariano test`) %>% 
  rowMeans() %>% matrix(.,nrow = 10)

DB.a.100 = sapply(1:50, function(i) MC.a.100[[i]]$`Diebold Mariano test`) %>% 
  rowMeans() %>% matrix(.,nrow = 10)

DB.b.50 = sapply(1:50, function(i) MC.b.50[[i]]$`Diebold Mariano test`) %>% 
  rowMeans() %>% matrix(.,nrow = 10)

DB.b.100 = sapply(1:50, function(i) MC.b.100[[i]]$`Diebold Mariano test`) %>% 
  rowMeans() %>% matrix(.,nrow = 10)

DB.grid = c("OLS" = 1, "Ridge" = 2,
            "Lasso" = 3, "ENET" = 4, "PCR"=5,
            "PLSR"=6, "GLM" = 7, "RandomF" = 8,
            "GBRT" = 9,
            "Oracle"=10)

colnames(DB.a.50) = names(DB.grid)
rownames(DB.a.50) = names(DB.grid)
colnames(DB.a.100) = names(DB.grid)
rownames(DB.a.100) = names(DB.grid)
colnames(DB.b.50) = names(DB.grid)
rownames(DB.b.50) = names(DB.grid)
colnames(DB.b.100) = names(DB.grid)
rownames(DB.b.100) = names(DB.grid)

print(xtable(DB.a.50, type = "latex"), file = "DB.a.50.tex")
print(xtable(DB.a.100, type = "latex"), file = "DB.a.100.tex")
print(xtable(DB.b.50, type = "latex"), file = "DB.b.50.tex")
print(xtable(DB.b.100, type = "latex"), file = "DB.b.100.tex")




#Simulation manipulatio

load("/Users/alexanderbech/Dropbox/Project/MCfinal.RData")

collect.MC.a.50 = sapply(1:13, 
                         function(j) sapply(1:50, function(i) MC.a.50[[i]][[j]][1:2]) %>% 
                           t()  %>% #ifelse(.>(-30), ., NA) %>%  
                           colMeans(na.rm = TRUE) )
collect.MC.b.50 = sapply(1:13, 
                         function(j) sapply(1:50, function(i) MC.b.50[[i]][[j]][1:2]) %>% 
                           t()  %>% #ifelse(.>(-30), ., NA) %>%  
                           colMeans(na.rm = TRUE) )
collect.MC.a.100 = sapply(1:13, 
                          function(j) sapply(1:50, function(i) MC.a.100[[i]][[j]][1:2]) %>% 
                            t()  %>% #ifelse(.>(-30), ., NA) %>%  
                            colMeans(na.rm = TRUE) )
collect.MC.b.100 = sapply(1:13, 
                          function(j) sapply(1:50, function(i) MC.b.100[[i]][[j]][1:2]) %>% 
                            t()  %>% #ifelse(.>(-30), ., NA) %>%  
                            colMeans(na.rm = TRUE) )

MC_table = data.frame(
  "$R^2(%)$" = names(MC.a.50[[1]][1:13]),
  "IS" = collect.MC.a.50[2,],
  "OOS" = collect.MC.a.50[1,],
  "IS" = collect.MC.a.100[2,],
  "OOS" = collect.MC.a.100[1,],
  "IS" = collect.MC.b.50[2,],
  "OOS" = collect.MC.b.50[1,],
  "IS" = collect.MC.b.100[2,],
  "OOS" = collect.MC.b.100[1,]
)

print(xtable(MC_table, type = "latex"), file = "simulation_table1.tex")


### For longer horizons #####
load("/Users/alexanderbech/Dropbox/Project/MC_quart2.RData")
load("/Users/alexanderbech/Dropbox/Project/MC_half.RData")
load("/Users/alexanderbech/Dropbox/Project/MC_annual.RData")

collect.MC.a.100.q = sapply(1:11, 
                          function(j) sapply(1:50, function(i) MC.a.100.q[[i]][[j]][1:2]) %>% 
                            t()  %>% #ifelse(.>(-30), ., NA) %>%  
                            colMeans(na.rm = TRUE) )
collect.MC.b.100.q = sapply(1:11, 
                          function(j) sapply(1:50, function(i) MC.b.100.q[[i]][[j]][1:2]) %>% 
                            t()  %>% #ifelse(.>(-30), ., NA) %>%  
                            colMeans(na.rm = TRUE) )

collect.MC.a.100.h = sapply(1:11, 
                          function(j) sapply(1:50, function(i) MC.a.100.h[[i]][[j]][1:2]) %>% 
                            t()  %>% #ifelse(.>(-30), ., NA) %>%  
                            colMeans(na.rm = TRUE) )
collect.MC.b.100.h = sapply(1:11, 
                          function(j) sapply(1:50, function(i) MC.b.100.h[[i]][[j]][1:2]) %>% 
                            t()  %>% #ifelse(.>(-30), ., NA) %>%  
                            colMeans(na.rm = TRUE) )

collect.MC.a.100.a = sapply(1:11, 
                          function(j) sapply(1:50, function(i) MC.a.100.a[[i]][[j]][1:2]) %>% 
                            t()  %>% #ifelse(.>(-30), ., NA) %>%  
                            colMeans(na.rm = TRUE) )
collect.MC.b.100.a = sapply(1:11, 
                          function(j) sapply(1:50, function(i) MC.b.100.a[[i]][[j]][1:2]) %>% 
                            t()  %>% #ifelse(.>(-30), ., NA) %>%  
                            colMeans(na.rm = TRUE) )


MC_table_horizon = data.frame(
  "$R^2(%)$" = names(MC.a.100.q[[1]][1:11]),
  "IS" = collect.MC.a.100.q[2,],
  "OOS" = collect.MC.a.100.q[1,],
  "IS" = collect.MC.a.100.h[2,],
  "OOS" = collect.MC.a.100.h[1,],
  "IS" = collect.MC.a.100.a[2,],
  "OOS" = collect.MC.a.100.a[1,],
  "IS" = collect.MC.b.100.q[2,],
  "OOS" = collect.MC.b.100.q[1,],
  "IS" = collect.MC.b.100.h[2,],
  "OOS" = collect.MC.b.100.h[1,],
  "IS" = collect.MC.b.100.a[2,],
  "OOS" = collect.MC.b.100.a[1,]
)

print(xtable(MC_table_horizon, type = "latex"), 
      file = "/Users/alexanderbech/Dropbox/Project/simulation_table_hor.tex")



#################################################################
############# Output from data ##########################
#################################################################
#test sample goes from 637 to 849

test.dates = seq(2000,2017)
load("/Users/alexanderbech/Dropbox/Project/data_output_2.RData")
lasso.compl = sapply(1:18, function(i) sum(output$models$lasso[[i]]$model!=0))
enet.compl = sapply(1:18, function(i) sum(output$models$enet[[i]]$model!=0))
pcr.compl = sapply(1:18, function(i) output$models$pcr[[i]]$ncomp)
plsr.compl = sapply(1:18, function(i) output$models$plsr[[i]]$ncomp)
RF.compl = sapply(1:18, function(i) output$models$RF[[i]]$model$min.node.size)
boost.compl = sapply(1:18, function(i) output$models$GBRT[[i]]$depth)

par(mfrow=c(3,2), mar=c(2,4,2,1))
plot(test.dates, lasso.compl,
     main="Lasso",
     ylab="df", xlab = "",
     type="l",
     col="blue")
grid()
plot(test.dates, enet.compl,
     main="ENet",
     ylab="df", xlab = "",
     type="l",
     col="blue")
grid()
plot(test.dates, pcr.compl,
           main="PCR",
           ylab="ncomp", xlab = "",
           type="l",
           col="blue")
grid()
plot(test.dates, plsr.compl,
           main="PLSR",
           ylab="ncomp", xlab = "",
           type="l",
           col="blue")
grid()
plot(test.dates, RF.compl,
     main="Random Forest",
     ylab="node size", xlab = "",
     type="l",
     col="blue")
grid()
plot(test.dates, boost.compl,
     main="GBRT",
     ylab="depth", xlab = "",
     type="l",
     col="blue")
grid()

Diebold = as.data.frame(output$Diebold)
DB.grid = c("OLS" = 1, "Ridge" = 2,
            "Lasso" = 3, "ENET" = 4, "PCR"=5,
            "PLSR"=6, "GLM" = 7, "RandomF" = 8,
            "GBRT" = 9)
colnames(Diebold) = names(DB.grid)
rownames(Diebold) = names(DB.grid)
dim(Diebold) #9x9
Diebold = Diebold[1:8, 2:9]


print(xtable(Diebold, type = "latex"), 
                 file = "/Users/alexanderbech/Dropbox/Project/Diebold.tex")

