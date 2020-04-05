x =arima.sim(n = 63, list(ar = c(0.8897, -0.4858), ma = c(-0.2279, 0.2488)),
             sd = sqrt(0.1796))
plot(x)
fit = lm(x ~ 1)

ACF = acf(x)
ACF$acf
VAR = var(x)

T = length(x)
L = round(0.75 * T^(1/3))

var=var(x)
lol <- embed(x,L+1)
cov=numeric(3)
for(h in 1:L){
  cov[h]=cov(lol[,4],lol[,4-h])
}
Sigmaxv=var+2*((1-1/L)*cov[1]+(1-2/L)*cov[2]+(1-3/L)*cov[3])
manualHAC=(1/(T))*Sigmaxv

(mine = (VAR + 2*( (2/3)*ACF$acf[2] + (1/3)*ACF$acf[3] ))/T)
(mine2 = (VAR + 2*( 0.5*ACF$acf[2] ))/T)
lrvar(x, type = "Newey-West", prewhite = FALSE) 
vcovHAC(fit)
manualHAC
