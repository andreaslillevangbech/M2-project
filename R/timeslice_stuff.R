createTimeSlices(1:15, 5, 3)
createTimeSlices(1:15, 5, 3, skip = 2)
createTimeSlices(1:25, 5, 3, skip = 3)

x = rep(1:10, 5)
t = sapply(1:5, function(i) rep(i, 10) )

df = tibble(x = x, t = c(t))
df$outcome = rnorm(50)

createTimeSlices(df$x, 30, 20)


lol = trainControl(method = "timeslice", 
                   fixedWindow = FALSE, 
                   initialWindow = 20, 
                   skip = 9)
train()
