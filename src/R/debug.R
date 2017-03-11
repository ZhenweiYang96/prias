surv = sapply(1:10, function(t){survivalFunc(t, 1)}, simplify = T)
qplot(x=1:10, y=surv, geom="line")
