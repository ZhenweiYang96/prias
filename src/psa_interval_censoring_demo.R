#1 generate random intercepts and slope
nObs = 10000

randEff = rnorm(nObs, 0, sqrt(10))

gender = ifelse(randEff>0, "m", "f")

rightCenstime = sapply(randEff, function(x){
    if(x>0){
        runif(1, 15,20)
    }else{
        runif(1, 1, 10)
    }
})

R = sapply(1:nObs, function(index){
    if(randEff[index]>0){
        runif()
    }else{
        Inf    
    }
})


L = sapply(1:nObs, function(index){
    if(randEff[index]>0){
        runif(1, 4, 15)
    }else{
        runif(1, 0, 8)
    }
})
    

modelAft = survreg(Surv(eventintervalStart, eventintervalEnd, type="interval2") ~ randEff[,2])
exp(modelAft$coefficients[2] * -1/modelAft$scale)
    

finaltime = pmin(eventtime, censtime)
eventInd = ifelse(eventtime<censtime, 1, 0)

summary(survreg(Surv(finaltime, eventInd) ~ randEff[,2]))

summary(coxph(Surv(finaltime, eventInd) ~ randEff[,2]))
