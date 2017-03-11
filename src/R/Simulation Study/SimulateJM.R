library(MASS)
library(splines)

nSub <- 1000 # number of subjects

fittedJointModel = mvJoint_psa_spline_pt1pt54_pt11_tdboth

generateLongtiudinalTimeBySchedule = function(){
  months = c(0, 3, 6, 9, 12, 15, 18, 21, 24, 30, 36, 42, 48, 54, 60, 66, 72, 78, 84, 
             90, 96, 102, 108, 114, 120)
  return(months/12)
}

getBetas = function(jointModel, weightedOnes = T){
  if(weightedOnes){
    jointModel$statistics$postwMeans$betas1  
  }else{
    jointModel$statistics$postMeans$betas1 
  }
}

getSigma = function(jointModel, weightedOnes = T){
  if(weightedOnes){
    jointModel$statistics$postwMeans$sigma1
  }else{
    jointModel$statistics$postMeans$sigma1
  }
}

getD = function(jointModel, weightedOnes = T){
  if(weightedOnes){
    jointModel$statistics$postwMeans$D
  }else{
    jointModel$statistics$postMeans$D
  }
}

getGamma = function(jointModel, weightedOnes = T){
  if(weightedOnes){
    jointModel$statistics$postwMeans$gammas
  }else{
    jointModel$statistics$postMeans$gammas
  }
}

getAlpha = function(jointModel, weightedOnes = T){
  if(weightedOnes){
    jointModel$statistics$postwMeans$alphas
  }else{
    jointModel$statistics$postMeans$alphas
  }
}

hazardFunc = function (s, i) {
  pdf_s = dweibull(s, shape = weibullShape, scale = weibullScale)
  survival_s = (1-pweibull(q = s,shape= weibullShape, scale = weibullScale))
  baselinehazard_s = pdf_s/survival_s
  
  df_s = data.frame(visitTimeYears = s, Age = Age[i])
  
  xi_s_val = model.matrix(fixedValueFormula, df_s)
  xi_s_slope = model.matrix(fixedSlopeFormula, df_s)
  
  zi_s_val = model.matrix(randomValueFormula, df_s)
  zi_s_slope = model.matrix(randomSlopeFormula, df_s)
  
  zib_val = zi_s_val %*% b[i, ]
  zib_slope = zi_s_slope %*% b[i, -1] #-1 to ignore intercept
  
  xBetaZb_s_value = xi_s_val %*% betas + zib_val
  xBetaZb_s_slope = xi_s_slope %*% betas[-c(1:3)] + zib_slope #-c(1:3) to ignore intercept, age and age^2
  
  y_Alpha = cbind(xBetaZb_s_value, xBetaZb_s_slope) %*% getAlpha(fittedJointModel)
  
  baselinehazard_s * exp(wGamma[i] + y_Alpha)
}

survivalFunc <- function (t, i) {
  exp(-integrate(hazardFunc, lower = 0, upper = t, i)$value)
}

invSurvival <- function (t, u, i) {
  integrate(hazardFunc, lower = 0, upper = t, i)$value + log(u)
}

###############################################
# parameters for the linear mixed effects model
###############################################
boundaryKnots <- c(0, 7)
fixedKnots <- c(0.1, 0.5, 4)
randomKnots <- fixedKnots[1]

fixedValueFormula = ~ 1 + I(Age - 70) + I((Age - 70)^2) + ns(visitTimeYears, knots = fixedKnots, Boundary.knots = boundaryKnots)
randomValueFormula = ~ 1 + ns(visitTimeYears, knots = randomKnots, Boundary.knots = boundaryKnots)

fixedSlopeFormula = ~ 0 + dns(visitTimeYears, knots = fixedKnots, Boundary.knots = boundaryKnots)
randomSlopeFormula = ~ 0 + dns(visitTimeYears, knots = randomKnots, Boundary.knots = boundaryKnots)

###############################################
# Design matrices for the longitudinal measurement model
###############################################
#For the age variable check the simulated and observed data distribution
longTimes <- c(replicate(nSub, generateLongtiudinalTimeBySchedule(), simplify = T)) # at which time points longitudinal measurements are supposed to be taken
timesPerSubject = length(longTimes) / nSub

qplot(x = c(prias.id$Age, rnorm(nrow(prias.id), mean=mean(prias.id$Age), sqrt(var(prias.id$Age)))), geom="density", color=c(rep("Obs", nrow(prias.id)), rep("Sim", nrow(prias.id))))
Age <- rnorm(n = nSub, mean=mean(prias.id$Age), sd=sqrt(var(prias.id$Age)))

subId <- rep(1:nSub)

simDs.id = data.frame(P_ID = subId, Age = Age)
simDs = data.frame(P_ID = rep(subId, each=timesPerSubject), 
                   Age = rep(Age, each=timesPerSubject),
                   visitTimeYears = longTimes)

X = model.matrix(fixedValueFormula, data = simDs)
Z = model.matrix(randomValueFormula, data = simDs)

betas <- getBetas(fittedJointModel)
sigma.y <- getSigma(fittedJointModel)
D <- getD(fittedJointModel)
b <- mvrnorm(nSub, mu = rep(0, nrow(D)), D)

Zb = unlist(lapply(1:nSub,function(i){
  Z[((i-1)*timesPerSubject + 1):(i*timesPerSubject),] %*% b[i, ]
}))

xBetaZb = X %*% betas + Zb
simDs$logpsa1 = rnorm(nSub * timesPerSubject, xBetaZb, sigma.y)

ggplot(data=simDs, aes(x=visitTimeYears, y=logpsa1)) + 
  geom_line(aes(group=P_ID)) + stat_smooth() + 
  ticksY(from=-2.5, max = 7.5, by = 1)

# design matrix for the survival model
survivalFormula = ~ 0 + I(Age - 70) + I((Age - 70)^2)
W <- model.matrix(survivalFormula, data = data.frame(Age))
gammas = getGamma(fittedJointModel)
wGamma <- as.vector(W %*% gammas)

weibullShape = 1
weibullScale = 7

u <- runif(nSub)
simDs.id$progression_time <- NA
for (i in 1:nSub) {
    Up <- 10
    tries <- 5
    Root <- try(uniroot(invSurvival, interval = c(1e-05, Up), u = u[i], i = i)$root, TRUE)
    while(inherits(Root, "try-error") && tries > 0) {
        tries <- tries - 1
        Up <- Up + 200
        Root <- try(uniroot(invSurvival, interval = c(1e-05, Up), u = u[i], i = i)$root, TRUE)
    }
    simDs.id$progression_time[i] <- if (!inherits(Root, "try-error")) Root else NA
}

pid_to_keep = simDs.id[!is.na(simDs.id$progression_time),]$P_ID
simDs = simDs[simDs$P_ID %in% pid_to_keep,]
simDs.id = simDs.id[simDs.id$P_ID %in% pid_to_keep,]
simDs$P_ID = droplevels(simDs$P_ID)
simDs.id$P_ID = droplevels(simDs.id$P_ID)

#Divide into training and test
trainingSize = round(nrow(simDs.id)/2)
trainingDs.id = simDs.id[1:trainingSize, ]
testDs.id = simDs.id[(trainingSize+1):nSub,]
trainingDs.id$P_ID = droplevels(trainingDs.id$P_ID)
testDs.id$P_ID = droplevels(testDs.id$P_ID)

trainingDs = simDs[simDs$P_ID %in% trainingDs.id$P_ID, ]
trainingDs$P_ID = droplevels(trainingDs$P_ID)
testDs = simDs[simDs$P_ID %in% testDs.id$P_ID, ]
testDs$P_ID = droplevels(testDs$P_ID)

# simulate censoring times from an exponential distribution for TRAINING DATA SET ONLY
# and calculate the observed event times, i.e., min(true event times, censoring times)
Ctimes <- rexp(trainingSize, 1/mean(prias.id[prias.id$progressed==0,]$progression_time))
trainingDs.id$progressed = trainingDs.id$progrssion_time > Ctimes
trainingDs.id$progression_time = min(trainingDs.id$progression_time, Ctimes)
trainingDs$progression_time = rep(trainingDs.id$progression_time, each=timesPerSubject)

# drop the longitudinal measurementsthat were taken after the observed event time for each subject.
trainingDs = trainingDs[trainingDs$visitTimeYears <= trainingDs$progression_time]

# delete all unused objects
rm(y, X, Z, id, n, na.ind, long.na.ind, ind, Ctimes, Time, event, W, 
    betas, sigma.y, gammas, alpha, eta.t, eta.y, phi, mean.Cens, t.max,
    trueTimes, u, Root, invSurvival, D, b, K, 
    times, group, i, tries, Up, boundaryKnots, fixedKnots, DF)
