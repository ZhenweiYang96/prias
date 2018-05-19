source("src/decision_analytic/load_lib.R")

#The simulations we will try are 
#1. expected failure time, choosing decreasing numbers
#2. expected failure time, within 6 months/1 year do biopsy
#3. annual schedule
#4. prias schedule
#5. 5, 10, 15, 20 risk
#6. F1 score
#7. Hybrid: F1 or expected
generatePSATimeBySchedule = function(){
  months = c(seq(0, 24, 3), seq(30, 10*12, 6))
  
  return(months/12)
}

survivalFunc <- function (t, patientId) {
  exp(-integrate(hazardFunc, lower = 0, upper = t, patientId)$value)
}

invSurvival <- function (t, u, patientId) {
  log(u) + integrate(hazardFunc, lower = 0, upper = t, patientId)$value
}

hazardFunc = function (s, patientId) {
  weibullShape = 13
  weibullScale = 7
  
  b_subject = b[patientId, ]
  Age = simDs.id$Age[patientId]  
  wGamma = simDs.id$wGamma[patientId]
  
  # pdf_s = dweibull(s, shape = weibullShape, scale = weibullScale)
  # survival_s = (1-pweibull(q = s, shape= weibullShape, scale = weibullScale))
  # baselinehazard_s = pdf_s/survival_s
  
  baselinehazard_s = (weibullShape/weibullScale)*(s/weibullScale)^(weibullShape-1)
  #baselinehazard_s = exp(splineDesign(mvJoint_dre_psa_dre_value$control$knots, s, 
  #                 ord = mvJoint_dre_psa_dre_value$control$ordSpline, 
  #                 outer.ok = T) %*% mvJoint_dre_psa_dre_value$statistics$postMeans$Bs_gammas)
  
  df_s = data.frame(visitTimeYears = s, Age = Age)
  
  xi_s_log2psa_val = model.matrix(fixedPSAFormula, df_s)
  xi_s_log2psa_slope = model.matrix(fixedPSASlopeFormula, df_s)
  
  zi_s_log2psa_val = model.matrix(randomPSAFormula, df_s)
  zi_s_log2psa_slope = model.matrix(randomPSASlopeFormula, df_s)
  
  zib_log2psa_val = zi_s_log2psa_val %*% b_subject
  zib_log2psa_slope = zi_s_log2psa_slope %*% b_subject[-1] #One less random effect to ignore intercept
  
  xBetaZb_s_log2psa_value = xi_s_log2psa_val %*% betas_psa + zib_log2psa_val
  xBetaZb_s_log2psa_slope = xi_s_log2psa_slope %*% betas_psa[-c(1:3)] + zib_log2psa_slope #-c(1:3) to ignore intercept, age and age^2
  
  y_Alpha = cbind(xBetaZb_s_log2psa_value, xBetaZb_s_log2psa_slope) %*% alphas
  #y_Alpha = cbind(xBetaZb_s_value) %*% alphas[1]
  
  baselinehazard_s * exp(wGamma + y_Alpha)
}

pSurvTime = function(survProb, patientId){
  Low = 1e-05
  Up <- 10
  Root <- try(uniroot(invSurvival, interval = c(Low, Up), 
    u = survProb, patientId)$root, TRUE)
  if(inherits(Root, "try-error")){
    return(NA)
  }else{
    return(Root)
  }
}

#####################################
#Code to generate the dataset begins here
#####################################
set.seed(2018)
nSub = 1500
nSubTraining = 750
nSubTest = 250
longTimes <- c(replicate(nSub, generatePSATimeBySchedule(), simplify = T)) 
timesPerSubject = length(longTimes) / nSub

mvJoint_psa_normaldist$mcmc$b = NULL
gammas = mvJoint_psa_normaldist$statistics$postMeans$gammas
alphas = mvJoint_psa_normaldist$statistics$postMeans$alphas
betas_psa = mvJoint_psa_normaldist$statistics$postMeans$betas1
sigma_psa = mvJoint_psa_normaldist$statistics$postMeans$sigma1
D = mvJoint_psa_normaldist$statistics$postMeans$D
b <- mvrnorm(nSub, mu = rep(0, nrow(D)), D)

###########################
# Generate the longitudinal dataset
###########################
Age <- rnorm(n = nSub, mean=mean(prias.id$Age), sd=sqrt(var(prias.id$Age)))
subId <- rep(1:nSub)

#ID dataset
simDs.id = data.frame(P_ID = subId, Age = Age, progression_time = NA)
survivalFormula = ~ 0 + I(Age - 70) + I((Age - 70)^2)
simDs.id$wGamma <- as.vector(model.matrix(survivalFormula, data = simDs.id) %*% gammas)

#Longitudinal dataset
simDs = data.frame(P_ID = rep(subId, each=timesPerSubject), 
 Age = rep(Age, each=timesPerSubject),
 visitTimeYears = longTimes,
 visitNumber = rep(1:timesPerSubject, nSub),
 log2psa = NA)

fixedPSAFormula = ~ 1 +I(Age - 70) +  I((Age - 70)^2) + ns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42))
randomPSAFormula = ~ 1 + ns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42))

fixedPSASlopeFormula = ~ 0 + dns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42))
randomPSASlopeFormula = ~ 0 + dns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42))

X_psa = model.matrix(fixedPSAFormula, data = simDs)
Z_psa = model.matrix(randomPSAFormula, data = simDs)

Zb_psa = unlist(lapply(1:nSub,function(i){
  Z_psa[((i-1)*timesPerSubject + 1):(i*timesPerSubject),] %*% b[i,]
}))

simDs$log2psa = rnorm(nSub * timesPerSubject, X_psa %*% betas_psa + Zb_psa, sigma_psa)
###################
#Now generate event times
###################
u <- runif(nSub)

ct = makeCluster(detectCores())
registerDoParallel(ct)
simDs.id$progression_time = foreach(i=1:nSub,.combine='c', 
               .export=c("pSurvTime", "invSurvival", "hazardFunc", 
                         "fixedPSAFormula", "randomPSAFormula", "fixedPSASlopeFormula",
                         "randomPSASlopeFormula",
                          "betas_psa", "alphas", "b", "simDs.id"),
               .packages = c("splines", "JMbayes")) %dopar%{
                 pSurvTime(u[i], i)
  }
percentageRejected = sum(is.na(simDs.id$progression_time))/nSub
print(paste("Percent reject", percentageRejected*100))
plot(density(simDs.id$progression_time, na.rm=T))
stopCluster(ct)

if(percentageRejected > 0.15){
  stop("Too many NA's sampled")
}

pid_to_keep = simDs.id[!is.na(simDs.id$progression_time),]$P_ID
pid_to_keep = pid_to_keep[1:(nSubTraining + nSubTest)]

simDs = simDs[simDs$P_ID %in% pid_to_keep,]
simDs.id = simDs.id[simDs.id$P_ID %in% pid_to_keep,]
b = b[pid_to_keep,]

simDs.id$P_ID = 1:nrow(simDs.id)
simDs$P_ID = rep(simDs.id$P_ID, each=timesPerSubject)

#Divide into training and test
trainingDs.id = simDs.id[1:nSubTraining, ]
testDs.id = simDs.id[(nSubTraining+1):nrow(simDs.id),]
trainingDs = simDs[simDs$P_ID %in% trainingDs.id$P_ID, ]
testDs = simDs[simDs$P_ID %in% testDs.id$P_ID, ]

#Ctimes <- rexp(nSubTraining, 1/mean(prias.id[prias.id$progressed==0,]$progression_time))
#Ctimes<-runif(nSubTraining, 0, 15)
Ctimes = rep(15, nSubTraining)

trainingDs.id$progressed = trainingDs.id$progression_time <= Ctimes
trainingDs.id$progression_time = pmin(trainingDs.id$progression_time, Ctimes)
trainingDs$progression_time = rep(trainingDs.id$progression_time, each=timesPerSubject)
trainingDs$progressed = rep(trainingDs.id$progressed, each=timesPerSubject)
testDs$progression_time = rep(testDs.id$progression_time, each=timesPerSubject)

# drop the longitudinal measurementsthat were taken after the observed event time for each subject.
trainingDs = trainingDs[trainingDs$visitTimeYears <= trainingDs$progression_time, ]

startTime_mvglmer_simDs = Sys.time()
mvglmer_psa_simDs = mvglmer(list(log2psa ~ I(Age - 70) +  I((Age - 70)^2) +
                                              ns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42)) + 
                                              (ns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42))|P_ID)),
                                   data=trainingDs, families = list(gaussian), engine="STAN")
endTime_mvglmer_simDs = Sys.time()

forms_simDs = list("log2psa" = "value",
                       "log2psa" = list(fixed = ~ 0 + dns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42)),
                                        random=~0 + dns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42)),
                                        indFixed = 4:7, indRandom=2:5, name = "slope"))

startTime_mvJoint_simDs = Sys.time()
survModel_simDs = coxph(Surv(progression_time, progressed) ~ I(Age - 70) +  I((Age - 70)^2),
                             data = trainingDs.id, x=T, model=T)
mvJoint_psa_simDs = mvJointModelBayes(mvglmer_psa_simDs, survModel_simDs, 
                                              timeVar = "visitTimeYears", Formulas = forms_simDs)
                                          #control = list(knots=seq(0, 10, 0.25), ObsTimes.knots=F, ordSpline=3))
endTime_mvJoint_simDs = Sys.time()

###########check baseline hazard
times = seq(0, 10, length.out = 500)

weibullShape = 13
weibullScale = 7
sim_hazard = log((weibullShape/weibullScale)*(times/weibullScale)^(weibullShape-1))

fit_hazard = splineDesign(mvJoint_psa_simDs$control$knots, times, 
                          ord = mvJoint_psa_simDs$control$ordSpline, outer.ok = T) %*% mvJoint_psa_simDs$statistics$postMeans$Bs_gammas
df = data.frame(hazard = c(sim_hazard, fit_hazard), type=rep(c("sim","fit"), each=500), time = rep(times,2))
ggplot(data=df[df$time>0.3 & df$time<9,]) + geom_line(aes(x=time, y=exp(hazard), color=type))

