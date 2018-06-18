#The simulations we will try are 
#1. expected failure time, choosing decreasing numbers
#2. expected failure time, within 6 months/1 year do biopsy
#3. annual schedule
#4. prias schedule
#5. 5, 10, 15, 20 risk
#6. F1 score
#7. Hybrid: F1 or expected
getNextSeed = function(lastSeed){
  lastSeed + 1
}

getTheoreticalHazard = function(timePoint, progression_speeds){
  
  getBaselineHazard = function(weibullScale, weibullShape, times){
    return((weibullShape/weibullScale)*(times/weibullScale)^(weibullShape-1))
  }
  
  theoreticalHazard = sapply(progression_speeds, function(progression_speed){
    getBaselineHazard(weibullScale = weibullScales[progression_speed], 
                      weibullShape = weibullShapes[progression_speed], times = timePoint)
  })
  
  unscaledWeights = sapply(progression_speeds, function(progression_speed){
    weibullScale = weibullScales[progression_speed]
    weibullShape = weibullShapes[progression_speed]
    return(exp(-(timePoint/weibullScale)^weibullShape))
  })
  weights = unscaledWeights/sum(unscaledWeights)
  
  return(sum(theoreticalHazard * weights))
}

generatePSATimeBySchedule = function(){
  months = c(seq(0, 24, 3), seq(30, MAX_FAIL_TIME*12, 6))
  
  return(months/12)
}

generateDRETimeBySchedule = function(){
  months = c(seq(0, MAX_FAIL_TIME*12,6))
  return(months/12)
}

generateTDistError = function(n){
  sigma_psa = mvJoint_dre_psa_dre_value_superlight$statistics$postwMeans$sigma2
  JMbayes:::rgt(n = n, mu = 0, sigma = sigma_psa, df = 3)
}

generateNormDistError = function(n){
  sigma_psa = mvJoint_dre_psa_dre_value_superlight$statistics$postwMeans$sigma2
  rnorm(n=n, mean = 0, sd = sigma_psa)
}

generateTruePSAProfile = function(Age, visitTimeYears, randomEff_psa){
  betas_psa = mvJoint_dre_psa_dre_value_superlight$statistics$postwMeans$betas2
  
  fixedPSAFormula = ~ 1 +I(Age - 70) +  I((Age - 70)^2) + ns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42))
  randomPSAFormula = ~ 1 + ns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42))
  
  df = data.frame(Age, visitTimeYears)
  model.matrix(fixedPSAFormula, df) %*% betas_psa + model.matrix(randomPSAFormula, df) %*% as.numeric(randomEff_psa)
}

generateTruePSASlope = function(visitTimeYears, randomEff_psa_slope){
  betas_psa_time = mvJoint_dre_psa_dre_value_superlight$statistics$postwMeans$betas2[4:7]
  
  fixedPSASlopeFormula = ~ 0 + dns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42))
  randomPSASlopeFormula = ~ 0 + dns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42))
  
  df = data.frame(visitTimeYears)
  model.matrix(fixedPSASlopeFormula, df) %*% betas_psa_time + model.matrix(randomPSASlopeFormula, df) %*% as.numeric(randomEff_psa_slope)
}

generateTrueDRELogOdds = function(Age, visitTimeYears, randomEff_dre){
  betas_dre = mvJoint_dre_psa_dre_value_superlight$statistics$postwMeans$betas1
  
  fixedDREFormula = ~ 1 + I(Age - 70) +  I((Age - 70)^2) + visitTimeYears
  randomDREFormula = ~ 1 + visitTimeYears
  
  df = data.frame(Age, visitTimeYears)
  
  model.matrix(fixedDREFormula, df) %*% betas_dre + model.matrix(randomDREFormula, df) %*% as.numeric(randomEff_dre)
}

generateHighDRE = function(logOdds){
  rbinom(length(logOdds), 1, plogis(logOdds))
}

getBaselineHazard = function(weibullScale, weibullShape, times){
  return((weibullShape/weibullScale)*(times/weibullScale)^(weibullShape-1))
}

generateSimulationData = function(seed, nSub, psaErrorDist = "t3",
                           progression_speeds, bNames){
  
  survivalFunc <- function (t, patientId) {
    exp(-integrate(hazardFunc, lower = 0, upper = t, patientId)$value)
  }
  
  invSurvival <- function (t, u, patientId) {
    log(u) + integrate(hazardFunc, lower = 0, upper = t, patientId)$value
  }
  
  hazardFunc = function (visitTimeYears, patientId) {
    alphas = mvJoint_dre_psa_dre_value_superlight$statistics$postwMeans$alphas
    
    b_subject = simDs.id[patientId, bNames]
    Age = simDs.id$Age[patientId]  
    wGamma = simDs.id$wGamma[patientId]
    
    # progression_speed = simDs.id$progression_speed[patientId]
    # baselinehazard = getBaselineHazard(weibullScale = weibullScales[progression_speed], 
    #                                      weibullShape = weibullShapes[progression_speed], s)
    
    baselinehazard = sapply(visitTimeYears, getTheoreticalHazard, progression_speeds=progression_speeds)
    
    truePSA = generateTruePSAProfile(Age, visitTimeYears, b_subject[3:7])
    truePSASlope = generateTruePSASlope(visitTimeYears, b_subject[4:7])
    trueDRELogOdds = generateTrueDRELogOdds(Age, visitTimeYears, b_subject[1:2])
    
    y_Alpha = cbind(trueDRELogOdds, truePSA, truePSASlope) %*% alphas
    
    baselinehazard * exp(wGamma + y_Alpha)
  }
  
  pSurvTime = function(survProb, patientId){
    Low = 1e-05
    Up <- MAX_FAIL_TIME
    Root <- try(uniroot(invSurvival, interval = c(Low, Up), 
                        u = survProb, patientId)$root, TRUE)
    if(inherits(Root, "try-error")){
      print(Root)
      return(NA)
    }else{
      return(Root)
    }
  }
  
  set.seed(seed)
  longTimes <- c(replicate(nSub, generatePSATimeBySchedule(), simplify = T)) 
  timesPerSubject = length(longTimes) / nSub
  
  ###########################
  # Generate the longitudinal dataset
  ###########################
  subId <- 1:nSub
  
  D = mvJoint_dre_psa_dre_value_superlight$statistics$postwMeans$D
  b <- mvrnorm(nSub, mu = rep(0, nrow(D)), D)
  colnames(b) = bNames
  
  Age <- rnorm(n = nSub, mean=70, sd=7)
  progression_speed = sample(progression_speeds, nSub, replace = T)
  
  #ID dataset
  simDs.id = data.frame(P_ID = subId, Age = Age, progression_speed=progression_speed, progression_time = NA)
  simDs.id$progression_speed = as.character(progression_speed)
  
  gammas = mvJoint_dre_psa_dre_value_superlight$statistics$postwMeans$gammas
  survivalFormula = ~ 0 + I(Age - 70) + I((Age - 70)^2)
  simDs.id$wGamma <- as.numeric(model.matrix(survivalFormula, data = simDs.id) %*% gammas)
  
  simDs.id = cbind(simDs.id, b)
  rm(b)
  
  #Longitudinal dataset
  simDs = data.frame(P_ID = rep(subId, each=timesPerSubject), 
                     Age = rep(Age, each=timesPerSubject),
                     progression_speed = rep(progression_speed, each=timesPerSubject),
                     visitTimeYears = longTimes,
                     visitNumber = rep(1:timesPerSubject, nSub),
                     log2psaplus1 = NA, high_dre=NA)
  simDs$progression_speed = as.character(simDs$progression_speed)
  
  simDs$trueLogOddsHighDRE = unlist(by(data=simDs, simDs$P_ID, function(data){
    pid = data$P_ID[1]
    randomEff = simDs.id[pid,bNames[1:2]]
    generateTrueDRELogOdds(Age = data$Age, visitTimeYears = data$visitTimeYears, randomEff_dre = randomEff)
  }))
  
  simDs$high_dre = generateHighDRE(logOdds = simDs$trueLogOddsHighDRE)
  simDs$high_dre[!(simDs$visitTimeYears %in% generateDRETimeBySchedule())] = NA
  
  simDs$trueLog2psaplus1Profile = unlist(by(data=simDs, simDs$P_ID, function(data){
    pid = data$P_ID[1]
    randomEff = simDs.id[pid,bNames[3:7]]
    generateTruePSAProfile(Age = data$Age, visitTimeYears = data$visitTimeYears, randomEff_psa = randomEff)
  }))
  
  if(psaErrorDist=="normal"){
    simDs$log2psaplus1 = simDs$trueLog2psaplus1Profile + generateNormDistError(nrow(simDs))
  }else{
    simDs$log2psaplus1 = simDs$trueLog2psaplus1Profile + generateTDistError(nrow(simDs))
  }
  
  ###################
  #Now generate event times
  ###################
  
  print("Done generating longitudinal data. Now trying to generate event times.")
  
  #set seed again
  set.seed(seed)
  simDs.id$survProbs <- runif(nSub)
  
  ct = makeCluster(max_cores)
  registerDoParallel(ct)
  simDs.id$progression_time = foreach(i=1:nSub,.combine='c', 
                                      .export=c("pSurvTime", "invSurvival", "hazardFunc", "MAX_FAIL_TIME","getBaselineHazard",
                                                "simDs.id", "weibullScales","weibullShapes", "getTheoreticalHazard",
                                                "mvJoint_dre_psa_dre_value_superlight",
                                                "generateTruePSASlope", "generateTrueDRELogOdds", "generateTruePSAProfile"),
                                      .packages = c("splines", "JMbayes")) %dopar%{
                                        pSurvTime(simDs.id$survProbs[i], i)
                                      }
  print(simDs.id$progression_time)
  percentageCensored = sum(is.na(simDs.id$progression_time))/nSub
  print(paste0("Percent of event times studyended-censored = ", percentageCensored*100, "%"))
  plot(density(simDs.id$progression_time, na.rm=T))
  stopCluster(ct)
  
  simDs.id$progressed = ifelse(is.na(simDs.id$progression_time), 0,1)
  simDs.id$studyend_censored = simDs.id$progressed
  simDs.id$progression_time[simDs.id$progressed==0] = MAX_FAIL_TIME
  
  # pid_to_keep = simDs.id[!is.na(simDs.id$progression_time),]$P_ID
  # pid_to_keep = pid_to_keep[1:(nSubTraining + nSubTest)]
  # 
  # simDs = simDs[simDs$P_ID %in% pid_to_keep,]
  # simDs.id = simDs.id[simDs.id$P_ID %in% pid_to_keep,]
  # b = b[pid_to_keep,]
  # 
  # simDs.id$P_ID = 1:nrow(simDs.id)
  # simDs$P_ID = rep(simDs.id$P_ID, each=timesPerSubject)
  
  return(list(simDs=simDs, simDs.id=simDs.id))
}

fitJointModelOnNewData = function(seed, simDs, simDs.id, nSubTraining, nSubTest, 
                                  mvglmer_iter=1000, 
                                  censStartTime=MAX_FAIL_TIME, censEndTime=MAX_FAIL_TIME,
                                  engine="STAN"){
  
  if(engine=="STAN" & mvglmer_iter > 10000){
    print("Too many iterations for STAN to handle. making it 1000")
    mvglmer_iter = 1000
  }
  
  if(engine=="JAGS" & mvglmer_iter < 28000){
    print("Too less iterations for JAGS. making it 28000")
    mvglmer_iter = 28000
  }
  
  #Divide into training and test
  trainingDs.id = simDs.id[1:nSubTraining, ]
  testDs.id = simDs.id[(nSubTraining+1):nrow(simDs.id),]
  trainingDs = simDs[simDs$P_ID %in% trainingDs.id$P_ID, ]
  testDs = simDs[simDs$P_ID %in% testDs.id$P_ID, ]
  
  #Dropout censoring for training patients
  Ctimes<-runif(nSubTraining, censStartTime, censEndTime)
  
  trainingDs.id$progressed = trainingDs.id$progression_time <= Ctimes
  trainingDs.id$progression_time = pmin(trainingDs.id$progression_time, Ctimes)
  trainingDs$progression_time = rep(trainingDs.id$progression_time, each=timesPerSubject)
  trainingDs$progressed = rep(trainingDs.id$progressed, each=timesPerSubject)
  testDs$progression_time = rep(testDs.id$progression_time, each=timesPerSubject)
  
  # drop the longitudinal measurementsthat were taken after the observed event time for each subject.
  trainingDs = trainingDs[trainingDs$visitTimeYears <= trainingDs$progression_time, ]
  
  print("Starting to fit long model with mvglmer")
  
  startTime_mvglmer_simDs = Sys.time()
  mvglmer_dre_psa_simDs = mvglmer(list(high_dre~I(Age - 70) +  I((Age - 70)^2) + visitTimeYears  +
                                         (visitTimeYears|P_ID),

                                       log2psaplus1 ~ I(Age - 70) +  I((Age - 70)^2) +
                                         ns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42)) +
                                         (ns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42))|P_ID)),

                                  data=trainingDs, families = list(binomial, gaussian), engine = engine,
                                  control = list(n.iter=mvglmer_iter, n.processors=max_cores))
  endTime_mvglmer_simDs = Sys.time()
  mvglmer_fitting_time = endTime_mvglmer_simDs - startTime_mvglmer_simDs
  
  print("Done fitting long model with mvglmer")
  print(mvglmer_fitting_time)
  
  forms_simDs = list("high_dre" = "value",
                     "log2psaplus1" = "value",
                     "log2psaplus1" = list(fixed = ~ 0 + dns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42)),
                                           random=~0 + dns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42)),
                                           indFixed = 4:7, indRandom=2:5, name = "slope"))
  
  survModel_simDs = coxph(Surv(progression_time, progressed) ~ I(Age - 70) +  I((Age - 70)^2),
                          data = trainingDs.id, x=T, model=T)
  
  print("Starting to fit joint model with mvJointModelBayes")
  startTime_mvJoint_simDs = Sys.time()
  mvJoint_dre_psa_simDs = mvJointModelBayes(mvglmer_dre_psa_simDs, survModel_simDs,
                                            timeVar = "visitTimeYears", Formulas = forms_simDs,
                                            control=list(n_cores=max_cores))
  endTime_mvJoint_simDs = Sys.time()
  mvjoint_fitting_time = endTime_mvJoint_simDs - startTime_mvJoint_simDs
  
  print("Done fitting joint model with mvJointModelBayes")
  print(mvjoint_fitting_time)
  
  out = list("trainingData"=list(trainingDs=trainingDs, trainingDs.id=trainingDs.id),
             "testData"=list(testDs=testDs, testDs.id=testDs.id),
             "censStartTime"=censStartTime, "censEndTime"=censEndTime, "mvglmer_iter"=mvglmer_iter,
             "seed"=seed, "weibullScales" = weibullScales, "weibullShapes"=weibullShapes,
             "survModel_simDs"=survModel_simDs, "progression_speeds"=progression_speeds,
             "mvglmer_dre_psa_simDs"=mvglmer_dre_psa_simDs, "mvglmer_fitting_time"=mvglmer_fitting_time,
             "mvJoint_dre_psa_simDs"=mvJoint_dre_psa_simDs, "mvjoint_fitting_time"=mvjoint_fitting_time)
  
  return(out)
}
#########################
#run younden, f1 calculation
#startTimeThresholds = Sys.time()
#mvJoint_dre_psa_simDs$mcmc$b = NULL
#thresholdsDt1= find_thresholds(mvJoint_dre_psa_simDs, trainingDs, idVar = "P_ID", Dt = 1, n_cores = 4)
#endTimeThresholds = Sys.time()
