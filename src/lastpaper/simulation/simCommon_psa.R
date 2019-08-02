fixed_psaFormula = ~ 1 + age + ns(I(year_visit-2)/2, knots=(c(0.5, 1.3, 3)-2)/2, Boundary.knots=(c(0, 6.3)-2)/2)
random_psaFormula = ~ 1 + ns(I(year_visit-2)/2, knots=(c(0.5, 1.3, 3)-2)/2, Boundary.knots=(c(0, 6.3)-2)/2)
fixed_random_psaSlopeFormula = ~ 0 + dns(I(year_visit-2)/2, knots=(c(0.5, 1.3, 3)-2)/2, Boundary.knots=(c(0, 6.3)-2)/2)

MAX_FOLLOW_UP = 10
PSA_CHECK_UP_TIME = c(seq(0, 2, 0.25), seq(2.5, MAX_FOLLOW_UP, 0.5))

generateTDistError = function(n){
  sigma_psa = mvJoint_psa_time_scaled$statistics$postMeans$sigma1
  JMbayes:::rgt(n = n, mu = 0, sigma = sigma_psa, df = 3)
}

generateNormDistError = function(n){
  sigma_psa = mvJoint_psa_time_scaled$statistics$postMeans$sigma1
  rnorm(n=n, mean = 0, sd = sigma_psa)
}

generateTruePSAProfile = function(age, year_visit, randomEff_psa){
  betas_psa = mvJoint_psa_time_scaled$statistics$postMeans$betas1
  
  df = data.frame(age, year_visit)
  model.matrix(fixed_psaFormula, df) %*% betas_psa + 
    model.matrix(random_psaFormula, df) %*% as.numeric(randomEff_psa)
}

generateTruePSASlope = function(year_visit, randomEff_psa_slope){
  betas_psa_time = mvJoint_psa_time_scaled$statistics$postMeans$betas1[3:6]
  
  df = data.frame(year_visit)
  model.matrix(fixed_random_psaSlopeFormula, df) %*% betas_psa_time + 
    model.matrix(fixed_random_psaSlopeFormula, df) %*% as.numeric(randomEff_psa_slope)
}

generateSimulationData = function(nSub, psaErrorDist = "t3", bNames){
  
  survivalFunc <- function (t, patientId) {
    exp(-integrate(hazardFunc, lower = 0, upper = t, patientId)$value)
  }
  
  invSurvival <- function (t, u, patientId) {
    log(u) + integrate(hazardFunc, lower = 0, upper = t, patientId)$value
  }
  
  hazardFunc = function (year_visit, patientId) {
    alphas = mvJoint_psa_time_scaled$statistics$postMeans$alphas
    
    b_subject = simDs.id[patientId, bNames]
    age = simDs.id$age[patientId]  
    wGamma = simDs.id$wGamma[patientId]
    
    baselinehazard = exp(splineDesign(mvJoint_psa_time_scaled$control$knots, year_visit, 
                                      ord = mvJoint_psa_time_scaled$control$ordSpline, outer.ok = T) %*% mvJoint_psa_time_scaled$statistics$postMeans$Bs_gammas)
    
    truePSA = generateTruePSAProfile(age, year_visit, b_subject)
    truePSASlope = generateTruePSASlope(year_visit, b_subject[-1])
    
    y_Alpha = cbind(truePSA, truePSASlope) %*% alphas
    
    baselinehazard * exp(wGamma + y_Alpha)
  }
  
  pSurvTime = function(survProb, patientId){
    Low = 1e-05
    Up <- MAX_FOLLOW_UP
    Root <- try(uniroot(invSurvival, interval = c(Low, Up), 
                        u = survProb, patientId)$root, TRUE)
    if(inherits(Root, "try-error")){
      return(NA)
    }else{
      return(Root)
    }
  }
  
  longTimes <- c(replicate(nSub, PSA_CHECK_UP_TIME, simplify = T)) 
  timesPerSubject = length(longTimes) / nSub
  
  ###########################
  # Generate the longitudinal dataset
  ###########################
  subId <- 1:nSub
  
  D = mvJoint_psa_time_scaled$statistics$postMeans$D
  b <- mvrnorm(nSub, mu = rep(0, nrow(D)), D)
  colnames(b) = bNames
  
  age <- rnorm(n = nSub, mean=65, sd=7)
  simDs.id = data.frame(P_ID = subId, age = age, progression_time = NA)
  
  gammas = mvJoint_psa_time_scaled$statistics$postMeans$gammas
  survivalFormula = ~ 0 + age
  simDs.id$wGamma <- as.numeric(model.matrix(survivalFormula, data = simDs.id) %*% gammas)
  
  simDs.id = cbind(simDs.id, b)
  rm(b)
  
  #Longitudinal dataset
  simDs = data.frame(P_ID = rep(subId, each=timesPerSubject), 
                     age = rep(age, each=timesPerSubject),
                     year_visit = longTimes,
                     visitNumber = rep(1:timesPerSubject, nSub),
                     log2psaplus1 = NA)
  
  simDs$trueLog2psaplus1Profile = unlist(by(data=simDs, simDs$P_ID, function(data){
    pid = data$P_ID[1]
    randomEff = simDs.id[pid,bNames]
    generateTruePSAProfile(age = data$age, year_visit = data$year_visit, randomEff_psa = randomEff)
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
  
  simDs.id$survProbs <- runif(nSub)
  
  simDs.id$progression_time = sapply(1:nSub, FUN = function(i){
    pSurvTime(simDs.id$survProbs[i], i)
  })
                                      
  print(simDs.id$progression_time)
  percentageCensored = sum(is.na(simDs.id$progression_time))/nSub
  print(paste0("Percent of event times studyended-censored = ", percentageCensored*100, "%"))
  plot(density(simDs.id$progression_time, na.rm=T))
  
  simDs.id$progressed = ifelse(is.na(simDs.id$progression_time), 0, 1)
  simDs.id$studyend_censored = simDs.id$progressed
  simDs.id$progression_time[simDs.id$progressed==0] = MAX_FOLLOW_UP
  
  return(list(simDs=simDs, simDs.id=simDs.id, timesPerSubject=timesPerSubject))
}

fitJointModelOnNewData = function(simDs, simDs.id, nSubTraining){
  
  nSub = nrow(simDs.id)
  timesPerSubject = nrow(simDs) / nSub
  
  #Divide into training and test
  trainingDs.id = simDs.id[1:nSubTraining, ]
  testDs.id = simDs.id[(nSubTraining+1):nSub,]
  trainingDs = simDs[simDs$P_ID %in% trainingDs.id$P_ID, ]
  testDs = simDs[simDs$P_ID %in% testDs.id$P_ID, ]
  
  trainingDs$progression_time = rep(trainingDs.id$progression_time, each=timesPerSubject)
  trainingDs$progressed = rep(trainingDs.id$progressed, each=timesPerSubject)
  testDs$progression_time = rep(testDs.id$progression_time, each=timesPerSubject)
  
  # drop the longitudinal measurementsthat were taken after the observed event time for each subject.
  trainingDs = trainingDs[trainingDs$year_visit <= trainingDs$progression_time, ]
  
  print("Starting to fit long model with mvglmer")
  
  startTime_mvglmer_simDs = Sys.time()
  mvglmer_psa_simDs = mvglmer(list(log2psaplus1 ~ age +
                                         ns(I((year_visit-2)/2), knots=(c(0.5, 1.3, 3)-2)/2, Boundary.knots=(c(0, 6.3)-2)/2) + 
                                         (ns(I((year_visit-2)/2), knots=(c(0.5, 1.3, 3)-2)/2, Boundary.knots=(c(0, 6.3)-2)/2)|P_ID)),
                                  data=trainingDs, families = list(gaussian))
  endTime_mvglmer_simDs = Sys.time()
  mvglmer_fitting_time = endTime_mvglmer_simDs - startTime_mvglmer_simDs
  
  print("Done fitting long model with mvglmer")
  print(mvglmer_fitting_time)
  
  forms_psa = list("log2psaplus1" = "value",
                   "log2psaplus1" = list(fixed = ~ 0 + dns(I((year_visit-2)/2), knots=(c(0.5, 1.3, 3)-2)/2, Boundary.knots=(c(0, 6.3)-2)/2),
                                         random=~0 + dns(I((year_visit-2)/2), knots=(c(0.5, 1.3, 3)-2)/2, Boundary.knots=(c(0, 6.3)-2)/2),
                                         indFixed = 3:6, indRandom=2:5, name = "slope"))
  
  survModel_simDs = coxph(Surv(progression_time, progressed) ~ age,
                          data = trainingDs.id, x=T, model=T)
  
  print("Starting to fit joint model with mvJointModelBayes")
  startTime_mvJoint_simDs = Sys.time()
  mvJoint_psa_simDs = mvJointModelBayes(mvglmer_psa_simDs, survModel_simDs,
                                            timeVar = "year_visit", Formulas = forms_psa)
                                            
  endTime_mvJoint_simDs = Sys.time()
  mvjoint_fitting_time = endTime_mvJoint_simDs - startTime_mvJoint_simDs
  
  print("Done fitting joint model with mvJointModelBayes")
  print(mvjoint_fitting_time)
  
  out = list("trainingData"=list(trainingDs=trainingDs, trainingDs.id=trainingDs.id),
             "testData"=list(testDs=testDs, testDs.id=testDs.id),
             "survModel_simDs"=survModel_simDs, "mvglmer_psa_simDs"=mvglmer_psa_simDs, "mvglmer_fitting_time"=mvglmer_fitting_time,
             "mvJoint_psa_simDs"=mvJoint_psa_simDs)
  
  return(out)
}
