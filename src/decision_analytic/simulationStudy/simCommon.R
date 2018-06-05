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

generatePSATimeBySchedule = function(){
  months = c(seq(0, 24, 3), seq(30, MAX_FAIL_TIME*12, 6))
  
  return(months/12)
}

generateDRETimeBySchedule = function(){
  months = c(seq(0, MAX_FAIL_TIME*12,6))
  return(months/12)
}

getBaselineHazard = function(weibullScale, weibullShape, times){
  return((weibullShape/weibullScale)*(times/weibullScale)^(weibullShape-1))
}

mvJoint_dre_psa_dre_value_light$mcmc$b = NULL
gammas = mvJoint_dre_psa_dre_value_light$statistics$postMeans$gammas
alphas = mvJoint_dre_psa_dre_value_light$statistics$postMeans$alphas
betas_dre = mvJoint_dre_psa_dre_value_light$statistics$postMeans$betas1
betas_psa = mvJoint_dre_psa_dre_value_light$statistics$postMeans$betas2
sigma_psa = mvJoint_dre_psa_dre_value_light$statistics$postMeans$sigma2
D = mvJoint_dre_psa_dre_value_light$statistics$postMeans$D

fixedPSAFormula = ~ 1 +I(Age - 70) +  I((Age - 70)^2) + ns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42))
randomPSAFormula = ~ 1 + ns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42))

fixedPSASlopeFormula = ~ 0 + dns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42))
randomPSASlopeFormula = ~ 0 + dns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42))

fixedDREFormula = ~ 1 + I(Age - 70) +  I((Age - 70)^2) + visitTimeYears
randomDREFormula = ~ 1 + visitTimeYears

survivalFormula = ~ 0 + I(Age - 70) + I((Age - 70)^2)

fitJointModelOnNewData = function(seed, nSub, nSubTraining, nSubTest, 
                                  mvglmer_iter=1000, 
                                  censStartTime=MAX_FAIL_TIME, censEndTime=MAX_FAIL_TIME,
                                  progression_type="Mixed", engine="STAN"){
  
  if(engine=="STAN" & mvglmer_iter > 10000){
    print("Too many iterations for STAN to handle. making it 1000")
    mvglmer_iter = 1000
  }
  
  if(engine=="JAGS" & mvglmer_iter < 28000){
    print("Too less iterations for JAGS. making it 28000")
    mvglmer_iter = 28000
  }
  
  weibullScales = c("Fast"=3, "Medium"=5, "Slow"=8)
  weibullShapes = c("Fast"=5, "Medium"=8, "Slow"=14)
  
  survivalFunc <- function (t, patientId) {
    exp(-integrate(hazardFunc, lower = 0, upper = t, patientId)$value)
  }
  
  invSurvival <- function (t, u, patientId) {
    log(u) + integrate(hazardFunc, lower = 0, upper = t, patientId)$value
  }
  
  hazardFunc = function (s, patientId) {
    
    b_subject = b[patientId, ]
    Age = simDs.id$Age[patientId]  
    wGamma = simDs.id$wGamma[patientId]
    progression_speed = simDs.id$progression_speed[patientId]
    
    baselinehazard_s = getBaselineHazard(weibullScale = weibullScales[progression_speed], 
                                         weibullShape = weibullShapes[progression_speed], s)
    
    df_s = data.frame(visitTimeYears = s, Age = Age)
    
    xi_s_logodds_high_dre_val = model.matrix(fixedDREFormula, df_s)
    xi_s_log2psaplus1_val = model.matrix(fixedPSAFormula, df_s)
    xi_s_log2psaplus1_slope = model.matrix(fixedPSASlopeFormula, df_s)
    
    zi_s_logodds_high_dre_val = model.matrix(randomDREFormula, df_s)
    zi_s_log2psaplus1_val = model.matrix(randomPSAFormula, df_s)
    zi_s_log2psaplus1_slope = model.matrix(randomPSASlopeFormula, df_s)
    
    #There are 7 random effects, 1,2 for DRE and 3,4,5,6,7 for PSA
    zib_logodds_high_dre_val = zi_s_logodds_high_dre_val %*% b_subject[1:2]
    zib_log2psaplus1_val = zi_s_log2psaplus1_val %*% b_subject[3:7]
    zib_log2psaplus1_slope = zi_s_log2psaplus1_slope %*% b_subject[4:7] #One less random effect to ignore intercept
    
    xBetaZb_s_logodds_high_dre_value = xi_s_logodds_high_dre_val  %*% betas_dre + zib_logodds_high_dre_val
    xBetaZb_s_log2psaplus1_value = xi_s_log2psaplus1_val %*% betas_psa + zib_log2psaplus1_val
    xBetaZb_s_log2psaplus1_slope = xi_s_log2psaplus1_slope %*% betas_psa[-c(1:3)] + zib_log2psaplus1_slope #-c(1:3) to ignore intercept, age and age^2
    
    y_Alpha = cbind(xBetaZb_s_logodds_high_dre_value, xBetaZb_s_log2psaplus1_value, xBetaZb_s_log2psaplus1_slope) %*% alphas
    #y_Alpha = cbind(xBetaZb_s_value) %*% alphas[1]
    
    baselinehazard_s * exp(wGamma + y_Alpha)
  }
  
  pSurvTime = function(survProb, patientId){
    Low = 1e-05
    Up <- MAX_FAIL_TIME
    Root <- try(uniroot(invSurvival, interval = c(Low, Up), 
                        u = survProb, patientId)$root, TRUE)
    if(inherits(Root, "try-error")){
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
  subId <- rep(1:nSub)
  b <- mvrnorm(nSub, mu = rep(0, nrow(D)), D)
  Age <- rnorm(n = nSub, mean=mean(prias.id$Age), sd=sd(prias.id$Age))
  progression_speed = switch(progression_type, 
                             "Fast"=rep("Fast", nSub), 
                             "Medium"=rep("Medium", nSub), 
                             "Slow"=rep("Slow", nSub), 
                             "Mixed"=sample(c("Slow", "Medium", "Fast"), nSub, replace = T))
  
  #ID dataset
  simDs.id = data.frame(P_ID = subId, Age = Age, progression_speed=progression_speed, progression_time = NA)
  simDs.id$progression_speed = as.character(progression_speed)
  simDs.id$wGamma <- as.vector(model.matrix(survivalFormula, data = simDs.id) %*% gammas)
  
  #Longitudinal dataset
  simDs = data.frame(P_ID = rep(subId, each=timesPerSubject), 
                     Age = rep(Age, each=timesPerSubject),
                     progression_speed = rep(progression_speed, each=timesPerSubject),
                     visitTimeYears = longTimes,
                     visitNumber = rep(1:timesPerSubject, nSub),
                     log2psaplus1 = NA, high_dre=NA)
  simDs$progression_speed = as.character(simDs$progression_speed)
  
  X_psa = model.matrix(fixedPSAFormula, data = simDs)
  Z_psa = model.matrix(randomPSAFormula, data = simDs)
  
  X_dre = model.matrix(fixedDREFormula, data = simDs)
  Z_dre = model.matrix(randomDREFormula, data = simDs)
  
  Zb_psa = unlist(lapply(1:nSub,function(i){
    Z_psa[((i-1)*timesPerSubject + 1):(i*timesPerSubject),] %*% b[i, 3:7]
  }))
  
  Zb_dre = unlist(lapply(1:nSub,function(i){
    Z_dre[((i-1)*timesPerSubject + 1):(i*timesPerSubject),] %*% b[i, 1:2]
  }))
  
  simDs$log2psaplus1 = JMbayes:::rgt(n = nSub * timesPerSubject, mu = X_psa %*% betas_psa + Zb_psa, sigma = sigma_psa, df = 3)
  simDs$high_dre = rbinom(nSub * timesPerSubject, 1, plogis(X_dre %*% betas_dre + Zb_dre))
  
  simDs$high_dre[!(simDs$visitTimeYears %in% generateDRETimeBySchedule())] = NA
  ###################
  #Now generate event times
  ###################
  
  print("Done generating longitudinal data. Now trying to generate event times.")
  
  u <- runif(nSub)
  
  ct = makeCluster(max_cores)
  registerDoParallel(ct)
  simDs.id$progression_time = foreach(i=1:nSub,.combine='c', 
                                      .export=c("pSurvTime", "invSurvival", "hazardFunc", "MAX_FAIL_TIME","getBaselineHazard",
                                                "fixedPSAFormula", "randomPSAFormula", "fixedPSASlopeFormula",
                                                "randomPSASlopeFormula", "fixedDREFormula", "randomDREFormula",
                                                "weibullScales","weibullShapes",
                                                "betas_dre", "betas_psa", "alphas", "b", "simDs.id"),
                                      .packages = c("splines", "JMbayes")) %dopar%{
                                        pSurvTime(u[i], i)
                                      }
  print(simDs.id$progression_time)
  percentageRejected = sum(is.na(simDs.id$progression_time))/nSub
  print(paste0("Percent of event times rejected = ", percentageRejected*100, "%"))
  plot(density(simDs.id$progression_time, na.rm=T))
  stopCluster(ct)
  
  if(percentageRejected > 0.1){
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
             "b"=b, "seed"=seed, "weibullScales" = weibullScales, "weibullShapes"=weibullShapes,
             "survModel_simDs"=survModel_simDs,
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
