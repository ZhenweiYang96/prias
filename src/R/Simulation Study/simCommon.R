plotBiopsy2DPlot = function(nbSummary, offsetSummary, nbMeasure = "Median", offsetMeasure = "Median"){
  methodNames = colnames(nbSummary)
  methodCategory = sapply(methodNames, substr, 1, 10)
  
  df = data.frame(nb = nbSummary[nbMeasure,], offset=offsetSummary[offsetMeasure,], 
                  method=colnames(nbSummary), methodCategory = methodCategory)
  p = ggplot(data=df, aes(x=nb, y=offset*12, label=method)) + 
    geom_label(aes(fill = methodCategory), colour = "white", fontface = "bold") + 
    xlab(paste(nbMeasure,"Number of biopsies")) + ylab(paste(offsetMeasure, "Offset (months)")) + 
    scale_y_continuous(breaks=seq(0, 100, 1)) + scale_x_continuous(breaks=seq(0, 100, 0.5))
  print(p)
}

plot2DBoxPlot = function(nbSummary, offsetSummary){
  methodNames = colnames(nbSummary)
  methodCategory = sapply(methodNames, substr, 1, 10)
  
  df = data.frame(nb_x1 = nbSummary["1st Qu.",], nb_x2 = nbSummary["3rd Qu.",],
                  offset_y1 = offsetSummary["1st Qu.",], offset_y2 = offsetSummary["3rd Qu.",],
                  method=colnames(nbSummary), methodCategory = methodCategory)
  
  p = ggplot(data=df) + 
    geom_rect(mapping=aes(xmin=nb_x1, xmax=nb_x2, ymin=offset_y1*12, ymax=offset_y2*12, fill=methodCategory), color="black", alpha=0.5) +
    xlab("Number of biopsies") + ylab("Offset (months)")
    
  print(p)
}

plotTrueSurvival = function(dsId, patientId){
  
  time = 1:15
  survProb = sapply(time, function(t){
    survivalFunc(t, dsId, patientId)
  })
  
  byYAxis = (max(survProb) - min(survProb))/10
  
  qplot(x=time, y = survProb, geom = "line", xlab = "Time (years)", ylab = "Probability") + 
    ticksX(from=0, max = 15, by=1) + ticksY(min(survProb), max(survProb), by=byYAxis) + 
    ggtitle(patientId)
}

plotTrueLongitudinal = function(dsId, patientId){
  simDs = simulatedDsList[[dsId]]$simDs
  ggplot(data=simDs[simDs$P_ID == patientId, ], aes(x=visitTimeYears, y=logpsa1)) + 
    geom_line() + geom_point(color="red") + ggtitle(patientId) + xlab("Time (years)") + 
    ylab("log(PSA + 1)")
}

invDynSurvival <- function (t, u, dsId, patientDs) {
  survProb = 1
  if(t>12){
    survFitObj = survfitJM(simulatedDsList[[dsId]]$models$simJointModel_replaced, patientDs, idVar="P_ID", survTimes = t)
    survProbHPDI = HPDinterval(as.mcmc(sapply(survFitObj$full.results, function(x){x[[1]][1]})))
    survProb = round(survProbHPDI[1],3)
  }else{
    survProb = survfitJM(simulatedDsList[[dsId]]$models$simJointModel_replaced, patientDs, idVar="P_ID", survTimes = t)$summaries[[1]][1, "Median"]
  }
  
  u - survProb
}

pDynSurvTime = function(survProb, dsId, patientDs){
  #Return the time at which the dynamic survival probability is say 90%
  
  if(is.na(survProb)){
    return(NA)
  }else if(survProb == 1.0){
    return(max(patientDs$visitTimeYears)) 
  }
  
  Low = max(patientDs$visitTimeYears) + 1e-05
  Up <- 10
  tries  = 0
  
  repeat{
    tries = tries + 1
    Root <- try(uniroot(invDynSurvival, interval = c(Low, Up), 
                        u = survProb, dsId = dsId, patientDs = patientDs)$root, TRUE)
    
    if(inherits(Root, "try-error")){
      if(tries >= 6){
        return(NA)
      }else{
        Up = Up + 1
      }
    }else{
      return(Root)
    }
  }
}

survivalFunc <- function (t, dsId, patientId) {
  exp(-integrate(hazardFunc, lower = 0, upper = t, dsId, patientId)$value)
}

invSurvival <- function (t, u, dsId, patientId) {
  log(u) + integrate(hazardFunc, lower = 0, upper = t, dsId, patientId)$value
}

hazardFunc = function (s, dsId, patientId) {
  weibullShape = simulatedDsList[[dsId]]$weibullShape
  weibullScale = simulatedDsList[[dsId]]$weibullScale
  Age = simulatedDsList[[dsId]]$Age[patientId]
  b = simulatedDsList[[dsId]]$b[patientId, ]
  wGamma = simulatedDsList[[dsId]]$wGamma[patientId]
  
  # pdf_s = dweibull(s, shape = weibullShape, scale = weibullScale)
  # survival_s = (1-pweibull(q = s, shape= weibullShape, scale = weibullScale))
  # baselinehazard_s = pdf_s/survival_s
  
  baselinehazard_s = (weibullShape/weibullScale)*(s/weibullScale)^(weibullShape-1)
  
  df_s = data.frame(visitTimeYears = s, Age = Age)
  
  xi_s_val = model.matrix(fixedValueFormula, df_s)
  xi_s_slope = model.matrix(fixedSlopeFormula, df_s)
  
  zi_s_val = model.matrix(randomValueFormula, df_s)
  zi_s_slope = model.matrix(randomSlopeFormula, df_s)
  
  zib_val = zi_s_val %*% b
  zib_slope = zi_s_slope %*% b[-1] #-1 to ignore intercept
  
  xBetaZb_s_value = xi_s_val %*% betas + zib_val
  xBetaZb_s_slope = xi_s_slope %*% betas[-c(1:3)] + zib_slope #-c(1:3) to ignore intercept, age and age^2
  
  y_Alpha = cbind(xBetaZb_s_value, xBetaZb_s_slope) %*% alphas
  #y_Alpha = cbind(xBetaZb_s_value) %*% alphas[1]
  
  baselinehazard_s * exp(wGamma + y_Alpha)
}

pSurvTime = function(survProb, dsId, patientId){
  Low = 1e-05
  Up <- 15
  tries  = 0
  
  repeat{
    tries = tries + 1
    Root <- try(uniroot(invSurvival, interval = c(Low, Up), 
                        u = survProb, dsId, patientId)$root, TRUE)
    
    if(inherits(Root, "try-error")){
      if(tries >= 5){
        return(NA)
      }else{
        Up = Up + 1    
      }
    }else{
      return(Root)
    }
  }
}


rLogPSA1 =  function(dsId, patientId, time){
  simDs.id = simulatedDsList[[dsId]]$simDs.id
  b = simulatedDsList[[dsId]]$b
  
  df_s = data.frame(visitTimeYears = time, Age = simDs.id$Age[patientId])
  
  xi_s_val = model.matrix(fixedValueFormula, df_s)
  zi_s_val = model.matrix(randomValueFormula, df_s)
  zib_val = zi_s_val %*% b[patientId, ]
  xBetaZb_s_value = xi_s_val %*% betas + zib_val
  
  sapply(xBetaZb_s_value, rnorm, n=1, sigma.y)
}

dynamicPredProb = function(futureTimes, dsId, patientDs, last.time=NULL){
  temp = survfitJM(simulatedDsList[[dsId]]$models$simJointModel_replaced, patientDs, 
                   idVar="P_ID", last.time=last.time, 
                   survTimes = futureTimes)$summaries[[1]][, "Median"]
  return(temp)
}

dynamicPredProbTimesTdiff = function(futureTimes, dsId, patientDs, last.time=NULL){
  temp = (futureTimes-last.time) * survfitJM(simulatedDsList[[dsId]]$models$simJointModel_replaced, patientDs, 
                   idVar="P_ID", last.time=last.time, 
                   survTimes = futureTimes)$summaries[[1]][, "Median"]
  return(temp)
}

expectedCondFailureTime = function(dsId, patientDs, last.time, upperLimitIntegral = 15){
  last.time + integrate(dynamicPredProb, last.time, upperLimitIntegral, dsId, patientDs, 
                            last.time = last.time,
                            abs.tol = 0.1)$value
}

varianceCondFailureTime = function(dsId, patientDs, last.time, upperLimitIntegral = 15){
  
  
  lhs = 2 * integrate(dynamicPredProbTimesTdiff, last.time, upperLimitIntegral, dsId, patientDs, 
                      last.time = last.time,
                      abs.tol = 0.1)$value
  
  rhs = (integrate(dynamicPredProb, last.time, upperLimitIntegral, dsId, patientDs, 
                   last.time = last.time,
                   abs.tol = 0.1)$value)^2
  
  return(lhs-rhs)
}


generateLongtiudinalTimeBySchedule = function(){
  #months = c(seq(0,24, 1), 27, 30, 33, 36, 39, 42, 48, 54, 60, 66, 72, 78, 84,
  #           90, 96, 102, 108, 114, 120)
  
  months = c(0, 3, 6, 9, 12, 15, 18, 21, 24, 30, 36, 42, 48, 54, 60, 66, 72, 78, 84,
              90, 96, 102, 108, 114, 120, 126, 132, 138, 144, 150, 156, 162, 168, 174, 180)
  
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

getBsGammas = function(jointModel, weightedOnes = T){
  if(weightedOnes){
    jointModel$statistics$postwMeans$Bs_gammas
  }else{
    jointModel$statistics$postMeans$Bs_gammas
  }
}

generateSimLongitudinalData = function(nSub){
  Age <- rnorm(n = nSub, mean=mean(prias.id$Age), sd=sqrt(var(prias.id$Age)))
  subId <- rep(1:nSub)
  simDs.id = data.frame(P_ID = subId, Age = Age)
  simDs = data.frame(P_ID = rep(subId, each=timesPerSubject), 
                     Age = rep(Age, each=timesPerSubject),
                     visitTimeYears = longTimes,
                     visitNumber = rep(1:timesPerSubject, nSub))
  
  X = model.matrix(fixedValueFormula, data = simDs)
  Z = model.matrix(randomValueFormula, data = simDs)
  
  b <- mvrnorm(nSub, mu = rep(0, nrow(D)), D)
  
  Zb = unlist(lapply(1:nSub,function(i){
    Z[((i-1)*timesPerSubject + 1):(i*timesPerSubject),] %*% b[i, ]
  }))
  
  xBetaZb = X %*% betas + Zb
  simDs$logpsa1 = rnorm(nSub * timesPerSubject, xBetaZb, sigma.y)
  
  W <- model.matrix(survivalFormula, data = simDs.id)
  wGamma <- as.vector(W %*% gammas)
  
  list(simDs = simDs, simDs.id = simDs.id, Age = Age, b=b, wGamma=wGamma)
}

generateSimJointData = function(dsId, nSub){
  u <- runif(nSub)
  
  simulatedDsList[[dsId]]$simDs.id$progression_time <- NA
  
  simulatedDsList[[dsId]]$simDs.id$progression_time = sapply(1:nSub, function(i){
                                        pSurvTime(u[i], dsId, i)
                                      })
  percentageRejected = sum(is.na(simulatedDsList[[dsId]]$simDs.id$progression_time[1:nSub]))/nSub
  
  if(percentageRejected > 0.2){
    stop("Too many NA's sampled")
  }
  
  pid_to_keep = simulatedDsList[[dsId]]$simDs.id[!is.na(simulatedDsList[[dsId]]$simDs.id$progression_time),]$P_ID
  
  simulatedDsList[[dsId]]$simDs = simulatedDsList[[dsId]]$simDs[simulatedDsList[[dsId]]$simDs$P_ID %in% pid_to_keep,]
  simulatedDsList[[dsId]]$simDs.id = simulatedDsList[[dsId]]$simDs.id[simulatedDsList[[dsId]]$simDs.id$P_ID %in% pid_to_keep,]
  simulatedDsList[[dsId]]$b = simulatedDsList[[dsId]]$b[pid_to_keep,]
  simulatedDsList[[dsId]]$wGamma = simulatedDsList[[dsId]]$wGamma[pid_to_keep]
  simulatedDsList[[dsId]]$Age = simulatedDsList[[dsId]]$Age[pid_to_keep]
  
  simulatedDsList[[dsId]]$simDs.id$P_ID = 1:nrow(simulatedDsList[[dsId]]$simDs.id)
  simulatedDsList[[dsId]]$simDs$P_ID = rep(simulatedDsList[[dsId]]$simDs.id$P_ID, each=timesPerSubject)
  
  #Divide into training and test
  trainingSize = round(nrow(simulatedDsList[[dsId]]$simDs.id)*0.75)
  trainingDs.id = simulatedDsList[[dsId]]$simDs.id[1:trainingSize, ]
  testDs.id = simulatedDsList[[dsId]]$simDs.id[(trainingSize+1):nrow(simulatedDsList[[dsId]]$simDs.id),]
  trainingDs = simulatedDsList[[dsId]]$simDs[simulatedDsList[[dsId]]$simDs$P_ID %in% trainingDs.id$P_ID, ]
  testDs = simulatedDsList[[dsId]]$simDs[simulatedDsList[[dsId]]$simDs$P_ID %in% testDs.id$P_ID, ]
  
  # simulate censoring times from an exponential distribution for TRAINING DATA SET ONLY
  # and calculate the observed event times, i.e., min(true event times, censoring times)
  
  #Ctimes <- rexp(trainingSize, 1/mean(prias.id[prias.id$progressed==0,]$progression_time))
  Ctimes<-runif(trainingSize, 3, 15)
  
  trainingDs.id$progressed = trainingDs.id$progression_time <= Ctimes
  trainingDs.id$progression_time = pmin(trainingDs.id$progression_time, Ctimes)
  trainingDs$progression_time = rep(trainingDs.id$progression_time, each=timesPerSubject)
  trainingDs$progressed = rep(trainingDs.id$progressed, each=timesPerSubject)
  testDs$progression_time = rep(testDs.id$progression_time, each=timesPerSubject)
  
  # drop the longitudinal measurementsthat were taken after the observed event time for each subject.
  trainingDs = trainingDs[trainingDs$visitTimeYears <= trainingDs$progression_time, ]
  
  list(simDs = simulatedDsList[[dsId]]$simDs, simDs.id = simulatedDsList[[dsId]]$simDs.id, 
       trainingDs = trainingDs, trainingDs.id = trainingDs.id, 
       testDs = testDs, testDs.id = testDs.id, 
       b=simulatedDsList[[dsId]]$b, wGamma=simulatedDsList[[dsId]]$wGamma, 
       weibullShape = simulatedDsList[[dsId]]$weibullShape, weibullScale = simulatedDsList[[dsId]]$weibullScale,
       Age = simulatedDsList[[dsId]]$Age, percentageRejected = percentageRejected)
}

fitJointModelSimDs = function(trainingDs.id, trainingDs){
  cox_Model_training = coxph(Surv(progression_time, progressed) ~ I(Age - 70) +  I((Age - 70)^2),
                             data = trainingDs.id, x=T, model=T)
  
  mvglmer_psa_training=mvglmer(list(logpsa1 ~  I(Age - 70) +  I((Age - 70)^2) + 
                                      ns(visitTimeYears, knots=c(0.5, 1.2, 2.5), Boundary.knots=c(0, 7)) + 
                                      (ns(visitTimeYears, knots=c(0.5), Boundary.knots=c(0, 7))|P_ID)), 
                               data = trainingDs, families = list(gaussian))
  
  forms_psa_training <- list("logpsa1" = "value",
                             "logpsa1" = list(fixed = ~ 0 + dns(visitTimeYears, knots=c(0.5, 1.2, 2.5), Boundary.knots=c(0, 7)),
                                              random=~0 + dns(visitTimeYears, knots=c(0.5), Boundary.knots=c(0, 7)),
                                              indFixed = 4:7, indRandom=2:3, name = "slope"))
  
  # mvJoint_psa_tdval_training=mvJointModelBayes(mvglmer_psa_training, cox_Model_training, timeVar = "visitTimeYears",
  #                                              priors = list(shrink_gammas = TRUE, shrink_alphas = TRUE))
  # 
  # mvJoint_psa_tdboth_training <- update(mvJoint_psa_tdval_training, Formulas = forms_psa_training,
  #                                       priors = list(shrink_gammas = TRUE, shrink_alphas = TRUE))
  #other way
  mvJoint_psa_tdval_training = NULL
  mvJoint_psa_tdboth_training <- mvJointModelBayes(mvglmer_psa_training, cox_Model_training, timeVar = "visitTimeYears",
                                                   Formulas = forms_psa_training,
                                                   priors = list(shrink_gammas = TRUE, shrink_alphas = TRUE))
  
  ###########################################
  #Also fit joint model using the jointModelAPI
  ###########################################
  lme_psa_training = lme(data=trainingDs,
                         fixed=logpsa1 ~  I(Age - 70) +  I((Age - 70)^2) + 
                           ns(visitTimeYears, knots=c(0.5, 1.2, 2.5), Boundary.knots=c(0, 7)),
                         random = ~ns(visitTimeYears, knots=c(0.5), Boundary.knots=c(0, 7))|P_ID,
                         control = lmeControl(opt = "optim"), method="ML")
  
  jmbayes_psa_tdval_training = jointModelBayes(lmeObject = lme_psa_training, survObject = cox_Model_training, 
                                               timeVar = "visitTimeYears", control = list(n.iter=1000))
  jmbayes_psa_tdboth_training = update(jmbayes_psa_tdval_training, param = "td-both", 
                                       extraForm = list(fixed = ~ 0 + dns(visitTimeYears, knots=c(0.5, 1.2, 2.5), Boundary.knots=c(0, 7)),
                                                        random = ~0 + dns(visitTimeYears, knots=c(0.5), Boundary.knots=c(0, 7)), 
                                                        indFixed = 4:7, indRandom=2:3))
  
  simJointModel_replaced = replaceMCMCContents(mvJoint_psa_tdboth_training, jmbayes_psa_tdboth_training)
  
  list(mvJoint_psa_tdval_training = mvJoint_psa_tdval_training, 
       mvJoint_psa_tdboth_training = mvJoint_psa_tdboth_training,
       cox_Model_training = cox_Model_training,
       mvglmer_psa_training = mvglmer_psa_training,
       simJointModel_replaced = simJointModel_replaced)
}

###############################################
# parameters for the linear mixed effects model
###############################################
fixedValueFormula = ~ 1 + I(Age - 70) + I((Age - 70)^2) + ns(visitTimeYears, knots = c(0.5, 1.2, 2.5), Boundary.knots = c(0, 7))
randomValueFormula = ~ 1 + ns(visitTimeYears, knots = c(0.5), Boundary.knots = c(0, 7))

fixedSlopeFormula = ~ 0 + dns(visitTimeYears, knots = c(0.5, 1.2, 2.5), Boundary.knots = c(0, 7))
randomSlopeFormula = ~ 0 + dns(visitTimeYears, knots = c(0.5), Boundary.knots = c(0, 7))

longTimes <- c(replicate(nSub, generateLongtiudinalTimeBySchedule(), simplify = T)) # at which time points longitudinal measurements are supposed to be taken
timesPerSubject = length(longTimes) / nSub

fittedJointModel = mvJoint_psa_spline_pt51pt22pt5_pt5_tdboth

betas <- getBetas(fittedJointModel, weightedOnes = F)
sigma.y <- getSigma(fittedJointModel, weightedOnes = F)
D <- getD(fittedJointModel, weightedOnes = F)
survivalFormula = ~ 0 + I(Age - 70) + I((Age - 70)^2)
gammas = getGamma(fittedJointModel, weightedOnes = F)
alphas = getAlpha(fittedJointModel, weightedOnes = F)
