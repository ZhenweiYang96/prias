
computeBiopsyTimes_temp = function(minVisits = 5, dsId, patientRowNum){
  
  cutoffValues = simulatedDsList[[dsId]]$cutoffValues
  
  patientDs_i = simulatedDsList[[dsId]]$biopsyTimes[[patientRowNum]]
  
  patientDs_i$survTimeAccuracy = NA
  patientDs_i$survTimeF1Score = NA
  
  #instead of times per subject I choose 28 because 28th time is 11 years.
  #expected time of failure is not available for last time point if faiure time for all subejccts is less than tht last time point
  #training max progression time is 11.216 and test is 10. something
  visitsOfInterest = minVisits:30
  res = foreach(j=visitsOfInterest, .packages = c("splines", "JMbayes", "coda"),
                .export=c("timesPerSubject", "dynamicCutOffTimes",
                          "expectedCondFailureTime", "dynamicPredProb",
                          "simulatedDsList", "pDynSurvTime", "invDynSurvival")) %dopar%{
                            persTestDs = patientDs_i[patientDs_i$visitNumber <= j,]
                            temp_lasttime = max(persTestDs$visitTimeYears)
                            nearest_time_index = which(abs(dynamicCutOffTimes-temp_lasttime)==min(abs(dynamicCutOffTimes-temp_lasttime)))
                            nearest_time_index = nearest_time_index[1]
                            
                            survTimeAccuracy = pDynSurvTime(survProb = cutoffValues[[nearest_time_index]]["accuracy"], dsId, persTestDs)
                            survTimeF1Score = pDynSurvTime(survProb =  cutoffValues[[nearest_time_index]]["f1score"], dsId, persTestDs)
                            
                            list(survTimeAccuracy, survTimeF1Score)
                          }
  
  patientDs_i$survTimeAccuracy[visitsOfInterest] = sapply(res, function(x){x[[1]]}, simplify = T)
  patientDs_i$survTimeF1Score[visitsOfInterest] = sapply(res, function(x){x[[2]]}, simplify = T)
  
  return(patientDs_i)
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
  
  #other way
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
  
  list(mvJoint_psa_tdboth_training = mvJoint_psa_tdboth_training,
       cox_Model_training = cox_Model_training,
       mvglmer_psa_training = mvglmer_psa_training,
       simJointModel_replaced = simJointModel_replaced)
}
