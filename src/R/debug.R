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

############################
ct= makeCluster(4)
registerDoParallel(ct)

testDs = simulatedDsList[[1]]$testDs
patientDsList = split(testDs, testDs$P_ID)

tStart = Sys.time()
tradRes = foreach(i=1:50, .packages = c("splines", "JMbayes", "coda"),
        .export=c("timesPerSubject", "dynamicCutOffTimes",
                  "expectedCondFailureTime", "dynamicPredProb",
                  "simulatedDsList", "pDynSurvTime", "invDynSurvival")) %dopar%{
                    patientDs_i = patientDsList[[i]][1:10,]
                    
                    list(expectedCondFailureTime(1, patientDs_i),
                    pDynSurvTime(survProb = 0.85, 1, patientDs_i),
                    pDynSurvTime(survProb = 0.93, 1, patientDs_i),
                    pDynSurvTime(survProb = 0.65, 1, patientDs_i),
                    pDynSurvTime(survProb = 0.34, 1, patientDs_i),
                    pDynSurvTime(survProb = 0.45, 1, patientDs_i),
                    pDynSurvTime(survProb = 0.27, 1, patientDs_i))
                  }

stopCluster(ct)
tEnd = Sys.time()
tEnd-tStart

#######################

ct= makeCluster(4)
registerDoParallel(ct)

testDs = simulatedDsList[[1]]$testDs
patientDsList = split(testDs, testDs$P_ID)

tStart = Sys.time()
forRes = foreach(i=1:50, .packages = c("splines", "JMbayes", "coda"),
        .export=c("timesPerSubject", "dynamicCutOffTimes",
                  "expectedCondFailureTime", "dynamicPredProb",
                  "simulatedDsList", "pDynSurvTime", "invDynSurvival")) %dopar%{
                    patientDs_i = patientDsList[[i]][1:10,]
                    
                    #upperTime = pDynSurvTime(survProb = 0.27, 1, patientDs_i)
                    
                    
                    tt = survfitJM(simulatedDsList[[1]]$models$simJointModel_replaced, 
                              newdata = patientDs_i, idVar = "P_ID", 
                              survTimes = seq(max(patientDs_i$visitTimeYears), 15, 0.05))
                    
                    times = tt$summaries[[1]][, 1]
                    probs = tt$summaries[[1]][, 3]
                    
                    list(expectedCondFailureTime(1, patientDs_i),
                    times[which(abs(probs-0.85)==min(abs(probs-0.85)))[1]],
                    times[which(abs(probs-0.93)==min(abs(probs-0.93)))[1]],
                    times[which(abs(probs-0.65)==min(abs(probs-0.65)))[1]],
                    times[which(abs(probs-0.34)==min(abs(probs-0.34)))[1]],
                    times[which(abs(probs-0.45)==min(abs(probs-0.45)))[1]],
                    times[which(abs(probs-0.27)==min(abs(probs-0.27)))[1]])
                  }

stopCluster(ct)
tEnd = Sys.time()


#######################

ct= makeCluster(4)
registerDoParallel(ct)

testDs = simulatedDsList[[1]]$testDs
patientDsList = split(testDs, testDs$P_ID)

tStart = Sys.time()
forAdaptiveRes = foreach(i=1:50, .packages = c("splines", "JMbayes", "coda"),
                 .export=c("timesPerSubject", "dynamicCutOffTimes",
                           "expectedCondFailureTime", "dynamicPredProb",
                           "simulatedDsList", "pDynSurvTime", "invDynSurvival")) %dopar%{
                             patientDs_i = patientDsList[[i]][1:10,]
                             lastVisitTime = max(patientDs_i$visitTimeYears)
                             #upperTime = pDynSurvTime(survProb = 0.27, 1, patientDs_i)
                             
                             tt = survfitJM(simulatedDsList[[1]]$models$simJointModel_replaced, 
                                            newdata = patientDs_i, idVar = "P_ID", 
                                            survTimes = seq(lastVisitTime, 15, 0.05))
                             
                             times = tt$summaries[[1]][, 1]
                             probs = tt$summaries[[1]][, 3]
                             
                             expectedFailureTime = lastVisitTime + integrate(function(evalTimes){
                               unlist(lapply(evalTimes, function(t){
                                probs[which(abs(times-t)==min(abs(times-t)))[1]]
                               }))
                             }, lastVisitTime, 15, abs.tol = 0.1)$value
                             
                             list(expectedFailureTime,
                                  times[which(abs(probs-0.85)==min(abs(probs-0.85)))[1]],
                                  times[which(abs(probs-0.93)==min(abs(probs-0.93)))[1]],
                                  times[which(abs(probs-0.65)==min(abs(probs-0.65)))[1]],
                                  times[which(abs(probs-0.34)==min(abs(probs-0.34)))[1]],
                                  times[which(abs(probs-0.45)==min(abs(probs-0.45)))[1]],
                                  times[which(abs(probs-0.27)==min(abs(probs-0.27)))[1]])
                           }

stopCluster(ct)
tEnd = Sys.time()

#######################
ct= makeCluster(4)
registerDoParallel(ct)

testDs = simulatedDsList[[1]]$testDs
patientDsList = split(testDs, testDs$P_ID)

tStart = Sys.time()
forAdaptiveRes2 = foreach(i=1:50, .packages = c("splines", "JMbayes", "coda"),
                         .export=c("timesPerSubject", "dynamicCutOffTimes",
                                   "expectedCondFailureTime", "dynamicPredProb",
                                   "simulatedDsList", "pDynSurvTime", "invDynSurvival")) %dopar%{
                                     patientDs_i = patientDsList[[i]][1:10,]
                                     lastVisitTime = max(patientDs_i$visitTimeYears)
                                     
                                     round1 = survfitJM(simulatedDsList[[1]]$models$simJointModel_replaced, 
                                                    newdata = patientDs_i, idVar = "P_ID", 
                                                    survTimes = seq(lastVisitTime, 15, 1))
                                     
                                     times = round1$summaries[[1]][, 1]
                                     probs = round1$summaries[[1]][, 3]
                                     
                                     probsOfInterest = c(0.85, 0.93, 0.65, 0.34, 0.45, 0.27)
                                     newTimes = do.call(c, lapply(probsOfInterest, function(p){
                                       nearest_time = times[which(abs(probs-p)==min(abs(probs-p)))[1]]  
                                       seq(nearest_time-1, nearest_time+1, 0.05)
                                     }))
                                     
                                     round2 = survfitJM(simulatedDsList[[1]]$models$simJointModel_replaced, 
                                                        newdata = patientDs_i, idVar = "P_ID", 
                                                        survTimes = unique(newTimes))
                                     
                                     times = round2$summaries[[1]][, 1]
                                     probs = round2$summaries[[1]][, 3]
                                     
                                     list(expectedCondFailureTime(1, patientDs_i),
                                          times[which(abs(probs-0.85)==min(abs(probs-0.85)))[1]],
                                          times[which(abs(probs-0.93)==min(abs(probs-0.93)))[1]],
                                          times[which(abs(probs-0.65)==min(abs(probs-0.65)))[1]],
                                          times[which(abs(probs-0.34)==min(abs(probs-0.34)))[1]],
                                          times[which(abs(probs-0.45)==min(abs(probs-0.45)))[1]],
                                          times[which(abs(probs-0.27)==min(abs(probs-0.27)))[1]])
                                   }

stopCluster(ct)
tEnd = Sys.time()



#######################
ct= makeCluster(4)
registerDoParallel(ct)

testDs = simulatedDsList[[1]]$testDs
patientDsList = split(testDs, testDs$P_ID)

tStart = Sys.time()
forAdaptiveRes3 = foreach(i=1:50, .packages = c("splines", "JMbayes", "coda"),
                          .export=c("timesPerSubject", "dynamicCutOffTimes",
                                    "expectedCondFailureTime", "dynamicPredProb",
                                    "simulatedDsList", "pDynSurvTime", "invDynSurvival")) %dopar%{
                                      patientDs_i = patientDsList[[i]][1:10,]
                                      lastVisitTime = max(patientDs_i$visitTimeYears)
                                      
                                      round1 = survfitJM(simulatedDsList[[1]]$models$simJointModel_replaced, 
                                                         newdata = patientDs_i, idVar = "P_ID", 
                                                         survTimes = seq(lastVisitTime, 15, 1))
                                      
                                      times = round1$summaries[[1]][, 1]
                                      probs = round1$summaries[[1]][, 3]
                                      
                                      probsOfInterest = c(0.85, 0.93, 0.65, 0.34, 0.45, 0.27)
                                      newTimes = do.call(c, lapply(probsOfInterest, function(p){
                                        nearest_time = times[which(abs(probs-p)==min(abs(probs-p)))[1]]  
                                        seq(nearest_time-1, nearest_time+1, 0.05)
                                      }))
                                      
                                      round2 = survfitJM(simulatedDsList[[1]]$models$simJointModel_replaced, 
                                                         newdata = patientDs_i, idVar = "P_ID", 
                                                         survTimes = unique(newTimes))
                                      
                                      times = round2$summaries[[1]][, 1]
                                      probs = round2$summaries[[1]][, 3]
                                      
                                      expectedFailureTime = lastVisitTime + integrate(function(evalTimes){
                                        print(evalTimes)
                                        survprobs = unlist(lapply(evalTimes, function(t){
                                          neartime = times[which(abs(times-t)==min(abs(times-t)))[1]]
                                          print(neartime)
                                          if(abs(neartime-t) <= 0.05){
                                            return(probs[neartime])
                                          }else{
                                            NA
                                          }
                                        }))
                                        
                                        naTimes = evalTimes[is.na(survprobs)]
                                        if(length(naTimes)>0){
                                          survprobs[is.na(survprobs)] = survfitJM(simulatedDsList[[1]]$models$simJointModel_replaced, 
                                                    patientDs_i, idVar="P_ID", 
                                                    survTimes = naTimes)$summaries[[1]][, "Median"]
                                        }
                                        
                                        return(survprobs)
                                      }, lastVisitTime, 15, abs.tol = 0.1)$value
                                      
                                      list(expectedFailureTime,
                                           times[which(abs(probs-0.85)==min(abs(probs-0.85)))[1]],
                                           times[which(abs(probs-0.93)==min(abs(probs-0.93)))[1]],
                                           times[which(abs(probs-0.65)==min(abs(probs-0.65)))[1]],
                                           times[which(abs(probs-0.34)==min(abs(probs-0.34)))[1]],
                                           times[which(abs(probs-0.45)==min(abs(probs-0.45)))[1]],
                                           times[which(abs(probs-0.27)==min(abs(probs-0.27)))[1]])
                                    }

stopCluster(ct)
tEnd = Sys.time()