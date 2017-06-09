dynamicPredProb = function(futureTimes, patientDs){
  temp = survfitJM(joint_psa_replaced, patientDs,
                   idVar="P_ID", survTimes = futureTimes)$summaries[[1]][, "Median"]
  return(temp)
}

dynamicPredProbTimesTdiff = function(futureTimes, patientDs){
  temp = (futureTimes- max(patientDs$visitTimeYears)) * survfitJM(joint_psa_replaced, patientDs,
                                                                  idVar="P_ID", survTimes = futureTimes)$summaries[[1]][, "Median"]
  return(temp)
}

# dynamicPredProb = function(futureTimes, patientDs){
#   survProbNearest = sapply(futureTimes, function(futureTime){
#     probs[which(abs(prob_times-futureTime)==min(abs(prob_times-futureTime)))[1]]
#   })
# 
#   return(survProbNearest)
# }
# 
# dynamicPredProbTimesTdiff = function(futureTimes, patientDs){
#   survProbNearest = sapply(futureTimes, function(futureTime){
#     probs[which(abs(prob_times-futureTime)==min(abs(prob_times-futureTime)))[1]]
#   })
#   temp = (futureTimes- max(patientDs$visitTimeYears)) * survProbNearest
#   return(temp)
# }


calcSurvProbs = function(patientDs){
  probs = c(1,survfitJM(joint_psa_replaced, patientDs, 
            idVar="P_ID", survTimes = prob_times)$summaries[[1]][, "Median"])
  return(probs)
}

expectedCondFailureTime = function(patientDs, upperLimitIntegral = 20){
  max(patientDs$visitTimeYears) + integrate(dynamicPredProb, lower=max(patientDs$visitTimeYears),
                                            upper=upperLimitIntegral, patientDs, rel.tol = 0.05)$value
}

varianceCondFailureTime = function(patientDs, upperLimitIntegral = 20){
  
  lhs = 2 * integrate(dynamicPredProbTimesTdiff, max(patientDs$visitTimeYears), 
                      upperLimitIntegral, patientDs, abs.tol = 0.1)$value
  
  rhs = (integrate(dynamicPredProb, max(patientDs$visitTimeYears), upperLimitIntegral, 
                   patientDs, abs.tol = 0.1)$value)^2
  
  return(lhs-rhs)
}

tStart = Sys.time()
#maxTime = max(ND$visitTimeYears)
#prob_times = seq(maxTime, 20, 0.1)
#probs = calcSurvProbs(ND)
eTime = expectedCondFailureTime(ND, upperLimitIntegral = 100)
#varTime = varianceCondFailureTime(ND)
tEnd=Sys.time()
print(tEnd-tStart)