expectedCondFailureTime = function(patientDs, upperLimitIntegral = 100){
  
  dynamicPredProb = function(futureTimes, patientDs){
    temp = survfitJM(simJointModel_replaced, patientDs, idVar="P_ID", survTimes = futureTimes)$summaries[[1]][, "Mean"]
    return(temp)
  }
  
  lastVisitTime = max(patientDs$visitTimeYears)
  lastVisitTime + integrate(dynamicPredProb, lastVisitTime, upperLimitIntegral, patientDs)$value
}

